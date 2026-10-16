"""Train the B9 neural student with a K/R-gain consistency term.

The model uses native residue, alternate residue, side, and distance. Rule
contrasts are exact linear functions of its 20 x 20 x 2 x 60 prediction table and
the report card's gene-by-native-residue weights. The pilot stage evaluates on
validation clusters, and the full stage evaluates on test clusters. Derived
AlphaMissense weights are not published.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import torch

import form_ab_lib as lib
from run_gate1 import read_topology

CELLS = 20 * 20 * 2 * 60
GRID = {"c1": (0.0, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0), "c2_protected": 100.0}
TRAINING = {"hidden": 32, "embedding": 8, "steps": 2500, "tune_steps": 1000, "lr": 3e-3, "tune_lr": 1e-3,
            "threads": 8}
OPERATING_FRACTION = 0.9
TARGET_RULE = "positive_charge_gain"
STAGES = {"pilot": dict(eval="validation", out=lib.ROOT / "results/form_ab_pilot_20261008/form_b9"),
          "full": dict(eval="test", out=lib.ROOT / "results/form_b_20261008/b9")}


def cell_ids(ref, alt, side, distance):
    return ((ref * 20 + alt) * 2 + side) * 60 + (distance - 1)


def aggregate(panel: dict) -> dict:
    rows, alts = np.nonzero(np.isfinite(panel["scores"]))
    index = cell_ids(panel["ref_aa"][rows], alts, panel["side"][rows], panel["distance"][rows])
    values = panel["scores"][rows, alts]
    return {"n": np.bincount(index, minlength=CELLS).astype(float),
            "sum": np.bincount(index, weights=values, minlength=CELLS),
            "sumsq": np.bincount(index, weights=values ** 2, minlength=CELLS)}


def rule_masks(tool) -> dict:
    """Numerator and comparator alternatives per native residue, exactly as the report card builds them."""
    out = {}
    for rule, (kind, group) in tool.RULES.items():
        num, den = np.zeros((20, 20), bool), np.zeros((20, 20), bool)
        for r, ref in enumerate(lib.ALPHABET):
            for a, alt in enumerate(lib.ALPHABET):
                if a == r:
                    continue
                if kind == "gain" and ref not in tool.CHARGED:
                    num[r, a] = alt in group
                    den[r, a] = alt not in tool.CHARGED and alt not in group
                elif kind == "loss" and ref in group:
                    num[r, a] = alt not in tool.CHARGED and alt not in group
                    den[r, a] = alt in group
        out[rule] = (num, den, num.any(1) & den.any(1))
    return out


def coefficients(window: dict, gene: np.ndarray, eligible_ref: np.ndarray) -> np.ndarray:
    """Weights A[ref, side, distance-1] with contrast = sum(A * penalty_table) over (gene, ref) strata."""
    keep = eligible_ref[window["ref_aa"]]
    ref, side, dist = window["ref_aa"][keep], window["side"][keep], window["distance"][keep]
    _, inverse = np.unique(gene[keep] * 20 + ref, return_inverse=True)
    size = inverse.max() + 1 if len(inverse) else 0
    n1, n0 = np.bincount(inverse[side == 0], minlength=size), np.bincount(inverse[side == 1], minlength=size)
    shared = (n1 > 0) & (n0 > 0)
    weight = np.where(shared, n1 * n0 / np.maximum(n1 + n0, 1), 0.0)
    total = weight.sum()
    contribution = np.where(side == 0, weight[inverse] / np.maximum(n1[inverse], 1),
                            -weight[inverse] / np.maximum(n0[inverse], 1)) / total
    flat = (ref * 2 + side) * 60 + (dist - 1)
    return np.bincount(flat, weights=contribution, minlength=20 * 2 * 60).reshape(20, 2, 60)


def penalty_table(table, num, den, eligible):
    top = (table * num[:, :, None, None]).sum(1) / num.sum(1).clamp(min=1)[:, None, None]
    bottom = (table * den[:, :, None, None]).sum(1) / den.sum(1).clamp(min=1)[:, None, None]
    return (top - bottom) * eligible[:, None, None]


class Student(torch.nn.Module):
    def __init__(self, hidden: int, embedding: int):
        super().__init__()
        self.ref, self.alt = torch.nn.Embedding(20, embedding), torch.nn.Embedding(20, embedding)
        self.net = torch.nn.Sequential(torch.nn.Linear(2 * embedding + 2, hidden), torch.nn.SiLU(),
                                       torch.nn.Linear(hidden, hidden), torch.nn.SiLU(),
                                       torch.nn.Linear(hidden, 1))
        r, a, s, d = np.meshgrid(np.arange(20), np.arange(20), np.arange(2), np.arange(60), indexing="ij")
        self.register_buffer("r", torch.tensor(r.ravel()))
        self.register_buffer("a", torch.tensor(a.ravel()))
        self.register_buffer("side", torch.tensor(s.ravel(), dtype=torch.float64))
        self.register_buffer("dist", torch.tensor(d.ravel() / 60.0))

    def forward(self) -> torch.Tensor:
        x = torch.cat([self.ref(self.r), self.alt(self.a), self.side[:, None], self.dist[:, None]], dim=1)
        return self.net(x).squeeze(1).view(20, 20, 2, 60)


class Problem:
    """Training-cluster cells, per-rule weights and the baseline contrasts a protected student must keep."""

    def __init__(self, tool, cells, coeffs, masks, sd):
        self.n = torch.tensor(cells["n"])
        self.mu = torch.tensor(np.divide(cells["sum"], cells["n"], out=np.zeros(CELLS), where=cells["n"] > 0))
        total = cells["n"].sum()
        self.variance = float(cells["sumsq"].sum() / total - (cells["sum"].sum() / total) ** 2)
        self.coeffs = {r: torch.tensor(c) for r, c in coeffs.items()}
        self.masks = {r: tuple(torch.tensor(x) if i < 2 else torch.tensor(x, dtype=torch.float64)
                               for i, x in enumerate((m[0].astype(float), m[1].astype(float), m[2])))
                      for r, m in masks.items()}
        self.sd, self.target = sd, -lib.CORRECTION["margin"] * sd
        self.rules = list(coeffs)

    def contrasts(self, table) -> dict:
        return {r: (self.coeffs[r] * penalty_table(table, *self.masks[r])).sum() for r in self.rules}

    def reducible_mse(self, table) -> torch.Tensor:
        return (self.n * (table.reshape(-1) - self.mu) ** 2).sum() / self.n.sum() / self.variance


def fit(problem, model, c1, c2, baseline, steps, lr):
    optimiser = torch.optim.Adam(model.parameters(), lr=lr)
    for _ in range(steps):
        optimiser.zero_grad()
        table = model()
        found = problem.contrasts(table)
        loss = problem.reducible_mse(table)
        if c1 > 0:
            loss = loss + c1 * torch.relu((found[TARGET_RULE] - problem.target) / problem.sd) ** 2
        if c2 > 0:
            loss = loss + c2 * sum(((found[r] - baseline[r]) / problem.sd) ** 2
                                   for r in problem.rules if r != TARGET_RULE)
        loss.backward()
        optimiser.step()
    return model


def snapshot(model) -> np.ndarray:
    with torch.no_grad():
        return model().numpy().copy()


def contrasts_numpy(table: np.ndarray, coeffs: dict, masks: dict) -> dict:
    out = {}
    for rule, (num, den, eligible) in masks.items():
        top = (table * num[:, :, None, None]).sum(1) / np.maximum(num.sum(1), 1)[:, None, None]
        bottom = (table * den[:, :, None, None]).sum(1) / np.maximum(den.sum(1), 1)[:, None, None]
        out[rule] = float((coeffs[rule] * (top - bottom) * eligible[:, None, None]).sum())
    return out


def student_scores(table: np.ndarray, panel: dict) -> np.ndarray:
    ref, side, dist = panel["ref_aa"], panel["side"], panel["distance"]
    out = np.full((len(ref), 20), np.nan)
    for a in range(20):
        out[:, a] = table[ref, a, side, dist - 1]
    out[np.arange(len(ref)), ref] = np.nan
    return out


def gene_codes(census: dict, panel: dict, sites: dict) -> np.ndarray:
    genes = np.array([census[str(a)].get("gene", str(a)) for a in panel["accessions"]])
    return np.unique(genes[sites["protein"]], return_inverse=True)[1]


def run(args) -> None:
    from scipy.stats import spearmanr
    if args.stage == "full":
        form_b = lib.ROOT / "results/form_b_20261008/criteria.json"
        form_b_passed = (form_b.exists()
                         and json.loads(form_b.read_text(encoding="utf-8"))
                         .get("all_gated_criteria_pass") is True)
        if not form_b_passed:
            raise SystemExit("B9 follows a clean B2 result: every gated Form B criterion must pass first")
    args.out.mkdir(parents=True, exist_ok=True)
    torch.manual_seed(lib.MASTER_SEED if args.stage == "full" else lib.PILOT_SEED)
    torch.set_num_threads(TRAINING["threads"])
    torch.set_default_dtype(torch.float64)
    tool, tool_sha = lib.report_card_tool()
    census = read_topology(lib.TOPOLOGY)
    split = lib.cluster_split(census)
    label = split["label"]
    panel = lib.P.load_panel(lib.AM_ARCHIVE, 60)
    train = lib.subset_panel(panel, lib.split_mask(panel, label, "train"))
    held = lib.subset_panel(panel, lib.split_mask(panel, label, args.eval))
    masks = rule_masks(tool)
    sd = lib.SD_FIXED["alphamissense"]
    coeffs = {}
    for name, part in (("train", train), ("held", held)):
        window = lib.window_panel(part)
        gene = gene_codes(census, part, window)
        coeffs[name] = {rule: coefficients(window, gene, masks[rule][2]) for rule in masks}
    problem = Problem(tool, aggregate(train), coeffs["train"], masks, sd)
    held_cells = aggregate(held)
    mu_held = torch.tensor(np.divide(held_cells["sum"], held_cells["n"], out=np.zeros(CELLS),
                                     where=held_cells["n"] > 0))

    plain = Student(TRAINING["hidden"], TRAINING["embedding"])
    fit(problem, plain, 0.0, 0.0, None, TRAINING["steps"], TRAINING["lr"])
    base_table = snapshot(plain)
    baseline = {r: float(v) for r, v in problem.contrasts(torch.tensor(base_table)).items()}
    baseline_t = {r: torch.tensor(v) for r, v in baseline.items()}
    held_rows = np.isfinite(held["scores"])
    ref_rows, alt_rows = np.nonzero(held_rows)
    cell_rows = cell_ids(held["ref_aa"][ref_rows], alt_rows, held["side"][ref_rows], held["distance"][ref_rows])
    truth = held["scores"][ref_rows, alt_rows]

    def evaluate(table):
        held_contrast = contrasts_numpy(table, coeffs["held"], masks)
        train_contrast = contrasts_numpy(table, coeffs["train"], masks)
        return {"train_contrast": train_contrast, "held_contrast": held_contrast,
                "train_reducible_mse": float(problem.reducible_mse(torch.tensor(table))),
                "held_reducible_mse": float(
                    (torch.tensor(held_cells["n"])
                     * (torch.tensor(table).reshape(-1) - mu_held) ** 2).sum()
                    / held_cells["n"].sum() / problem.variance),
                "held_spearman_vs_alphamissense": float(spearmanr(table.reshape(-1)[cell_rows], truth).statistic)}

    curve = []
    models = {"plain": evaluate(base_table)}
    tables = {"plain": base_table}
    for protected, c2 in (("protected", GRID["c2_protected"]), ("unprotected", 0.0)):
        for c1 in GRID["c1"]:
            if c1 == 0:
                continue
            model = Student(TRAINING["hidden"], TRAINING["embedding"])
            model.load_state_dict(plain.state_dict())
            fit(problem, model, c1, c2, baseline_t, TRAINING["tune_steps"], TRAINING["tune_lr"])
            table = snapshot(model)
            result = evaluate(table)
            others = [r for r in problem.rules if r != TARGET_RULE]
            curve.append({"series": protected, "c1": c1, "c2": c2,
                          "train_target_sd_units": result["train_contrast"][TARGET_RULE] / sd,
                          "held_target_sd_units": result["held_contrast"][TARGET_RULE] / sd,
                          "train_reducible_mse": result["train_reducible_mse"],
                          "held_spearman": result["held_spearman_vs_alphamissense"],
                          "spearman_loss_vs_plain": models["plain"]["held_spearman_vs_alphamissense"]
                          - result["held_spearman_vs_alphamissense"],
                          "max_collateral_sd_units": max(
                              abs(result["held_contrast"][r]
                                  - models["plain"]["held_contrast"][r])
                              for r in others) / sd})
            tables[(protected, c1)] = table
            models[(protected, c1)] = result
    chosen = next((row for row in curve if row["series"] == "protected"
                   and row["train_target_sd_units"] <= OPERATING_FRACTION * (-lib.CORRECTION["margin"])), None)

    cards, report = {}, {}
    if chosen is not None:
        picks = {"plain": tables["plain"], "protected": tables[("protected", chosen["c1"])],
                 "unprotected": tables[("unprotected", chosen["c1"])]}
        for name, table in picks.items():
            out = args.out / f"report_card_{name}"
            out.mkdir(parents=True, exist_ok=True)
            cards[name] = lib.run_report_card(tool, held, student_scores(table, held), out)
        rules = {n: c["rule_audit"]["within_10_residues"]["rules"] for n, c in cards.items()}
        functional = contrasts_numpy(picks["plain"], coeffs["held"], masks)
        report = {
            "chosen_c1": chosen["c1"],
            "functional_vs_tool_max_abs_error_plain": max(abs(functional[r] - rules["plain"][r]["estimate"])
                                                          for r in functional),
            "target_after": {n: rules[n][TARGET_RULE] for n in rules},
            "verdict_changes_vs_plain": {
                n: {r: [rules["plain"][r]["verdict"], rules[n][r]["verdict"]] for r in rules[n]
                    if rules[n][r]["verdict"] != rules["plain"][r]["verdict"]} for n in ("protected", "unprotected")},
            "max_collateral_tool_sd_units": {
                n: max(abs(rules[n][r]["estimate"] - rules["plain"][r]["estimate"])
                       for r in rules[n] if r != TARGET_RULE) / sd for n in ("protected", "unprotected")}}
        target = rules["protected"][TARGET_RULE]
        spearman_loss = chosen["spearman_loss_vs_plain"]
        report["criterion_b10_5"] = gate_b10_5(target, sd, spearman_loss,
                                               report["verdict_changes_vs_plain"]["protected"],
                                               report["functional_vs_tool_max_abs_error_plain"])
    summary = {"stage": args.stage, "as_of_utc": datetime.now(timezone.utc).isoformat(), "evaluated_on": args.eval,
               "grid": GRID, "training": TRAINING, "operating_rule": f"smallest protected c1 with train contrast <= "
               f"{OPERATING_FRACTION} x target", "plain": models["plain"], "baseline_train_contrasts": baseline,
               "tradeoff_curve": curve, "chosen": chosen, "report": report,
               "report_card_tool_sha256": tool_sha, "packages": lib.package_versions(),
               "implementation_sha256": {"form_ab_lib": lib.P.sha256(Path(lib.__file__)),
                                         "form_b9": lib.P.sha256(Path(__file__))},
               "weights_published": False, "clinical_labels_used": False}
    (args.out / "b9_results.json").write_text(json.dumps(summary, indent=2, default=float) + "\n", encoding="utf-8")
    print(json.dumps({"chosen": chosen, "criterion_b10_5": report.get("criterion_b10_5")}, indent=2, default=float))


def gate_b10_5(target: dict, sd: float, spearman_loss: float, verdict_changes: dict,
               functional_error: float) -> dict:
    """Gate G11 for the protected student. The target rule's own verdict is meant to change; no other rule's may."""
    others = [rule for rule in verdict_changes if rule != TARGET_RULE]
    gate = {
        "target_met": bool(target["estimate"] <= -0.03 * sd and target["ci95"][1] < 0),
        "fidelity_loss_at_most_0.01": bool(spearman_loss <= 0.01),
        "no_other_rule_changes_verdict": not others,
        "functional_matches_tool": bool(functional_error <= 1e-9),
    }
    gate["passed"] = all(gate.values())
    return gate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("pilot", "full"), required=True)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    args.eval = STAGES[args.stage]["eval"]
    args.out = args.out or STAGES[args.stage]["out"]
    run(args)


if __name__ == "__main__":
    main()

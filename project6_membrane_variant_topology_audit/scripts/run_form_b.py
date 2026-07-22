"""Run the Form B correction and evaluate it on held-out clusters.

The pilot stage fits on training clusters and evaluates on validation clusters.
The full stage evaluates on test clusters. Corrected scores are temporary audit
inputs and are deleted after use. No clinical labels are read.
"""
from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import form_ab_lib as lib
from run_gate1 import read_topology

STAGES = {
    "pilot": dict(eval="validation", out=lib.ROOT / "results/form_ab_pilot_20261008/form_b"),
    "full": dict(eval="test", out=lib.ROOT / "results/form_b_20261008"),
}
SCORERS = {"alphamissense": lib.AM_ARCHIVE, "esm1b": lib.ESM_ARCHIVE}
TARGET_RULE = "positive_charge_gain"


def contrast(window: dict, scores: np.ndarray, draws: int = lib.PANEL_DRAWS) -> dict:
    result = lib.audit_matrix(window, scores, draws=draws)[0]
    return {"estimate": result["estimate_raw"], "ci95": result["ci95_raw"], "penalty_sd": result["penalty_sd"],
            "n_sites": result["n_sites"], "n_strata": result["n_strata"], "n_proteins": result["n_proteins"]}


def correct_and_audit(train_w, held, held_w, held_in_window, sd, margin, length, declared, draws):
    """Fit beta on the training window, apply it to the held-out clusters, audit the held-out window."""
    base_train = lib.audit_matrix(train_w, train_w["scores"], draws=2)[0]["estimate_raw"]
    k = lib.correction_weight(train_w, length)
    beta, target = lib.fit_beta(base_train, k, sd, margin, declared)
    corrected = lib.apply_correction(held, beta, length, declared)
    after = contrast(held_w, corrected[held_in_window], draws)
    return {"delta0_train": base_train, "K_train": k, "beta": beta, "target": target, "after": after}, corrected


def evaluate(scorer, panel, label, key_of, kind, args, tool, primary: bool) -> dict:
    sd = lib.SD_FIXED[scorer]
    train = lib.subset_panel(panel, lib.split_mask(panel, label, "train"))
    held = lib.subset_panel(panel, lib.split_mask(panel, label, args.eval))
    train_w, held_w = lib.window_panel(train), lib.window_panel(held)
    held_in_window = held["distance"] <= 10
    margin, length = lib.CORRECTION["margin"], lib.CORRECTION["decay_length"]
    out = {"scorer": scorer, "split": kind, "evaluated_on": args.eval, "sd_fixed": sd,
           "n_sites": {"train": len(train["protein"]), "held": len(held["protein"]),
                       "train_window": len(train_w["protein"]), "held_window": len(held_w["protein"])}}
    before = contrast(held_w, held_w["scores"])
    out["baseline"] = {"delta0_held": before, "delta0_held_sd_units": before["estimate"] / sd,
                       "penalty_sd_train": contrast(train_w, train_w["scores"], 2)["penalty_sd"],
                       "penalty_sd_held": before["penalty_sd"]}
    fit, corrected = correct_and_audit(train_w, held, held_w, held_in_window, sd, margin, length, -1,
                                       lib.PANEL_DRAWS)
    out["correction"] = {**fit, "margin": margin, "decay_length": length,
                         "after_sd_units": fit["after"]["estimate"] / sd,
                         "target_sd_units": fit["target"] / sd}
    if fit["beta"] > 0:
        train_scores = lib.apply_correction(train_w, fit["beta"], length)
        again = lib.audit_matrix(train_w, train_scores, draws=2)[0]["estimate_raw"]
        out["unit_test"] = {"train_after": again, "error_vs_target": again - fit["target"],
                            "passed": abs(again - fit["target"]) < 1e-9}
    else:
        out["unit_test"] = {"passed": None, "note": "beta = 0, scores unchanged"}
    parts = lib.protein_parts(held_w, length)
    cluster_of = np.array([key_of[str(panel["accessions"][p])] for p in parts["protein"]])
    rebuilt = (parts["D"].sum() - fit["beta"] * parts["K"].sum()) / parts["W"].sum()
    intervals = {"before": lib.cluster_interval(parts["D"], parts["W"], cluster_of),
                 "after": lib.cluster_interval(parts["D"] - fit["beta"] * parts["K"], parts["W"], cluster_of)}
    out["cluster_intervals"] = {**intervals, "unit": "clusters of the split (sequence clusters or symbol families)",
                                "n_clusters": int(len(np.unique(cluster_of))),
                                "parts_vs_tool_error": rebuilt - fit["after"]["estimate"]}
    ok = (fit["after"]["estimate"] <= -0.03 * sd) and intervals["after"][1] < 0
    out["criterion_b10_1"] = {"target_met": bool(ok), "threshold_sd_units": -0.03,
                              "after_sd_units": fit["after"]["estimate"] / sd,
                              "cluster_ci95_upper": intervals["after"][1],
                              "tool_protein_ci95_upper": fit["after"]["ci95"][1]}
    out["fidelity"] = {"all_rows": lib.fidelity(held["scores"], corrected),
                       "window_rows": lib.fidelity(held["scores"][held_in_window], corrected[held_in_window])}
    held_scores_w = held_w["scores"]
    after_scores_w = corrected[held_in_window]
    for tag, scores in (("before", held_scores_w), ("after", after_scores_w)):
        per = lib.protein_contrasts({**held_w, "scores": scores})
        out.setdefault("per_protein", {})[tag] = {
            "n_proteins": int(len(per["protein"])), "fraction_negative": float((per["contrast"] < 0).mean()),
            "median": float(np.median(per["contrast"]))}
        if primary:
            args.out.mkdir(parents=True, exist_ok=True)
            with (args.out / f"per_protein_{scorer}_{tag}.csv").open("w", newline="", encoding="utf-8") as stream:
                writer = csv.writer(stream)
                writer.writerow(["protein_index", "accession", "weight", "n_strata", "contrast"])
                for p, w, n, c in zip(per["protein"], per["weight"], per["n_strata"], per["contrast"]):
                    writer.writerow([int(p), str(panel["accessions"][p]), float(w), int(n), float(c)])
    if not primary:
        return out

    rng = np.random.default_rng(lib.MASTER_SEED if args.stage == "full" else lib.PILOT_SEED)
    controls = {}
    if fit["beta"] > 0:
        changes = []
        for _ in range(args.shuffles):
            side = lib.shuffle_side(held, rng)
            scores = lib.apply_correction(held, fit["beta"], length, side=side)[held_in_window]
            changes.append(lib.audit_matrix(held_w, scores, draws=2)[0]["estimate_raw"] - before["estimate"])
        real = fit["after"]["estimate"] - before["estimate"]
        changes = np.array(changes)
        controls["shuffled_at_application"] = {
            "n_shuffles": int(args.shuffles), "real_change": real, "mean_change": float(changes.mean()),
            "sd_change": float(changes.std(ddof=1)), "max_abs_change": float(np.abs(changes).max()),
            "mean_ratio_to_real": float(abs(changes.mean()) / abs(real)),
            "max_ratio_to_real": float(np.abs(changes).max() / abs(real)),
            "passed": bool(abs(changes.mean()) <= 0.1 * abs(real) and np.abs(changes).max() <= 0.5 * abs(real))}
        shuffled_train = lib.shuffle_side(train_w, rng)
        shuffled_est = lib.audit_matrix({**train_w, "side": shuffled_train}, train_w["scores"], draws=2)[0]
        k_shuffled = lib.correction_weight({**train_w, "side": shuffled_train}, length)
        beta_s, _ = lib.fit_beta(shuffled_est["estimate_raw"], k_shuffled, sd, margin, -1)
        moved = lib.audit_matrix(held_w, lib.apply_correction(held, beta_s, length)[held_in_window], draws=2)[0]
        controls["shuffled_at_fit_descriptive"] = {
            "delta0_shuffled_train": shuffled_est["estimate_raw"], "beta_shuffled": beta_s, "beta_real": fit["beta"],
            "held_contrast_after": moved["estimate_raw"], "held_contrast_after_sd_units": moved["estimate_raw"] / sd,
            "note": "beta is anchored to the configured target, so a shuffled fit is not expected to give zero"}
    wrong, _ = correct_and_audit(train_w, held, held_w, held_in_window, sd, margin, length, +1, lib.PANEL_DRAWS)
    plus_ok = wrong["beta"] > 0 and wrong["after"]["estimate"] >= 0.03 * sd and wrong["after"]["ci95"][0] > 0
    controls["declared_plus_one"] = {**wrong, "after_sd_units": wrong["after"]["estimate"] / sd,
                                     "moved_to_declared_side": bool(plus_ok)}
    out["controls"] = controls

    # Refit the correction across the configured sensitivity grid.
    grid = []
    for sens_length in (4.0, 8.0, 16.0, None):
        for sens_margin in (0.025, 0.05, 0.075):
            row, _ = correct_and_audit(train_w, held, held_w, held_in_window, sd, sens_margin, sens_length, -1, 2)
            grid.append({"decay_length": sens_length, "margin": sens_margin, "beta": row["beta"],
                         "K_train": row["K_train"], "held_contrast_sd_units": row["after"]["estimate"] / sd})
    out["sensitivity"] = grid

    if args.report_card:
        card_dir = args.out / "report_card"
        (card_dir / f"{scorer}_before").mkdir(parents=True, exist_ok=True)
        cards = {"before": lib.run_report_card(tool, held, held["scores"], card_dir / f"{scorer}_before")}
        if fit["beta"] > 0:
            (card_dir / f"{scorer}_after").mkdir(parents=True, exist_ok=True)
            cards["after"] = lib.run_report_card(tool, held, corrected, card_dir / f"{scorer}_after")
            diff = lib.card_diff(tool, cards["before"], cards["after"], TARGET_RULE)
            diff["target_rule_before_after"] = {
                scope: {"before": cards["before"]["rule_audit"][scope]["rules"][TARGET_RULE],
                        "after": cards["after"]["rule_audit"][scope]["rules"][TARGET_RULE]}
                for scope in cards["before"]["rule_audit"]}
            out["report_card_diff"] = diff
        else:
            out["report_card_diff"] = {"note": "beta = 0: corrected scores are identical to the originals",
                                       "n_unexpected": 0}
    return out


def criteria(results: dict, report_card_validated: bool) -> dict:
    am, esm = results["alphamissense"]["sequence"], results["esm1b"]["sequence"]
    unchanged = (am.get("report_card_diff", {}).get("n_unexpected") == 0) if "report_card_diff" in am else None
    out = {
        "report_card_validated_by_a11": report_card_validated,
        "B10_1_alphamissense_sequence_split": am["criterion_b10_1"]["target_met"],
        "B10_1_alphamissense_symbol_split": results["alphamissense"]["symbol"]["criterion_b10_1"]["target_met"],
        "B10_2_other_rules_unchanged": unchanged if report_card_validated else None,
        "B10_3_fidelity_all_rows": am["fidelity"]["all_rows"]["spearman"] >= 0.999,
        "B10_3_fidelity_window_rows": am["fidelity"]["window_rows"]["spearman"] >= 0.999,
        "B10_4a_shuffled_at_application": am["controls"].get("shuffled_at_application", {}).get("passed"),
        "B10_4b_esm1b_beta_zero": esm["correction"]["beta"] == 0.0,
        "B10_4c_esm1b_declared_plus_one_moves": esm["controls"]["declared_plus_one"]["moved_to_declared_side"],
        "unit_test_alphamissense": am["unit_test"]["passed"],
    }
    gated = {k: v for k, v in out.items() if k not in ("B10_1_alphamissense_symbol_split",
                                                       "report_card_validated_by_a11")}
    out["all_gated_criteria_pass"] = bool(report_card_validated) and all(v is True for v in gated.values())
    out["family_wording_allowed"] = bool(out["B10_1_alphamissense_symbol_split"])
    return out


def run(args) -> None:
    if args.stage == "full":
        form_a = lib.ROOT / "results/form_a_20261008/criteria.json"
        if not form_a.exists() or not json.loads(form_a.read_text(encoding="utf-8")).get("form_b_may_start"):
            raise SystemExit("Form A has not passed; Form B does not start")
        a11 = lib.ROOT / "results/form_a_20261008/a11/a11_criteria.json"
        validated = a11.exists() and json.loads(a11.read_text(encoding="utf-8")).get("a11_passed") is True
    else:
        validated = True
    args.out.mkdir(parents=True, exist_ok=True)
    census = read_topology(lib.TOPOLOGY)
    tool, tool_sha = lib.report_card_tool()
    splits = {"sequence": lib.cluster_split(census), "symbol": lib.cluster_split(census, kind="symbol")}
    manifest = {
        "stage": args.stage, "as_of_utc": datetime.now(timezone.utc).isoformat(), "evaluated_on": args.eval,
        "splits": {k: {x: v[x] for x in ("kind", "seed", "counts", "assignment_sha256")}
                   for k, v in splits.items()},
        "correction": lib.CORRECTION, "sd_fixed": lib.SD_FIXED, "panel_draws": lib.PANEL_DRAWS,
        "panel_seed": lib.PANEL_SEED, "report_card_settings": lib.REPORT_CARD, "shuffles": args.shuffles,
        "validated_tool": lib.validated_tool(), "report_card_tool_sha256": tool_sha,
        "packages": lib.package_versions(),
        "archive_sha256": {"alphamissense": lib.P.sha256(lib.AM_ARCHIVE), "esm1b": lib.P.sha256(lib.ESM_ARCHIVE)},
        "topology_sha256": lib.P.sha256(lib.TOPOLOGY),
        "implementation_sha256": {"form_ab_lib": lib.P.sha256(Path(lib.__file__)),
                                  "run_form_b": lib.P.sha256(Path(__file__))},
        "clinical_labels_used": False, "corrected_scores_published": False}
    (args.out / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    results = {}
    for scorer, path in SCORERS.items():
        print(f"loading {scorer}", flush=True)
        panel = lib.P.load_panel(path, 60)
        results[scorer] = {}
        for kind in ("sequence", "symbol"):
            print(f"  {scorer} / {kind} split", flush=True)
            results[scorer][kind] = evaluate(scorer, panel, splits[kind]["label"], splits[kind]["key_of"], kind,
                                             args, tool, primary=(kind == "sequence"))
            (args.out / f"{scorer}_{kind}.json").write_text(
                json.dumps(results[scorer][kind], indent=2) + "\n", encoding="utf-8")
    verdict = criteria(results, validated)
    (args.out / "criteria.json").write_text(json.dumps(verdict, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(verdict, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("pilot", "full"), required=True)
    parser.add_argument("--shuffles", type=int, default=100)
    parser.add_argument("--no-report-card", dest="report_card", action="store_false")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    defaults = STAGES[args.stage]
    args.eval = defaults["eval"]
    args.out = args.out or defaults["out"]
    run(args)


if __name__ == "__main__":
    main()

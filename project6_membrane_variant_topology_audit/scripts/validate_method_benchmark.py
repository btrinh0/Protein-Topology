"""Benchmark the retained K/R-gain estimator against known-answer generators and explicit
baselines.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from membrane_rule_audit.panel import (  # noqa: E402
    ALPHABET, NEUTRAL, gain_contrast, sha256, stratum_parts, summarize_parts,
)

PLAN = ROOT / "docs" / "METHOD_COMPLETION_PLAN_20261004.json"
OUT = ROOT / "results" / "method_completion_20261004"
METHODS = (
    "naive_charge_arm_only", "naive_unstratified_paired",
    "packaged_within_site_stratified", "standard_fixed_effects",
)


def wilson(k: int, n: int) -> list[float]:
    z = 1.959963984540054
    p = k / n
    den = 1 + z * z / n
    center = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [center - half, center + half]


def binomial_two_sided(k: int, n: int, p: float) -> float:
    log_mass = np.array([
        math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
        + i * math.log(p) + (n - i) * math.log1p(-p)
        for i in range(n + 1)
    ])
    return float(min(1.0, np.exp(log_mass[log_mass <= log_mass[k] + 1e-12]).sum()))


def holm(values: list[float]) -> list[float]:
    order = np.argsort(values)
    adjusted = np.empty(len(values))
    previous = 0.0
    for rank, index in enumerate(order):
        previous = max(previous, min(1.0, (len(values) - rank) * values[index]))
        adjusted[index] = previous
    return adjusted.tolist()


def layout(spec: dict, cfg: dict) -> dict:
    rows = []
    refs = cfg["reference_residues"]
    for protein in range(cfg["independent_proteins"]):
        position = 0
        for j, ref in enumerate(refs):
            for side in (0, 1):
                n = cfg["balanced_sites_per_protein_reference_side"]
                if spec["layout"] == "reference_imbalance":
                    n = 8 if ((ref in "NQ") == (side == 0)) else 2
                elif spec["layout"] == "uneven_proteins":
                    n *= 1 + protein % 4
                for k in range(n):
                    position += 1
                    distance = 1 + (3 * k + protein + 2 * j + 5 * side) % 10
                    rows.append((protein, position, side, distance, ALPHABET.index(ref)))
    array = np.asarray(rows, dtype=np.int64)
    return {
        "protein": array[:, 0], "position": array[:, 1], "side": array[:, 2],
        "distance": array[:, 3], "ref_aa": array[:, 4],
        "tmd_index": np.zeros(len(array), dtype=np.int64),
    }


def generate(panel: dict, spec: dict, cfg: dict, rng: np.random.Generator,
             deterministic: bool = False, uniform_charge_side: float = 0.0) -> dict:
    protein = panel["protein"]
    n_proteins = cfg["independent_proteins"]
    cyto = (panel["side"] == 0).astype(float)
    noise = 0.0 if deterministic else 1.0
    constraint = noise * rng.normal(0, cfg["protein_constraint_sd"], n_proteins)[protein]
    constraint += spec["shared_cytosolic_constraint"] * cyto
    constraint += cfg["shared_distance_constraint_per_residue"] * panel["distance"]
    constraint += noise * rng.normal(0, cfg["shared_site_noise_sd"], len(protein))
    ref_offset = np.array([
        cfg["reference_charge_offsets"][ALPHABET[int(x)]] for x in panel["ref_aa"]
    ])
    charge_intercept = noise * rng.normal(
        0, cfg["protein_charge_intercept_sd"], n_proteins)[protein]
    charge_side_sd = spec.get("protein_by_side_charge_sd", cfg["default_protein_by_side_charge_sd"])
    charge_side = noise * rng.normal(0, charge_side_sd, n_proteins)[protein]
    charge_side += uniform_charge_side
    damping = (np.exp(-panel["distance"] / spec["decay_lambda"])
               if "decay_lambda" in spec else np.ones(len(protein)))
    charge = ref_offset + charge_intercept + cyto * (spec["theta"] * damping + charge_side)
    kr_mask = np.array([a in "KR" for a in ALPHABET], dtype=float)
    chemistry = np.linspace(-0.15, 0.15, 20)
    scores = constraint[:, None] + chemistry[None, :] + charge[:, None] * kr_mask
    if not deterministic:
        scores += rng.normal(0, cfg["variant_noise_sd"], scores.shape)
    scores[np.arange(len(protein)), panel["ref_aa"]] = np.nan
    return dict(panel, scores=scores)


def independent_oracle(panel: dict, spec: dict) -> float:
    """Expected target from declared cell weights, independent of score helpers."""
    numerator = denominator = 0.0
    for protein in np.unique(panel["protein"]):
        for ref in np.unique(panel["ref_aa"]):
            cell = (panel["protein"] == protein) & (panel["ref_aa"] == ref)
            cytosolic = cell & (panel["side"] == 0)
            outer = cell & (panel["side"] == 1)
            nc, no = int(cytosolic.sum()), int(outer.sum())
            if not nc or not no:
                continue
            damping = (np.exp(-panel["distance"][cytosolic] / spec["decay_lambda"])
                       if "decay_lambda" in spec else np.ones(nc))
            weight = nc * no / (nc + no)
            numerator += weight * spec["theta"] * float(damping.mean())
            denominator += weight
    return numerator / denominator


def independent_fixed_effects(panel: dict, penalty: np.ndarray,
                              protein_weights: np.ndarray | None = None) -> float:
    """Demean site outcomes and the side regressor, then ordinary least squares."""
    groups = np.stack([panel["protein"], panel["ref_aa"]], axis=1)
    _, inverse = np.unique(groups, axis=0, return_inverse=True)
    n = np.bincount(inverse)
    x = (panel["side"] == 0).astype(float)
    x_center = x - (np.bincount(inverse, weights=x) / n)[inverse]
    y_center = penalty - (np.bincount(inverse, weights=penalty) / n)[inverse]
    weight = (np.ones(len(penalty)) if protein_weights is None
              else protein_weights[panel["protein"]])
    return float(np.dot(weight * x_center, y_center) / np.dot(weight * x_center, x_center))


def ingredients(panel: dict) -> dict:
    paired = gain_contrast(panel)
    primary = stratum_parts(panel, paired["penalty"], paired["eligible"])
    kr = np.mean(panel["scores"][:, [ALPHABET.index("K"), ALPHABET.index("R")]], axis=1)
    arm = stratum_parts(panel, kr, paired["eligible"])
    n = len(primary["denominator"])
    cyto = panel["side"] == 0
    p = panel["protein"]
    cyto_sum = np.bincount(p[cyto], weights=paired["penalty"][cyto], minlength=n)
    cyto_n = np.bincount(p[cyto], minlength=n)
    outer_sum = np.bincount(p[~cyto], weights=paired["penalty"][~cyto], minlength=n)
    outer_n = np.bincount(p[~cyto], minlength=n)
    primary_point = float(primary["numerator"].sum() / primary["denominator"].sum())
    fe = independent_fixed_effects(panel, paired["penalty"])
    points = np.array([
        arm["numerator"].sum() / arm["denominator"].sum(),
        cyto_sum.sum() / cyto_n.sum() - outer_sum.sum() / outer_n.sum(),
        primary_point, fe,
    ])
    columns = np.column_stack([
        primary["numerator"], primary["denominator"], arm["numerator"], arm["denominator"],
        cyto_sum, cyto_n, outer_sum, outer_n,
    ])
    return {"points": points, "columns": columns,
            "fixed_effects_error": abs(primary_point - fe), "parts": primary}


def interval_samples(columns: np.ndarray, counts: np.ndarray) -> np.ndarray:
    summed = counts @ columns
    primary = summed[:, 0] / summed[:, 1]
    return np.column_stack([
        summed[:, 2] / summed[:, 3],
        summed[:, 4] / summed[:, 5] - summed[:, 6] / summed[:, 7],
        primary, primary,
    ])


def literal_arm_fixture() -> dict:
    """Independent numeric answer, including coverage and exclusion edge cases."""
    scores = np.full((7, 20), 2.0)
    scores[:, ALPHABET.index("K")] = 3.0
    scores[:, ALPHABET.index("R")] = 7.0
    for aa in "DEH":
        scores[:, ALPHABET.index(aa)] = 10000.0
    refs = np.array([ALPHABET.index(a) for a in "AAAAKHD"])
    scores[np.arange(len(refs)), refs] = np.nan
    scores[1, ALPHABET.index("K")] = np.nan
    scores[2] = np.nan
    scores[2, ALPHABET.index("K")] = 3.0
    scores[2, ALPHABET.index("R")] = 7.0
    scores[2, ALPHABET.index("Q")] = 4.0
    scores[3] = np.nan
    scores[3, ALPHABET.index("K")] = 3.0
    scores[3, ALPHABET.index("R")] = 7.0
    answer = gain_contrast({"scores": scores, "ref_aa": refs})
    expected = np.array([3.0, 5.0, 1.0, np.nan, 5.0, 3.0, 3.0])
    expected_eligibility = np.array([True, True, True, False, False, False, False])
    error = float(np.nanmax(np.abs(answer["penalty"] - expected)))
    passed = (error < 1e-10 and np.isnan(answer["penalty"][3])
              and np.array_equal(answer["eligible"], expected_eligibility))
    return {"name": "literal_charge_and_comparator_arms", "expected_penalties": [3, 5, 1, None, 5, 3, 3],
            "expected_eligibility": expected_eligibility.tolist(), "max_error": error,
            "native_reference": "A excluded from its own comparator; K/H/D references ineligible",
            "excluded_alternatives": "D/E/H are scored 10000 and must not enter the neutral arm",
            "passed": bool(passed)}


def run(plan_path: Path, out: Path) -> dict:
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    cfg = plan["benchmark"]
    reps, draws, n_proteins = (cfg["replicates_per_scenario"], cfg["bootstrap_draws"],
                               cfg["independent_proteins"])
    cells, fixtures, replicate_rows = [], [], []
    arm_fixture = literal_arm_fixture()
    if not arm_fixture["passed"]:
        raise ValueError("Literal arm fixture failed before benchmark execution")
    max_fe_error = 0.0
    for scenario_index, spec in enumerate(cfg["scenarios"]):
        panel = layout(spec, cfg)
        truth = independent_oracle(panel, spec)
        rng = np.random.default_rng(np.random.SeedSequence([plan["seed"], scenario_index, 0]))
        deterministic = generate(panel, spec, cfg, rng, deterministic=True)
        check = ingredients(deterministic)
        shifted = generate(panel, spec, cfg, rng, deterministic=True, uniform_charge_side=0.25)
        shared = dict(deterministic, scores=deterministic["scores"] +
                      7.5 * (panel["side"] == 0)[:, None])
        base_penalty = gain_contrast(deterministic)["penalty"]
        shared_error = float(np.nanmax(np.abs(gain_contrast(shared)["penalty"] - base_penalty)))
        injection_error = abs((ingredients(shifted)["points"][2] - check["points"][2]) - 0.25)
        oracle_error = abs(check["points"][2] - truth)
        bridge_panel = generate(panel, spec, cfg, np.random.default_rng(
            np.random.SeedSequence([plan["seed"], scenario_index, 10000])))
        bridge = ingredients(bridge_panel)
        bridge_seed = plan["seed"] + scenario_index + 70000
        production_summary = summarize_parts(bridge["parts"], draws, bridge_seed)
        bridge_counts = np.random.default_rng(bridge_seed).multinomial(
            n_proteins, np.full(n_proteins, 1 / n_proteins), draws)
        bridge_samples = interval_samples(bridge["columns"], bridge_counts)
        bridge_point_error = abs(production_summary["estimate_raw"] - bridge["points"][2])
        bridge_interval_error = float(np.max(np.abs(np.percentile(bridge_samples[:, 2],
            [2.5, 97.5]) - np.asarray(production_summary["ci95_raw"]))))
        bridge_penalty = gain_contrast(bridge_panel)["penalty"]
        bootstrap_fe_error = max(abs(independent_fixed_effects(bridge_panel, bridge_penalty,
            bridge_counts[i]) - bridge_samples[i, 2]) for i in range(25))
        fixtures.append({
            "scenario": spec["name"], "oracle": truth, "n_sites": len(panel["protein"]),
            "oracle_error": oracle_error, "fixed_effects_error": check["fixed_effects_error"],
            "shared_constraint_cancellation_error": shared_error,
            "uniform_protein_by_side_injection_error": injection_error,
            "production_point_bridge_error": bridge_point_error,
            "production_interval_bridge_error": bridge_interval_error,
            "weighted_bootstrap_fixed_effects_error": bootstrap_fe_error,
            "weighted_bootstrap_fixed_effects_draws_checked": 25,
            "passed": bool(max(oracle_error, shared_error, injection_error,
                               check["fixed_effects_error"], bridge_point_error,
                               bridge_interval_error, bootstrap_fe_error) < 1e-10),
        })
        print(f"Benchmark {scenario_index + 1}/{len(cfg['scenarios'])}: {spec['name']}; "
              f"truth={truth:+.6f}; {reps} replicates", flush=True)
        points = np.empty((reps, len(METHODS)))
        intervals = np.empty((reps, 2, len(METHODS)))
        for replicate in range(reps):
            simulation_rng = np.random.default_rng(np.random.SeedSequence(
                [plan["seed"], scenario_index, replicate, 1]))
            bootstrap_rng = np.random.default_rng(np.random.SeedSequence(
                [plan["seed"], scenario_index, replicate, 2]))
            sample = generate(panel, spec, cfg, simulation_rng)
            part = ingredients(sample)
            points[replicate] = part["points"]
            max_fe_error = max(max_fe_error, part["fixed_effects_error"])
            counts = bootstrap_rng.multinomial(n_proteins, np.full(n_proteins, 1 / n_proteins), draws)
            intervals[replicate] = np.percentile(interval_samples(part["columns"], counts),
                                                [2.5, 97.5], axis=0)
        for m, method in enumerate(METHODS):
            estimates = points[:, m]
            lo, hi = intervals[:, 0, m], intervals[:, 1, m]
            covered = int(np.count_nonzero((lo <= truth) & (hi >= truth)))
            rejects = int(np.count_nonzero((lo > 0) | (hi < 0)))
            correct = int(np.count_nonzero(lo > 0 if truth > 0 else hi < 0)) if truth else None
            bias = float(estimates.mean() - truth)
            mcse = float(estimates.std(ddof=1) / math.sqrt(reps))
            squared_error = (estimates - truth) ** 2
            rmse = float(np.sqrt(squared_error.mean()))
            rmse_mcse = float(squared_error.std(ddof=1) / math.sqrt(reps) / (2 * rmse)) if rmse else 0.0
            direction_n = int(np.count_nonzero(estimates * truth > 0)) if truth else None
            cell = {
                "scenario": spec["name"], "method": method, "target": truth,
                "replicates": reps, "n_proteins": n_proteins, "n_sites": len(panel["protein"]),
                "mean_estimate": float(estimates.mean()), "bias": bias, "bias_mcse": mcse,
                "bias_mc95": [bias - 1.959963984540054 * mcse, bias + 1.959963984540054 * mcse],
                "rmse": rmse, "rmse_mcse_delta_method": rmse_mcse,
                "rmse_mc95_delta_method": [max(0.0, rmse - 1.959963984540054 * rmse_mcse),
                                          rmse + 1.959963984540054 * rmse_mcse],
                "coverage": covered / reps, "coverage_mc95": wilson(covered, reps),
                "coverage_exact_binomial_p": binomial_two_sided(reps - covered, reps, 0.05),
                "zero_rejection": rejects / reps, "zero_rejection_mc95": wilson(rejects, reps),
                "zero_rejection_kind": "false_claim_rate_for_charge_specific_target" if truth == 0 else "power",
                "point_correct_direction_fraction": direction_n / reps if direction_n is not None else None,
                "point_correct_direction_mc95": wilson(direction_n, reps) if direction_n is not None else None,
                "correct_direction_rejection": correct / reps if correct is not None else None,
                "correct_direction_rejection_mc95": wilson(correct, reps) if correct is not None else None,
                "interval": "Percentile whole-protein resampling; same draws across methods within each replicate",
                "fixed_effects_interval_note": "Fixed-effects point coefficients and sampled coefficients agree algebraically; reuse primary bootstrap draws rather than claim independent uncertainty" if method == "standard_fixed_effects" else None,
            }
            cells.append(cell)
            for replicate in range(reps):
                replicate_rows.append({"scenario": spec["name"], "method": method,
                    "replicate": replicate, "target": truth, "estimate": float(estimates[replicate]),
                    "ci_low": float(lo[replicate]), "ci_high": float(hi[replicate])})
        print(f"Completed {spec['name']}", flush=True)
    primary = [c for c in cells if c["method"] == "packaged_within_site_stratified"]
    adjusted = holm([c["coverage_exact_binomial_p"] for c in primary])
    for cell, p in zip(primary, adjusted):
        cell["coverage_holm_adjusted_p"] = p
        cell["bias_check_passed"] = abs(cell["bias"]) < max(0.005, 4 * cell["bias_mcse"])
    calibration = all(p >= 0.05 for p in adjusted)
    fixtures_pass = arm_fixture["passed"] and all(f["passed"] for f in fixtures) and max_fe_error < 1e-10
    bias_pass = all(c["bias_check_passed"] for c in primary)
    summary = {
        "date": plan["date"], "scope": "Known-answer validation for the predeclared additive-score generators only",
        "plan_sha256": sha256(plan_path), "script_sha256": sha256(Path(__file__)),
        "production_sha256": sha256(ROOT / "tools/membrane_rule_audit/panel.py"),
        "python": platform.python_version(), "numpy": np.__version__, "seed": plan["seed"],
        "replicates_per_scenario": reps, "bootstrap_draws": draws,
        "scenario_count": len(cfg["scenarios"]), "method_count": len(METHODS),
        "literal_arm_fixture": arm_fixture, "fixtures": fixtures,
        "max_simulated_fixed_effects_error": max_fe_error,
        "cells": cells, "oracle_checks_passed": fixtures_pass,
        "primary_bias_checks_passed": bias_pass, "limited_calibration_check_passed": calibration,
        "limited_benchmark_passed": bool(fixtures_pass and bias_pass and calibration),
        "limitations": [
            "Validates explicit additive generators, not every nonlinear predictor or biological mechanism",
            "Independent simulated proteins do not establish homology-disjoint generalization",
            "Naive baselines intentionally omit one control and target different observed quantities; their charge-specific errors are not numerical bugs",
            "Standard fixed effects is established prior art and agrees with the packaged point estimate",
            "Interval calibration cannot be transferred to small-gene Phase C or to clinical outcomes",
        ],
    }
    out.mkdir(parents=True, exist_ok=True)
    with (out / "benchmark_replicates.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(replicate_rows[0]))
        writer.writeheader()
        writer.writerows(replicate_rows)
    (out / "benchmark_results.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: summary[k] for k in ("limited_benchmark_passed", "oracle_checks_passed",
          "primary_bias_checks_passed", "limited_calibration_check_passed",
          "max_simulated_fixed_effects_error")}, indent=2), flush=True)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, default=PLAN)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()
    result = run(args.plan, args.out)
    if not result["limited_benchmark_passed"]:
        raise SystemExit("Benchmark did not meet its predeclared limited-validation checks; retain all results.")

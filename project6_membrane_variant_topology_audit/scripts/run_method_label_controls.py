"""Run side-label destruction controls and count declared relative-score candidates."""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from membrane_rule_audit.panel import (  # noqa: E402
    ALPHABET, NEUTRAL, arm_mean, gain_contrast, load_panel, sha256, stratum_parts,
)

PLAN = ROOT / "config" / "METHOD_COMPLETION_PLAN_20261004.json"
OUT = ROOT / "results" / "method_completion_20261004"
SOURCES = {
    "alphamissense": ROOT / "data/processed/gate0/am_full_flank_scores.npz",
    "esm1b": ROOT / "data/processed/gate0/esm1b_flank_scores.npz",
}


def candidate_counts(panel: dict, penalty: np.ndarray, eligible: np.ndarray) -> dict:
    outer = eligible & (panel["side"] == 1)
    neutral, n_neutral = arm_mean(panel["scores"], NEUTRAL)
    variants = {}
    any_candidate = np.zeros(len(penalty), dtype=bool)
    for aa in "KR":
        value = panel["scores"][:, ALPHABET.index(aa)]
        covered = outer & np.isfinite(value)
        candidate = covered & (value < neutral)
        any_candidate |= candidate
        variants[aa] = {"eligible_variants": int(covered.sum()), "candidate_variants": int(candidate.sum()),
            "unique_proteins_with_candidates": int(len(np.unique(panel["protein"][candidate])))}
    pooled = outer & (penalty < 0)
    full = outer & (n_neutral == len(NEUTRAL) - 1)
    full &= np.isfinite(panel["scores"][:, ALPHABET.index("K")])
    full &= np.isfinite(panel["scores"][:, ALPHABET.index("R")])
    return {
        "definition": "Outer K/R gain whose score is below the measured same-site uncharged-alternative mean",
        "threshold": 0.0, "window_residues": 10,
        "eligible_outer_sites": int(outer.sum()),
        "eligible_outer_proteins": int(len(np.unique(panel["protein"][outer]))),
        "fully_covered_outer_sites": int(full.sum()), "by_introduced_residue": variants,
        "count_population": "All eligible outer sites, not only sites in two-sided matched strata",
        "introduced_arm_coverage": {
            str(n): int(np.count_nonzero(outer & (np.isfinite(panel["scores"][:,
                [ALPHABET.index("K"), ALPHABET.index("R")]]).sum(axis=1) == n))) for n in (1, 2)
        },
        "neutral_arm_coverage": {str(int(n)): int(np.count_nonzero(outer & (n_neutral == n)))
                                 for n in np.unique(n_neutral[outer])},
        "candidate_variants_K_plus_R": sum(v["candidate_variants"] for v in variants.values()),
        "sites_with_at_least_one_candidate_variant": int(any_candidate.sum()),
        "proteins_with_at_least_one_candidate_variant": int(len(np.unique(panel["protein"][any_candidate]))),
        "pooled_KR_negative_sites": int(pooled.sum()),
        "pooled_KR_negative_proteins": int(len(np.unique(panel["protein"][pooled]))),
        "interpretation": "Counts of relative score-pattern candidates, not confirmed insertion defects, reversed internal rules, clinical errors or patients",
    }


def run_model(model: str, path: Path, plan: dict, model_index: int) -> dict:
    cfg = plan["shuffled_label_control"]
    full_panel = load_panel(path, cfg["window"])
    score_direction = cfg["stored_score_directions"][model]
    if score_direction == "higher-is-tolerated":
        full_panel["scores"] = -full_panel["scores"]
    elif score_direction != "higher-is-damaging":
        raise ValueError("Declare the stored score orientation")
    contrast = gain_contrast(full_panel)
    impact = candidate_counts(full_panel, contrast["penalty"], contrast["eligible"])
    group = full_panel["protein"].astype(int) * 20 + full_panel["ref_aa"].astype(int)
    eligible = contrast["eligible"] & np.isfinite(contrast["penalty"])
    cyto = full_panel["side"] == 0
    nc = np.bincount(group[eligible & cyto], minlength=int(group.max()) + 1)
    no = np.bincount(group[eligible & ~cyto], minlength=len(nc))
    shared = (nc > 0) & (no > 0)
    keep = eligible & shared[group]
    panel = {k: full_panel[k][keep].copy() for k in ("protein", "position", "side", "ref_aa", "distance", "tmd_index")}
    penalty = contrast["penalty"][keep]
    _, inverse = np.unique(group[keep], return_inverse=True)
    order = np.argsort(inverse, kind="stable")
    panel = {k: v[order] for k, v in panel.items()}
    penalty, inverse = penalty[order], inverse[order]
    n_group = np.bincount(inverse)
    n_cyto = np.bincount(inverse, weights=panel["side"] == 0).astype(int)
    starts = np.r_[0, np.cumsum(n_group)[:-1]]
    ranks = np.arange(len(penalty)) - np.repeat(starts, n_group)
    cyto_template = ranks < np.repeat(n_cyto, n_group)
    weights = n_cyto * (n_group - n_cyto) / n_group
    denominator = float(weights.sum())
    group_sum = np.bincount(inverse, weights=penalty)
    residual = penalty - (group_sum / n_group)[inverse]
    group_sample_variance = np.bincount(inverse, weights=residual ** 2) / (n_group - 1)
    exact_variance = float(np.sum(weights * group_sample_variance) / denominator ** 2)
    closed_offset = float(np.sum(n_cyto / n_group * group_sum))
    original = stratum_parts(panel, penalty, np.ones(len(penalty), dtype=bool))
    observed = float(original["numerator"].sum() / original["denominator"].sum())
    rng = np.random.default_rng(np.random.SeedSequence([plan["seed"], model_index, 3]))
    points = []
    max_agreement_error = 0.0
    count_failures = 0
    original_side = panel["side"].copy()
    print(f"{model}: {len(penalty):,} contributing sites; {cfg['shuffles_per_model']} label shuffles", flush=True)
    for replicate in range(cfg["shuffles_per_model"]):
        shuffled_order = np.argsort(inverse + rng.random(len(inverse)), kind="stable")
        shuffled_cyto = np.empty(len(penalty), dtype=bool)
        shuffled_cyto[shuffled_order] = cyto_template
        check_counts = np.bincount(inverse, weights=shuffled_cyto).astype(int)
        count_failures += int(not np.array_equal(check_counts, n_cyto))
        panel["side"] = (~shuffled_cyto).astype(np.int8)
        part = stratum_parts(panel, penalty, np.ones(len(penalty), dtype=bool))
        point = float(part["numerator"].sum() / part["denominator"].sum())
        independent = float((penalty[shuffled_cyto].sum() - closed_offset) / denominator)
        max_agreement_error = max(max_agreement_error, abs(point - independent))
        points.append(point)
        if (replicate + 1) % 250 == 0:
            print(f"{model}: {replicate + 1} shuffles completed", flush=True)
    panel["side"] = original_side
    values = np.array(points)
    mcse = float(values.std(ddof=1) / math.sqrt(len(values)))
    exact_mean_mcse = math.sqrt(exact_variance / len(values))
    empirical_variance = float(values.var(ddof=1))
    approximate_variance_mcse = exact_variance * math.sqrt(2 / (len(values) - 1))
    quantiles = np.percentile(values, [2.5, 50, 97.5]).tolist()
    checks = {
        "central_95_contains_zero": quantiles[0] <= 0 <= quantiles[2],
        "mean_centered_within_four_mcse": abs(float(values.mean())) < max(1e-10, 4 * exact_mean_mcse),
        "all_stratum_counts_preserved": count_failures == 0,
        "production_and_independent_statistic_agree": max_agreement_error < 1e-10,
        "variance_agrees_with_exact_oracle_within_approx_four_mcse":
            abs(empirical_variance - exact_variance) <= 4 * approximate_variance_mcse,
    }
    return {
        "model": model, "input_sha256": sha256(path), "stored_score_direction": score_direction,
        "output_score_direction": "higher-is-damaging", "observed_fixed_comparator_raw": observed,
        "contributing_sites": len(penalty), "contributing_strata": len(n_group),
        "contributing_proteins": original["numerator"].size,
        "shuffles": len(values), "shuffled_mean_raw": float(values.mean()),
        "shuffled_sd_raw": float(values.std(ddof=1)), "mean_mcse_raw": mcse,
        "exact_shuffled_mean": 0.0, "exact_shuffled_variance_raw": exact_variance,
        "exact_mean_mcse_raw": exact_mean_mcse, "empirical_shuffled_variance_raw": empirical_variance,
        "empirical_to_exact_variance_ratio": empirical_variance / exact_variance,
        "variance_mcse_raw_approximate": approximate_variance_mcse,
        "variance_check_note": "Gaussian large-sample Monte Carlo approximation; the finite-population variance oracle itself is exact",
        "shuffled_central_95_and_median_raw": quantiles,
        "mean_mc95_raw": [float(values.mean()) - 1.959963984540054 * mcse,
                          float(values.mean()) + 1.959963984540054 * mcse],
        "original_penalty_sd": original["penalty_sd"],
        "max_production_statistic_agreement_error": max_agreement_error,
        "count_preservation_failures": count_failures, "checks": checks,
        "pipeline_negative_control_passed": all(checks.values()),
        "original_data_randomization_p_value": None,
        "interpretation": "Conditional label-destruction pipeline control; no exchangeability-based claim about original observational scores",
        "shuffled_estimates_raw": points, "candidate_impact_counts": impact,
    }


def run(plan_path: Path, out: Path) -> dict:
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    models = {name: run_model(name, SOURCES[name], plan, i)
              for i, name in enumerate(plan["shuffled_label_control"]["models"])}
    summary = {
        "date": plan["date"], "plan_sha256": sha256(plan_path),
        "script_sha256": sha256(Path(__file__)),
        "production_sha256": sha256(ROOT / "tools/membrane_rule_audit/panel.py"),
        "seed": plan["seed"], "models": models,
        "all_pipeline_controls_passed": all(m["pipeline_negative_control_passed"] for m in models.values()),
        "clinical_labels_used": False, "comparator_or_scale_search_performed": False,
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "label_controls_and_counts.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"all_pipeline_controls_passed": summary["all_pipeline_controls_passed"],
                      "models": {k: {key: v[key] for key in ("shuffled_mean_raw", "mean_mcse_raw", "pipeline_negative_control_passed")}
                                 for k, v in models.items()}}, indent=2), flush=True)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, default=PLAN)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()
    result = run(args.plan, args.out)
    if not result["all_pipeline_controls_passed"]:
        raise SystemExit("A pipeline control did not meet its declared checks; keep all results.")

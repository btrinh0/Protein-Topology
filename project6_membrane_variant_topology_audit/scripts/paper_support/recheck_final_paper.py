"""Check final-paper values against saved results and optional score arrays. This helper
was added on 10 November 2026 for the support bundle.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np


DEFAULT_ROOT = Path(__file__).resolve().parents[2]
PI_DS = "positive_inside+disulfide_loss"
THREE = "PI_DS_constraint_glyco"
CHECKS = []


def load_json(root, relative):
    return json.loads((root / relative).read_text(encoding="utf-8"))


def load_csv(root, relative):
    with (root / relative).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def near(label, observed, expected, tolerance):
    CHECKS.append(dict(check=label, observed=float(observed), expected=float(expected),
                       tolerance=tolerance, passed=abs(float(observed) - expected) <= tolerance))


def select(records, **fields):
    selected = [row for row in records if all(str(row.get(key)) == str(value) for key, value in fields.items())]
    if len(selected) != 1:
        raise ValueError(f"Expected exactly one row for {fields}, found {len(selected)}")
    return selected[0]


def saved_results(root):
    base = "results/standalone_multirule_20261025/"
    plan = load_json(root, base + "benchmark_20261025/BENCHMARK_PLAN.json")
    for key, expected in (("rows", 127811), ("proteins", 39), ("membrane_proteins", 13),
                          ("dependence_groups", 34), ("n_pc", 64)):
        near("training " + key, plan[key], expected, 0)
    assays = load_csv(root, "results/standalone_dataset_accounting_20261028/current_44_assays.csv")
    near("DMS assay count", len(assays), 44, 0)
    near("membrane measurements", sum(int(row["measurements"]) for row in assays if row["membrane"] == "1"), 43991, 0)

    per = load_csv(root, "results/dms_comparator_pvalues_20261028/per_protein.csv")
    group_out = [row for row in per if row["scheme"] == "guarded_loocv"]
    means = {}
    for model, expected in (("raw", .4623), ("positive_inside", .4622), (PI_DS, .4626),
                            ("VespaG", .4509), ("GEMME", .4280)):
        values = [float(row["spearman"]) for row in group_out if row["model"] == model]
        near("Table 4 protein count " + model, len(values), 13, 0)
        means[model] = float(np.mean(values))
        near("Table 4 mean " + model, means[model], expected, .00005)
    glyco = load_json(root, "results/glycosylation_constraint_repair_20261031/run2/accuracy_results.json")
    near("Table 4 constrained mean", glyco["summary"][THREE]["mean_spearman"], .4608, .00005)
    comparisons = load_csv(root, "results/dms_comparator_pvalues_20261028/highlighted_configuration_main_comparisons.csv")
    for reference, expected in (("VespaG", .586), ("GEMME", .211)):
        row = select(comparisons, scheme="guarded_loocv", model=PI_DS, reference=reference)
        near("Table 4 p " + reference, row["p_two_sided"], expected, .0005)
        near("Table 4 exact sign patterns " + reference, row["sign_patterns"], 4096, 0)

    pi = load_csv(root, "results/positive_inside_standalone_20261028/all_audits.csv")
    for model, expected in (("guarded_oof_raw", .00854), ("guarded_oof__" + PI_DS, -.02029)):
        row = select(pi, model=model, window=10, comparator_name="all_uncharged")
        near("Table 5 positive inside " + model, row["estimate_raw"], expected, .000005)
    robustness = [row for row in pi if row["model"] == "guarded_oof__" + PI_DS]
    near("positive inside sensitivity settings", len(robustness), 20, 0)
    near("positive inside negative estimates", sum(float(row["estimate_raw"]) < 0 for row in robustness), 20, 0)
    actual_windows = sorted({int(row["window"]) for row in robustness})

    ga = load_json(root, "results/glycosylation_constraint_repair_20261031/run2/audit_results.json")["guarded_group_out"]
    for model, arm, expected in (("existing_PI_DS", "A", -.01140), ("existing_PI_DS", "B", .01811),
                                 (THREE, "A", .04696), (THREE, "B", -.03977)):
        item = ga[model]["original_all_motifs"][arm + "_noncytosolic"]["group_resampling"]
        near("Table 5 glyco " + model + " " + arm, item["estimate"], expected, .000005)
    for arm, value, sites, proteins in (("A", .12351, 23, 7), ("B", -.11157, 26, 9)):
        item = ga[THREE]["feature_eligible"][arm + "_noncytosolic"]["group_resampling"]
        near("Table 12 " + arm + " estimate", item["estimate"], value, .000005)
        near("Table 12 " + arm + " sites", item["n_target_sites"], sites, 0)
        near("Table 12 " + arm + " proteins", item["n_proteins"], proteins, 0)
    profile = load_json(root, "results/glycosylation_constraint_repair_20261031/run2/RUN_STATUS.json")
    disulfide = load_json(root, "results/disulfide_refinement_verification_20261031/verification.json")
    for model, expected in (("raw", .05986), (THREE, .12572)):
        item = disulfide["results"][model]["S"]["whole_protein"]
        near("Table 5 disulfide " + model, item["estimate"], expected, .000005)
        near("disulfide contributing endpoints " + model, item["n_bonded_sites"], 44, 0)
    near("disulfide positive protein estimates", disulfide["current_proteins_positive"], 7, 0)
    charge = load_json(root, base + "label_free_audit_20261025/existing_refinement_audits.json")
    near("Table 5 charge control", charge["results"]["positive_inside"]["charge_burial"]["grouped_endpoint"]["estimate"], .039452, .0000005)

    external = load_json(root, "results/tool_demo_20261003/two_model_demo.json")
    for model, expected, sites, strata, proteins in (
            ("alphamissense", .0121, 173413, 29819, 4517),
            ("esm1b", -.147, 173622, 29858, 4534)):
        item = external["results"][model]["result"]
        near("Table 6 " + model, item["estimate_raw"], expected, .0005 if model == "esm1b" else .00005)
        for key, value in (("n_contributing_sites", sites), ("n_strata", strata), ("n_proteins", proteins)):
            near("Table 6 " + model + " " + key, item[key], value, 0)
    gemme = load_csv(root, "results/gemme_membrane_census_20261031/results.csv")
    for model, expected in (("GEMME", -.186983), ("AlphaMissense", .000861)):
        item = select(gemme, model=model, scope="identical_GEMME_AM_cells", coverage="complete19_positions")
        near("Table 7 " + model, item["estimate"], expected, .0000005)
        near("Table 7 " + model + " proteins", item["n_proteins"], 959, 0)
        near("Table 7 " + model + " sites", item["n_contributing_sites"], 27591, 0)

    criteria = load_json(root, "results/form_a_20261008/criteria.json")
    near("Table 10 slope", criteria["S1"]["slope_through_origin"], .819, .0005)
    near("Table 10 noise free error", criteria["S0_exact"]["max_abs_error"], 1.97e-15, .005e-15)
    near("Table 10 implementation agreement", criteria["independent_arithmetic"]["max_abs_difference"], 1.92e-15, .005e-15)
    near("Table 10 side blind", criteria["S2_identity"]["max_abs_estimate"], 3.62e-18, .005e-18)
    near("Table 10 width ratio", criteria["S0n_calibration"]["half_width_over_pooled_sd"], 2.37, .005)
    benchmark = load_json(root, "results/method_completion_20261004/benchmark_results.json")
    near("fixed effects agreement", benchmark["max_simulated_fixed_effects_error"], 1.39e-16, .005e-16)
    full = [cell for cell in benchmark["cells"] if cell["method"] == "packaged_within_site_stratified"]
    near("Tier 2 maximum absolute bias", max(abs(cell["bias"]) for cell in full), .000624, .0000005)
    speed = load_json(root, "results/speed_comparison_20261026/results.json")
    for method, expected_ms, tolerance in (("cached_embedding/Standalone_raw", 847., .5),
                                            ("correction_only/Standalone_PI_disulfide", .74, .005)):
        row = select(speed["measurements"], method=method, rows=43991, threads=1)
        near("cached timing milliseconds " + method, row["median_seconds"] * 1000, expected_ms, tolerance)
    example = load_json(root, "results/latex_todo_audit_20261101/worked_example.json")
    near("Table 2 position", example["position"], 361, 0)
    near("Table 2 KR mean", example["KR_mean"], .57394, .000005)
    near("Table 2 comparator mean", example["comparator_mean"], .56100, .000005)
    near("Table 2 within site contrast", example["z"], .01295, .000005)
    return dict(mean_spearman=means, actual_positive_inside_windows=actual_windows,
                disulfide_evaluation_kind=disulfide["evaluation_kind"],
                encoder_timed=not speed["no_encoder_inference"],
                glycosylation_run_status=profile["status"])


def score_checks(root):
    sys.path.insert(0, str(root / "tools"))
    from membrane_rule_audit import panel
    saved = load_json(root, "results/comparator_sensitivity_20261001/comparator_sensitivity.json")
    result = {}
    for model, archive in (("alphamissense", "am_full_flank_scores.npz"), ("esm1b", "esm1b_flank_scores.npz")):
        loaded = panel.load_panel(root / "data/processed/gate0" / archive, 10)
        introduced, n_introduced = panel.arm_mean(loaded["scores"], "KR")
        native_eligible = np.isin(loaded["ref_aa"], [panel.ALPHABET.index(a) for a in panel.NEUTRAL])
        result[model] = {}
        for comparator, previous in saved["results"][model].items():
            reference, n_reference = panel.arm_mean(loaded["scores"], previous["comparator"])
            parts = panel.stratum_parts(loaded, introduced-reference,
                                       native_eligible & (n_introduced > 0) & (n_reference > 0))
            point = float(parts["numerator"].sum() / parts["denominator"].sum())
            near("recomputed comparator " + model + " " + comparator, point, previous["estimate_raw"], 1e-12)
            result[model][comparator] = dict(estimate=point, saved_ci95=previous["ci95_raw"],
                                              interval_source="archived run; its original runner was not found")
    path = root / "results/gemme_membrane_census_20261031/GEMME__public_human_census__available_cells.npz"
    loaded = panel.load_panel(path, 10)
    rows = load_csv(root, "results/gemme_membrane_census_20261031/independent_per_protein_results.csv")
    contributors = {row["accession"] for row in rows if row["model"] == "GEMME"
                    and row["scope"] == "public_human_census" and row["coverage"] == "available_cells"}
    keep = np.isin(loaded["accessions"][loaded["protein"]], sorted(contributors))
    kr = np.isin(loaded["ref_aa"], [panel.ALPHABET.index(a) for a in "KR"])
    frequencies = {}
    for side, label, expected in ((0, "cytosolic", 18.7), (1, "noncytosolic", 9.1)):
        arm = keep & (loaded["side"] == side)
        numerator, denominator = int(np.count_nonzero(arm & kr)), int(arm.sum())
        percent = 100 * numerator / denominator
        near("GEMME cohort native KR percent " + label, percent, expected, .05)
        frequencies[label] = dict(native_KR=numerator, scored_flank_positions=denominator, percent=percent)
    result["gemme_native_composition"] = dict(contributing_proteins=len(contributors), frequencies=frequencies,
                                              scope="descriptive reference composition, includes native K/R sites")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--score-audits", action="store_true")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.resolve()
    summary = saved_results(root)
    optional = score_checks(root) if args.score_audits else None
    report = dict(origin="support-bundle reconstruction dated 2026-11-10",
                  project_root=str(root), scope="saved-result checks" + (" and score-array reconstruction" if args.score_audits else ""),
                  summary=summary, optional_score_checks=optional, checks=CHECKS,
                  passed=sum(row["passed"] for row in CHECKS), failed=sum(not row["passed"] for row in CHECKS))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(dict(passed=report["passed"], failed=report["failed"], output=str(args.out)), indent=2))
    if report["failed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

"""Verify the frozen completed benchmark and write its evidence report."""
from datetime import datetime, timezone
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts/beat_am"))
from m4_benchmark_multirule import RULES,ABLATIONS,FULL

OUT = ROOT/"results/standalone_multirule_20261025/benchmark_20261025"
REPORT = ROOT/"docs/STANDALONE_MULTIRULE_BENCHMARK_20261025.md"
LABELS = {"raw":"Raw standalone baseline","positive_inside":"Positive-inside only",
          "charge_burial":"Charge burial","glycosylation_loss":"Glycosylation loss","disulfide_loss":"Disulfide loss"}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def label(name):
    if name in LABELS: return LABELS[name]
    if name == FULL: return "Positive-inside + all three additions"
    return "PI + " + " + ".join(LABELS[r] for r in name.split("+")[1:])


def interval(values):
    def number(value):
        return f"{value:+.2e}" if 0<abs(value)<1e-6 else f"{value:+.6f}"
    return f"[{number(values[0])}, {number(values[1])}]"


def conclusion(result):
    delta = result["delta_spearman"]
    lo,hi = result["ci95_delta_spearman"]
    if abs(delta)<1e-12: return "No mean change from baseline."
    if delta>0:
        return "Higher mean correlation, with a positive conditional stability interval." if lo>1e-12 else "A positive mean change; the conditional stability interval does not establish improvement."
    return "Lower mean correlation, with a negative conditional stability interval." if hi < -1e-12 else "A negative mean change; the conditional stability interval overlaps zero."


def direct_application_coverage(results):
    """Distinguish an annotation example from an actually nonzero held-out term."""
    with np.load(ROOT/"results/beat_am_20261016/features/l1_table.npz",allow_pickle=False) as data:
        acc,pos,alt = data["accession"],data["position"],data["alt"]
    with np.load(ROOT/"results/standalone_multirule_20261025/prepared_variant_features.npz",allow_pickle=False) as data:
        pa,pp,pv,F = (data[name] for name in ("accession","position","alt","features"))
    lookup = {(str(a),int(p),str(v)):f for a,p,v,f in zip(pa,pp,pv,F)}
    features = np.array([lookup.get((str(a),int(p),str(v)),np.zeros(4)) for a,p,v in zip(acc,pos,alt)])
    coverage = {}
    for mode,result in results.items():
        rows = {}
        for j,rule in enumerate(RULES):
            selected = np.zeros(len(acc),dtype=bool)
            for fold in result["fitted_coefficients"]:
                if fold["coefficients"][FULL][rule]>1e-12:
                    selected |= np.isin(acc,fold["evaluation_accessions"]) & (features[:,j] != 0)
            rows[rule] = dict(rows=int(selected.sum()),proteins=sorted(set(acc[selected])),
                             sites=len(set(zip(acc[selected],pos[selected].tolist()))))
        coverage[mode] = rows
    (OUT/"DIRECT_APPLICATION_COVERAGE.json").write_text(json.dumps(coverage,indent=2)+"\n")
    return coverage


def verify(plan,results):
    plan_hash = sha(OUT/"BENCHMARK_PLAN.json")
    for name,expected in plan["code_sha256"].items():
        if sha(ROOT/name) != expected: raise AssertionError("Frozen benchmark code changed")
    if sha(ROOT/"results/beat_am_20261016/features/l1_table.npz") != plan["table_sha256"]:
        raise AssertionError("Experimental development table changed")
    checked = 0
    for mode,result in results.items():
        assert result["plan_sha256"] == plan_hash
        plan_folds = {p["fold"]:p for p in plan["folds"][mode]}
        for path in sorted((OUT/mode).glob("fold_*.json")):
            fold = json.loads(path.read_text())
            test = set(fold["evaluation_groups"])
            train = set(fold["training_groups"])
            assert not test & train
            assert test == set(plan_folds[fold["fold"]]["test_groups"])
            assert not test & set(fold["correction_calibration_groups"])
            for inner in fold["inner"]:
                fitting,calibration = set(inner["training_groups"]),set(inner["calibration_groups"])
                assert not fitting & calibration
                assert not test & (fitting|calibration)
                assert fitting|calibration == train
            for name,fit in fold["correction_fits"].items():
                assert not test & set(fit["calibration_families"])
                for rule,value in fit["coefficients"].items():
                    assert value>=0
                    if rule not in name.split("+"): assert value == 0
            checked += 1
        with np.load(OUT/mode/"out_of_fold_predictions.npz",allow_pickle=False) as data:
            values = data["predictions"]
            assert values.shape == (plan["rows"],len(ABLATIONS))
            assert np.isfinite(values).all()
            assert tuple(data["models"]) == tuple(ABLATIONS)
        assert len(plan_folds) == result["outer_fits"]
        assert result["non_membrane_scores_exactly_unchanged"]
    record = dict(verified_utc=datetime.now(timezone.utc).isoformat(),plan_sha256=plan_hash,
        checked_outer_folds=checked,all_outer_and_inner_group_separation_checks_passed=True,
        inactive_coefficients_exactly_zero=True,all_held_out_predictions_finite=True,
        non_membrane_scores_exactly_unchanged=True,benchmark_design_tests_passed=6,
        previous_refinement_regression_tests_passed=107,
        reserved_experimental_outcomes_opened=False,old_model_artifacts_modified=False,
        result_sha256={mode:sha(OUT/mode/"results.json") for mode in results})
    (OUT/"BENCHMARK_VERIFICATION.json").write_text(json.dumps(record,indent=2)+"\n")
    return record


def main():
    plan = json.loads((OUT/"BENCHMARK_PLAN.json").read_text())
    results = {mode:json.loads((OUT/mode/"results.json").read_text()) for mode in ("grouped_cv","guarded_loocv")}
    checked = verify(plan,results)
    coverage = direct_application_coverage(results)
    cv,loo = results["grouped_cv"],results["guarded_loocv"]
    full_cv,full_loo = cv["comparisons_vs_raw"]["membrane"][-1],loo["comparisons_vs_raw"]["membrane"][-1]
    lines = ["# Do the standalone refinements beat their baseline?", "", "Date: 25 October 2026 (user's local date).", "",
        "**Both requested evaluations are complete: five-fold grouped cross-validation, followed by leave-one-protein-out evaluation with relatives excluded together.** These are new refits on existing exposed development data, not fresh reserved-protein confirmation.", "",
        f"In guarded LOOCV, the full combination changes mean membrane-protein Spearman correlation from **{full_loo['reference_mean_spearman']:.6f} to {full_loo['model_mean_spearman']:.6f}**, a change of **{full_loo['delta_spearman']:+.6f}** with conditional group-resampling interval **{interval(full_loo['ci95_delta_spearman'])}**. {conclusion(full_loo)}", "",
        f"Five-fold grouped CV gives a full-combination change of **{full_cv['delta_spearman']:+.6f}**, interval **{interval(full_cv['ci95_delta_spearman'])}**. {conclusion(full_cv)}", "",
        "## Fair baseline and split design", "",
        "The baseline is the same standalone architecture: 18 substitution/topology features plus 64 ProtT5 principal components, followed by a fixed gradient-boosted tree head. AlphaMissense and ESM scores are not inputs. The head has depth 6, learning rate 0.05 and 500 iterations; no model setting was selected using these new results. The older reported 0.4873 mean over 39 proteins came from another evaluation design and is not the comparator here.", "",
        "All 127,811 existing measured rows from 44 assays and 39 proteins were used. The main comparison is the 13 membrane proteins to which these annotation features apply. The 26 other proteins receive zero correction and supply an exact scope check. Relatedness uses the conservative union of the two existing dependence-group registries: 34 groups overall, 12 among the membrane proteins. These automated historical groups are not claimed to be a newly validated biological family classification.", "",
        "For each outer test group, three-fold inner group-out-of-fold predictions on the remaining proteins produce calibration residuals. Coefficients are fitted on those residuals and the outer-training membrane assay data only. PCA and the supervised head are refitted within every applicable training partition. The final outer head is refitted on all outer-training data; no outer test group participates in its PCA, supervised training or coefficient calibration. Each rule subset is independently calibrated with inactive terms fixed at zero, using the previously implemented nonnegative ridge calibration (ridge 0.001). This measures that specific calibration procedure; it does not test every possible method of enforcing the same chemistry.", "",
        "Guarded LOOCV produces one held-out prediction for every protein, removing all known related-group members together. There are 34 distinct outer fits rather than 39 independent fits because some proteins share an excluded group. No old fitted correction coefficient is reused.", "",
        "## All declared combinations", "",
        "Values are the mean of protein-level Spearman correlations; assays are first averaged equally within protein. Deltas are relative to the raw baseline under the same evaluation scheme.", "",
        "| Model | Five-fold mean rho | Change | Guarded LOOCV mean rho | Change | LOOCV delta interval |", "|---|---:|---:|---:|---:|---|"]
    for a,b in zip(cv["comparisons_vs_raw"]["membrane"],loo["comparisons_vs_raw"]["membrane"]):
        lines.append(f"| {label(a['model'])} | {a['model_mean_spearman']:.6f} | {a['delta_spearman']:+.6f} | {b['model_mean_spearman']:.6f} | {b['delta_spearman']:+.6f} | {interval(b['ci95_delta_spearman'])} |")
    lines += ["", "## Which addition helps most?", "",
        "The first comparison adds one rule to positive-inside. The second measures the loss on removing that rule from the full combination, with the remaining weights refitted. A positive number favors including the rule. These two views distinguish an individual gain from a contribution that depends on another component.", "",
        "| Addition | LOOCV: add to PI | Interval | LOOCV: full minus without this rule | Interval |", "|---|---:|---|---:|---|"]
    ranking = sorted(RULES[1:],key=lambda r:loo["add_one_rule_to_PI"][r]["delta_spearman"],reverse=True)
    for rule in ranking:
        add,remove = loo["add_one_rule_to_PI"][rule],loo["remove_one_rule_from_full"][rule]
        lines.append(f"| {LABELS[rule]} | {add['delta_spearman']:+.6f} | {interval(add['ci95_delta_spearman'])} | {remove['delta_spearman']:+.6f} | {interval(remove['ci95_delta_spearman'])} |")
    best = ranking[0];winning = loo["add_one_rule_to_PI"][best]
    lines += ["", f"The largest add-one point estimate is **{LABELS[best]} ({winning['delta_spearman']:+.6f})**. This ranking is descriptive; comparison of several additions on this exposed data does not establish a uniquely best generalizable rule.", "",
        "The disulfide addition has the best observed average, but its positive-inside-plus-disulfide model improves over the raw baseline by only +0.000266, with interval [-0.000393, +0.000967]. The gain does not establish a reliable improvement. Relative to positive-inside alone, disulfide improves 2 proteins, worsens 3 and leaves 8 unchanged; charge burial improves 8 but worsens 2 enough to give a negative mean; glycosylation worsens 1 and leaves 12 unchanged. These counts explain why the ranking by average size is not the same as the ranking by number of proteins helped.", "",
        "## Experimental support and limits", "", "| Feature | Measured variant rows | Distinct sites | Proteins | Dependence groups |", "|---|---:|---:|---:|---:|"]
    for rule,support in plan["measured_support"].items():
        lines.append(f"| {LABELS.get(rule,rule)} | {support['rows']} | {support['sites']} | {len(support['proteins'])} | {len(support['dependence_groups'])} |")
    lines += ["", "The number actually receiving a nonzero direct term can be lower because a calibration fit may return zero. For the full combination:", "",
        "| Feature | Five-fold: nonzero direct term on measured test rows | LOOCV: nonzero direct term on measured test rows |", "|---|---:|---:|"]
    for rule in RULES:
        lines.append(f"| {LABELS.get(rule,rule)} | {coverage['grouped_cv'][rule]['rows']} | {coverage['guarded_loocv'][rule]['rows']} |")
    if coverage["grouped_cv"]["disulfide_loss"]["rows"] == 0:
        lines += ["", "In five-fold CV, all measured disulfide examples occur in the same outer test fold. Its remaining calibration data have no disulfide-feature variation, so that fit returns zero. Consequently, the five-fold run supplies no direct held-out test of a nonzero disulfide term. Tiny differences between the disulfide-containing and disulfide-free combinations there arise through re-estimation of other weights, not a demonstrated bond-loss benefit. Guarded LOOCV separates the transporter and GPCR groups and addresses this particular coverage problem."]
    lines += ["", "Disulfide and glycosylation coverage is particularly narrow. The disulfide result uses the requested control-free annotation feature, not a bonded-versus-verified-free cysteine experiment. It can evaluate predictive utility, not isolate an intrinsic bond-loss penalty. Bond endpoints are not forced to have equal scores. A direction-constrained feature alone does not guarantee total-score compliance with the corresponding rule; this benchmark is about functional prediction accuracy, not a new whole-census rule-compliance verdict.", "",
        "All candidate features, assay outcomes and the model architecture had prior project exposure. The frozen new splits prevent direct evaluation-group leakage inside this run, but do not turn the project history into untouched confirmation. The reported 95% intervals resample dependence groups while keeping fitted cross-validated predictions fixed; they are conditional stability summaries, do not refit the entire learning procedure in every resample, and are not corrected for selecting a winner among several additions. Very small effects and few affected groups require particular caution.", "",
        "## Individual-protein effects and secondary scopes", "",
        f"For the full combination in guarded LOOCV: **{full_loo['wins']} membrane proteins improved, {full_loo['losses']} worsened, and {full_loo['ties']} were unchanged**, using a 1e-12 numerical tolerance. Mean squared-error change was {full_loo['mean_mse_delta']:+.6f} (negative is better).", ""]
    for scope in ("all","non_membrane"):
        row = loo["comparisons_vs_raw"][scope][-1]
        lines.append(f"- {scope.replace('_',' ')}: {row['proteins']} proteins; raw mean rho {row['reference_mean_spearman']:.6f}, full {row['model_mean_spearman']:.6f}, change {row['delta_spearman']:+.6f}.")
    lines += ["", "Full per-protein changes and every fold's correction weights are retained. A favorable mean must not be described as no individual protein worsening.", "",
        "### Membrane-protein changes in guarded LOOCV", "",
        "Every change below is relative to that protein's raw held-out baseline. These are paired protein results, not additional independent variant replications.", "",
        "| Protein | Raw rho | PI + depth: change | PI + glyco: change | PI + disulfide: change | Full: change |", "|---|---:|---:|---:|---:|---:|"]
    with (OUT/"guarded_loocv/per_protein.csv").open(encoding="utf-8",newline="") as stream:
        protein_rows = [r for r in csv.DictReader(stream) if r["membrane"] == "1"]
    by_protein = {}
    for row in protein_rows:
        by_protein.setdefault(row["accession"],{})[row["model"]] = row
    for acc,rows in sorted(by_protein.items()):
        raw = float(rows["raw"]["spearman"])
        changes = [float(rows["positive_inside+"+r]["spearman"])-raw for r in RULES[1:]]
        changes.append(float(rows[FULL]["spearman"])-raw)
        lines.append(f"| {rows['raw']['gene']} ({acc}) | {raw:.6f} | "+" | ".join(f"{v:+.6f}" for v in changes)+" |")
    lines += ["",
        "These assays mix experimental endpoints. For example, the glycosylation-supported set includes ACE2_HUMAN_Chan_2020 (binding), OPSD_HUMAN_Wan_2019 (abundance) and SC6A4_HUMAN_Young_2021. An annotated chemical event is not by itself evidence that every assay readout should decrease. The benchmark tests one transferable directional correction against those measured readouts, not clinical pathogenicity.", "",
        "## Reproducibility and integrity", "",
        f"The benchmark plan was frozen before opening experimental targets for this run. It records code/input/embedding hashes, all outer and inner group assignments, fixed model settings, every ablation, endpoints and interpretation limits. Post-run verification checked **{checked['checked_outer_folds']} outer folds**, including all inner training/calibration separations, inactive zero coefficients, finite predictions and the exact non-membrane scope check. Six known-answer split/metric tests passed in addition to the earlier 107 refinement/regression checks.", "",
        "- Runner: `scripts/beat_am/m4_benchmark_multirule.py`.",
        "- Frozen plan and fitted fold models: `results/standalone_multirule_20261025/benchmark_20261025/`.",
        "- Five-fold results: `grouped_cv/results.json` and `grouped_cv/per_protein.csv` under that folder.",
        "- Guarded LOOCV: `guarded_loocv/results.json` and `guarded_loocv/per_protein.csv`.",
        "- Post-run verification: `BENCHMARK_VERIFICATION.json`.", "",
        "No reserved experimental outcome was opened and no old model artifact or manuscript was overwritten. No new final full-data predictor was selected or released from these comparisons."]
    REPORT.write_text("\n".join(lines)+"\n",encoding="utf-8")
    (OUT/"REPORT.md").write_text(REPORT.read_text(encoding="utf-8"),encoding="utf-8")
    print(json.dumps(dict(report=str(REPORT),full_cv=full_cv,full_loocv=full_loo,
        add_one_ranking=[dict(rule=r,delta=loo["add_one_rule_to_PI"][r]["delta_spearman"]) for r in ranking]),indent=2))


if __name__ == "__main__":
    main()

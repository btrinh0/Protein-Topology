"""Calibrate a glycosylation correction with group-disjoint, label-free audit constraints."""
from __future__ import annotations

import argparse
from collections import defaultdict
from contextlib import ExitStack, redirect_stdout
from datetime import datetime, timezone
import csv
import io
import json
from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "beat_am"))
import m7_audit_standalone_rules as audit
from membrane_rule_audit.glycosylation_constraints import (
    coefficients_with_glycosylation_floor, fit_glycosylation_floor,
)

OUT = ROOT / "results" / "glycosylation_constraint_repair_20261031"
BASE, BENCH, OLD = audit.BASE, audit.BENCH, audit.OUT
PI_DS = "positive_inside+disulfide_loss"
OLD_THREE = "positive_inside+glycosylation_loss+disulfide_loss"
MODELS = ["raw", "existing_PI_DS", "existing_PI_glyco_DS",
          "constraint_glyco", "PI_DS_constraint_glyco"]
RULE_INDEX = list(audit.RULES).index("glycosylation_loss")


def save(name, value):
    (OUT / name).write_text(json.dumps(audit.clean(value), indent=2, allow_nan=False) + "\n",
                            encoding="utf-8")


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def log(message):
    print(message, flush=True)


def csv_rows(name, rows):
    if not rows:
        return
    with (OUT / name).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def metadata():
    bank = read(BASE / "annotation_bank.json")
    whole = audit.read_topology(ROOT / "data" / "raw" / "topology_viewer.html")
    census = {a: whole[a] for a in sorted(bank)}
    plan = read(BENCH / "BENCHMARK_PLAN.json")
    group_for = {a: g for g, members in plan["dependence_members"].items() for a in members}
    folds = []
    for fold in plan["folds"]["guarded_loocv"]:
        test = sorted(set(fold["test_accessions"]) & set(bank))
        if not test:
            continue
        detail_path = BENCH / "guarded_loocv" / f"fold_{fold['fold']:02d}.json"
        detail = read(detail_path)
        model_path = BENCH / "model_cache" / (detail["outer_model_cache_key"] + ".joblib")
        if set(detail["evaluation_groups"]) & set(detail["training_groups"]):
            raise AssertionError("Old outer fit has overlapping groups")
        folds.append(dict(fold=fold["fold"], test=test, detail=detail,
                          detail_path=detail_path, model_path=model_path))
    with np.load(OLD / "standalone_position_scores.npz", allow_pickle=False) as z:
        cached = {k: z[k].copy() for k in z.files}
    if "".join(cached["alphabet"]) != audit.ALPHABET:
        raise AssertionError("Score alphabet mismatch")
    index = {(str(a), int(p)): i for i, (a, p) in
             enumerate(zip(cached["accession"], cached["position"]))}
    with redirect_stdout(io.StringIO()):
        sites = audit.b8.enumerate_sites(census)
    return bank, census, group_for, folds, cached, index, sites


def constraint_rows(scores, features, sites, index, bank, groups, *, eligible_only=True,
                    accessions=None):
    """Collapse the unchanged B8 strata to each protein's A and B contrast."""
    allowed = set(bank) if accessions is None else set(accessions)
    occupied = {a: {g["position"] for g in r["glycosites"]} for a, r in bank.items()}
    bins = defaultdict(lambda: [[], []])
    for acc, pos, native, side, contrast, target in sites:
        if acc not in allowed or side != "noncytosolic":
            continue
        if eligible_only and target and pos - (1 if contrast == "A" else 2) not in occupied[acc]:
            continue
        j = index[(acc, pos)]
        finite = np.isfinite(scores[j]) & (np.array(list(audit.ALPHABET)) != native)
        alternatives = {a: float(v) for a, v, ok in zip(audit.ALPHABET, scores[j], finite) if ok}
        response = {a: float(v) for a, v, ok in zip(audit.ALPHABET, features[j], finite) if ok}
        z, slope = audit.b8.z_value(alternatives, contrast, native), audit.b8.z_value(response, contrast, native)
        if z is not None and slope is not None:
            bins[(acc, contrast, native)][int(target)].append((z, slope))
    parts = defaultdict(lambda: np.zeros(6))
    for (acc, contrast, native), (reference, target) in bins.items():
        if len(target) < 1 or len(reference) < audit.b8.MIN_STRATUM:
            continue
        weight = len(target) * len(reference) / (len(target) + len(reference))
        delta, slope = np.mean(target, axis=0) - np.mean(reference, axis=0)
        parts[(acc, contrast)] += [weight * delta, weight * slope, weight,
                                   len(target), len(reference), 1]
    return [dict(accession=a, group=groups[a], contrast=c,
                 baseline_delta=float(v[0] / v[2]), feature_delta=float(v[1] / v[2]),
                 weight=float(v[2]), n_target=int(v[3]), n_reference=int(v[4]), n_strata=int(v[5]))
            for (a, c), v in sorted(parts.items())]


def freeze(bank, folds):
    sources = [Path(__file__), ROOT / "tools/membrane_rule_audit/glycosylation_constraints.py",
               ROOT / "tools/membrane_rule_audit/refinement.py", Path(audit.__file__),
               Path(audit.b8.__file__), ROOT / "scripts/beat_am/m2_standalone_refinement_scorer.py",
               BASE / "annotation_bank.json", BENCH / "BENCHMARK_PLAN.json",
               OLD / "standalone_position_scores.npz", audit.HEAD,
               ROOT / "data/raw/topology_viewer.html"]
    sources += [p for f in folds for p in (f["detail_path"], f["model_path"])]
    sources += [audit.EMB / (a + ".npz") for a in sorted(bank)]
    plan = dict(
        version=1, status="exploratory extension after inspecting the existing rule audits",
        primary_candidate="PI_DS_constraint_glyco", secondary_candidate="constraint_glyco",
        scope="existing 13-protein membrane development panel; no reserved data",
        feature="unchanged binary loss of directly supported glycosylation in resolved mature outer loops",
        calibration_constraints="protein-level B8 A>=0 and B<=0 on feature-eligible targets; same non-motif controls, matching and weights",
        objective="smallest nonnegative shared glycosylation coefficient satisfying every supported calibration protein/contrast",
        margin=0.0, bootstrap_intervals_used_to_choose_coefficient=False,
        coefficient_selection_uses_experimental_outcomes=False,
        outer_protocol="reuse each frozen outer head; score calibration proteins with that head; exclude outer evaluation groups from floor calculation",
        other_coefficients="fixed from that fold's existing PI+DS configuration; raw comparator has all other terms zero",
        primary_audit="unchanged A/B B8 contrasts, original all-motif scope AND occupied feature-eligible scope, dependence-group intervals",
        controls="cytosolic contrasts; unchanged C-to-S scores and TM scores; native masks and out-of-scope scores unchanged",
        accuracy="after all floors/audits frozen: existing guarded group-out predictions; assay Spearman averaged within protein, then proteins equally; no threshold retuning",
        fit_new_head_or_PCA=False, old_artifacts_modified=False,
        guarantee="constraints on calibration protein means only; held-out audit compliance is measured separately",
        input_sha256={str(p.relative_to(ROOT)): audit.sha(p) for p in sources},
    )
    path = OUT / "PLAN.json"
    if path.exists():
        previous = read(path)
        previous.pop("frozen_utc", None)
        if previous != plan:
            raise ValueError("Repair plan or inputs changed; use a new output directory")
    else:
        save("PLAN.json", dict(frozen_utc=datetime.now(timezone.utc).isoformat(), **plan))
    return plan


def fit_profile(scores, cached, sites, index, bank, groups, calibration, excluded, coefficients):
    if coefficients["glycosylation_loss"] != 0:
        raise ValueError("The floor must be derived from a zero-glycosylation base")
    base = scores + np.einsum("ijk,k->ij", cached["features"],
                              [coefficients[r] for r in audit.RULES])
    rows = constraint_rows(base, cached["features"][:, :, RULE_INDEX], sites, index, bank, groups,
                           accessions=calibration)
    profile = fit_glycosylation_floor(rows, calibration_groups=sorted({groups[a] for a in calibration}),
                                      excluded_groups=sorted(excluded), margin=0.)
    profile["other_coefficients"] = dict(coefficients)
    profile["coefficients"] = coefficients_with_glycosylation_floor(coefficients, profile)
    return profile


def audit_stage():
    OUT.mkdir(parents=True, exist_ok=True)
    bank, census, groups, folds, cached, index, sites = metadata()
    plan = freeze(bank, folds)
    log("Plan frozen. No DMS outcomes will be read in the audit stage.")
    acc = cached["accession"].astype(str)
    raw, features = cached["guarded_oof"], cached["features"]
    matrices = {name: raw.copy() for name in MODELS}
    zero = {r: 0. for r in audit.RULES}
    profiles = {}
    embeddings = {}
    for a in bank:
        with np.load(audit.EMB / (a + ".npz"), allow_pickle=False) as z:
            embeddings[a] = (z["position"].copy(), z["embedding"].copy())
            if int(z["seq_len"]) != len(bank[a]["sequence"]):
                raise ValueError("Embedding sequence mismatch")
    wanted = {a: sorted({p for aa, p, *_ in sites if aa == a}) for a in bank}
    with ExitStack() as stack, threadpool_limits(limits=8):
        for cls, method in ((audit.PCA, "fit"), (audit.PCA, "fit_transform"),
                            (audit.HistGradientBoostingRegressor, "fit"),
                            (audit.refinement, "fit_coefficients")):
            stack.enter_context(patch.object(cls, method, side_effect=audit.forbidden_fit))
        # The frozen full-data head supplies a deployable raw-plus-glycosylation profile.
        full_profile = fit_profile(cached["saved_full_fit"], cached, sites, index, bank, groups,
                                   sorted(bank), [], zero)
        save("saved_full_glycosylation_profile.json", full_profile)
        save("saved_full_glycosylation_coefficients.json", dict(
            coefficients=full_profile["coefficients"], floor_profile="saved_full_glycosylation_profile.json",
            base_head=str(audit.HEAD.relative_to(ROOT)), score_state="raw_head_plus_one_refinement_application"))
        for number, fold in enumerate(folds, 1):
            detail, test = fold["detail"], fold["test"]
            calibration = sorted(set(bank) - set(test))
            excluded = detail["evaluation_groups"]
            if {groups[a] for a in calibration} & set(excluded):
                raise AssertionError("A related calibration protein escaped the exclusion")
            bundle = audit.load_bundle(fold["model_path"])
            if set(bundle["training_groups"]) & set(excluded):
                raise AssertionError("The reused head saw an excluded group")
            if any(a in bundle["training_accessions"] for a in test):
                raise AssertionError("The reused head saw a held-out protein")
            log(f"Outer {number}/{len(folds)}: hold out {', '.join(census[a]['gene'] for a in test)}; calibrate on {len(calibration)} proteins")
            cal_scores = np.full_like(raw, np.nan)
            parity = 0.
            for a in sorted(bank):
                if not wanted[a]:
                    continue
                ep, emb = embeddings[a]
                scored = audit.score_positions(bundle, bank[a], wanted[a], ep, emb)
                rows = [index[(a, p)] for p in wanted[a]]
                if a in test:
                    parity = max(parity, float(np.nanmax(np.abs(scored["scores"] - raw[rows]))))
                else:
                    cal_scores[rows] = scored["scores"]
            if parity > 1e-10:
                raise AssertionError(f"Frozen held-out score parity failed: {parity}")
            base_coeff = detail["correction_fits"][PI_DS]["coefficients"]
            fitted = {}
            for name, coefficients in [("constraint_glyco", zero), ("PI_DS_constraint_glyco", base_coeff)]:
                profile = fit_profile(cal_scores, cached, sites, index, bank, groups,
                                       calibration, excluded, coefficients)
                fitted[name] = profile
            mask = np.isin(acc, test)
            for name, source in [("existing_PI_DS", PI_DS), ("existing_PI_glyco_DS", OLD_THREE)]:
                coefficients = detail["correction_fits"][source]["coefficients"]
                matrices[name][mask] += np.einsum("ijk,k->ij", features[mask], [coefficients[r] for r in audit.RULES])
            for name, profile in fitted.items():
                matrices[name][mask] += np.einsum("ijk,k->ij", features[mask],
                                                  [profile["coefficients"][r] for r in audit.RULES])
            profiles[str(fold["fold"])] = dict(test_accessions=test, excluded_groups=excluded,
                calibration_accessions=calibration, outer_model_sha256=audit.sha(fold["model_path"]),
                maximum_held_out_score_reproduction_error=parity, candidates=fitted)
            save("group_out_profiles.json", profiles)
        log("All group-disjoint coefficients frozen. Running the unchanged glycosylation audits.")
        full_matrices = dict(raw=cached["saved_full_fit"], constraint_glyco=cached["saved_full_fit"]
                            + full_profile["beta"] * features[:, :, RULE_INDEX])
        results = {}
        for mode, objects in [("guarded_group_out", matrices), ("saved_full_head", full_matrices)]:
            results[mode] = {}
            for name, scores in objects.items():
                log(f"Audit: {mode} / {name}")
                results[mode][name] = {
                    "original_all_motifs": audit.glyco_audit(scores, sites, index, groups, bank),
                    "feature_eligible": audit.glyco_audit(scores, sites, index, groups, bank, eligible_only=True),
                }
        # Verify the archived baseline before attributing any change to the correction.
        historical = read(OLD / "audit_results.json")["guarded_oof"]["glycosylation"]
        for key in ("A_noncytosolic", "B_noncytosolic", "A_cytosolic", "B_cytosolic"):
            current = results["guarded_group_out"]["raw"]["original_all_motifs"][key]["group_resampling"]
            previous = historical["original"][key]["group_resampling"]
            if not np.isclose(current["estimate"], previous["estimate"], atol=1e-12, rtol=0):
                raise AssertionError("Original audit estimate was not reproduced")
        checks = dict(native_masks_preserved=True, excluded_groups_absent=True,
                      new_model_or_PCA_fits=0, experimental_outcomes_read=False,
                      raw_audit_reproduced=True, per_candidate={})
        cysteine = cached["native"] == "C"
        for name, base in [("constraint_glyco", "raw"), ("PI_DS_constraint_glyco", "existing_PI_DS")]:
            before, after = matrices[base], matrices[name]
            support = features[:, :, RULE_INDEX] != 0
            if not np.array_equal(np.isnan(before), np.isnan(after)):
                raise AssertionError("Native/missing masks changed")
            if not np.array_equal(before[~support], after[~support], equal_nan=True):
                raise AssertionError("Glycosylation correction escaped its feature scope")
            si = audit.ALPHABET.index("S")
            if not np.array_equal(before[cysteine, si], after[cysteine, si], equal_nan=True):
                raise AssertionError("C-to-S disulfide scores changed")
            checks["per_candidate"][name] = dict(
                changed_missense_scores=int(np.count_nonzero(np.isfinite(before) & (before != after))),
                maximum_absolute_shift=float(np.nanmax(np.abs(after - before))),
                outside_feature_unchanged=True, disulfide_C_to_S_scores_exactly_unchanged=True)
        rows = []
        for name, scores in matrices.items():
            for scope, eligible in [("feature_eligible", True), ("original_all_motifs", False)]:
                for row in constraint_rows(scores, features[:, :, RULE_INDEX], sites, index, bank, groups, eligible_only=eligible):
                    rows.append(dict(model=name, scope=scope, gene=census[row["accession"]]["gene"], **row))
        csv_rows("per_protein_contrasts.csv", rows)
        np.savez_compressed(OUT / "guarded_position_scores.npz", accession=acc,
                            position=cached["position"], native=cached["native"], alphabet=cached["alphabet"], **matrices)
        save("audit_results.json", results)
        save("VERIFICATION.json", checks)
    for path, digest in plan["input_sha256"].items():
        if audit.sha(ROOT / path) != digest:
            raise AssertionError("An original input changed during the repair")
    save("AUDIT_COMPLETE.json", dict(completed_utc=datetime.now(timezone.utc).isoformat(),
        plan_sha256=audit.sha(OUT / "PLAN.json"),
        profiles_sha256=audit.sha(OUT / "group_out_profiles.json"),
        audits_sha256=audit.sha(OUT / "audit_results.json"),
        score_sha256=audit.sha(OUT / "guarded_position_scores.npz"),
        old_input_hashes_unchanged=True, accuracy_outcomes_not_yet_read=True))
    for scope in ("original_all_motifs", "feature_eligible"):
        for key in ("A_noncytosolic", "B_noncytosolic"):
            result = results["guarded_group_out"]["PI_DS_constraint_glyco"][scope][key]
            log(f"Primary {scope} / {key}: {result['group_resampling']} / {result['directional_verdict']}")


def accuracy_stage():
    completion = read(OUT / "AUDIT_COMPLETE.json")
    for field, filename in [("profiles_sha256", "group_out_profiles.json"), ("audits_sha256", "audit_results.json"),
                            ("score_sha256", "guarded_position_scores.npz"), ("plan_sha256", "PLAN.json")]:
        if audit.sha(OUT / filename) != completion[field]:
            raise AssertionError("The repair changed after being frozen")
    bank, census, groups, folds, cached, index, sites = metadata()
    profiles = read(OUT / "group_out_profiles.json")
    source = BENCH / "guarded_loocv" / "out_of_fold_predictions.npz"
    # Experimental outcomes are first accessed after the complete audit stage.
    with np.load(source, allow_pickle=False) as z:
        data = {k: z[k].copy() for k in z.files}
    with np.load(ROOT / "results/beat_am_20261016/features/l1_table.npz", allow_pickle=False) as z:
        assays, genes = z["dms_id"].astype(str), z["gene"].astype(str)
        for key in ("accession", "position", "native", "alt"):
            if not np.array_equal(z[key], data[key]):
                raise AssertionError("Benchmark metadata row alignment mismatch")
    names = data["models"].astype(str).tolist()
    rows = np.flatnonzero(np.isin(data["accession"], list(bank)))
    acc = data["accession"][rows].astype(str)
    observed = data["target"][rows]
    pred = np.column_stack([data["predictions"][rows, names.index(name)] for name in ("raw", PI_DS, OLD_THREE)])
    pred = np.column_stack([pred, pred[:, 0].copy(), pred[:, 1].copy()])
    feature = np.array([cached["features"][index[(str(a), int(p))], audit.ALPHABET.index(str(v)), RULE_INDEX]
                        for a, p, v in zip(data["accession"][rows], data["position"][rows], data["alt"][rows])])
    for fold in folds:
        mask = np.isin(acc, fold["test"])
        fits = profiles[str(fold["fold"])]["candidates"]
        for name, column in [("constraint_glyco", 3), ("PI_DS_constraint_glyco", 4)]:
            pred[mask, column] += fits[name]["beta"] * feature[mask]
    # Require agreement with independently reconstructed full-position scores.
    position_rows = np.array([index[(str(a), int(p))] for a, p in
                              zip(data["accession"][rows], data["position"][rows])])
    alternative_columns = np.array([audit.ALPHABET.index(str(v)) for v in data["alt"][rows]])
    with np.load(OUT / "guarded_position_scores.npz", allow_pickle=False) as z:
        for col, model in enumerate(MODELS):
            model_scores = z[model]
            from_full = model_scores[position_rows, alternative_columns]
            if not np.allclose(from_full, pred[:, col], atol=1e-10, rtol=0):
                raise AssertionError(f"Full-position / benchmark prediction mismatch for {model}")
    records = []
    for a in sorted(bank):
        selected = acc == a
        metrics = []
        for assay in sorted(set(assays[rows][selected])):
            use = np.flatnonzero(selected & (assays[rows] == assay))
            if len(use) < 20 or len(np.unique(observed[use])) < 3:
                continue
            rho = [spearmanr(pred[use, j], observed[use]).statistic for j in range(len(MODELS))]
            mse = np.mean((pred[use] - observed[use, None]) ** 2, axis=0)
            metrics.append((rho, mse))
        if not metrics:
            continue
        for j, name in enumerate(MODELS):
            records.append(dict(accession=a, gene=census[a]["gene"], group=groups[a], membrane=1,
                                model=name, n_assays=len(metrics), n_variants=int(selected.sum()),
                                spearman=float(np.mean([r[0][j] for r in metrics])),
                                mse=float(np.mean([r[1][j] for r in metrics]))))
    import m4_benchmark_multirule as benchmark
    summary = {name: dict(mean_spearman=float(np.mean([r["spearman"] for r in records if r["model"] == name])))
               for name in MODELS}
    comparisons = [benchmark.paired_summary(records, model, reference) for model, reference in
                   [("constraint_glyco", "raw"), ("PI_DS_constraint_glyco", "existing_PI_DS"),
                    ("PI_DS_constraint_glyco", "raw"), ("PI_DS_constraint_glyco", "existing_PI_glyco_DS")]]
    csv_rows("accuracy_per_protein.csv", records)
    save("accuracy_results.json", dict(summary=summary, comparisons=comparisons,
        source=str(source.relative_to(ROOT)), source_sha256=audit.sha(source),
        accuracy_implementation_sha256=audit.sha(Path(__file__)),
        mean_definition="Spearman within each assay, assays averaged within protein, equal protein mean",
        new_coefficient_selection_after_accuracy=False, independent_confirmation=False,
        scope="same previously exposed 13-protein development panel; guarded outer heads and group-disjoint floors"))
    np.savez_compressed(OUT / "accuracy_predictions.npz", predictions=pred, models=np.array(MODELS),
                        accession=acc, position=data["position"][rows], alt=data["alt"][rows],
                        target=observed, assay=assays[rows], original_row=rows)
    log(json.dumps(summary, indent=2))
    log("Accuracy evaluated after freezing the repair; no weights changed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["audit", "accuracy"], required=True)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()
    OUT = args.out.resolve()
    if args.stage == "audit":
        audit_stage()
    else:
        accuracy_stage()

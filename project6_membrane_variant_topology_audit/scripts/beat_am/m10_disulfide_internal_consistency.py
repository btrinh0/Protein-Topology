"""Check circular agreement between corrected disulfide scores and annotation-defined
targets.
"""
from contextlib import ExitStack,redirect_stdout
from datetime import datetime,timezone
import json
from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingRegressor

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts/beat_am"))
import m7_audit_standalone_rules as audit

OUT=audit.BASE/"disulfide_internal_consistency_20261025"
REPORT=ROOT/"docs/DISULFIDE_INTERNAL_CONSISTENCY_20261025.md"


def criterion(result):
    if not result: return "not_evaluable"
    lo,hi=result["ci95"]
    return "positive_annotation_association" if lo>0 else "opposite_annotation_association" if hi<0 else "inconclusive_annotation_association"


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    bank=json.loads((audit.BASE/"annotation_bank.json").read_text(encoding="utf-8"))
    plan=json.loads((audit.BENCH/"BENCHMARK_PLAN.json").read_text())
    paths=[];pi={};joint={}
    for fold in plan["folds"]["guarded_loocv"]:
        members=[a for a in fold["test_accessions"] if a in bank]
        if not members: continue
        path=audit.BENCH/"guarded_loocv"/f"fold_{fold['fold']:02d}.json"
        paths.append(path);metadata=json.loads(path.read_text())
        for a in members:
            pi[a]=metadata["correction_fits"]["positive_inside"]["coefficients"]
            joint[a]=metadata["correction_fits"]["positive_inside+disulfide_loss"]["coefficients"]
    source_paths=[audit.OUT/"standalone_position_scores.npz",audit.BASE/"annotation_bank.json",audit.b10.UNIPROT_TSV,*paths]
    source_hashes={str(p.relative_to(ROOT)):audit.sha(p) for p in source_paths}
    with np.load(audit.OUT/"standalone_position_scores.npz",allow_pickle=False) as source:
        acc,pos,F,raw,full=(source[k] for k in ("accession","position","features","guarded_oof","saved_full_fit"))
    beta0=np.array([[pi[str(a)][r] for r in audit.RULES] for a in acc])
    beta1=np.array([[joint[str(a)][r] for r in audit.RULES] for a in acc])
    reference=raw+np.sum(F*beta0[:,None,:],axis=2)
    refined=raw+np.sum(F*beta1[:,None,:],axis=2)
    direct=F[:,:,3]*beta1[:,None,3]
    other=np.sum(F[:,:,:3]*(beta1-beta0)[:,None,:3],axis=2)
    finite=np.isfinite(raw)
    if not np.allclose((refined-reference)[finite],(direct+other)[finite],atol=1e-12):
        raise AssertionError("Correction decomposition failed")
    eligible=(F[:,:,3] != 0)&finite
    if np.any(direct[finite]<0) or np.any(direct[finite&~eligible] != 0):
        raise AssertionError("Direct bond-loss component escaped its annotation mask")
    if not np.array_equal(np.isnan(refined),np.isnan(raw)):
        raise AssertionError("Missing/native identity cells changed")
    census0=audit.read_topology(ROOT/"data/raw/topology_viewer.html")
    census={a:census0[a] for a in bank}
    index={(str(a),int(p)):i for i,(a,p) in enumerate(zip(acc,pos))}
    results={}
    with ExitStack() as stack:
        for cls,method in ((PCA,"fit"),(PCA,"fit_transform"),(HistGradientBoostingRegressor,"fit"),(audit.refinement,"fit_coefficients")):
            stack.enter_context(patch.object(cls,method,side_effect=audit.forbidden_fit))
        with (OUT/"legacy_enumeration.log").open("w",encoding="utf-8") as stream,redirect_stdout(stream):
            uni=audit.b10.load_uniprot();uni={a:uni[a] for a in bank}
            sites=audit.b10.enumerate_sites(census,uni)
        targets={(a,b[k]) for a,r in bank.items() for b in r["disulfides"] for k in ("start","end")}
        target=np.array([(str(a),int(p)) in targets for a,p in zip(sites["acc"],sites["pos"])])
        keep=target|(sites["bonded"]==0)
        strata=audit.b10.stratum_ids(sites,("comp",))
        for name,matrix in (("saved_full_fit_raw",full),("guarded_oof_raw",raw),
                            ("guarded_oof_PI",reference),("guarded_oof_PI_plus_disulfide",refined)):
            values=audit.score_rows(matrix,index,sites["acc"],sites["pos"])
            result={}
            for probe in ("S","A","all19_mean"):
                v=np.nanmean(values,axis=1) if probe=="all19_mean" else values[:,audit.ALPHABET.index(probe)]
                row=audit.b10.matched_gap(sites["prot"][keep],strata[keep],target[keep].astype(int),v[keep],
                    len(sites["prot_names"]),np.random.default_rng(audit.b10.SEED))
                # Preserve the legacy schema without treating missing annotations as evidence of a free cysteine.
                result[probe]=dict(measurement=row,internal_status=criterion(row),
                    background="unannotated cysteines, not verified native free thiols")
            results[name]=result
    before=results["guarded_oof_PI"]["S"]["measurement"]
    after=results["guarded_oof_PI_plus_disulfide"]["S"]["measurement"]
    if not np.isclose(before["estimate"],.05985714789839747,atol=1e-12) or not np.isclose(after["estimate"],.12571905739302744,atol=1e-12):
        raise AssertionError("Existing audit reproduction failed")
    if not all(audit.sha(ROOT/p)==digest for p,digest in source_hashes.items()):
        raise AssertionError("Source changed")
    output=dict(created_utc=datetime.now(timezone.utc).isoformat(),evaluation_kind="circular_internal_annotation_consistency",
        circularity="The same bond annotations activate the positive offset and label the target sites in the metric",
        criterion="positive matched annotated-minus-unannotated gap; intervals are internal resampling summaries, not independent validation",
        protocol_status="post-hoc interpretation of already evaluated outputs; no new biological preregistration",
        controls_are_verified_free=False,results=results,
        direct_component_checks=dict(eligible_missense_cells=int(eligible.sum()),
            directly_shifted_cells=int(np.count_nonzero(direct[finite])),negative_direct_penalties=0,
            noneligible_direct_changes=0,missing_native_cells_preserved=True,decomposition_verified=True),
        new_fits=0,new_experimental_outcome_reads=0,old_artifacts_unchanged=True,source_sha256=source_hashes,
        independent_biological_validation=False,independent_predictive_accuracy_improvement=False)
    (OUT/"internal_consistency.json").write_text(json.dumps(audit.clean(output),indent=2,allow_nan=False)+"\n",encoding="utf-8")
    lines=["# Disulfide: explicitly circular internal-consistency evaluation", "", "Date: 25 October 2026 (user's local date).", "",
        "At the user's request, disulfide is evaluated against its own annotation-defined target. The same annotations that trigger the correction also identify the target sites in the evaluation. This is an internal constraint/association check; it is not independent biological validation or a new predictive-accuracy result.", "",
        "## Criterion and comparison", "",
        "At the same protein and annotated compartment, compare the score for the identical replacement at eligible annotated bond endpoints with other unannotated Cys positions. The internal target is a positive annotated-minus-unannotated gap. Unannotated sites are background for this circular evaluation, not verified free cysteines. The original estimator, minimum reference-arm size and 1,999 whole-protein resamples were retained.", "",
        "| Existing output | C→S gap [internal 95% interval] | Internal result |", "|---|---|---|"]
    for name,result in results.items():
        row=result["S"]["measurement"]
        lines.append(f"| {name} | {row['estimate']:+.6f} [{row['ci95'][0]:+.6f}, {row['ci95'][1]:+.6f}] | {result['S']['internal_status']} |")
    lines += ["", f"For the guarded outputs, the existing PI-plus-disulfide model increases the gap from **{before['estimate']:+.6f} to {after['estimate']:+.6f}**, an increase of **{after['estimate']-before['estimate']:+.6f} ranking-score units**. This matched comparison contains {after['n_bonded_sites']} annotated endpoint sites and {after['n_free_sites']} unannotated background sites across {after['n_proteins']} proteins. The legacy field named `n_free_sites` denotes unannotated background only.", "",
        "The primary saved standalone model already has a positive annotation-defined gap before this correction. Its increase therefore does not establish that a missing biological rule was repaired. The weaker guarded baseline interval crosses zero, and the corrected interval is positive under this deliberately circular criterion.", "",
        "## Substitution sensitivity", "", "| Replacement summary | Guarded PI gap | Guarded PI + disulfide gap |", "|---|---:|---:|"]
    for probe in ("S","A","all19_mean"):
        left,right=results["guarded_oof_PI"][probe]["measurement"],results["guarded_oof_PI_plus_disulfide"][probe]["measurement"]
        lines.append(f"| {probe} | {left['estimate']:+.6f} | {right['estimate']:+.6f} |")
    checks=output["direct_component_checks"]
    lines += ["", "## Scope and arithmetic checks", "",
        f"The direct disulfide component shifted {checks['directly_shifted_cells']} eligible missense scores nonnegatively, out of {checks['eligible_missense_cells']} eligible cells. It changed zero noneligible cells and preserved all missing/native-identity cells. Any extra changes in other substitutions between the two complete model variants are accounted for by the separately recalibrated positive-inside coefficient; the component decomposition was verified to 1e-12 tolerance.", "",
        "A positive offset on the annotation-selected target sites mechanically increases this metric. That is precisely why this check is useful for verifying enforcement but cannot establish that the magnitude is biologically correct. The disulfide pairing statistic remains unchanged by a common endpoint offset, so it is not used as a target for this magnitude-only component.", "",
        "## Supported wording", "",
        "> The disulfide correction increased the annotation-defined bond-loss score contrast and met the internal consistency criterion. Because the same annotations define both the correction and its evaluation, this result demonstrates enforcement of the encoded preference; it does not independently validate the biological effect size.", "",
        "The earlier overall Spearman comparison, experimental targeted-error follow-up, and primary controlled-audit conclusions are unchanged. No new model or coefficient was fitted, and no experimental label table was read. An increase in this contrast is not a percentage improvement in predictive accuracy.", "",
        "Reproduction: `scripts/beat_am/m10_disulfide_internal_consistency.py`. Results and input hashes: `results/standalone_multirule_20261025/disulfide_internal_consistency_20261025/internal_consistency.json`."]
    REPORT.write_text("\n".join(lines)+"\n",encoding="utf-8")
    (OUT/"REPORT.md").write_text(REPORT.read_text(encoding="utf-8"),encoding="utf-8")
    print(json.dumps(dict(report=str(REPORT),before=before["estimate"],after=after["estimate"],
        interpretation="circular internal consistency only",scope_checks=checks,new_fits=0),indent=2))


if __name__=="__main__": main()

"""Check rule-contrast shifts using previously fitted leave-one-group-out coefficients."""
from contextlib import ExitStack,redirect_stdout
from datetime import datetime,timezone
import json
from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingRegressor
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts/beat_am"))
import m7_audit_standalone_rules as audit


def main():
    OUT=audit.OUT;BASE=audit.BASE;BENCH=audit.BENCH
    bank=json.loads((BASE/"annotation_bank.json").read_text(encoding="utf-8"))
    all_census=audit.read_topology(ROOT/"data/raw/topology_viewer.html")
    census={a:all_census[a] for a in bank}
    bplan=json.loads((BENCH/"BENCHMARK_PLAN.json").read_text())
    group_for={a:g for g,members in bplan["dependence_members"].items() for a in members}
    model_names=["positive_inside",*["positive_inside+"+r for r in audit.RULES[1:]]]
    paths=[];fits={name:{} for name in model_names}
    for fold in bplan["folds"]["guarded_loocv"]:
        members=[a for a in fold["test_accessions"] if a in bank]
        if not members: continue
        path=BENCH/"guarded_loocv"/f"fold_{fold['fold']:02d}.json"
        paths.append(path)
        data=json.loads(path.read_text())
        for name in model_names:
            coefficients=data["correction_fits"][name]["coefficients"]
            for acc in members: fits[name][acc]=coefficients
    audit.save("EXISTING_REFINEMENT_PLAN.json",dict(created_utc=datetime.now(timezone.utc).isoformat(),
        purpose="movement of the same label-free statistics under existing calibrated refinements",
        new_fit=False,experimental_outcomes_read=False,coefficient_source="frozen fold metadata only",
        comparisons=model_names,coefficient_metadata_hashes={str(p.relative_to(ROOT)):audit.sha(p) for p in paths},
        inference="audit movement and accuracy movement are separate; no code defect inferred solely from zero accuracy delta"))
    with np.load(OUT/"standalone_position_scores.npz",allow_pickle=False) as source:
        acc,pos,raw,F=(source[k] for k in ("accession","position","guarded_oof","features"))
    index={(str(a),int(p)):i for i,(a,p) in enumerate(zip(acc,pos))}
    matrices={}
    for name in model_names:
        beta=np.array([[fits[name][str(a)][rule] for rule in audit.RULES] for a in acc])
        matrices[name]=raw+np.sum(F*beta[:,None,:],axis=2)
    with ExitStack() as stack,threadpool_limits(limits=8):
        for cls,method in ((PCA,"fit"),(PCA,"fit_transform"),(HistGradientBoostingRegressor,"fit"),(audit.refinement,"fit_coefficients")):
            stack.enter_context(patch.object(cls,method,side_effect=audit.forbidden_fit))
        with (OUT/"existing_refinement_audit.log").open("w",encoding="utf-8") as stream,redirect_stdout(stream):
            gs=audit.b8.enumerate_sites(census)
            cs=audit.b9.enumerate_sites(census)
            uni=audit.b10.load_uniprot();uni={a:uni[a] for a in bank}
            ds=audit.b10.enumerate_sites(census,uni)
            results={}
            for name,rule in [("positive_inside",None),*[('positive_inside+'+r,r) for r in audit.RULES[1:]]]:
                matrix=matrices[name];result={}
                if rule in (None,"glycosylation_loss"):
                    result["glycosylation"]=dict(original=audit.glyco_audit(matrix,gs,index,group_for,bank),
                        feature_eligible=audit.glyco_audit(matrix,gs,index,group_for,bank,eligible_only=True))
                if rule in (None,"charge_burial"):
                    result["charge_burial"]=audit.charge_audit(audit.score_rows(matrix,index,cs["acc"],cs["pos"]),cs,group_for)
                if rule in (None,"disulfide_loss"):
                    result["disulfide"]=audit.disulfide_audit(audit.score_rows(matrix,index,ds["acc"],ds["pos"]),ds,bank,uni)
                results[name]=result
    audit.save("existing_refinement_audits.json",dict(new_fit=False,labels_read=False,results=results))
    print("Existing correction coefficients audited without fitting or outcome reads",flush=True)


if __name__=="__main__": main()

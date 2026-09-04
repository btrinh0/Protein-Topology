"""Audit frozen standalone predictions against sequon, charge-burial and disulfide
criteria.
"""
from __future__ import annotations

from collections import defaultdict
from contextlib import ExitStack,redirect_stdout
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import sys
from unittest.mock import patch

import joblib
import numpy as np
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingRegressor
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
sys.path.insert(0,str(ROOT/"scripts/beat_am"))
sys.path.insert(0,str(ROOT/"tools"))
import b8_sequon_rule as b8
import b9_charge_burial as b9
import b10_disulfide_rule as b10
import b10d_pairing_domain_control as domains
import membrane_rule_audit.refinement as refinement
from membrane_rule_audit.refinement import ALPHABET,RULES,variant_features
from m2_standalone_refinement_scorer import load_bundle,score_positions
from run_gate1 import read_topology

BASE=ROOT/"results/standalone_multirule_20261025"
OUT=BASE/"label_free_audit_20261025"
HEAD=ROOT/"results/beat_am_20261016/speed_benchmark/predictor_free_model.joblib"
EMB=ROOT/"results/beat_am_20261016/_from_instance_gsq37aoc/results/beat_am_20261016/features/embeddings_prot_t5_all"
BENCH=BASE/"benchmark_20261025"
PROBE=.05


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(1024*1024),b""): h.update(block)
    return h.hexdigest()


def clean(value):
    if isinstance(value,dict): return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)): return [clean(v) for v in value]
    if isinstance(value,np.ndarray): return clean(value.tolist())
    if isinstance(value,(np.floating,float)): return float(value) if np.isfinite(value) else None
    if isinstance(value,np.integer): return int(value)
    if isinstance(value,np.bool_): return bool(value)
    return value


def save(name,value):
    (OUT/name).write_text(json.dumps(clean(value),indent=2,allow_nan=False)+"\n",encoding="utf-8")


def log(value):
    print(value,flush=True)


def forbidden_fit(*args,**kwargs):
    raise AssertionError("This audit forbids every model/PCA/correction fit")


def matched_cluster(rows,group_for):
    """B8's fixed estimate, with whole-dependence-group uncertainty sensitivity."""
    buckets=defaultdict(lambda:[[],[]])
    for acc,native,target,value in rows: buckets[(acc,native)][0 if target else 1].append(value)
    parts=defaultdict(lambda:[0.,0.]);sites=[0,0];strata=0;proteins=set()
    for (acc,native),(target,reference) in buckets.items():
        if len(target)<1 or len(reference)<b8.MIN_STRATUM: continue
        w=len(target)*len(reference)/(len(target)+len(reference))
        parts[group_for[acc]][0]+=w*(np.mean(target)-np.mean(reference));parts[group_for[acc]][1]+=w
        sites[0]+=len(target);sites[1]+=len(reference);strata+=1;proteins.add(acc)
    if not parts: return None
    units=sorted(parts);num=np.array([parts[g][0] for g in units]);den=np.array([parts[g][1] for g in units])
    point=float(num.sum()/den.sum());interval=None
    if len(units)>=2:
        counts=np.random.default_rng(b8.SEED).multinomial(len(units),np.full(len(units),1/len(units)),size=b8.N_BOOT)
        interval=np.quantile(counts@num/(counts@den),[.025,.975]).tolist()
    return dict(estimate=point,ci95=interval,n_groups=len(units),n_proteins=len(proteins),n_strata=strata,
                n_target_sites=sites[0],n_reference_sites=sites[1])


def matched_verdict(result,direction):
    if not result or result.get("ci95") is None: return "insufficient_data"
    lo,hi=result["ci95"]
    if direction>0:
        return "expected_direction_supported" if lo>0 else "opposite_direction_supported" if hi<0 else "inconclusive"
    return "expected_direction_supported" if hi<0 else "opposite_direction_supported" if lo>0 else "inconclusive"


def glyco_audit(score,sites,index,group_for,bank,*,eligible_only=False):
    output={}
    occupied={a:{g["position"] for g in r["glycosites"]} for a,r in bank.items()}
    for contrast in ("A","B"):
        for side in ("noncytosolic","cytosolic"):
            rows=[]
            for acc,pos,native,s,c,is_seq in sites:
                if s!=side or c!=contrast: continue
                if eligible_only and is_seq:
                    start=pos-(1 if c=="A" else 2)
                    if start not in occupied[acc]: continue
                values=score[index[(acc,pos)]]
                alternatives={a:float(v) for a,v in zip(ALPHABET,values) if a!=native and np.isfinite(v)}
                z=b8.z_value(alternatives,contrast,native)
                if z is not None: rows.append((acc,native,is_seq,z))
            key=contrast+"_"+side
            original=b8.stratified_delta(rows,np.random.default_rng(b8.SEED))
            grouped=matched_cluster(rows,group_for)
            if original and grouped and not np.isclose(original["estimate"],grouped["estimate"],atol=1e-12):
                raise AssertionError("B8 estimator changed in grouping sensitivity")
            direction=1 if contrast=="A" else -1
            output[key]=dict(legacy_whole_protein=original,group_resampling=grouped,
                directional_verdict=matched_verdict(grouped,direction) if side=="noncytosolic" else "placebo_check_not_equivalence_test")
    return output


def charge_groups(score,sites,group_for):
    c=b9.contrasts(score,sites["nat"].astype(int));selected=c["usable"]
    prot,nat,grp,z=(sites["prot"][selected],sites["nat"][selected],sites["grp"][selected],c["z"][selected])
    N,M,protein,n=b9.stratum_table(prot,nat,grp,z,len(sites["prot_names"]),b9.N_GROUPS)
    num,den,ns,nt,nr=b9.contrast_arm(N,M,protein,n,4,0)
    groups=sorted({group_for[a] for a in sites["prot_names"]})
    gn=np.zeros(len(groups));gd=np.zeros(len(groups))
    for i,a in enumerate(sites["prot_names"]):
        j=groups.index(group_for[a]);gn[j]+=num[i];gd[j]+=den[i]
    live=gd>0;gn,gd=gn[live],gd[live]
    if not len(gn): return None
    interval=None
    if len(gn)>1:
        counts=np.random.default_rng(b9.SEED).multinomial(len(gn),np.full(len(gn),1/len(gn)),size=b9.N_BOOT)
        interval=np.quantile(counts@gn/(counts@gd),[.025,.975]).tolist()
    return dict(estimate=float(gn.sum()/gd.sum()),ci95=interval,n_groups=len(gn),n_strata=ns,
                n_centre_sites=nt,n_interface_sites=nr)


def charge_audit(score,sites,group_for):
    result=b9.run_scorer("standalone",score,sites)
    grouped=charge_groups(score,sites,group_for)
    if grouped and not np.isclose(grouped["estimate"],result["z"]["endpoint_gap"]["paired_stratified_delta"],atol=1e-12):
        raise AssertionError("B9 estimator changed in grouping sensitivity")
    return dict(legacy_results=result,legacy_verdicts=b9.verdicts(result),grouped_endpoint=grouped,
                grouped_endpoint_verdict=matched_verdict(grouped,1))


def eligible_pair_sites(sites,bank):
    result=dict(sites)
    result["bonds_by_prot"]={a:[(b["start"],b["end"],True) for b in r["disulfides"]] for a,r in bank.items() if r["disulfides"]}
    return result


def association(score,sites,bank):
    eligible={(a,b[k]) for a,r in bank.items() for b in r["disulfides"] for k in ("start","end")}
    is_eligible=np.array([(a,int(p)) in eligible for a,p in zip(sites["acc"],sites["pos"])])
    keep=(sites["bonded"]==0)|is_eligible
    stratum=b10.stratum_ids(sites,("comp",))
    return b10.matched_gap(sites["prot"][keep],stratum[keep],is_eligible[keep].astype(int),
        score[keep,ALPHABET.index("S")],len(sites["prot_names"]),np.random.default_rng(b10.SEED))


def disulfide_audit(score,sites,bank,uni):
    original=b10.run_scorer("standalone",score,sites)
    v=score[:,ALPHABET.index("S")]
    per_prot,counts=domains.collect(sites,domains.domain_intervals(uni),v)
    domain_result=dict(counts=counts,n_usable_proteins=len(per_prot))
    if per_prot:
        domain_result.update(subset_pairing=domains.subset_p2(per_prot,v),
            primary=domains.domain_matched(per_prot,v,True),
            position_order_sensitivity=domains.domain_matched(per_prot,v,False))
    eligible=eligible_pair_sites(sites,bank)
    return dict(legacy_results=original,legacy_verdicts=b10.verdicts(original),domain_control=domain_result,
        refinement_eligible_pairing=b10.pairing_test(eligible,v,np.random.default_rng(b10.SEED)),
        refinement_eligible_magnitude_association=association(score,sites,bank),
        magnitude_rule_status="not_identifiable_with_verified_free_controls; legacy free means unannotated",
        pairing_is_not_loss_magnitude=True)


def score_rows(matrix,index,accessions,positions):
    return matrix[[index[(str(a),int(p))] for a,p in zip(accessions,positions)]]


def freeze_and_metadata():
    bank=json.loads((BASE/"annotation_bank.json").read_text(encoding="utf-8"))
    all_census=read_topology(ROOT/"data/raw/topology_viewer.html")
    census={a:all_census[a] for a in sorted(bank)}
    bplan=json.loads((BENCH/"BENCHMARK_PLAN.json").read_text())
    bsha=sha(BENCH/"BENCHMARK_PLAN.json")
    group_for={a:g for g,members in bplan["dependence_members"].items() for a in members}
    paths={};models={}
    for fold in bplan["folds"]["guarded_loocv"]:
        key=hashlib.sha256(json.dumps(dict(groups=sorted(fold["train_groups"]),plan=bsha),sort_keys=True).encode()).hexdigest()[:20]
        for acc in fold["test_accessions"]:
            if acc in bank: models[acc]=BENCH/"model_cache"/(key+".joblib")
    for path in [HEAD,BASE/"annotation_bank.json",ROOT/"data/raw/topology_viewer.html",b10.UNIPROT_TSV,
                 Path(__file__),Path(b8.__file__),Path(b9.__file__),Path(b10.__file__),Path(domains.__file__),
                 ROOT/"scripts/beat_am/m2_standalone_refinement_scorer.py",ROOT/"tools/membrane_rule_audit/refinement.py",
                 *set(models.values()),*(EMB/(a+".npz") for a in bank)]:
        paths[str(path.relative_to(ROOT))]=sha(path)
    plan=dict(audit="label-free standalone score behavior; no new fit or experimental outcome",
        accessions=sorted(bank),n_proteins=len(bank),n_positions=sum(len(r["sequence"]) for r in bank.values()),
        primary_head=str(HEAD.relative_to(ROOT)),secondary_heads="already fitted guarded-LOOCV heads, excluding each evaluated dependence group",
        controls="B8 matched motif controls and cytosolic placebo; B9 same-protein/native depth; B10 pairing, proximity and domain controls",
        unchanged_legacy_audits=["B8","B9","B10","B10D"],
        disulfide_caveat="P2 tests pairing, not loss magnitude; legacy unannotated cysteines are not verified-free controls",
        refinement_probes=dict(coefficient=PROBE,one_at_a_time=True,fitted=False,meaning="feature/statistic alignment, never accuracy or trained refinement"),
        additional_declared_sensitivities=["B8 direct occupied sites versus non-motif controls","B9 length>=17 versus original>=10","B10 refinement-eligible pairing endpoints","same estimator with dependence-group resampling for B8 and B9 endpoint"],
        interpretation="a baseline audit failure plus zero accuracy gain does not alone imply a code bug; check scope, metric sensitivity and calibration",
        no_target_or_DMS_table_read=True,no_predictor_score_values_read=True,no_network=True,source_sha256=paths)
    path=OUT/"AUDIT_PLAN.json"
    if path.exists():
        old=json.loads(path.read_text());old.pop("frozen_utc",None)
        if old!=plan: raise ValueError("Frozen audit inputs/design changed")
    else: save("AUDIT_PLAN.json",dict(frozen_utc=datetime.now(timezone.utc).isoformat(),**plan))
    return bank,census,group_for,models,plan


def make_scores(bank,models):
    accessions=[];positions=[];natives=[];full=[];held=[];phis=[]
    primary=load_bundle(HEAD)
    for acc,record in sorted(bank.items()):
        log("Scoring existing heads: "+acc)
        with np.load(EMB/(acc+".npz"),allow_pickle=False) as z:
            ep,emb,seq_len=z["position"],z["embedding"],int(z["seq_len"])
        if seq_len!=len(record["sequence"]): raise ValueError("Embedding length mismatch")
        pos=np.arange(1,seq_len+1)
        raw=score_positions(primary,record,pos,ep,emb)
        secondary=load_bundle(models[acc])
        if acc in secondary.get("training_accessions",[]): raise AssertionError("Scored held-out protein was fitted")
        oof=score_positions(secondary,record,pos,ep,emb)
        native=np.array(list(record["sequence"]))
        vp=np.repeat(pos,19);vn=np.repeat(native,19)
        va=np.concatenate([np.array([a for a in ALPHABET if a!=n]) for n in native])
        feature=variant_features(record,vp,vn,va)
        tensor=np.zeros((len(pos),20,4))
        for i,n in enumerate(native):
            tensor[i,[j for j,a in enumerate(ALPHABET) if a!=n]]=feature[i*19:(i+1)*19]
        accessions.extend([acc]*len(pos));positions.extend(pos);natives.extend(native)
        full.append(raw["scores"]);held.append(oof["scores"]);phis.append(tensor)
    result=dict(accession=np.array(accessions),position=np.array(positions),native=np.array(natives),
        alphabet=np.array(list(ALPHABET)),saved_full_fit=np.concatenate(full),guarded_oof=np.concatenate(held),features=np.concatenate(phis))
    for name in ("saved_full_fit","guarded_oof"):
        if not np.all(np.isfinite(result[name]).sum(axis=1)==19): raise AssertionError("Missense/native masking failed")
    np.savez_compressed(OUT/"standalone_position_scores.npz",**result)
    return result


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    bank,census,group_for,models,plan=freeze_and_metadata()
    log("Audit plan frozen; original input hashes recorded")
    with ExitStack() as stack,threadpool_limits(limits=8):
        for cls,method in ((PCA,"fit"),(PCA,"fit_transform"),(HistGradientBoostingRegressor,"fit"),(refinement,"fit_coefficients")):
            stack.enter_context(patch.object(cls,method,side_effect=forbidden_fit))
        scores=make_scores(bank,models)
        index={(str(a),int(p)):i for i,(a,p) in enumerate(zip(scores["accession"],scores["position"]))}
        with (OUT/"legacy_audit.log").open("w",encoding="utf-8") as stream,redirect_stdout(stream):
            glyco=b8.enumerate_sites(census)
            charge=b9.enumerate_sites(census)
            uni=b10.load_uniprot();uni={a:uni[a] for a in bank}
            cysteine=b10.enumerate_sites(census,uni)
        save("site_inventory.json",dict(glyco_site_rows=len(glyco),charge_sites=len(charge["pos"]),
            cysteine_sites=len(cysteine["pos"]),cysteine_inventory=cysteine["counts"],
            cys_background_is_verified_free=False,complete_scored_missense_cells=19*len(scores["position"])))
        results={}
        for name in ("saved_full_fit","guarded_oof"):
            log("Auditing raw standalone scores: "+name)
            matrix=scores[name]
            probe_matrices={r:matrix+PROBE*scores["features"][:,:,j] for j,r in enumerate(RULES) if r!="positive_inside"}
            with (OUT/"legacy_audit.log").open("a",encoding="utf-8") as stream,redirect_stdout(stream):
                g=glyco_audit(matrix,glyco,index,group_for,bank)
                ge=glyco_audit(matrix,glyco,index,group_for,bank,eligible_only=True)
                gp=glyco_audit(probe_matrices["glycosylation_loss"],glyco,index,group_for,bank)
                gep=glyco_audit(probe_matrices["glycosylation_loss"],glyco,index,group_for,bank,eligible_only=True)
                cm=score_rows(matrix,index,charge["acc"],charge["pos"])
                cp=score_rows(probe_matrices["charge_burial"],index,charge["acc"],charge["pos"])
                c=charge_audit(cm,charge,group_for)
                cprobe=charge_audit(cp,charge,group_for)
                # Filter short TM segments while preserving the original loop-control pool.
                length={a:{p:int(t["end"])-int(t["start"])+1 for t in census[a]["tmds"] for p in range(int(t["start"]),int(t["end"])+1)} for a in bank}
                keep=np.array([grp==b9.LOOP_GROUP or length[str(a)].get(int(p),0)>=17 for a,p,grp in zip(charge["acc"],charge["pos"],charge["grp"])])
                cs={k:(v[keep] if isinstance(v,np.ndarray) and len(v)==len(keep) else v) for k,v in charge.items()}
                c17=charge_audit(cm[keep],cs,group_for)
                ss=score_rows(matrix,index,cysteine["acc"],cysteine["pos"])
                ssp=score_rows(probe_matrices["disulfide_loss"],index,cysteine["acc"],cysteine["pos"])
                d=disulfide_audit(ss,cysteine,bank,uni)
                dprobe=dict(magnitude_association=association(ssp,cysteine,bank),
                    eligible_pairing=b10.pairing_test(eligible_pair_sites(cysteine,bank),ssp[:,ALPHABET.index("S")],np.random.default_rng(b10.SEED)))
            results[name]=dict(glycosylation=dict(original=g,feature_eligible=ge,fixed_probe=gp,feature_eligible_probe=gep),
                charge_burial=dict(original=c,length17_sensitivity=c17,fixed_probe=cprobe),
                disulfide=dict(original=d,fixed_probe=dprobe))
            save("audit_results.json",results)
        unchanged=all(sha(ROOT/path)==digest for path,digest in plan["source_sha256"].items())
        if not unchanged: raise AssertionError("An input changed during audit")
        save("VERIFICATION.json",dict(completed_utc=datetime.now(timezone.utc).isoformat(),
            new_model_or_pca_fit=False,new_correction_fit=False,fit_methods_blocked_during_execution=True,
            experimental_labels_read=False,new_network_data=False,original_inputs_unchanged=unchanged,
            score_positions=len(scores["position"]),scored_missense_per_head=19*len(scores["position"]),
            audit_plan_sha256=sha(OUT/"AUDIT_PLAN.json"),results_sha256=sha(OUT/"audit_results.json")))
    log("Completed standalone label-free audit; no fits or outcome reads")


if __name__=="__main__": main()

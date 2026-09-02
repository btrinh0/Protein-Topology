"""Benchmark standalone refinements with grouped cross-validation and guarded
leave-one-group-out evaluation.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time

import joblib
import numpy as np
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import GroupKFold
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(Path(__file__).resolve().parent))
from membrane_rule_audit.refinement import RULES,fit_coefficients,apply_refinement,ablation_coefficients,write_json
from m2_standalone_refinement_scorer import FEATURES

TABLE = ROOT/"results/beat_am_20261016/features/l1_table.npz"
MANIFEST = ROOT/"results/beat_am_20261016/inventory_frozen/group_manifest_frozen.csv"
PREP = ROOT/"results/standalone_multirule_20261025"
EMB = ROOT/"results/beat_am_20261016/_from_instance_gsq37aoc/results/beat_am_20261016/features/embeddings_prot_t5_all"
SEED = 20261025
BASE_SEED = 20261016
CONFIG = dict(max_depth=6,learning_rate=.05,max_iter=500,early_stopping=False,random_state=BASE_SEED)
ABLATIONS = list(ablation_coefficients({r:1.0 for r in RULES}))
FULL = "+".join(RULES)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def log(text):
    print(f"[{datetime.now(timezone.utc).strftime('%H:%M:%S')} UTC] {text}",flush=True)


def conservative_groups(accessions, old, frozen):
    """Union the two pre-existing dependence groupings, without score selection."""
    unique = sorted(set(accessions))
    parent = {a:a for a in unique}
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    def join(a,b):
        a,b = find(a),find(b)
        if a != b: parent[max(a,b)] = min(a,b)
    for values in (old,frozen):
        buckets = defaultdict(list)
        for acc,value in zip(accessions,values): buckets[str(value)].append(str(acc))
        for members in buckets.values():
            for member in members[1:]: join(members[0],member)
    return np.array(["DG_"+find(str(a)) for a in accessions])


def fold_plans(groups, mode):
    indices = np.arange(len(groups))
    if mode == "grouped_cv":
        outer = list(GroupKFold(n_splits=5).split(indices,groups=groups))
    elif mode == "guarded_loocv":
        outer = [(np.flatnonzero(groups != g),np.flatnonzero(groups == g)) for g in sorted(set(groups))]
    else:
        raise ValueError("Unknown evaluation scheme")
    plans = []
    for fold,(train,test) in enumerate(outer,1):
        inner = [(train[tr],train[va]) for tr,va in GroupKFold(n_splits=3).split(train,groups=groups[train])]
        te_groups = set(groups[test])
        if te_groups & set(groups[train]): raise AssertionError("Outer group leakage")
        for tr,va in inner:
            if set(groups[tr]) & set(groups[va]) or te_groups & set(groups[np.r_[tr,va]]):
                raise AssertionError("Inner calibration/outer evaluation leakage")
        plans.append(dict(fold=fold,train=train,test=test,inner=inner))
    coverage = np.concatenate([p["test"] for p in plans])
    if not np.array_equal(np.sort(coverage),indices): raise AssertionError("Outer coverage must be exactly once")
    return plans


def load_metadata():
    names = ["accession","position","native","alt","group_id","gene","dms_id","is_membrane","endpoint_type"]
    with np.load(TABLE,allow_pickle=False) as table:
        data = {name:table[name] for name in names}
        data["X"] = np.column_stack([table[name] for name in FEATURES]).astype(np.float32)
    acc = data["accession"].astype(str)
    entries = {r["accession"]:r for r in csv.DictReader(MANIFEST.open(encoding="utf-8",newline=""))}
    if any("confirmation" in entries[a]["role"] for a in set(acc)):
        raise ValueError("Reserved confirmation protein in development table")
    frozen = np.array([entries[a]["frozen_group"] for a in acc])
    data["groups"] = conservative_groups(acc,data["group_id"].astype(str),frozen)
    bank = json.loads((PREP/"annotation_bank.json").read_text(encoding="utf-8"))
    with np.load(PREP/"prepared_variant_features.npz",allow_pickle=False) as prepared:
        if tuple(prepared["feature_names"]) != RULES: raise ValueError("Feature order mismatch")
        pa,pp,pn,pv,features = (prepared[k] for k in ("accession","position","native","alt","features"))
    lookup = {(str(a),int(p),str(v)):(str(n),f) for a,p,n,v,f in zip(pa,pp,pn,pv,features)}
    if len(lookup) != len(pa): raise ValueError("Duplicate prepared feature key")
    F = np.zeros((len(acc),4),dtype=float)
    for i in np.flatnonzero(data["is_membrane"] == 1):
        key = (acc[i],int(data["position"][i]),str(data["alt"][i]))
        native,row = lookup[key]
        if native != data["native"][i]: raise ValueError("Prepared native mismatch")
        F[i] = row
    data["F"] = F
    data["bank_accessions"] = sorted(bank)
    data["role_entries"] = entries
    # Align embeddings by explicit canonical position and reject missing positions.
    E = np.empty((len(acc),1024),dtype=np.float32)
    embedding_sources = []
    for a in sorted(set(acc)):
        path = EMB/f"{a}.npz"
        with np.load(path,allow_pickle=False) as source:
            pos,values = source["position"],source["embedding"]
            seq_len = int(source["seq_len"])
        if values.shape != (len(pos),1024) or len(set(pos)) != len(pos):
            raise ValueError("Embedding shape or duplicate coordinate")
        positions = {int(p):j for j,p in enumerate(pos)}
        selected = np.flatnonzero(acc == a)
        requested = data["position"][selected].astype(int)
        if any(int(p) not in positions or not 1 <= p <= seq_len for p in requested):
            raise ValueError("Missing canonical embedding")
        if a in bank:
            sequence = bank[a]["sequence"]
            if seq_len != len(sequence) or any(sequence[int(p)-1] != n for p,n in zip(requested,data["native"][selected])):
                raise ValueError("Embedding/reference native mismatch")
        E[selected] = values[[positions[int(p)] for p in requested]].astype(np.float32)
        embedding_sources.append(dict(accession=a,path=str(path.relative_to(ROOT)),sha256=sha(path)))
    if not np.isfinite(E).all() or not np.isfinite(data["X"]).all():
        raise ValueError("Nonfinite baseline input")
    data["E"],data["embedding_sources"] = E,embedding_sources
    return data


def freeze_plan(out,data,plans):
    acc,groups = data["accession"].astype(str),data["groups"]
    members = {g:sorted(set(acc[groups == g])) for g in sorted(set(groups))}
    plan = dict(design_version=1,study="development-only standalone refinement comparison",
        rows=len(acc),proteins=len(set(acc)),membrane_proteins=len(set(acc[data["is_membrane"]==1])),
        dependence_groups=len(members),dependence_members=members,
        grouping="conservative union of existing training and frozen dependence groupings; not a new certified homology classification",
        features=FEATURES,n_pc=64,pca_solver="auto; refit on each supervised training partition only",
        baseline_backend="sklearn HistGradientBoostingRegressor",baseline_config=CONFIG,
        baseline_seed=BASE_SEED,bootstrap_seed=SEED,bootstrap_draws=9999,
        rule_features=list(RULES),ridge=.001,minimum_calibration_families=3,
        inner_calibration="3-fold group-out-of-fold residuals on outer-training membrane proteins only",
        correction_fit="nonnegative, within-assay centered, equal family/protein/assay weights; no hyperparameter tuning",
        ablations=ABLATIONS,ablation_fit="each subset independently calibrated by locking inactive coefficients at zero",
        primary_scope="13 existing membrane proteins",secondary_scopes=["all 39 proteins","26 non-membrane proteins"],
        primary_metric="mean protein Spearman; each protein averages qualifying assays (>=20 rows and >=3 target values)",
        secondary_metric="mean protein assay-averaged squared error on within-assay rank damage",
        primary_comparison="all four components versus raw standalone baseline",
        component_comparisons="PI plus each addition versus PI alone; full combination versus each leave-one-addition-out combination",
        uncertainty="paired dependence-group percentile resampling of fixed cross-validated protein metrics; conditional stability, not fresh confirmation",
        no_result_selected_configuration=True,no_reserved_outcome_files_opened=True,
        old_saved_head_reused=False,old_positive_inside_beta_reused=False,
        previous_exposure="All 39 proteins/experimental assays were already exposed during project development; this is not independent new-data confirmation",
        table_sha256=sha(TABLE),manifest_sha256=sha(MANIFEST),
        prepared_features_sha256=sha(PREP/"prepared_variant_features.npz"),
        annotation_bank_sha256=sha(PREP/"annotation_bank.json"),embedding_sources=data["embedding_sources"],
        code_sha256={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__),ROOT/"tools/membrane_rule_audit/refinement.py",ROOT/"scripts/beat_am/m2_standalone_refinement_scorer.py")},
        measured_support={r:dict(rows=int(np.count_nonzero(data["F"][:,j])),
            proteins=sorted(set(acc[data["F"][:,j] != 0])),
            dependence_groups=sorted(set(groups[data["F"][:,j] != 0])),
            sites=len(set(zip(acc[data["F"][:,j] != 0],data["position"][data["F"][:,j] != 0].tolist())))) for j,r in enumerate(RULES)},
        folds={mode:[dict(fold=p["fold"],test_accessions=sorted(set(acc[p["test"]])),
            test_groups=sorted(set(groups[p["test"]])),train_groups=sorted(set(groups[p["train"]])),
            inner=[dict(train_groups=sorted(set(groups[tr])),calibration_groups=sorted(set(groups[va]))) for tr,va in p["inner"]]) for p in values] for mode,values in plans.items()})
    path = out/"BENCHMARK_PLAN.json"
    if path.exists():
        previous = json.loads(path.read_text())
        previous.pop("frozen_utc",None)
        if previous != plan: raise ValueError("Frozen benchmark plan changed; use a new output directory")
    else:
        write_json(path,dict(frozen_utc=datetime.now(timezone.utc).isoformat(),**plan))
    return sha(path)


def model_key(train_groups,plan_hash):
    return hashlib.sha256(json.dumps(dict(groups=sorted(train_groups),plan=plan_hash),sort_keys=True).encode()).hexdigest()[:20]


def fit_or_load(data,train,out,plan_hash):
    train_groups = sorted(set(data["groups"][train]))
    key = model_key(train_groups,plan_hash)
    path = out/"model_cache"/f"{key}.joblib"
    if path.exists():
        bundle = joblib.load(path)
        if bundle["training_groups"] != train_groups or bundle["plan_sha256"] != plan_hash:
            raise ValueError("Model cache provenance mismatch")
        return bundle,key
    start = time.perf_counter()
    pca = PCA(n_components=64,random_state=BASE_SEED).fit(data["E"][train])
    X = np.hstack([data["X"][train],pca.transform(data["E"][train])]).astype(np.float32)
    model = HistGradientBoostingRegressor(**CONFIG).fit(X,data["y"][train])
    bundle = dict(model=model,pca=pca,training_groups=train_groups,training_accessions=sorted(set(data["accession"][train])),
                  plan_sha256=plan_hash,features=FEATURES,n_pc=64,fit_seconds=time.perf_counter()-start)
    joblib.dump(bundle,path)
    return bundle,key


def predict(bundle,data,indices):
    X = np.hstack([data["X"][indices],bundle["pca"].transform(data["E"][indices])]).astype(np.float32)
    predicted = bundle["model"].predict(X)
    if not np.isfinite(predicted).all(): raise ValueError("Nonfinite baseline prediction")
    return predicted


def protein_metrics(pred,data):
    result = []
    acc = data["accession"].astype(str)
    assays = data["dms_id"].astype(str)
    for a in sorted(set(acc)):
        selected = acc == a
        by_assay = []
        for assay in sorted(set(assays[selected])):
            idx = np.flatnonzero(selected & (assays == assay))
            if len(idx)<20 or len(np.unique(data["y"][idx]))<3: continue
            values = pred[idx]
            rho = np.array([spearmanr(values[:,j],data["y"][idx]).statistic for j in range(values.shape[1])])
            if not np.isfinite(rho).all(): raise ValueError("Nonfinite protein correlation")
            mse = np.mean((values-data["y"][idx,None])**2,axis=0)
            by_assay.append((rho,mse))
        if by_assay:
            first = np.flatnonzero(selected)[0]
            for j,name in enumerate(ABLATIONS):
                result.append(dict(accession=a,gene=str(data["gene"][first]),group=str(data["groups"][first]),
                    membrane=int(data["is_membrane"][first]),model=name,n_assays=len(by_assay),n_variants=int(selected.sum()),
                    spearman=float(np.mean([row[0][j] for row in by_assay])),mse=float(np.mean([row[1][j] for row in by_assay]))))
    return result


def paired_summary(records,model,reference,*,scope="membrane",draws=9999,seed=SEED):
    keep = [r for r in records if scope == "all" or r["membrane"] == (1 if scope == "membrane" else 0)]
    left = {r["accession"]:r for r in keep if r["model"] == model}
    right = {r["accession"]:r for r in keep if r["model"] == reference}
    accessions = sorted(set(left)&set(right))
    if not accessions: return None
    delta = np.array([left[a]["spearman"]-right[a]["spearman"] for a in accessions])
    mse_delta = np.array([left[a]["mse"]-right[a]["mse"] for a in accessions])
    cluster = np.array([left[a]["group"] for a in accessions])
    groups = sorted(set(cluster))
    numerator = np.array([delta[cluster == g].sum() for g in groups])
    denominator = np.array([(cluster == g).sum() for g in groups])
    rng = np.random.default_rng(seed)
    counts = rng.multinomial(len(groups),np.full(len(groups),1/len(groups)),size=draws)
    samples = counts@numerator/(counts@denominator)
    lo,hi = np.quantile(samples,[.025,.975])
    tolerance = 1e-12
    return dict(model=model,reference=reference,scope=scope,proteins=len(accessions),groups=len(groups),
        model_mean_spearman=float(np.mean([left[a]["spearman"] for a in accessions])),
        reference_mean_spearman=float(np.mean([right[a]["spearman"] for a in accessions])),
        delta_spearman=float(delta.mean()),ci95_delta_spearman=[float(lo),float(hi)],
        equal_group_mean_delta=float(np.mean(numerator/denominator)),
        wins=int((delta>tolerance).sum()),losses=int((delta < -tolerance).sum()),ties=int((np.abs(delta)<=tolerance).sum()),
        mean_mse_delta=float(mse_delta.mean()),
        leave_one_group_out_delta=[dict(omitted=g,delta=float(delta[cluster != g].mean())) for g in groups] if len(groups)>1 else [],
        interpretation="conditional_cross_validation_stability_not_independent_confirmation")


def evaluate_scheme(mode,data,plans,out,plan_hash):
    target = out/mode
    target.mkdir(exist_ok=True)
    all_predictions = np.full((len(data["y"]),len(ABLATIONS)),np.nan)
    fold_records = []
    for plan in plans:
        fold = plan["fold"]
        cache = target/f"fold_{fold:02d}.npz"
        detail_path = target/f"fold_{fold:02d}.json"
        train,test = plan["train"],plan["test"]
        if cache.exists() and detail_path.exists():
            detail = json.loads(detail_path.read_text())
            if detail["plan_sha256"] != plan_hash: raise ValueError("Fold cache plan mismatch")
            with np.load(cache) as z:
                if not np.array_equal(z["test_indices"],test): raise ValueError("Fold cache index mismatch")
                all_predictions[test] = z["predictions"]
            fold_records.append(detail)
            continue
        log(f"{mode}: outer {fold}/{len(plans)}; held out {', '.join(sorted(set(data['gene'][test])))}")
        calibration = np.full(len(data["y"]),np.nan)
        inner_provenance = []
        for inner,(tr,va) in enumerate(plan["inner"],1):
            bundle,key = fit_or_load(data,tr,out,plan_hash)
            if set(bundle["training_groups"]) & set(data["groups"][np.r_[va,test]]):
                raise AssertionError("Inner model or PCA used held-out group")
            calibration[va] = predict(bundle,data,va)
            inner_provenance.append(dict(inner=inner,model_cache_key=key,training_groups=bundle["training_groups"],
                calibration_groups=sorted(set(data["groups"][va]))))
        cal = train[data["is_membrane"][train] == 1]
        if not np.isfinite(calibration[cal]).all(): raise ValueError("Missing calibration OOF prediction")
        bundle,key = fit_or_load(data,train,out,plan_hash)
        if set(bundle["training_groups"]) & set(data["groups"][test]): raise AssertionError("Outer model leakage")
        baseline = predict(bundle,data,test)
        predictions = np.empty((len(test),len(ABLATIONS)))
        fitted = {}
        for j,name in enumerate(ABLATIONS):
            if name == "raw":
                predictions[:,j] = baseline
                continue
            active = set(name.split("+"))
            fit = fit_coefficients(data["F"][cal],calibration[cal],data["y"][cal],
                accession=data["accession"][cal],assay=data["dms_id"][cal],family=data["groups"][cal],
                role=np.repeat("calibration",len(cal)),target_scale="within_assay_rank_damage",
                label_provenance="existing exposed 44-assay development table; inner group-OOF residuals only",
                ridge=.001,minimum_families=3,locked_zero_rules=set(RULES)-active)
            predictions[:,j] = apply_refinement(baseline,data["F"][test],fit["coefficients"])["scores"]
            fitted[name] = fit
        all_predictions[test] = predictions
        detail = dict(fold=fold,plan_sha256=plan_hash,outer_model_cache_key=key,
            training_groups=bundle["training_groups"],evaluation_groups=sorted(set(data["groups"][test])),
            evaluation_accessions=sorted(set(data["accession"][test])),inner=inner_provenance,
            correction_calibration_groups=sorted(set(data["groups"][cal])),correction_fits=fitted,
            calibrated_on_outer_evaluation=False,experimental_confirmation_opened=False)
        write_json(detail_path,detail)
        np.savez_compressed(cache,test_indices=test,predictions=predictions,models=np.array(ABLATIONS))
        fold_records.append(detail)
    if not np.isfinite(all_predictions).all(): raise ValueError("Incomplete held-out predictions")
    # Non-membrane proteins receive zero correction as an exact scope control.
    non = data["is_membrane"] == 0
    if not np.array_equal(all_predictions[non],np.repeat(all_predictions[non,:1],len(ABLATIONS),axis=1)):
        raise AssertionError("Correction escaped its declared membrane annotation scope")
    np.savez_compressed(target/"out_of_fold_predictions.npz",predictions=all_predictions,
        models=np.array(ABLATIONS),accession=data["accession"],position=data["position"],native=data["native"],alt=data["alt"],
        groups=data["groups"],target=data["y"])
    records = protein_metrics(all_predictions,data)
    with (target/"per_protein.csv").open("w",encoding="utf-8",newline="") as stream:
        writer = csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    comparisons = {scope:[paired_summary(records,name,"raw",scope=scope) for name in ABLATIONS] for scope in ("membrane","all","non_membrane")}
    additions,removals = {},{}
    for rule in RULES[1:]:
        additions[rule] = paired_summary(records,"positive_inside+"+rule,"positive_inside")
        reduced = "+".join(r for r in RULES if r != rule)
        removals[rule] = paired_summary(records,FULL,reduced)
    result = dict(scheme=mode,plan_sha256=plan_hash,outer_fits=len(plans),
        baseline_refitted_per_outer_partition=True,preprocessing_refitted_per_inner_partition=True,
        outer_evaluation_groups_excluded_from_all_training_and_calibration=True,
        comparisons_vs_raw=comparisons,add_one_rule_to_PI=additions,remove_one_rule_from_full=removals,
        per_protein_file=str((target/"per_protein.csv").relative_to(ROOT)),
        old_model_artifacts_modified=False,new_experimental_confirmation_opened=False,
        scope="previously exposed development data only",non_membrane_scores_exactly_unchanged=True,
        fitted_coefficients=[dict(fold=f["fold"],evaluation_accessions=f["evaluation_accessions"],
            coefficients={name:fit["coefficients"] for name,fit in f["correction_fits"].items()}) for f in fold_records])
    write_json(target/"results.json",result)
    full = comparisons["membrane"][-1]
    log(f"{mode} finished: full vs raw mean protein Spearman change {full['delta_spearman']:+.6f}; interval {full['ci95_delta_spearman']}")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=PREP/"benchmark_20261025")
    parser.add_argument("--stage",choices=["all","grouped_cv","guarded_loocv"],default="all")
    parser.add_argument("--threads",type=int,default=8)
    args = parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True);(args.out/"model_cache").mkdir(exist_ok=True)
    log("Loading existing development metadata and embeddings; experimental targets remain unopened until the plan is frozen")
    data = load_metadata()
    plans = {mode:fold_plans(data["groups"],mode) for mode in ("grouped_cv","guarded_loocv")}
    plan_hash = freeze_plan(args.out,data,plans)
    with np.load(TABLE,allow_pickle=False) as table: data["y"] = table["target"].astype(float)
    if not np.isfinite(data["y"]).all(): raise ValueError("Nonfinite experimental target")
    log(f"Plan frozen: {len(data['y']):,} variants, 39 proteins, {len(set(data['groups']))} conservative dependence groups")
    with threadpool_limits(limits=args.threads):
        for mode in ("grouped_cv","guarded_loocv"):
            if args.stage in ("all",mode): evaluate_scheme(mode,data,plans[mode],args.out,plan_hash)
    write_json(args.out/"RUN_STATUS.json",dict(completed_utc=datetime.now(timezone.utc).isoformat(),
        plan_sha256=plan_hash,completed_schemes=[mode for mode in plans if (args.out/mode/"results.json").exists()],
        new_reserved_outcomes_opened=False,old_model_artifacts_modified=False))


if __name__ == "__main__":
    main()

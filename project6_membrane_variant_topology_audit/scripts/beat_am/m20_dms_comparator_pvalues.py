"""Compute paired group-level sign-flip tests for the membrane cross-validation
comparisons.
"""
from __future__ import annotations

from datetime import datetime, timezone
import itertools
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'scripts/beat_am'),str(ROOT/'tools')]
from m4_benchmark_multirule import ABLATIONS,sha,TABLE
from membrane_rule_audit.refinement import write_json

OUT=ROOT/'results/dms_comparator_pvalues_20261028'
BENCH=ROOT/'results/standalone_multirule_20261025/benchmark_20261025'
EXISTING=ROOT/'results/external_predictor_comparison_20261026'
MAIN=('GEMME','VESPA','VespaG')
SUPPLEMENT=('VESPAl',)
COMPARATORS=MAIN+SUPPLEMENT
CANDIDATE='positive_inside+disulfide_loss'
PAPER_CANDIDATE='positive_inside+glycosylation_loss+disulfide_loss'
FULL='positive_inside+charge_burial+glycosylation_loss+disulfide_loss'
SEED,DRAWS=20261025,9999
SCHEMES=('guarded_loocv','grouped_cv')


def read_json(path):return json.loads(path.read_text(encoding='utf-8'))


def holm(pvalues):
    """Holm step-down adjusted p in original hypothesis order."""
    p=np.asarray(pvalues,dtype=float)
    if p.ndim!=1 or not np.isfinite(p).all() or np.any((p<0)|(p>1)):raise ValueError('Invalid p-values')
    order=np.argsort(p,kind='stable');out=np.empty(len(p));previous=0.
    for rank,index in enumerate(order):
        previous=max(previous,min(1.,(len(p)-rank)*p[index]));out[index]=previous
    return out


def paired_test(delta,groups):
    """Equal-protein mean statistic, with signs shared within each group."""
    d=np.asarray(delta,dtype=float);g=np.asarray(groups).astype(str)
    if d.ndim!=1 or g.shape!=d.shape or not len(d) or not np.isfinite(d).all():raise ValueError('Invalid paired data')
    units=sorted(set(g))
    if len(units)>20:raise ValueError('Exact enumeration limited to 20 groups')
    numerator=np.array([d[g==unit].sum() for unit in units])
    denominator=np.array([(g==unit).sum() for unit in units])
    signs=np.array(list(itertools.product([-1.,1.],repeat=len(units))))
    null=signs@numerator/len(d);observed=float(d.mean())
    p=float(np.mean(np.abs(null)>=abs(observed)-1e-14))
    counts=np.random.default_rng(SEED).multinomial(len(units),np.full(len(units),1/len(units)),size=DRAWS)
    sample=counts@numerator/(counts@denominator)
    interval=np.quantile(sample,[.025,.975]).astype(float).tolist()
    return dict(delta_spearman=observed,ci95_delta_spearman=interval,p_two_sided=p,
        proteins=len(d),dependence_groups=len(units),sign_patterns=len(signs),
        nonzero_group_sums=int((np.abs(numerator)>1e-12).sum()),
        wins=int((d>1e-12).sum()),losses=int((d< -1e-12).sum()),ties=int((np.abs(d)<=1e-12).sum()))


def source_paths():
    paths=[TABLE,BENCH/'BENCHMARK_PLAN.json',EXISTING/'COMPARISON_PLAN.json',
        EXISTING/'aligned_external_scores.npz',Path(__file__)]
    for mode in SCHEMES:
        paths.extend([BENCH/mode/'out_of_fold_predictions.npz',
            EXISTING/mode/'requested_common_proteins.csv',EXISTING/mode/'requested_common_assays.csv',
            EXISTING/mode/'paired_to_VespaG_proteins.csv',EXISTING/mode/'paired_to_VespaG_assays.csv'])
    return paths


def freeze_plan():
    OUT.mkdir(parents=True,exist_ok=True)
    plan=dict(study='paired p-values on the original membrane cross-validated scores',
        status='exploratory, conditional on existing predictions and previously selected configurations',
        original_scores_reused=True,new_model_fit=False,new_calibration=False,split_shuffling=False,
        primary_presentation_scheme='guarded_loocv',sensitivity_scheme='grouped_cv',
        primary_panel='13 membrane proteins, 15 assays, 43991 measurements, 12 dependence groups',
        displayed_candidate=CANDIDATE,other_previously_reported_paper_candidate=PAPER_CANDIDATE,
        standalone_configurations=ABLATIONS,main_comparators=list(MAIN),supplementary_comparators=list(SUPPLEMENT),
        excluded_comparators=['AlphaMissense','ESM1b'],
        statistic='equal protein mean of model-minus-comparator within-protein assay-mean Spearman differences',
        test='two-sided exhaustive block sign-flip, all 4096 patterns for 12 dependence groups',
        null_assumptions='independent/exchangeable group differences with sign symmetry; shared CV training induces dependence not modeled by this fixed-score calculation',
        multiplicity=['Holm across three main comparators for each fixed configuration and scheme',
            'Holm across all 27 main-comparator/configuration contrasts per scheme',
            'Holm across all 36 comparator/configuration contrasts per scheme, including VESPA-light'],
        uncertainty='paired dependence-group percentile bootstrap; fixed scores, 9999 draws, seed 20261025',
        interpretation='not an independent confirmation and not uncertainty over retraining or model/configuration selection',
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in source_paths()})
    path=OUT/'ANALYSIS_PLAN.json'
    if path.exists():
        old=read_json(path);old.pop('created_utc')
        if old!=plan:raise ValueError('Analysis specification changed; preserve the original')
    else:write_json(path,dict(created_utc=datetime.now(timezone.utc).isoformat(),**plan))
    return plan


def rebuild_metrics(mode):
    """Recalculate scores from original aligned rows, reading requested models only."""
    with np.load(TABLE,allow_pickle=False) as z:
        frame=pd.DataFrame({name:z[name] for name in ['accession','gene','dms_id','position','native','alt','is_membrane','target']})
    with np.load(BENCH/mode/'out_of_fold_predictions.npz',allow_pickle=False) as z:
        own=z['predictions'];names=z['models'].astype(str).tolist();frame['group']=z['groups'].astype(str)
        assert names==ABLATIONS and np.array_equal(z['target'],frame.target)
        for key in ['accession','position','native','alt']:assert np.array_equal(z[key],frame[key])
    with np.load(EXISTING/'aligned_external_scores.npz',allow_pickle=False) as z:
        for key in ['accession','position','native','alt','dms_id']:assert np.array_equal(z[key],frame[key])
        external_names=z['models'].astype(str).tolist()
        external=z['scores'][:,[external_names.index(name) for name in COMPARATORS]]
    values=np.column_stack([own,external]);names=names+list(COMPARATORS)
    mask=frame.is_membrane.eq(1).to_numpy() & np.isfinite(values).all(axis=1)
    assert mask.sum()==43991 and frame.loc[mask,'accession'].nunique()==13
    assert frame.loc[mask,'dms_id'].nunique()==15 and frame.loc[mask,'group'].nunique()==12
    assay_records=[];max_rank_difference=0.
    for assay,part in frame[mask].groupby('dms_id',sort=True):
        assert len(part)>=20 and part.target.nunique()>=3
        for j,name in enumerate(names):
            x=values[part.index,j];y=part.target.to_numpy()
            rho=float(spearmanr(x,y).statistic)
            independent=float(np.corrcoef(rankdata(x),rankdata(y))[0,1])
            max_rank_difference=max(max_rank_difference,abs(rho-independent))
            assert np.isfinite(rho) and abs(rho-independent)<1e-12
            assay_records.append(dict(scheme=mode,assay=assay,accession=part.accession.iloc[0],
                gene=part.gene.iloc[0],group=part.group.iloc[0],model=name,rows=len(part),spearman=rho))
    adf=pd.DataFrame(assay_records);proteins=[]
    for (acc,name),part in adf.groupby(['accession','model'],sort=True):
        proteins.append(dict(scheme=mode,accession=acc,gene=part.gene.iloc[0],group=part.group.iloc[0],
            model=name,assays=len(part),rows=int(part.rows.sum()),spearman=float(part.spearman.mean())))
    pdf=pd.DataFrame(proteins)
    # Reproduce the archived per-protein metrics as well as their mean.
    common=pd.read_csv(EXISTING/mode/'requested_common_proteins.csv')
    vg=pd.read_csv(EXISTING/mode/'paired_to_VespaG_proteins.csv')
    for row in proteins:
        source=vg if row['model']=='VespaG' else common
        match=source[source.scope.eq('membrane') & source.accession.eq(row['accession']) & source.model.eq(row['model'])]
        assert len(match)==1 and int(match.iloc[0]['rows'])==row['rows']
        assert abs(float(match.iloc[0]['spearman'])-row['spearman'])<1e-12
    return pdf,adf,dict(reconstructed_proteins=13,reconstructed_assays=15,reconstructed_rows=43991,
        reconstructed_dependence_groups=12,independent_rank_correlation_max_error=max_rank_difference,
        all_requested_metrics_reproduce_existing_results=True,split_shuffle=False,new_fit=False)


def contrast(proteins,model,reference):
    left=proteins[proteins.model.eq(model)].set_index('accession').sort_index()
    right=proteins[proteins.model.eq(reference)].set_index('accession').sort_index()
    if not left.index.equals(right.index):raise ValueError('Unmatched protein sets')
    if not np.array_equal(left.group,right.group) or not np.array_equal(left.rows,right.rows):raise ValueError('Unmatched groups/rows')
    delta=(left.spearman-right.spearman).to_numpy()
    return dict(scheme=left.scheme.iloc[0],model=model,reference=reference,
        model_mean_spearman=float(left.spearman.mean()),reference_mean_spearman=float(right.spearman.mean()),
        **paired_test(delta,left.group.to_numpy()))


def analyze(mode,proteins):
    rows=[contrast(proteins,name,reference) for name in ABLATIONS for reference in COMPARATORS]
    for name in ABLATIONS:
        selected=[r for r in rows if r['model']==name and r['reference'] in MAIN]
        for row,adjusted in zip(selected,holm([r['p_two_sided'] for r in selected])):
            row['p_holm_three_main_for_configuration']=float(adjusted)
    main=[r for r in rows if r['reference'] in MAIN]
    for row,adjusted in zip(main,holm([r['p_two_sided'] for r in main])):
        row['p_holm_27_main_across_configurations']=float(adjusted)
    for row,adjusted in zip(rows,holm([r['p_two_sided'] for r in rows])):
        row['p_holm_36_all_across_configurations']=float(adjusted)
    own=[contrast(proteins,name,'raw') for name in ABLATIONS if name!='raw']
    for row,adjusted in zip(own,holm([r['p_two_sided'] for r in own])):row['p_holm_eight_refinements']=float(adjusted)
    # Check exact null tails with an independent subset-sum enumeration.
    for row in rows+own:
        left=proteins[proteins.model.eq(row['model'])].set_index('accession').sort_index()
        right=proteins[proteins.model.eq(row['reference'])].set_index('accession').sort_index()
        delta=left.spearman-right.spearman
        blocks=np.array([delta[left.group.eq(g)].sum() for g in sorted(set(left.group))])
        absolute=abs(float(delta.mean()));hits=0
        for bits in range(1<<len(blocks)):
            flipped=sum(v for j,v in enumerate(blocks) if bits&(1<<j))
            statistic=(float(blocks.sum())-2*flipped)/len(delta)
            hits+=abs(statistic)>=absolute-1e-14
        if hits/(1<<len(blocks))!=row['p_two_sided']:raise AssertionError('Independent exact p mismatch')
    summaries=[dict(scheme=mode,model=name,mean_spearman=float(part.spearman.mean()),
        proteins=len(part),groups=part.group.nunique(),assays=int(part.assays.sum()),rows=int(part.rows.sum()))
        for name,part in proteins.groupby('model',sort=True)]
    return dict(scheme=mode,summaries=summaries,comparisons=rows,refinements_vs_raw=own,
        independent_exact_p_check_passed=True)


def render(result):
    primary=result['guarded_loocv'];by_name={r['model']:r for r in primary['summaries']}
    lines=['# Original DMS cross-validation: VESPA, VespaG and GEMME p-values','',
        'Computed 28 October 2026 (America/Los_Angeles). This report excludes AlphaMissense and ESM1b as requested. It reuses the original cross-validated scores without fitting a model, changing a coefficient or constructing a new split.','',
        '**Evaluation naming:** the primary scores are guarded leave-group-out predictions. The original five-fold grouped predictions are a sensitivity analysis. Both procedures used unshuffled group splits; calling these predictions "shuffled scoring" would be inaccurate. Sign flips below are statistical null calculations, not new train/test splits.','',
        'The primary panel contains 13 membrane proteins, 15 assays and 43,991 measurements across 12 dependence groups. Correlations are averaged equally across assays within each protein, then equally across proteins. All displayed models use the same rows.','',
        '## Score table','','| Model | Mean Spearman |','|---|---:|']
    for name in ['raw',CANDIDATE,PAPER_CANDIDATE,FULL,*MAIN,*SUPPLEMENT]:
        lines.append(f"| {name} | {by_name[name]['mean_spearman']:.6f} |")
    for model,title in [(CANDIDATE,'Previously highlighted PI plus disulfide configuration'),('raw','Unchanged standalone predictor'),(PAPER_CANDIDATE,'PI plus glycosylation plus disulfide configuration'),(FULL,'All four corrections')]:
        lines.extend(['',f'## {title}','','| Comparator | Difference in mean Spearman | Conditional 95% interval | Two-sided exact p | Holm p: three main comparisons | Holm p: all 27 main/configuration comparisons |','|---|---:|---|---:|---:|---:|'])
        for ref in MAIN:
            row=next(r for r in primary['comparisons'] if r['model']==model and r['reference']==ref)
            lo,hi=row['ci95_delta_spearman']
            lines.append(f"| {ref} | {row['delta_spearman']:+.6f} | [{lo:+.6f}, {hi:+.6f}] | {row['p_two_sided']:.6f} | {row['p_holm_three_main_for_configuration']:.6f} | {row['p_holm_27_main_across_configurations']:.6f} |")
        light=next(r for r in primary['comparisons'] if r['model']==model and r['reference']=='VESPAl')
        lines.extend(['',f"Supplementary VESPA-light contrast: difference {light['delta_spearman']:+.6f}, exact p = {light['p_two_sided']:.6f}, Holm p across all 36 configuration/comparator contrasts = {light['p_holm_36_all_across_configurations']:.6f}."])
    lines.extend(['','## Does refinement improve the unchanged standalone?','','| Configuration | Difference | Exact p | Holm p: eight refinements |','|---|---:|---:|---:|'])
    for row in primary['refinements_vs_raw']:
        lines.append(f"| {row['model']} | {row['delta_spearman']:+.6f} | {row['p_two_sided']:.6f} | {row['p_holm_eight_refinements']:.6f} |")
    lines.extend(['','## Existing five-fold sensitivity','','| Comparator | PI+D mean | Comparator mean | Difference | Exact p | Holm p: three comparisons |','|---|---:|---:|---:|---:|---:|'])
    for ref in MAIN:
        row=next(r for r in result['grouped_cv']['comparisons'] if r['model']==CANDIDATE and r['reference']==ref)
        lines.append(f"| {ref} | {row['model_mean_spearman']:.6f} | {row['reference_mean_spearman']:.6f} | {row['delta_spearman']:+.6f} | {row['p_two_sided']:.6f} | {row['p_holm_three_main_for_configuration']:.6f} |")
    lines.extend(['','## Statistical interpretation','',
        'The two-sided exact block sign-flip test enumerates all 4,096 sign assignments across 12 dependence groups. OPSD and CCR5 share a group and always receive the same sign. The statistic preserves the equal-protein mean; groups are not incorrectly reweighted as one equal observation when they contain different numbers of proteins. The analysis does not treat 43,991 mutations as independent replicates.','',
        'The null calculation assumes sign-symmetric/exchangeable group differences. Because cross-validated heads share training data, this is a conditional exploratory comparison of existing predictions, not an independent repeated-training significance result. The percentile-bootstrap confidence interval and exact test are different statistical constructions and need not agree at a 0.05 threshold. Report the actual test result rather than choosing the method that is favorable.','',
        'Holm adjustment across three comparisons applies to a fixed displayed configuration. It does not correct selecting that configuration after viewing earlier results. More conservative columns adjust all 27 main-comparator/configuration contrasts, and the CSV/JSON also include adjustment over 36 contrasts when VESPA-light is included. Both schemes are saved in full. No claim of model-selection-adjusted or fresh confirmatory significance follows from a nominal p-value.','',
        'The external methods have different training objectives and supervision histories. The chosen DMS comparators are relevant to this outcome; omitting other predictors does not establish matched pretraining or supervision across the remaining methods.','',
        'Fixed hold-out experiments were archived at the user\'s request with byte-for-byte verification. Their historical outcome exposure remains in archive/fixed_holdout_20261028/. No archived hold-out scores were used to fit or select anything in this calculation.','',
        '## Reproduce','','```powershell','.\\.venv_ml\\Scripts\\python.exe scripts\\beat_am\\m20_dms_comparator_pvalues.py','```','',
        'Full per-assay/per-protein metrics, every configuration comparison, analysis specification, source hashes and independent exact-p verification are saved in results/dms_comparator_pvalues_20261028/.'])
    text='\n'.join(lines)+'\n';(OUT/'REPORT.md').write_text(text,encoding='utf-8')


def main():
    plan=freeze_plan();result={};verification={};all_proteins=[];all_assays=[]
    for mode in SCHEMES:
        proteins,assays,check=rebuild_metrics(mode)
        result[mode]=analyze(mode,proteins);verification[mode]=check
        all_proteins.extend(proteins.to_dict('records'));all_assays.extend(assays.to_dict('records'))
    result['plan_sha256']=sha(OUT/'ANALYSIS_PLAN.json')
    write_json(OUT/'results.json',result)
    pd.DataFrame(all_proteins).to_csv(OUT/'per_protein.csv',index=False)
    pd.DataFrame(all_assays).to_csv(OUT/'per_assay.csv',index=False)
    pd.DataFrame([r for mode in SCHEMES for r in result[mode]['comparisons']]).to_csv(OUT/'all_comparisons.csv',index=False)
    pd.DataFrame([r for mode in SCHEMES for r in result[mode]['refinements_vs_raw']]).to_csv(OUT/'refinements_vs_raw.csv',index=False)
    verification['all_original_source_hashes_unchanged']=all(sha(ROOT/p)==h for p,h in plan['source_sha256'].items())
    assert verification['all_original_source_hashes_unchanged']
    verification['excluded_models_absent_from_outputs']=not any(r['model'] in ['AlphaMissense','ESM1b'] for r in all_proteins)
    assert verification['excluded_models_absent_from_outputs']
    write_json(OUT/'VERIFICATION.json',verification);render(result)
    for row in result['guarded_loocv']['comparisons']:
        if row['model']==CANDIDATE:print(json.dumps(row),flush=True)


if __name__=='__main__':main()

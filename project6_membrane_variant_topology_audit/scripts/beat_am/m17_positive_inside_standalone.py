"""Audit positive-inside contrasts in frozen standalone outputs without refitting or
reading outcomes.
"""
from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timezone
import csv
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
from membrane_rule_audit import panel as estimator
from membrane_rule_audit.refinement import ALPHABET, RULES, topology_map

BASE=ROOT/'results/standalone_multirule_20261025'
BENCH=BASE/'benchmark_20261025'
CACHE=BASE/'label_free_audit_20261025'
OUT=ROOT/'results/positive_inside_standalone_20261028'
WINDOWS=[5,10,20,60]
COMPARATORS={'all_uncharged':'ACFGILMNPQSTVWY','polar_neutral':'NQST',
             'hydrophobic':'AFILMV','small':'AGS','aromatic':'FWY'}
DRAWS=1999
SEED=20261028
SOURCES={}


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def read(path):
    path=Path(path)
    SOURCES[str(path.relative_to(ROOT))]=sha(path)
    return json.loads(path.read_text(encoding='utf-8'))


def write(name,value):
    (OUT/name).write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')


def independent_parts(protein,native,side,matrix,comparator):
    """Plain Python implementation, separate from packaged vectorized strata."""
    rows=defaultdict(lambda:[[],[]]);per_site=[]
    eligible=[]
    for i,(a,n,s,scores) in enumerate(zip(protein,native,side,matrix)):
        if n not in estimator.NEUTRAL:continue
        charged=[scores[ALPHABET.index(k)] for k in 'KR']
        reference=[scores[ALPHABET.index(k)] for k in comparator if k!=n]
        if not reference or not np.isfinite(charged).all() or not np.isfinite(reference).all():continue
        z=float(np.mean(charged)-np.mean(reference))
        rows[(str(a),str(n))][int(s)].append(z)
        per_site.append(z);eligible.append(i)
    per_protein=defaultdict(lambda:[0.,0.]);strata=[]
    for (a,n),(cyto,outer) in sorted(rows.items()):
        if not cyto or not outer:continue
        nc,no=len(cyto),len(outer)
        weight=nc*no/(nc+no)
        delta=float(np.mean(cyto)-np.mean(outer))
        per_protein[a][0]+=weight*delta
        per_protein[a][1]+=weight
        strata.append(dict(accession=a,native=n,cytosolic_sites=nc,outer_sites=no,
                           weight=weight,delta=delta,numerator=weight*delta))
    return per_protein,strata,np.asarray(per_site),np.asarray(eligible,dtype=int)


def resample_parts(parts):
    keys=sorted(parts)
    numerator=np.array([parts[k][0] for k in keys],dtype=float)
    denominator=np.array([parts[k][1] for k in keys],dtype=float)
    if not keys:return None,None
    point=float(numerator.sum()/denominator.sum())
    if len(keys)<2:return point,None
    rng=np.random.default_rng(SEED)
    counts=rng.multinomial(len(keys),np.full(len(keys),1/len(keys)),size=DRAWS)
    values=(counts@numerator)/(counts@denominator)
    return point,np.percentile(values,[2.5,97.5]).tolist()


def group_parts(per_protein,group_for):
    grouped=defaultdict(lambda:[0.,0.])
    for acc,(num,den) in per_protein.items():
        grouped[group_for[acc]][0]+=num;grouped[group_for[acc]][1]+=den
    return dict(grouped)


def verdict(interval):
    if interval is None:return 'not_evaluable'
    if interval[1]<0:return 'expected_direction_supported'
    if interval[0]>0:return 'opposite_direction_supported'
    return 'inconclusive'


def known_answers():
    acc=np.repeat(['A','B','C'],4)
    native=np.tile(['A','A','Q','Q'],3)
    side=np.tile([0,1,0,1],3)
    base=np.arange(12,dtype=float)*.17
    matrix=np.repeat(base[:,None],20,axis=1)
    kr=[ALPHABET.index(k) for k in 'KR']
    matrix[:,kr]+=np.where(side==0,-.02,.02)[:,None]
    matrix[np.arange(12),[ALPHABET.index(n) for n in native]]=np.nan
    def delta(m,s=side):
        pp,_,_,_=independent_parts(acc,native,s,m,estimator.NEUTRAL)
        return resample_parts(pp)[0]
    point=delta(matrix)
    checks=dict(known_negative_direction=np.isclose(point,-.04,atol=1e-12),
                side_reversal=np.isclose(delta(matrix,1-side),-point,atol=1e-12),
                same_site_offset_cancels=np.isclose(delta(matrix+np.arange(12)[:,None]),point,atol=1e-12))
    null=matrix.copy();null[:,kr]-=np.where(side==0,-.02,.02)[:,None]
    checks['null_reads_zero']=np.isclose(delta(null),0,atol=1e-12)
    panel=dict(protein=np.repeat(np.arange(3),4),ref_aa=np.array([ALPHABET.index(n) for n in native]),
               side=side,scores=matrix)
    contrast=estimator.gain_contrast(panel)
    original=estimator.stratum_parts(panel,contrast['penalty'],contrast['eligible'])
    checks['known_answer_matches_packaged']=np.isclose(original['numerator'].sum()/original['denominator'].sum(),point,atol=1e-12)
    if not all(checks.values()):raise AssertionError(checks)
    return {k:bool(v) for k,v in checks.items()}


def freeze():
    OUT.mkdir(parents=True,exist_ok=True)
    bank=read(BASE/'annotation_bank.json')
    plan=read(BENCH/'BENCHMARK_PLAN.json')
    oldplan=read(CACHE/'AUDIT_PLAN.json')
    cached_provenance={}
    for name,expected in oldplan['source_sha256'].items():
        actual=sha(ROOT/name)
        cached_provenance[name]=(actual==expected)
        if actual!=expected:raise AssertionError(f'Cached generation source changed: {name}')
    for p in [CACHE/'standalone_position_scores.npz',Path(estimator.__file__),
              ROOT/'tools/membrane_rule_audit/refinement.py',Path(__file__)]:
        SOURCES[str(p.relative_to(ROOT))]=sha(p)
    folds={}
    for fold in plan['folds']['guarded_loocv']:
        if not set(fold['test_accessions']) & set(bank):continue
        p=BENCH/'guarded_loocv'/f"fold_{fold['fold']:02d}.json"
        folds[fold['fold']]=read(p)
    runplan=dict(created_utc=datetime.now(timezone.utc).isoformat(),
        question='Do standalone predictions express the positive-inside score contrast?',
        scope='Existing exposed 13-membrane-protein panel; score behavior, not prospective accuracy confirmation',
        primary=dict(window=10,comparator=COMPARATORS['all_uncharged'],
                     native_eligibility=estimator.NEUTRAL,charged_alternatives='KR',
                     sign='cytosolic minus non-cytosolic within-position K/R penalty; expected negative'),
        matching='protein x native residue; harmonic arm-size weights',
        primary_uncertainty='percentile dependence-group resampling, fixed predictions',
        secondary_uncertainty='whole-protein resampling',draws=DRAWS,seed=SEED,
        windows=WINDOWS,comparators=COMPARATORS,
        outputs=['saved_full_fit_raw','guarded_oof_raw',*['guarded_oof__'+m for m in plan['ablations'] if m!='raw']],
        coefficient_policy='Existing held-out fold coefficients only; no calibration, target margin fitting or selection',
        protocol_status='Recorded before this run, on historically exposed data; not new preregistered confirmation',
        source_sha256=SOURCES.copy(),cached_generation_sources_match=cached_provenance,
        no_DMS_target_table_read=True,no_external_predictor_scores_read=True,no_model_fit=True,no_network=True)
    file=OUT/'AUDIT_PLAN.json'
    if file.exists():
        old=json.loads(file.read_text(encoding='utf-8'))
        old.pop('created_utc',None);new=runplan.copy();new.pop('created_utc',None)
        if old!=new:raise ValueError('Frozen positive-inside plan changed; use a new output folder')
    else:write('AUDIT_PLAN.json',runplan)
    return bank,plan,folds,runplan


def main():
    bank,bplan,folds,runplan=freeze()
    print('Positive-inside design frozen; loading the existing score archive',flush=True)
    checks=known_answers()
    group_for={a:g for g,members in bplan['dependence_members'].items() for a in members}
    with np.load(CACHE/'standalone_position_scores.npz',allow_pickle=False) as z:
        acc,pos,native=z['accession'],z['position'],z['native']
        if ''.join(z['alphabet'])!=ALPHABET:raise AssertionError('Alphabet mismatch')
        full,raw,F=z['saved_full_fit'],z['guarded_oof'],z['features']
    assert len(pos)==10185 and len(set(acc))==13
    native_index=np.array([ALPHABET.index(n) for n in native])
    for matrix in [full,raw]:
        if not np.isnan(matrix[np.arange(len(pos)),native_index]).all():raise AssertionError('Native cells not masked')
        if not np.all(np.isfinite(matrix).sum(axis=1)==19):raise AssertionError('Missing missense scores')
    contexts={a:topology_map(r['topology']) for a,r in bank.items()}
    geo=[contexts[str(a)].get(int(p),dict(side=-1,distance=-1)) for a,p in zip(acc,pos)]
    source_side=np.array([g['side'] for g in geo]);distance=np.array([g['distance'] for g in geo])
    sides=np.where(source_side==1,0,np.where(source_side==2,1,-1))
    for a,p,n in zip(acc,pos,native):
        if bank[str(a)]['sequence'][int(p)-1]!=n:raise AssertionError('Canonical sequence mismatch')
    matrices={'saved_full_fit_raw':full,'guarded_oof_raw':raw}
    for name in bplan['ablations']:
        if name=='raw':continue
        fits={}
        for metadata in folds.values():
            if set(metadata['evaluation_groups']) & set(metadata['correction_calibration_groups']):
                raise AssertionError('Calibration includes scored evaluation group')
            for a in metadata['evaluation_accessions']:
                if a in bank:fits[a]=metadata['correction_fits'][name]['coefficients']
        beta=np.array([[fits[str(a)][r] for r in RULES] for a in acc])
        if np.any(beta<0):raise AssertionError('Negative coefficient')
        matrices['guarded_oof__'+name]=raw+np.sum(F*beta[:,None,:],axis=2)
    results=[];contributions=[];primary_parts={}
    accessions=sorted(bank);protein_ids={a:i for i,a in enumerate(accessions)}
    for model,matrix in matrices.items():
        for window in WINDOWS:
            keep=(sides>=0)&(distance>=1)&(distance<=window)
            ap,nat,s,mat=acc[keep],native[keep],sides[keep],matrix[keep]
            for comparator_name,comparator in COMPARATORS.items():
                pp,strata,z,eligible=independent_parts(ap,nat,s,mat,comparator)
                estimate,protein_ci=resample_parts(pp)
                gp=group_parts(pp,group_for)
                group_estimate,group_ci=resample_parts(gp)
                if estimate is not None and not np.isclose(estimate,group_estimate,atol=1e-12):raise AssertionError('Grouping changed point estimate')
                sd=float(np.std(z,ddof=1)) if len(z)>1 else None
                packaged=None
                if comparator_name=='all_uncharged':
                    panel=dict(protein=np.array([protein_ids[str(a)] for a in ap]),
                               ref_aa=np.array([ALPHABET.index(n) for n in nat]),side=s,scores=mat)
                    c=estimator.gain_contrast(panel)
                    parts=estimator.stratum_parts(panel,c['penalty'],c['eligible'])
                    packaged=float(parts['numerator'].sum()/parts['denominator'].sum()) if len(parts['denominator']) else None
                    if estimate is not None and not np.isclose(packaged,estimate,atol=1e-12):raise AssertionError('Independent estimator disagrees')
                    if parts['n_contributing_sites']!=sum(r['cytosolic_sites']+r['outer_sites'] for r in strata):raise AssertionError('Count disagreement')
                record=dict(model=model,window=window,comparator_name=comparator_name,comparator=comparator,
                    estimate_raw=estimate,ci95_group=group_ci,ci95_protein=protein_ci,
                    estimate_standardized_within_scorer=estimate/sd if sd and estimate is not None else None,
                    penalty_sd=sd,n_proteins=len(pp),n_groups=len(gp),n_strata=len(strata),
                    n_eligible_sites=len(eligible),n_contributing_sites=sum(r['cytosolic_sites']+r['outer_sites'] for r in strata),
                    n_contributing_cytosolic=sum(r['cytosolic_sites'] for r in strata),
                    n_contributing_outer=sum(r['outer_sites'] for r in strata),
                    contributing_accessions=sorted(pp),primary_group_verdict=verdict(group_ci),
                    packaged_point=packaged)
                results.append(record)
                if window==10 and comparator_name=='all_uncharged':
                    primary_parts[model]=gp
                    for row in strata:contributions.append(dict(model=model,group=group_for[row['accession']],**row))
        print('Audited '+model,flush=True)
    primary=[r for r in results if r['window']==10 and r['comparator_name']=='all_uncharged']
    changes=[]
    reference=primary_parts['guarded_oof_raw']
    for model,parts in primary_parts.items():
        if not model.startswith('guarded_oof__'):continue
        if set(reference)!=set(parts):raise AssertionError('Paired comparison uses different groups')
        difference={}
        for group in reference:
            if not np.isclose(reference[group][1],parts[group][1],atol=1e-12):raise AssertionError('Paired weights differ')
            difference[group]=[parts[group][0]-reference[group][0],parts[group][1]]
        point,ci=resample_parts(difference)
        changes.append(dict(model=model,reference='guarded_oof_raw',change_raw=point,ci95_group=ci))
    pi=next(r for r in changes if r['model']=='guarded_oof__positive_inside')
    checks['PI_only_moves_in_expected_direction']=bool(pi['change_raw']<=1e-12)
    if not all(checks.values()):raise AssertionError(checks)
    write('results.json',dict(primary=primary,paired_changes=changes,all_sensitivities=results,
         score_behavior_only=True,new_fit=False,experimental_outcomes_read=False,comparisons_selected_by_result=False))
    with (OUT/'primary_stratum_contributions.csv').open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(contributions[0]));writer.writeheader();writer.writerows(contributions)
    flat=[]
    for rec in results:
        row={k:v for k,v in rec.items() if k not in ['ci95_group','ci95_protein','contributing_accessions']}
        for kind in ['group','protein']:
            ci=rec['ci95_'+kind]
            row[kind+'_ci_low']=ci[0] if ci else None;row[kind+'_ci_high']=ci[1] if ci else None
        flat.append(row)
    with (OUT/'all_audits.csv').open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(flat[0]));writer.writeheader();writer.writerows(flat)
    unchanged=all(sha(ROOT/name)==value for name,value in SOURCES.items())
    if not unchanged:raise AssertionError('Input changed')
    write('VERIFICATION.json',dict(completed_utc=datetime.now(timezone.utc).isoformat(),checks=checks,
         independent_primary_points_match_packaged=True,source_files_unchanged=True,
         cached_generation_source_hashes_match=True,native_cells_masked=True,canonical_sequence_matches=True,
         explicit_side_conversion='standalone codes 1/2 become audit codes 0/1; unknown/TM excluded',
         total_scores_per_raw_head=193515,n_results=len(results),new_model_fits=0,new_correction_fits=0,
         experimental_target_reads=0,plan_sha256=sha(OUT/'AUDIT_PLAN.json'),results_sha256=sha(OUT/'results.json')))
    report(primary,changes,results)
    print(json.dumps(primary,indent=2),flush=True)


def report(primary,changes,results):
    names={'saved_full_fit_raw':'Saved full-fit raw','guarded_oof_raw':'Group-out raw'}
    def label(model):return names.get(model,model.replace('guarded_oof__','Group-out: ').replace('positive_inside','PI').replace('charge_burial','charge').replace('glycosylation_loss','glyco').replace('disulfide_loss','disulfide'))
    def number(value):return f'{value:+.6f}' if value is not None else 'NA'
    def ci(value):return f'[{number(value[0])}, {number(value[1])}]' if value else 'not estimable'
    lines=['# Positive-inside audit of the standalone predictor','',
        '28 October 2026. Completed on existing standalone predictions, with no new fit, new coefficient calibration, network data, or experimental-target reads.','',
        '## Primary result','',
        'The primary statistic is the cytosolic-minus-outer difference in within-position K/R-introduction penalty, matched within protein and native residue, at distances 1–10 residues from a membrane boundary. Negative supports the specified positive-inside expectation. The comparator is the original 15 uncharged residues, excluding the native identity cell.','',
        '| Output | Raw contrast | 95% group-resampling interval | Result |',
        '|---|---:|---|---|']
    for r in primary:lines.append(f"| {label(r['model'])} | {number(r['estimate_raw'])} | {ci(r['ci95_group'])} | {r['primary_group_verdict'].replace('_',' ')} |")
    r=primary[0]
    lines.extend(['',f"The primary comparison contains {r['n_contributing_sites']} contributing positions ({r['n_contributing_cytosolic']} cytosolic and {r['n_contributing_outer']} outer), {r['n_strata']} matched strata, {r['n_proteins']} proteins and {r['n_groups']} dependence groups. The broader score archive contains 10,185 positions and 193,515 missense scores per raw head across 13 membrane proteins; that is not the independent sample size.",'',
        'The saved full-fit head describes behavior on the fitted panel. Existing group-out heads exclude each scored dependence group from supervised fitting and correction calibration. Both use historically exposed project data; neither is new prospective confirmation.','',
        '## Correction movement','',
        '| Existing correction | Paired change versus group-out raw | 95% group-resampling interval |',
        '|---|---:|---|'])
    for r in changes:lines.append(f"| {label(r['model'])} | {number(r['change_raw'])} | {ci(r['ci95_group'])} |")
    lines.extend(['','The coefficients were read unchanged from the existing guarded-fold metadata. A negative movement is movement in the declared audit direction; it is not an experimental accuracy gain.','',
                  '## Comparator sensitivity at the primary window','',
                  '| Output | Comparator | Contrast | 95% group-resampling interval | Result |',
                  '|---|---|---:|---|---|'])
    selected=['saved_full_fit_raw','guarded_oof_raw','guarded_oof__positive_inside']
    for r in results:
        if r['window']==10 and r['model'] in selected:
            lines.append(f"| {label(r['model'])} | {r['comparator_name']} | {number(r['estimate_raw'])} | {ci(r['ci95_group'])} | {r['primary_group_verdict'].replace('_',' ')} |")
    lines.extend(['','## Window sensitivity with the original comparator','',
                  '| Output | Window | Contrast | 95% group-resampling interval | Result |',
                  '|---|---:|---:|---|---|'])
    for r in results:
        if r['comparator_name']=='all_uncharged' and r['model'] in selected:
            lines.append(f"| {label(r['model'])} | {r['window']} | {number(r['estimate_raw'])} | {ci(r['ci95_group'])} | {r['primary_group_verdict'].replace('_',' ')} |")
    lines.extend(['','## Scope and reproducibility','',
        'The primary interval resamples dependence groups 1,999 times with predictions fixed. Whole-protein intervals and within-scorer standardization are secondary columns in `all_audits.csv`; standardized magnitudes must not be used as cross-model effect-size rankings. Comparator and window sensitivities were recorded before this run and are all retained. A confidence interval spanning zero is inconclusive, not proof of no biological preference.','',
        'Canonical sequences, all native-cell masks, explicit side-code conversion and cached-generation source hashes were checked. A plain-Python estimator independently reproduces every primary-comparator point from the packaged audit. Known-answer direction, null, offset-cancellation and side-reversal checks passed. The original inputs and model/coefficient artifacts are unchanged.','',
        'Runner: `scripts/beat_am/m17_positive_inside_standalone.py`. Results: `results/positive_inside_standalone_20261028/`. The plan, full sensitivity table, per-stratum contributions and verification record are retained there.'])
    text='\n'.join(lines)+'\n'
    (OUT/'REPORT.md').write_text(text,encoding='utf-8')


if __name__=='__main__':main()

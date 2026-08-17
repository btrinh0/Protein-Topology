"""Conditional Pro-gain audit in independently defined native structural contexts.

Native Pro loss is excluded. A TM-context contrast is not by itself proof of a
membrane-specific physical mechanism or clinical pathogenicity. Geometry and
matching are fixed before score access; comparator/scale sensitivity is retained.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import hashlib
from importlib import resources
import json
import math
from pathlib import Path
import numpy as np

ALPHABET='ACDEFGHIKLMNPQRSTVWY'
NATIVES='ACFILMNQSTVWY'
COMPARATORS={'aliphatic':'AILMV','ser_thr':'ST','glycine_diagnostic':'G'}
PRIMARY='aliphatic'


def dihedral(a,b,c,d):
    """Signed backbone dihedral in degrees, rejecting degenerate coordinates."""
    a,b,c,d=(np.asarray(x,dtype=float) for x in (a,b,c,d))
    u=-(b-a);axis=c-b;v=d-c
    norm=np.linalg.norm(axis)
    if norm<1e-10: return None
    axis=axis/norm
    u=u-np.dot(u,axis)*axis;v=v-np.dot(v,axis)*axis
    if min(np.linalg.norm(u),np.linalg.norm(v))<1e-10: return None
    return float(np.degrees(np.arctan2(np.dot(np.cross(axis,u),v),np.dot(u,v))))


def helix_angles(phi,psi):
    return phi is not None and psi is not None and -100<=phi<=-30 and -80<=psi<=-5


def beta_like(phi,psi):
    return phi is not None and psi is not None and -180<=phi<=-70 and 50<=psi<=180


def axis_bend(ca):
    """Angle between oriented principal axes of overlapping five-CA windows.

    This reproducible local straightness filter is not a universal helix-axis
    definition. The nine-coordinate window and threshold are reported explicitly.
    """
    ca=np.asarray(ca,dtype=float)
    if ca.shape!=(9,3) or not np.isfinite(ca).all(): return None
    axes=[]
    for block in (ca[:5],ca[4:]):
        _,s,v=np.linalg.svd(block-block.mean(axis=0),full_matrices=False)
        if s[0]<1e-8: return None
        axis=v[0]
        if np.dot(axis,block[-1]-block[0])<0: axis=-axis
        axes.append(axis)
    return float(np.degrees(np.arccos(np.clip(np.dot(*axes),-1,1))))


def classify_geometry(sequence,pos,atoms,tm_segments,helical_positions,sheet_positions,
                      packing_count,bend_max=20):
    """Classify a canonically mapped native site, without predictor outcomes."""
    if not 1<=pos<=len(sequence) or sequence[pos-1] not in NATIVES:
        return None,'native_not_eligible'
    window=range(pos-4,pos+5)
    if min(window)<1 or max(window)>len(sequence): return None,'terminal_window'
    if any(sequence[p-1] in 'PG' for p in window): return None,'native_pro_or_gly_nearby'
    containing=[(s,e) for s,e in tm_segments if s<=pos<=e]
    full_window=all(p in atoms and 'CA' in atoms[p] for p in window)
    if containing and not full_window: return None,'unresolved_or_non_native_window'
    if any(p not in atoms or 'CA' not in atoms[p] for p in range(pos-1,pos+2)):
        return None,'unresolved_or_non_native_local_control'
    torsions={}
    for p in range(pos-2,pos+3):
        try:
            phi=dihedral(atoms[p-1]['C'],atoms[p]['N'],atoms[p]['CA'],atoms[p]['C'])
            psi=dihedral(atoms[p]['N'],atoms[p]['CA'],atoms[p]['C'],atoms[p+1]['N'])
        except KeyError:
            if p==pos or containing: return None,'missing_backbone_atoms'
            continue
        if phi is None or psi is None:
            if p==pos or containing: return None,'degenerate_backbone_geometry'
            continue
        torsions[p]=(phi,psi)
    regular=len(torsions)==5 and all(helix_angles(*torsions[p]) for p in torsions)
    bend=axis_bend([atoms[p]['CA'] for p in window]) if full_window else None
    regular=regular and bend is not None and bend<=bend_max
    intersect=[(s,e) for s,e in tm_segments if any(s<=p<=e for p in window)]
    boundaries=[x for segment in tm_segments for x in segment]
    distance=min(abs(pos-x) for x in boundaries) if boundaries else None
    context=None
    if containing:
        s,e=containing[0]
        if e-s+1<17 or pos-s<4 or e-pos<4: return None,'tm_helix_end_or_short_segment'
        if not regular: return None,'tm_window_not_geometrically_regular'
        context='tm_helix'
    elif not intersect and regular:
        context='non_tm_helix'
    elif not intersect and distance is not None and 5<=distance<=60:
        phi,psi=torsions[pos]
        if pos in helical_positions or pos in sheet_positions or helix_angles(phi,psi) or beta_like(phi,psi) or phi>0:
            return None,'control_not_coil_like'
        context='coil'
    if context is None: return None,'outside_declared_contexts'
    count=int(packing_count.get(pos,0))
    return {'position':pos,'native':sequence[pos-1],'context':context,
            'packing_count':count,'packing_bin':0 if count<=3 else 1 if count<=8 else 2,
            'phi':torsions[pos][0],'psi':torsions[pos][1],
            'local_axis_bend_degrees':bend,'distance_to_tm_boundary':distance},None


def position_contrasts(scores,natives):
    """Complete-case Pro gain contrasts; native/synonymous values are excluded."""
    scores=np.asarray(scores,dtype=float).copy()
    natives=np.asarray(natives).astype(str)
    if scores.shape!=(len(natives),20): raise ValueError('Scores must be sites x the declared 20-aa alphabet')
    if not set(natives)<=set(NATIVES): raise ValueError('Native Pro/Gly/charged sites are not Pro-gain eligibility')
    for i,native in enumerate(natives): scores[i,ALPHABET.index(native)]=np.nan
    means={};complete=np.isfinite(scores[:,ALPHABET.index('P')])
    for name,pool in COMPARATORS.items():
        arm=np.full(len(natives),np.nan)
        for native in np.unique(natives):
            sel=natives==native
            cols=[ALPHABET.index(a) for a in pool if a!=native]
            data=scores[sel][:,cols]
            good=np.isfinite(data).all(axis=1)
            values=np.full(len(data),np.nan)
            values[good]=data[good].mean(axis=1)
            arm[sel]=values
        means[name]=arm;complete &= np.isfinite(arm)
    all_measured=np.isfinite(scores).sum(axis=1)==19
    spread=np.full(len(scores),np.nan)
    spread[all_measured]=np.nanstd(scores[all_measured],axis=1,ddof=1)
    scale=np.max(np.where(np.isfinite(scores),np.abs(scores),0),axis=1)
    stable=all_measured & (spread>64*np.finfo(float).eps*np.maximum(scale,np.finfo(float).tiny))
    result={}
    for name,mean in means.items():
        z=scores[:,ALPHABET.index('P')]-mean
        z[~complete]=np.nan
        normalized=np.full(len(z),np.nan)
        use=complete & stable
        normalized[use]=z[use]/spread[use]
        result[name]={'raw':z,'position_scaled':normalized}
    return result,{'n_requested':len(scores),'n_complete_common_positions':int(complete.sum()),
                   'n_complete_position_scaled':int((complete&stable).sum()),
                   'n_missing_or_incomplete':int((~complete).sum()),
                   'n_incomplete_or_zero_spread':int((~stable).sum())}


def context_effect(sites,values,target='tm_helix',reference='coil',group_field='family',
                   draws=1999,seed=20261025,min_groups=6):
    """Protein/native/packing matched contrast with whole-group resampling."""
    if len(sites)!=len(values) or draws<2 or min_groups<2: raise ValueError('Invalid audit dimensions or uncertainty configuration')
    buckets=defaultdict(lambda:defaultdict(list))
    group_for={}
    for site,value in zip(sites,values):
        if not np.isfinite(value) or site['context'] not in {target,reference}: continue
        key=(site['accession'],site['native'],site['packing_bin'])
        group=site[group_field]
        if site['accession'] in group_for and group_for[site['accession']]!=group:
            raise ValueError('One protein has inconsistent uncertainty-group assignments')
        group_for[site['accession']]=group
        buckets[key][site['context']].append(float(value))
    parts=defaultdict(lambda:[0.0,0.0]);strata=[];proteins=set()
    for key,cells in sorted(buckets.items()):
        if not cells[target] or not cells[reference]: continue
        nt,nr=len(cells[target]),len(cells[reference])
        weight=nt*nr/(nt+nr)
        difference=float(np.mean(cells[target])-np.mean(cells[reference]))
        group=group_for[key[0]];parts[group][0]+=weight*difference;parts[group][1]+=weight
        proteins.add(key[0]);strata.append({'accession':key[0],'native':key[1],'packing_bin':key[2],
            'n_target':nt,'n_reference':nr,'difference':difference,'weight':weight,'group':group})
    keys=sorted(parts);n=len(keys)
    result={'target_context':target,'reference_context':reference,'group_field':group_field,
            'n_strata':len(strata),'n_proteins':len(proteins),'n_groups':n,
            'n_target_sites':sum(r['n_target'] for r in strata),
            'n_reference_sites':sum(r['n_reference'] for r in strata),
            'estimate':None,'ci95':None,'verdict':'insufficient_replication',
            'strata':strata,'group_effects':{},'draws':draws,'seed':seed}
    if not n: return result
    numerator=np.array([parts[k][0] for k in keys]);denominator=np.array([parts[k][1] for k in keys])
    estimate=float(numerator.sum()/denominator.sum())
    result.update(estimate=estimate,group_effects={k:float(parts[k][0]/parts[k][1]) for k in keys},
                  equal_group_mean=float(np.mean(numerator/denominator)))
    if n<2: return result
    rng=np.random.default_rng(seed);samples=[]
    for start in range(0,draws,128):
        counts=rng.multinomial(n,np.full(n,1/n),size=min(128,draws-start))
        samples.extend((counts@numerator/(counts@denominator)).tolist())
    lo,hi=np.percentile(samples,[2.5,97.5]).tolist()
    magnitude=max((abs(x) for v in buckets.values() for cell in v.values() for x in cell),default=0)
    tolerance=float(64*np.finfo(float).eps*magnitude)
    verdict=('supported_conditional_direction' if lo>tolerance else
             'opposite_conditional_direction' if hi < -tolerance else 'indeterminate')
    result.update(ci95=[lo,hi],numerical_tolerance=tolerance,
                  verdict=verdict if n>=min_groups else 'insufficient_replication')
    return result


def audit_scores(sites,scores,*,score_direction,draws=1999,seed=20261025,min_groups=6):
    if score_direction not in {'higher-is-damaging','higher-is-tolerated'}:
        raise ValueError('Score direction must be declared')
    keys=[(r['accession'],r['position']) for r in sites]
    if len(set(keys))!=len(keys): raise ValueError('One scored row per canonical native position is required')
    factor=1 if score_direction=='higher-is-damaging' else -1
    contrasts,coverage=position_contrasts(factor*np.asarray(scores),[r['native'] for r in sites])
    results={}
    for name,scales in contrasts.items():
        results[name]={}
        for scale,values in scales.items():
            results[name][scale]={}
            for field in ('accession','family'):
                results[name][scale][field]=context_effect(sites,values,group_field=field,
                    draws=draws,seed=seed,min_groups=min_groups)
            results[name][scale]['tm_minus_non_tm_helix_diagnostic']=context_effect(
                sites,values,reference='non_tm_helix',group_field='family',
                draws=draws,seed=seed,min_groups=min_groups)
            results[name][scale]['tm_minus_non_tm_helix_diagnostic']['direction_prespecified']=False
            results[name][scale]['tm_minus_non_tm_helix_diagnostic']['verdict']='descriptive_no_direction_prespecified'
    main=results[PRIMARY]['raw']['family']
    sensitivity=results['ser_thr']['raw']['family']
    scaled=results[PRIMARY]['position_scaled']['family']
    robust=all(r['verdict']=='supported_conditional_direction' for r in (main,sensitivity,scaled))
    opposite=all(r['verdict']=='opposite_conditional_direction' for r in (main,sensitivity,scaled))
    status=('conditional_compliance_robust' if robust else
            'candidate_opposite_direction_requires_experimental_validation' if opposite else
            'mixed_or_insufficient_evidence')
    return {'rule':'proline_introduction_in_regular_tm_helices','input_score_direction':score_direction,
            'output_score_direction':'higher-is-damaging','coverage':coverage,'results':results,
            'overall_status':status,'new_biological_failure_demonstrated':False,
            'correction_policy':'zero_for_a_supported_rule; no nonzero correction without independent endpoint evidence',
            'clinical_interpretation':False}


def average_tie_ranks(values):
    """NumPy-only within-assay average ranks, with finite values required."""
    values=np.asarray(values,dtype=float)
    if not np.isfinite(values).all(): raise ValueError('Rank only finite measured outcomes')
    order=np.argsort(values,kind='stable');sorted_values=values[order]
    boundaries=np.r_[0,np.flatnonzero(np.diff(sorted_values)!=0)+1,len(values)]
    result=np.empty(len(values),dtype=float)
    for a,b in zip(boundaries[:-1],boundaries[1:]): result[order[a:b]]=(a+1+b)/2
    return result


def load_score_archive(path,sites):
    """Align a named, native-identity-checked NPZ matrix to the structural bank."""
    with np.load(path,allow_pickle=False) as data:
        required={'scores','accession','position','native','alphabet'}
        if not required<=set(data.files): raise ValueError('Score archive lacks required named coordinate/native metadata')
        if ''.join(data['alphabet'].astype(str).tolist())!=ALPHABET:
            raise ValueError('Unexpected score alphabet')
        values=data['scores'];accession=data['accession'].astype(str)
        positions=data['position'];native=data['native'].astype(str)
        if values.shape!=(len(accession),20) or len(positions)!=len(accession) or len(native)!=len(accession):
            raise ValueError('Conflicting score archive dimensions')
        index={}
        for i,(acc,pos) in enumerate(zip(accession,positions)):
            if int(pos)!=pos or int(pos)<1: raise ValueError('Invalid canonical score position')
            key=(acc,int(pos))
            if key in index: raise ValueError('Duplicate canonical score row')
            index[key]=i
        result=np.full((len(sites),20),np.nan)
        for i,site in enumerate(sites):
            row=index.get((site['accession'],int(site['position'])))
            if row is None: continue
            if native[row]!=site['native']: raise ValueError('Score archive native identity disagrees with structural bank')
            result[i]=values[row]
    return result


def run_proline_audit(sites_path,scores_path,out_dir,*,score_direction,draws=1999,
                      seed=20261025,min_groups=6,model='supplied_predictor'):
    sites=json.loads(Path(sites_path).read_text(encoding='utf-8'))
    scores=load_score_archive(scores_path,sites)
    result=audit_scores(sites,scores,score_direction=score_direction,draws=draws,seed=seed,min_groups=min_groups)
    def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    result.update(model=model,sites_sha256=digest(sites_path),scores_sha256=digest(scores_path),
                  implementation_sha256=hashlib.sha256(resources.files(__package__).joinpath('proline.py').read_bytes()).hexdigest(),
                  context_provenance='Caller-supplied structural site bank; this command audits scores, not experimental annotation truth')
    out_dir=Path(out_dir);out_dir.mkdir(parents=True,exist_ok=True)
    (out_dir/'proline_audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sites',type=Path,required=True)
    parser.add_argument('--scores',type=Path,required=True)
    parser.add_argument('--score-direction',choices=('higher-is-damaging','higher-is-tolerated'),required=True)
    parser.add_argument('--model',default='supplied_predictor')
    parser.add_argument('--draws',type=int,default=1999)
    parser.add_argument('--seed',type=int,default=20261025)
    parser.add_argument('--min-groups',type=int,default=6)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    result=run_proline_audit(args.sites,args.scores,args.out,score_direction=args.score_direction,
        draws=args.draws,seed=args.seed,min_groups=args.min_groups,model=args.model)
    primary=result['results'][PRIMARY]['raw']['family']
    print(json.dumps({'model':args.model,'status':result['overall_status'],
        'primary':{k:primary[k] for k in ('estimate','ci95','n_proteins','n_groups','n_target_sites','n_reference_sites','verdict')}},indent=2))
    return 0


if __name__=='__main__': raise SystemExit(main())

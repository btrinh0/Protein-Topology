"""Audit structural disulfide-loss scores with evidence-qualified matched controls."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np

VERSION = '1.0.0'
MATCH_FIELDS = ('accession', 'domain', 'compartment', 'secondary_structure')
UNKNOWN = {'', 'unknown', 'unreviewed', 'unclassified', 'na', 'n/a', 'none'}
PRIMARY_ALT = 'S'
SENSITIVITY_ALT = 'A'


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as fh:
        for chunk in iter(lambda: fh.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def sequence_sha256(sequence):
    return hashlib.sha256(sequence.encode('ascii')).hexdigest()


def _true(value):
    return value is True or str(value).strip().lower() in {'true', '1', 'yes'}


def _false(value):
    return value is False or str(value).strip().lower() in {'false', '0', 'no'}


def _known(value):
    return str(value or '').strip().lower() not in UNKNOWN


def _number(value):
    try:
        x = float(value)
    except (ValueError, TypeError):
        return None
    return x if math.isfinite(x) else None


def eligible_sites(rows, sequences):
    """Return validated eligible sites and explicit exclusion reasons.

    sequences maps accession to the actual canonical amino-acid sequence.
    Wrong native identity, duplicate/conflicting sites, stale sequence hashes,
    inconsistent family roles and invalid bond endpoints are hard errors.
    Incomplete evidence/context instead produces an auditable exclusion.
    """
    sites, excluded = [], []
    all_annotations = {}
    seen, families, protein_families = set(), {}, {}
    for raw in rows:
        row = dict(raw)
        acc = str(row.get('accession', '')).strip()
        try:
            pos = int(row.get('position', ''))
        except (TypeError, ValueError):
            raise ValueError(f'Invalid residue coordinate for {acc}') from None
        key = (acc, pos)
        if key in seen:
            raise ValueError(f'Duplicate site annotation: {key}')
        seen.add(key)
        sequence = sequences.get(acc)
        if not sequence or not 1 <= pos <= len(sequence):
            raise ValueError(f'Missing canonical sequence or out-of-range position: {key}')
        if sequence[pos-1] != 'C' or row.get('native') != 'C':
            raise ValueError(f'Native cysteine mismatch: {key}')
        if row.get('sequence_sha256') != sequence_sha256(sequence):
            raise ValueError(f'Stale or missing canonical-sequence hash: {key}')
        row['position'] = pos
        row['site_id'] = f'{acc}:{pos}'
        all_annotations[key] = row
        family = str(row.get('protein_family', '')).strip()
        role = str(row.get('role', '')).strip()
        if _known(family):
            old = protein_families.setdefault(acc, family)
            if old != family:
                raise ValueError(f'One protein assigned to multiple families: {acc}')
            old_role = families.setdefault(family, role)
            if old_role != role:
                raise ValueError(f'One family crosses data roles: {family}')
        reasons = []
        if not _known(family):
            reasons.append('unknown_family')
        if role not in {'development', 'calibration', 'confirmation'}:
            reasons.append('unknown_data_role')
        if role == 'confirmation' and row.get('exposure_status') != 'unseen_and_homology_checked':
            reasons.append('confirmation_exposure_not_cleared')
        if not _true(row.get('mapping_verified')):
            reasons.append('sequence_mapping_not_verified')
        if row.get('evidence_type') != 'experimental' or not _known(row.get('evidence_source')):
            reasons.append('no_direct_experimental_bond_state_evidence')
        if not _true(row.get('native_context_verified')):
            reasons.append('native_chemical_context_not_verified')
        state = row.get('bond_state')
        if state == 'bonded':
            try:
                partner = int(row.get('bond_partner', ''))
            except (TypeError, ValueError):
                raise ValueError(f'Missing bonded-Cys partner: {key}') from None
            if partner == pos or not 1 <= partner <= len(sequence) or sequence[partner-1] != 'C':
                raise ValueError(f'Invalid canonical bond partner: {key}, {partner}')
            row['bond_partner'] = partner
            row['bond_id'] = f'{acc}:{min(pos, partner)}-{max(pos, partner)}'
            if row.get('bond_role') != 'structural':
                reasons.append('bond_role_not_verified_structural')
        elif state == 'free':
            if _known(row.get('bond_partner')) and str(row.get('bond_partner')) != '0':
                raise ValueError(f'Free cysteine also assigned a bond partner: {key}')
            if not _true(row.get('control_reviewed')):
                reasons.append('free_control_not_reviewed')
        else:
            reasons.append('unknown_bond_state_is_not_free')
        for flag in ('metal_ligand', 'signal_peptide', 'other_functional_cysteine', 'redox_or_catalytic'):
            if not _false(row.get(flag)):
                reasons.append(f'{flag}_present_or_unreviewed')
        for field in MATCH_FIELDS[1:]:
            if not _known(row.get(field)):
                reasons.append(f'unknown_{field}')
        if row.get('secondary_structure') not in {'H', 'E', 'C'}:
            reasons.append('invalid_secondary_structure')
        rsa = _number(row.get('rsa'))
        if rsa is None or not 0 <= rsa <= 1:
            reasons.append('unknown_or_invalid_rsa')
        else:
            row['rsa'] = rsa
        if not _known(row.get('context_source')):
            reasons.append('missing_structural_context_source')
        if reasons:
            excluded.append({'site_id': row['site_id'], 'reasons': sorted(set(reasons))})
        else:
            row['_eligibility_verified'] = True
            sites.append(row)
    for row in all_annotations.values():
        if row.get('bond_state') != 'bonded':
            continue
        other = all_annotations.get((row['accession'],row['bond_partner']))
        if other and (other.get('bond_state') == 'free' or
                      (other.get('bond_state') == 'bonded' and
                       other['bond_partner'] != row['position'])):
            raise ValueError(f'Conflicting bond states/partners at {row["site_id"]}')
    return sites, excluded


def match_sites(sites, max_rsa_difference=0.15):
    """Annotation-only deterministic 1:1 matching without replacement."""
    if not 0 <= max_rsa_difference <= 1:
        raise ValueError('RSA caliper must be between zero and one')
    strata = defaultdict(list)
    for row in sites:
        strata[tuple(row[f] for f in MATCH_FIELDS)].append(row)
    pairs, unmatched = [], []
    for group in sorted(strata):
        rows = strata[group]
        bonded = [r for r in rows if r['bond_state'] == 'bonded']
        free = [r for r in rows if r['bond_state'] == 'free']
        candidates = []
        for b in bonded:
            for f in free:
                distance = abs(b['rsa']-f['rsa'])
                if distance <= max_rsa_difference + 1e-12:
                    candidates.append((distance, b['position'], f['position'], b, f))
        used_b, used_f = set(), set()
        for distance, bp, fp, b, f in sorted(candidates, key=lambda x:x[:3]):
            if bp in used_b or fp in used_f:
                continue
            used_b.add(bp)
            used_f.add(fp)
            pairs.append({'pair_id':f'{b["site_id"]}|{f["site_id"]}',
                          'bonded':b, 'free':f, 'rsa_difference':float(distance),
                          'protein_family':b['protein_family'], 'role':b['role']})
        unmatched.extend({'site_id':r['site_id'], 'bond_state':r['bond_state']}
                         for r in rows if r['position'] not in
                         (used_b if r['bond_state']=='bonded' else used_f))
    return pairs, unmatched


def qualification(pairs, partition, min_families):
    if partition not in {'development', 'calibration', 'confirmation', 'all'}:
        raise ValueError('Declare an existing data partition')
    if min_families < 2:
        raise ValueError('At least two families required; frozen confirmatory default is six')
    selected = [p for p in pairs if partition == 'all' or p['role']==partition]
    families = sorted({p['protein_family'] for p in selected})
    reasons = []
    if len(families) < min_families:
        reasons.append(f'Only {len(families)} matched families; need {min_families}')
    # Recheck the exposure gate immediately before authorizing score access.
    if partition == 'confirmation':
        if any(r.get('exposure_status') != 'unseen_and_homology_checked'
               for p in selected for r in (p['bonded'],p['free'])):
            reasons.append('Confirmation exposure/homology clearance is incomplete')
    return {'status':'ready' if not reasons else 'blocked_insufficient_evidence',
            'partition':partition,'minimum_families':min_families,
            'n_pairs':len(selected),'n_families':len(families),
            'n_proteins':len({p['bonded']['accession'] for p in selected}),
            'families':families,'blocking_reasons':reasons}, selected


def read_scores(path, requested_sites, alternates=('S','A')):
    """Read explicit canonical variant-score CSV/TSV, retaining only requested Cys.

    Input columns: accession, position, reference, alternate, score. Duplicate
    variants are allowed only when numerically identical. Synonymous entries,
    nonfinite values and conflicting native identities are rejected for selected
    sites. AlphaMissense's raw released schema needs an explicit extraction step.
    """
    path=Path(path)
    opener=gzip.open if path.suffix=='.gz' else open
    delimiter='\t' if '.tsv' in path.name else ','
    table={}
    with opener(path,'rt',encoding='utf-8-sig',newline='') as fh:
        reader=csv.DictReader(fh,delimiter=delimiter)
        required={'accession','position','reference','alternate','score'}
        if not required <= set(reader.fieldnames or []):
            raise ValueError('Scores need canonical accession,position,reference,alternate,score columns')
        for row in reader:
            try:
                position=int(row['position'])
            except (ValueError,TypeError):
                raise ValueError('Invalid canonical score coordinate') from None
            site=(row['accession'],position)
            if site not in requested_sites or row['alternate'] not in alternates:
                continue
            if row['reference'] != 'C':
                raise ValueError(f'Conflicting native residue in score table: {site}')
            value=_number(row['score'])
            if value is None:
                raise ValueError(f'Nonfinite selected variant score: {site}')
            key=(*site,row['alternate'])
            if key in table and table[key] != value:
                raise ValueError(f'Conflicting duplicate variant score: {key}')
            table[key]=value
    return table


def audit_pairs(pairs, scores, alternate='S', score_direction=None,
                draws=1999, seed=20261025, min_families=6):
    if alternate not in {'S','A'}:
        raise ValueError('The frozen protocol specifies C>S primary and C>A sensitivity')
    if score_direction not in {'higher-is-damaging','higher-is-tolerated'}:
        raise ValueError('Declare the score direction explicitly')
    if draws < 2 or min_families < 2:
        raise ValueError('At least two bootstrap draws and families required')
    sign=1 if score_direction=='higher-is-damaging' else -1
    per_family=defaultdict(list)
    measured=[]
    missing=[]
    for pair in pairs:
        b,f=pair['bonded'],pair['free']
        bk=(b['accession'],b['position'],alternate)
        fk=(f['accession'],f['position'],alternate)
        if bk not in scores or fk not in scores:
            missing.append(pair['pair_id'])
            continue
        bv,fv=_number(scores[bk]),_number(scores[fk])
        if bv is None or fv is None:
            raise ValueError(f'Nonfinite paired value: {pair["pair_id"]}')
        difference=sign*(bv-fv)
        per_family[pair['protein_family']].append(difference)
        measured.append({'pair_id':pair['pair_id'],'protein_family':pair['protein_family'],
                         'bonded_score_oriented':sign*bv,'free_score_oriented':sign*fv,
                         'difference':difference})
    means={k:float(np.mean(v)) for k,v in sorted(per_family.items())}
    n=len(means)
    result={'rule':'verified_structural_disulfide_loss','version':VERSION,
            'substitution':'C>'+alternate,'expected_direction':'bonded_minus_free > 0',
            'estimand':'equal-family mean of annotation-matched identical-substitution differences',
            'input_score_direction':score_direction,'output_score_direction':'higher-is-damaging',
            'n_pairs_requested':len(pairs),'n_pairs_measured':len(measured),
            'missing_pair_ids':missing,'n_families':n,'family_effects':means,
            'estimate_raw':None,'ci95_raw':None,'draws':draws,'seed':seed,
            'paired_values':measured,'verdict':'blocked_insufficient_evidence',
            'clinical_pathogenicity_verdict':False}
    if n < min_families:
        result['blocking_reasons']=[f'Only {n} measured families; need {min_families}']
        return result
    values=np.array(list(means.values()),dtype=float)
    rng=np.random.default_rng(seed)
    samples=[]
    for start in range(0,draws,256):
        counts=rng.multinomial(n,np.full(n,1/n),size=min(256,draws-start))
        samples.extend((counts@values/n).tolist())
    lo,hi=np.percentile(samples,[2.5,97.5]).tolist()
    scale=max((abs(r[k]) for r in measured for k in
               ('bonded_score_oriented','free_score_oriented')),default=0.0)
    tolerance=float(64*np.finfo(float).eps*scale)
    result.update(estimate_raw=float(values.mean()),ci95_raw=[float(lo),float(hi)],
                  numerical_tolerance=tolerance,
                  interval_unit='homology family; all bond endpoints remain together',
                  verdict=('supported_conditional_direction' if lo>tolerance else
                           'opposite_conditional_direction' if hi < -tolerance else 'indeterminate'))
    return result


def fit_correction(calibration_pairs, scores, validation_result, *, alternate='S',
                   score_direction=None, target_margin=None, draws=1999,
                   seed=20261025, min_families=3):
    """Fit one nonnegative offset on CALIBRATION families only.

    The caller must provide an endpoint-aligned experimental validation result
    from development/calibration data, never confirmation data. Experimental
    measurements are independent of predictor outputs; calibration families may
    be shared, as clarified before outcome access in protocol amendment 01.
    This prerequisite validates the expectation, not correction benefit.
    """
    if not calibration_pairs or any(p['role']!='calibration' for p in calibration_pairs):
        raise ValueError('Only explicitly assigned calibration pairs may fit a correction')
    margin=_number(target_margin)
    if margin is None or margin <= 0:
        raise ValueError('Declare a finite positive target margin on this score scale')
    if validation_result.get('kind')!='experimental_validation' or \
            validation_result.get('verdict')!='supported_conditional_direction' or \
            validation_result.get('substitution')!='C>'+alternate or \
            not _known(validation_result.get('assay_endpoint')) or \
            not _known(validation_result.get('source_sha256')) or \
            validation_result.get('n_families',0)<min_families or \
            validation_result.get('partition') not in {'development','calibration'} or \
            not _true(validation_result.get('assay_context_reviewed')):
        raise ValueError('Independent endpoint-aligned experimental validation is required')
    ci=validation_result.get('ci95_raw')
    if not isinstance(ci,list) or len(ci)!=2 or any(_number(x) is None for x in ci) or \
            not 0 < ci[0] <= ci[1]:
        raise ValueError('Experimental validation has no positive finite interval')
    families={p['protein_family'] for p in calibration_pairs}
    result=audit_pairs(calibration_pairs,scores,alternate,score_direction,
                       draws,seed,min_families)
    if result['ci95_raw'] is None:
        raise ValueError('Insufficient measured calibration families')
    beta=max(0.0,margin-result['ci95_raw'][0])
    return {'rule':'verified_structural_disulfide_loss','version':VERSION,
            'alternate':alternate,'beta_oriented':float(beta),'target_margin':margin,
            'calibration_families':sorted(families),'calibration_audit':result,
            'experimental_validation':validation_result,
            'scope':'verified eligible structural-disulfide endpoints, declared replacement only',
            'status':'zero_correction_already_passes' if beta==0 else 'calibrated_pending_held_out_evaluation',
            'held_out_accuracy_improvement_established':False}


def apply_correction(variants, sites, fitted, score_direction):
    """Return copies with only eligible Cys-loss scores shifted on the raw scale."""
    if score_direction not in {'higher-is-damaging','higher-is-tolerated'}:
        raise ValueError('Declare the score direction explicitly')
    beta=_number(fitted.get('beta_oriented'))
    alt=fitted.get('alternate')
    if beta is None or beta < 0 or alt not in {'S','A'}:
        raise ValueError('Invalid disulfide correction artifact')
    if fitted.get('rule') != 'verified_structural_disulfide_loss' or \
            fitted.get('version') != VERSION or \
            fitted.get('experimental_validation',{}).get('kind') != 'experimental_validation':
        raise ValueError('Use a validated fit_correction artifact')
    if any(r.get('_eligibility_verified') is not True for r in sites):
        raise ValueError('Apply corrections only to the output of eligible_sites')
    eligible={(r['accession'],int(r['position'])) for r in sites if r['bond_state']=='bonded'}
    shift=beta if score_direction=='higher-is-damaging' else -beta
    result=[]
    for raw in variants:
        row=dict(raw)
        key=(row['accession'],int(row['position']))
        if key in eligible and row['reference']=='C' and row['alternate']==alt:
            value=_number(row['score'])
            if value is None:
                raise ValueError('Cannot correct a nonfinite selected score')
            row['score']=value+shift
        result.append(row)
    return result


def prepare_inventory(sites_path, sequences_path, partition='confirmation',
                      min_families=6, max_rsa_difference=0.15):
    sequences=json.loads(Path(sequences_path).read_text(encoding='utf-8'))
    with Path(sites_path).open(encoding='utf-8-sig',newline='') as fh:
        raw=list(csv.DictReader(fh))
    sites,excluded=eligible_sites(raw,sequences)
    pairs,unmatched=match_sites(sites,max_rsa_difference)
    gate,selected=qualification(pairs,partition,min_families)
    return {'tool':'disulfide-rule-audit','version':VERSION,
            'sites_sha256':sha256(sites_path),'sequences_sha256':sha256(sequences_path),
            'n_input_sites':len(raw),'n_eligible_sites':len(sites),
            'input_data_kinds':sorted({r.get('dataset_kind','unspecified') for r in raw}),
            'n_eligible_bonded':sum(r['bond_state']=='bonded' for r in sites),
            'n_eligible_free':sum(r['bond_state']=='free' for r in sites),
            'exclusions':excluded,'unmatched_sites':unmatched,
            'gate':gate,'score_files_opened':False,
            'matching':{'exact':list(MATCH_FIELDS),'rsa_caliper':max_rsa_difference,
                        'without_replacement':True,'uses_score_values':False}}, sites, selected


def run_disulfide_audit(sites_path, sequences_path, scores_path, out_dir, *,
                        score_direction, partition='confirmation', min_families=6,
                        max_rsa_difference=0.15, draws=1999, seed=20261025,
                        kind='predictor_audit', assay_endpoint=None, assay_context_reviewed=False):
    if score_direction not in {'higher-is-damaging','higher-is-tolerated'}:
        raise ValueError('Declare the score direction explicitly')
    if draws < 2:
        raise ValueError('At least two bootstrap draws required')
    if kind not in {'predictor_audit','experimental_validation'}:
        raise ValueError('Declare the evidence kind')
    if kind=='experimental_validation' and not _known(assay_endpoint):
        raise ValueError('Experimental validation requires the measured assay endpoint')
    inventory,sites,pairs=prepare_inventory(sites_path,sequences_path,partition,
                                           min_families,max_rsa_difference)
    out_dir=Path(out_dir)
    out_dir.mkdir(parents=True,exist_ok=True)
    (out_dir/'eligibility.json').write_text(json.dumps(inventory,indent=2)+'\n',encoding='utf-8')
    (out_dir/'matched_pairs.json').write_text(json.dumps(pairs,indent=2)+'\n',encoding='utf-8')
    if inventory['gate']['status']!='ready':
        result={'status':'blocked_insufficient_evidence','score_files_opened':False,
                'gate':inventory['gate'],'new_model_failure_demonstrated':False}
    else:
        requested={(r['accession'],r['position']) for p in pairs for r in (p['bonded'],p['free'])}
        scores=read_scores(scores_path,requested)
        results={}
        for alt in ('S','A'):
            entry=audit_pairs(pairs,scores,alt,score_direction,draws,seed,min_families)
            entry.update(kind=kind,assay_endpoint=assay_endpoint,
                         partition=partition,
                         assay_context_reviewed=bool(assay_context_reviewed),
                         source_sha256=sha256(scores_path))
            results['C>'+alt]=entry
        status=('blocked_insufficient_measured_coverage'
                if results['C>S']['verdict']=='blocked_insufficient_evidence' else 'audited')
        result={'status':status,'score_files_opened':True,'partition':partition,
                'input_data_kinds':inventory['input_data_kinds'],
                'kind':kind,'source_sha256':sha256(scores_path),'results':results,
                'new_model_failure_demonstrated':False,
                'interpretation':'Conditional score direction; experimental endpoint validation is required before calling a predictor discrepancy a biological error.'}
    (out_dir/'disulfide_audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sites',type=Path,required=True)
    parser.add_argument('--sequences',type=Path,required=True)
    parser.add_argument('--scores',type=Path,required=True)
    parser.add_argument('--score-direction',required=True,
                        choices=('higher-is-damaging','higher-is-tolerated'))
    parser.add_argument('--partition',choices=('development','calibration','confirmation','all'),
                        default='confirmation')
    parser.add_argument('--min-families',type=int,default=6)
    parser.add_argument('--draws',type=int,default=1999)
    parser.add_argument('--seed',type=int,default=20261025)
    parser.add_argument('--rsa-caliper',type=float,default=0.15)
    parser.add_argument('--kind',choices=('predictor_audit','experimental_validation'),default='predictor_audit')
    parser.add_argument('--assay-endpoint')
    parser.add_argument('--assay-context-reviewed',action='store_true',
                        help='Affirm the assay measures an appropriate endpoint in a context permitting these native bonds')
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    result=run_disulfide_audit(args.sites,args.sequences,args.scores,args.out,
        score_direction=args.score_direction,partition=args.partition,min_families=args.min_families,
        max_rsa_difference=args.rsa_caliper,draws=args.draws,seed=args.seed,
        kind=args.kind,assay_endpoint=args.assay_endpoint,
        assay_context_reviewed=args.assay_context_reviewed)
    summary={k:result[k] for k in ('status','score_files_opened','partition','input_data_kinds','gate') if k in result}
    if 'results' in result:
        summary['results']={name:{k:entry[k] for k in
            ('substitution','verdict','n_families','n_pairs_measured','estimate_raw','ci95_raw')}
            for name,entry in result['results'].items()}
    print(json.dumps(summary,indent=2))
    return 2 if result['status'].startswith('blocked_') else 0


if __name__=='__main__':
    raise SystemExit(main())

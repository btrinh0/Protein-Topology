"""Test disulfide evidence gates, matched controls and correction invariants."""
from copy import deepcopy
import csv
from contextlib import contextmanager
import json
from pathlib import Path
import random
import sys
import shutil
import unittest
import uuid
from unittest.mock import patch

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from membrane_rule_audit import disulfide as d
from membrane_rule_audit import panel


@contextmanager
def test_directory():
    # Use inherited-permission workspace folders because the Windows sandbox rejects private tempfile ACLs.
    parent=ROOT/'results/disulfide_rule_20261025/software_tests'
    root=parent/f'test-{uuid.uuid4().hex}'
    root.mkdir(parents=True)
    try:
        yield root
    finally:
        resolved=root.resolve()
        if resolved.parent != parent.resolve():
            raise RuntimeError('Refusing cleanup outside the dedicated test directory')
        shutil.rmtree(resolved)


def fixture(n=9, role='confirmation', offset=0.3):
    sequences={f'P{i}':'ACACACACAC' for i in range(n)}
    rows=[]
    scores={}
    for i,(acc,seq) in enumerate(sequences.items()):
        common={'accession':acc,'native':'C','sequence_sha256':d.sequence_sha256(seq),
                'sequence_version':'1','protein_family':f'F{i}','role':role,
                'exposure_status':'unseen_and_homology_checked',
                'bond_role':'structural','evidence_type':'experimental',
                'evidence_source':'synthetic known-answer experimental record',
                'native_context_verified':True,'mapping_verified':True,
                'control_reviewed':True,'metal_ligand':False,'signal_peptide':False,
                'other_functional_cysteine':False,'redox_or_catalytic':False,
                'domain':'D1','compartment':'extracellular','secondary_structure':'C',
                'context_source':'synthetic reference structure'}
        for pos,state,rsa in [(2,'bonded',0.2),(4,'bonded',0.3),(6,'free',0.22),(8,'free',0.32)]:
            rows.append({**common,'position':pos,'bond_state':state,'rsa':rsa,
                         'bond_partner':4 if pos==2 else 2 if pos==4 else ''})
            baseline=i*10+0.1
            for alt in ('S','A'):
                scores[(acc,pos,alt)]=baseline+(offset if state=='bonded' else 0)
    return sequences,rows,scores


def bank(n=9,role='confirmation',offset=0.3):
    seq,raw,scores=fixture(n,role,offset)
    sites,excluded=d.eligible_sites(raw,seq)
    assert not excluded
    pairs,_=d.match_sites(sites)
    return seq,sites,pairs,scores


def experiment():
    _,_,pairs,scores=bank(6,offset=0.4)
    # Keep evaluation evidence families disjoint from correction fitting.
    for pair in pairs:
        pair['protein_family']='EX_'+pair['protein_family']
    result=d.audit_pairs(pairs,scores,score_direction='higher-is-damaging',draws=199)
    return {**result,'kind':'experimental_validation','assay_endpoint':'secretion',
            'partition':'development',
            'assay_context_reviewed':True,'source_sha256':'synthetic-experimental-source'}


class EvidenceGates(unittest.TestCase):
    def test_unknown_bond_state_never_becomes_free(self):
        seq,rows,_=fixture(1)
        rows[2].update(bond_state='unknown',evidence_type='unknown',evidence_source='')
        sites,excluded=d.eligible_sites(rows,seq)
        self.assertEqual(len(sites),3)
        self.assertIn('unknown_bond_state_is_not_free',excluded[0]['reasons'])

    def test_native_identity_and_sequence_hash_are_verified(self):
        for change in ({'position':1},{'sequence_sha256':'old-version'}):
            seq,rows,_=fixture(1);rows[0].update(change)
            with self.assertRaises(ValueError): d.eligible_sites(rows,seq)

    def test_unreviewed_free_context_excluded(self):
        for change in ({'control_reviewed':False},{'native_context_verified':False},
                       {'evidence_type':'inferred'},{'signal_peptide':True},
                       {'metal_ligand':'unknown'},{'other_functional_cysteine':True}):
            seq,rows,_=fixture(1);rows[2].update(change)
            _,excluded=d.eligible_sites(rows,seq)
            self.assertTrue(excluded,change)

    def test_unknown_geometry_is_not_a_matching_category(self):
        seq,rows,_=fixture(1);rows[2].update(secondary_structure='unknown',rsa='')
        sites,excluded=d.eligible_sites(rows,seq)
        self.assertEqual(len(sites),3)
        self.assertIn('unknown_or_invalid_rsa',excluded[0]['reasons'])

    def test_unknown_bond_role_is_not_structural(self):
        seq,rows,_=fixture(1);rows[0]['bond_role']='unclassified'
        _,excluded=d.eligible_sites(rows,seq)
        self.assertIn('bond_role_not_verified_structural',excluded[0]['reasons'])

    def test_invalid_partner_and_conflicting_states_rejected(self):
        for change in ({'bond_partner':1},{'bond_state':'free','bond_partner':''}):
            seq,rows,_=fixture(1);rows[0].update(change)
            with self.assertRaises(ValueError): d.eligible_sites(rows,seq)

    def test_duplicate_sites_and_family_role_leakage_rejected(self):
        seq,rows,_=fixture(1)
        with self.assertRaises(ValueError): d.eligible_sites(rows+[dict(rows[0])],seq)
        rows[2]['role']='calibration'
        with self.assertRaises(ValueError): d.eligible_sites(rows,seq)

    def test_exposed_confirmation_is_ineligible(self):
        seq,rows,_=fixture(1)
        for r in rows: r['exposure_status']='previously_exposed'
        sites,excluded=d.eligible_sites(rows,seq)
        self.assertFalse(sites)
        self.assertEqual(len(excluded),4)


class Matching(unittest.TestCase):
    def test_matching_ignores_input_order_and_never_reuses_controls(self):
        _,sites,expected,_=bank(3)
        random.Random(7).shuffle(sites)
        actual,_=d.match_sites(sites)
        self.assertEqual([p['pair_id'] for p in expected],[p['pair_id'] for p in actual])
        ids=[p['free']['site_id'] for p in actual]
        self.assertEqual(len(ids),len(set(ids)))

    def test_exact_context_and_rsa_caliper_restrict_matches(self):
        _,sites,_,_=bank(1)
        sites[2]['compartment']='cytosolic'
        sites[3]['rsa']=0.95
        pairs,unmatched=d.match_sites(sites)
        self.assertFalse(pairs)
        self.assertEqual(len(unmatched),4)

    def test_pairing_does_not_use_predictor_scores(self):
        _,sites,pairs,_=bank(2)
        self.assertEqual(len(pairs),4)
        self.assertFalse(any('score' in k for p in pairs for k in p))


class Estimation(unittest.TestCase):
    def test_identical_substitution_recovers_effect_despite_protein_baselines(self):
        _,_,pairs,scores=bank(9,offset=0.3)
        result=d.audit_pairs(pairs,scores,score_direction='higher-is-damaging',draws=299)
        self.assertAlmostEqual(result['estimate_raw'],0.3)
        self.assertEqual(result['verdict'],'supported_conditional_direction')

    def test_direction_is_explicit_and_flips_the_effect(self):
        _,_,pairs,scores=bank()
        with self.assertRaises(ValueError): d.audit_pairs(pairs,scores)
        result=d.audit_pairs(pairs,scores,score_direction='higher-is-tolerated',draws=99)
        self.assertAlmostEqual(result['estimate_raw'],-0.3)
        self.assertEqual(result['verdict'],'opposite_conditional_direction')

    def test_family_not_variant_count_is_the_replication_unit(self):
        _,_,pairs,scores=bank(6)
        for i,p in enumerate(pairs):
            scores[(p['bonded']['accession'],p['bonded']['position'],'S')]+=i/20
        one=d.audit_pairs(pairs,scores,score_direction='higher-is-damaging',draws=299)
        doubled=d.audit_pairs(pairs*20,scores,score_direction='higher-is-damaging',draws=299)
        self.assertEqual(one['n_families'],6)
        np.testing.assert_allclose(one['ci95_raw'],doubled['ci95_raw'])

    def test_many_pairs_in_one_family_cannot_pass_the_gate(self):
        _,_,pairs,scores=bank(6)
        for p in pairs: p['protein_family']='ONE'
        result=d.audit_pairs(pairs,scores,score_direction='higher-is-damaging',draws=99)
        self.assertEqual(result['verdict'],'blocked_insufficient_evidence')
        self.assertIsNone(result['estimate_raw'])

    def test_missing_paired_scores_are_counted_and_cannot_be_imputed(self):
        _,_,pairs,scores=bank(6)
        scores={k:v for k,v in scores.items() if k[0]!='P0'}
        result=d.audit_pairs(pairs,scores,score_direction='higher-is-damaging',draws=99)
        self.assertEqual(result['n_families'],5)
        self.assertEqual(len(result['missing_pair_ids']),2)
        self.assertEqual(result['verdict'],'blocked_insufficient_evidence')

    def test_null_is_indeterminate(self):
        _,_,pairs,scores=bank(6,offset=0)
        result=d.audit_pairs(pairs,scores,score_direction='higher-is-damaging',draws=99)
        self.assertEqual(result['verdict'],'indeterminate')
        self.assertEqual(result['ci95_raw'],[0,0])

    def test_roundoff_does_not_manufacture_a_directional_verdict(self):
        _,_,pairs,scores=bank(6)
        scores={k:(1.0+1e-15 if k[1] in (2,4) else 1.0) for k in scores}
        result=d.audit_pairs(pairs,scores,score_direction='higher-is-damaging',draws=99)
        self.assertEqual(result['verdict'],'indeterminate')

    def test_numerical_tolerance_scales_with_the_score_units(self):
        _,_,pairs,scores=bank(6)
        scores={k:(1.1e-15 if k[1] in (2,4) else 1e-15) for k in scores}
        result=d.audit_pairs(pairs,scores,score_direction='higher-is-damaging',draws=99)
        self.assertEqual(result['verdict'],'supported_conditional_direction')


class Corrections(unittest.TestCase):
    def test_already_passed_rule_gets_zero_correction(self):
        _,_,pairs,scores=bank(3,'calibration',offset=0.4)
        fitted=d.fit_correction(pairs,scores,experiment(),score_direction='higher-is-damaging',
                                target_margin=0.1,draws=99)
        self.assertEqual(fitted['beta_oriented'],0)

    def test_calibration_only_and_no_confirmation_outcome_gating(self):
        _,_,pairs,scores=bank(3,'confirmation',offset=-0.1)
        with self.assertRaises(ValueError):
            d.fit_correction(pairs,scores,experiment(),score_direction='higher-is-damaging',target_margin=0.1)
        _,_,pairs,scores=bank(3,'calibration',offset=-0.1)
        validation=experiment();validation['partition']='confirmation'
        with self.assertRaises(ValueError):
            d.fit_correction(pairs,scores,validation,score_direction='higher-is-damaging',target_margin=0.1)
        validation=experiment();validation['assay_context_reviewed']=False
        with self.assertRaises(ValueError):
            d.fit_correction(pairs,scores,validation,score_direction='higher-is-damaging',target_margin=0.1)

    def test_correction_changes_only_verified_scope_and_respects_direction(self):
        _,sites,pairs,scores=bank(3,'calibration',offset=-0.1)
        fitted=d.fit_correction(pairs,scores,experiment(),score_direction='higher-is-damaging',
                                target_margin=0.1,draws=99)
        self.assertAlmostEqual(fitted['beta_oriented'],0.2)
        variants=[{'accession':acc,'position':p,'reference':'C','alternate':a,'score':v}
                  for (acc,p,a),v in scores.items()]
        corrected=d.apply_correction(variants,sites,fitted,'higher-is-damaging')
        changed=0
        for a,b in zip(variants,corrected):
            is_target=a['position'] in (2,4) and a['alternate']=='S'
            self.assertAlmostEqual(b['score']-a['score'],0.2 if is_target else 0)
            changed+=is_target
        self.assertEqual(changed,6)
        lower=d.apply_correction(variants,sites,fitted,'higher-is-tolerated')
        self.assertAlmostEqual(lower[0]['score']-variants[0]['score'],-0.2)
        unverified=deepcopy(sites);unverified[0].pop('_eligibility_verified')
        with self.assertRaises(ValueError): d.apply_correction(variants,unverified,fitted,'higher-is-damaging')

    def test_disulfide_shift_can_change_the_existing_charge_comparator(self):
        # Validate jointly because shifting an uncharged comparator changes the K/R penalty even when K/R scores are unchanged.
        matrix=np.zeros((1,20),dtype=float)
        matrix[0,panel.ALPHABET.index('C')]=np.nan
        old=panel.gain_contrast({'scores':matrix,'ref_aa':np.array([1])})['penalty'][0]
        matrix[0,panel.ALPHABET.index('S')]=0.3
        new=panel.gain_contrast({'scores':matrix,'ref_aa':np.array([1])})['penalty'][0]
        self.assertLess(new,old)


class FileInterfaces(unittest.TestCase):
    def write_bank(self,root,n=6):
        seq,rows,scores=fixture(n)
        sp=root/'sites.csv';qp=root/'sequences.json';yp=root/'scores.csv'
        with sp.open('w',newline='',encoding='utf-8') as fh:
            writer=csv.DictWriter(fh,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        qp.write_text(json.dumps(seq),encoding='utf-8')
        with yp.open('w',newline='',encoding='utf-8') as fh:
            writer=csv.DictWriter(fh,fieldnames=['accession','position','reference','alternate','score'])
            writer.writeheader();writer.writerows({'accession':a,'position':p,'reference':'C','alternate':alt,'score':v}
                for (a,p,alt),v in scores.items())
        return sp,qp,yp

    def test_blocked_inventory_never_opens_scores(self):
        with test_directory() as td:
            root=Path(td).resolve()
            self.assertTrue(root.is_relative_to(ROOT.resolve()))
            sp,qp,yp=self.write_bank(root,n=1)
            with patch.object(d,'read_scores',side_effect=AssertionError('Outcomes opened')):
                result=d.run_disulfide_audit(sp,qp,root/'does_not_exist.csv',root/'out',
                                            score_direction='higher-is-damaging')
            self.assertEqual(result['status'],'blocked_insufficient_evidence')
            self.assertFalse(result['score_files_opened'])

    def test_complete_known_answer_file_path_audits_both_substitutions(self):
        with test_directory() as td:
            root=Path(td).resolve()
            self.assertTrue(root.is_relative_to(ROOT.resolve()))
            sp,qp,yp=self.write_bank(root)
            result=d.run_disulfide_audit(sp,qp,yp,root/'out',score_direction='higher-is-damaging',draws=99)
            self.assertEqual(result['status'],'audited')
            for a in ('S','A'):
                self.assertAlmostEqual(result['results']['C>'+a]['estimate_raw'],0.3)
            self.assertFalse(result['new_model_failure_demonstrated'])

    def test_conflicting_score_duplicates_are_errors(self):
        with test_directory() as td:
            root=Path(td).resolve()
            self.assertTrue(root.is_relative_to(ROOT.resolve()))
            _,_,yp=self.write_bank(root)
            with yp.open('a',encoding='utf-8') as fh: fh.write('P0,2,C,S,999\n')
            with self.assertRaises(ValueError): d.read_scores(yp,{('P0',2)})


if __name__=='__main__':
    unittest.main()

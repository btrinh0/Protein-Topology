"""Test context-specific proline contrasts and their confound controls."""
from copy import deepcopy
from contextlib import contextmanager
import json
from pathlib import Path
import sys
import shutil
import unittest
import uuid
from unittest.mock import patch
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from membrane_rule_audit import proline as p


@contextmanager
def workspace_case():
    parent=ROOT/'results/proline_rule_20261025/software_tests'
    root=parent/f'test-{uuid.uuid4().hex}';root.mkdir(parents=True)
    try: yield root
    finally:
        if root.resolve().parent!=parent.resolve(): raise RuntimeError('Unsafe test cleanup path')
        shutil.rmtree(root)


def fixture(n=8,delta=0.3,chemistry=True):
    sites=[];scores=[]
    chemistry=np.linspace(-0.2,0.2,20) if chemistry else np.zeros(20)
    for i in range(n):
        for j,context in enumerate(['tm_helix','coil','non_tm_helix']):
            for k,native in enumerate('AL'):
                site={'accession':f'P{i}','position':j*10+k+1,'native':native,
                      'context':context,'packing_bin':0,'family':f'F{i}'}
                values=chemistry.copy()+10*i+3*j
                values[p.ALPHABET.index('P')]+=0.5+(delta if context=='tm_helix' else 0)
                values[p.ALPHABET.index(native)]=np.nan
                sites.append(site);scores.append(values)
    return sites,np.asarray(scores)


class NativeGeometry(unittest.TestCase):
    def test_signed_dihedral_and_degenerate_geometry(self):
        self.assertAlmostEqual(p.dihedral([1,0,0],[0,0,0],[0,1,0],[0,1,1]),-90)
        self.assertAlmostEqual(p.dihedral([1,0,0],[0,0,0],[0,1,0],[0,1,-1]),90)
        self.assertIsNone(p.dihedral([0,0,0],[0,0,0],[0,0,0],[1,0,0]))

    def test_local_axis_distinguishes_regular_and_bent_helices(self):
        t=np.arange(9)*np.radians(100)
        ca=np.column_stack([2.3*np.cos(t),2.3*np.sin(t),1.5*np.arange(9)])
        self.assertLess(p.axis_bend(ca),20)
        rotation=np.array([[0,0,1],[0,1,0],[-1,0,0]])
        bent=ca.copy();bent[5:]=(ca[5:]-ca[4])@rotation.T+ca[4]
        self.assertGreater(p.axis_bend(bent),30)

    def atom_fixture(self):
        return {i:{'CA':np.array([0.,0.,i*1.5]),'N':np.array([1.,0.,i*1.5]),
                   'C':np.array([0.,1.,i*1.5])} for i in range(1,41)}

    def test_native_proline_and_neighbor_proline_are_excluded(self):
        atoms=self.atom_fixture();sequence='A'*40
        sequence=sequence[:14]+'P'+sequence[15:]
        site,reason=p.classify_geometry(sequence,15,atoms,[(8,30)],set(),set(),{})
        self.assertIsNone(site);self.assertEqual(reason,'native_not_eligible')
        site,reason=p.classify_geometry(sequence,12,atoms,[(8,30)],set(),set(),{})
        self.assertIsNone(site);self.assertEqual(reason,'native_pro_or_gly_nearby')

    def test_regular_internal_helix_is_eligible_but_ends_are_not(self):
        with patch.object(p,'dihedral',side_effect=[-60.,-40.]*5):
            site,reason=p.classify_geometry('A'*40,15,self.atom_fixture(),[(8,30)],set(),set(),{})
        self.assertEqual(site['context'],'tm_helix')
        with patch.object(p,'dihedral',side_effect=[-60.,-40.]*5):
            site,reason=p.classify_geometry('A'*40,9,self.atom_fixture(),[(8,30)],set(),set(),{})
        self.assertEqual(reason,'tm_helix_end_or_short_segment')

    def test_unresolved_native_window_is_not_guessed_regular(self):
        atoms=self.atom_fixture();atoms.pop(13)
        site,reason=p.classify_geometry('A'*40,15,atoms,[(8,30)],set(),set(),{})
        self.assertIsNone(site);self.assertEqual(reason,'unresolved_or_non_native_window')

    def test_local_coil_does_not_require_a_helix_bend_window(self):
        atoms=self.atom_fixture();atoms.pop(31)
        with patch.object(p,'dihedral',side_effect=[-65.,120.]*5):
            site,reason=p.classify_geometry('A'*40,35,atoms,[(8,30)],set(),set(),{})
        self.assertEqual(site['context'],'coil')
        self.assertIsNone(site['local_axis_bend_degrees'])


class WithinPosition(unittest.TestCase):
    def test_planted_target_effect_survives_additive_site_and_protein_baselines(self):
        sites,scores=fixture()
        contrasts,coverage=p.position_contrasts(scores,[r['native'] for r in sites])
        effect=p.context_effect(sites,contrasts['aliphatic']['raw'],draws=199)
        self.assertAlmostEqual(effect['estimate'],0.3)
        self.assertEqual(effect['verdict'],'supported_conditional_direction')

    def test_native_synonymous_value_never_enters_a_comparator(self):
        sites,scores=fixture()
        old,_=p.position_contrasts(scores,[r['native'] for r in sites])
        for i,s in enumerate(sites): scores[i,p.ALPHABET.index(s['native'])]=999
        new,_=p.position_contrasts(scores,[r['native'] for r in sites])
        np.testing.assert_allclose(old['aliphatic']['raw'],new['aliphatic']['raw'])

    def test_missing_alternative_removes_same_position_from_all_comparators(self):
        sites,scores=fixture();scores[0,p.ALPHABET.index('S')]=np.nan
        contrasts,coverage=p.position_contrasts(scores,[r['native'] for r in sites])
        self.assertEqual(coverage['n_missing_or_incomplete'],1)
        for arm in contrasts.values(): self.assertTrue(np.isnan(arm['raw'][0]))

    def test_native_pro_loss_is_not_a_pro_gain(self):
        with self.assertRaises(ValueError): p.position_contrasts(np.zeros((1,20)),['P'])

    def test_duplicate_canonical_sites_are_rejected(self):
        sites,scores=fixture()
        with self.assertRaises(ValueError):
            p.audit_scores(sites+[sites[0]],np.vstack([scores,scores[0]]),score_direction='higher-is-damaging')

    def test_constant_rows_have_no_scale_diagnostic(self):
        sites,scores=fixture();scores[:]=0.1
        contrasts,coverage=p.position_contrasts(scores,[r['native'] for r in sites])
        self.assertEqual(coverage['n_complete_position_scaled'],0)
        self.assertTrue(np.isnan(contrasts['aliphatic']['position_scaled']).all())


class ConfoundsAndUncertainty(unittest.TestCase):
    def test_native_composition_imbalance_does_not_manufacture_context_effect(self):
        sites=[];values=[]
        for i in range(8):
            for native in 'AL':
                for context in ('tm_helix','coil'):
                    n=10 if (native=='L')==(context=='tm_helix') else 1
                    for j in range(n):
                        sites.append({'accession':f'P{i}','family':f'F{i}','native':native,
                                      'packing_bin':0,'context':context})
                        values.append(1.0 if native=='L' else 0.0)
        raw=np.mean([v for s,v in zip(sites,values) if s['context']=='tm_helix'])-np.mean([v for s,v in zip(sites,values) if s['context']=='coil'])
        self.assertGreater(raw,0.5)
        result=p.context_effect(sites,values,draws=199)
        self.assertEqual(result['estimate'],0)
        self.assertEqual(result['verdict'],'indeterminate')

    def test_no_shared_packing_stratum_produces_no_estimate(self):
        sites,scores=fixture()
        for site in sites:
            if site['context']=='coil': site['packing_bin']=2
        values=p.position_contrasts(scores,[r['native'] for r in sites])[0]['aliphatic']['raw']
        result=p.context_effect(sites,values,draws=99)
        self.assertIsNone(result['estimate'])

    def test_homologous_proteins_do_not_count_as_independent_groups(self):
        sites,scores=fixture()
        for site in sites: site['family']='ONE_FAMILY'
        result=p.audit_scores(sites,scores,score_direction='higher-is-damaging',draws=99)
        self.assertEqual(result['results']['aliphatic']['raw']['family']['n_groups'],1)
        self.assertEqual(result['overall_status'],'mixed_or_insufficient_evidence')

    def test_score_orientation_is_required_and_reverses_the_contrast(self):
        sites,scores=fixture()
        with self.assertRaises(ValueError): p.audit_scores(sites,scores,score_direction='unspecified')
        result=p.audit_scores(sites,scores,score_direction='higher-is-tolerated',draws=99)
        self.assertAlmostEqual(result['results']['aliphatic']['raw']['family']['estimate'],-0.3)

    def test_scale_diagnostic_removes_a_pure_multiplicative_context_effect(self):
        sites,scores=fixture()
        base=np.linspace(0,1,20)
        base[p.ALPHABET.index('P')]=1.5
        for i,site in enumerate(sites):
            scores[i]=base*(2 if site['context']=='tm_helix' else 1)
        result=p.audit_scores(sites,scores,score_direction='higher-is-damaging',draws=199)
        self.assertGreater(result['results']['aliphatic']['raw']['family']['estimate'],0)
        self.assertAlmostEqual(result['results']['aliphatic']['position_scaled']['family']['estimate'],0)
        self.assertFalse(result['new_biological_failure_demonstrated'])

    def test_normalization_can_also_remove_a_real_isolated_pro_effect(self):
        sites,scores=fixture(chemistry=False)
        result=p.audit_scores(sites,scores,score_direction='higher-is-damaging',draws=99)
        self.assertAlmostEqual(result['results']['aliphatic']['raw']['family']['estimate'],0.3)
        self.assertAlmostEqual(result['results']['aliphatic']['position_scaled']['family']['estimate'],0)
        # A failed scale diagnostic alone does not establish a biological-rule violation.
        self.assertFalse(result['new_biological_failure_demonstrated'])

    def test_non_tm_helix_comparison_has_no_invented_directional_verdict(self):
        sites,scores=fixture()
        result=p.audit_scores(sites,scores,score_direction='higher-is-damaging',draws=99)
        diag=result['results']['aliphatic']['raw']['tm_minus_non_tm_helix_diagnostic']
        self.assertFalse(diag['direction_prespecified'])
        self.assertEqual(diag['verdict'],'descriptive_no_direction_prespecified')


class RankOrientation(unittest.TestCase):
    def test_tied_measured_scores_get_average_ranks(self):
        np.testing.assert_allclose(p.average_tie_ranks([4,1,1,2]),[4,1.5,1.5,3])

    def test_nonfinite_measured_outcomes_are_not_ranked(self):
        with self.assertRaises(ValueError): p.average_tie_ranks([1,np.nan])


class ScoreArchive(unittest.TestCase):
    def archive(self,root,sites,scores):
        path=root/'scores.npz'
        np.savez(path,scores=scores,accession=np.array([s['accession'] for s in sites]),
                 position=np.array([s['position'] for s in sites]),
                 native=np.array([s['native'] for s in sites]),alphabet=np.array(list(p.ALPHABET)))
        return path

    def test_named_rows_align_independently_of_file_order(self):
        sites,scores=fixture()
        with workspace_case() as root:
            path=self.archive(root,sites[::-1],scores[::-1])
            np.testing.assert_allclose(p.load_score_archive(path,sites),scores,equal_nan=True)

    def test_native_mismatch_and_duplicate_archive_positions_are_rejected(self):
        sites,scores=fixture()
        with workspace_case() as root:
            bad=deepcopy(sites);bad[0]['native']='C'
            path=self.archive(root,bad,scores)
            with self.assertRaises(ValueError): p.load_score_archive(path,sites)
            path=self.archive(root,sites+[sites[0]],np.vstack([scores,scores[0]]))
            with self.assertRaises(ValueError): p.load_score_archive(path,sites)

    def test_missing_archive_site_stays_unmeasured(self):
        sites,scores=fixture()
        with workspace_case() as root:
            path=self.archive(root,sites[1:],scores[1:])
            result=p.load_score_archive(path,sites)
            self.assertTrue(np.isnan(result[0]).all())

    def test_public_file_interface_recovers_known_answer(self):
        sites,scores=fixture()
        with workspace_case() as root:
            path=self.archive(root,sites,scores);source=root/'sites.json'
            source.write_text(json.dumps(sites),encoding='utf-8')
            result=p.run_proline_audit(source,path,root/'out',score_direction='higher-is-damaging',draws=99)
            self.assertAlmostEqual(result['results']['aliphatic']['raw']['family']['estimate'],0.3)
            self.assertTrue((root/'out/proline_audit.json').exists())


if __name__=='__main__': unittest.main()

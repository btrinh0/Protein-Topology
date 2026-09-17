"""Test grouped exact-tail calculations and multiple-comparison handling."""
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts/beat_am'))
from m20_dms_comparator_pvalues import paired_test,holm


class TestPairedPvalues(unittest.TestCase):
    def test_exact_floor_uses_groups_not_protein_count(self):
        result=paired_test([.1]*13,['a','a']+list('bcdefghijkl'))
        self.assertEqual(result['proteins'],13)
        self.assertEqual(result['dependence_groups'],12)
        self.assertEqual(result['sign_patterns'],4096)
        self.assertEqual(result['p_two_sided'],2/4096)
        np.testing.assert_allclose(result['ci95_delta_spearman'],[.1,.1])

    def test_unequal_groups_keep_equal_protein_statistic(self):
        result=paired_test([.1,.1,-.3],['a','a','b'])
        self.assertAlmostEqual(result['delta_spearman'],-.1/3)
        self.assertEqual(result['sign_patterns'],4)
        self.assertEqual(result['p_two_sided'],1.)

    def test_zero_differences_have_p_one(self):
        result=paired_test(np.zeros(12),np.arange(12))
        self.assertEqual(result['p_two_sided'],1.)
        self.assertEqual(result['ci95_delta_spearman'],[0.,0.])

    def test_two_sided_test_is_direction_symmetric(self):
        a=paired_test([.2,-.1,.4,.05],list('abcd'))
        b=paired_test([-.2,.1,-.4,-.05],list('abcd'))
        self.assertEqual(a['p_two_sided'],b['p_two_sided'])
        self.assertAlmostEqual(a['delta_spearman'],-b['delta_spearman'])

    def test_holm_monotonicity_original_order_and_ties(self):
        np.testing.assert_allclose(holm([.01,.04,.03]),[.03,.06,.06])
        np.testing.assert_allclose(holm([.01,.01,.3]),[.03,.03,.3])
        np.testing.assert_allclose(holm([1.,1.,1.]),[1.,1.,1.])

    def test_invalid_inputs_fail(self):
        for d,g in [([1,np.nan],['a','b']),([1,2],['a']),([],[])]:
            with self.assertRaises(ValueError):paired_test(d,g)
        with self.assertRaises(ValueError):holm([-.1,.2])


if __name__=='__main__':unittest.main()

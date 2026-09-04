"""Test label-free audit adaptation with synthetic predictions and annotations."""
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts/beat_am"))
import m7_audit_standalone_rules as audit


class TestStandaloneAudits(unittest.TestCase):
    def test_estimator_point_preserved_while_grouping_bootstrap_units(self):
        rows=[]
        for a in ("a","b","c"):
            rows += [(a,"A",1,.6),(a,"A",1,.8),(a,"A",0,.1),(a,"A",0,.3)]
        original=audit.b8.stratified_delta(rows,np.random.default_rng(1))
        grouped=audit.matched_cluster(rows,{"a":"g1","b":"g1","c":"g2"})
        self.assertAlmostEqual(original["estimate"],.5)
        self.assertAlmostEqual(grouped["estimate"],.5)
        self.assertEqual(grouped["n_groups"],2)
        self.assertEqual(grouped["n_proteins"],3)

    def test_no_controls_is_not_a_pass(self):
        out=audit.matched_cluster([("a","C",1,.8)],{"a":"g"})
        self.assertIsNone(out)
        self.assertEqual(audit.matched_verdict(out,1),"insufficient_data")

    def test_uniform_bond_loss_offset_cannot_change_pairing_statistic(self):
        sites=dict(acc=np.array(["a"]*4+["b"]*4),pos=np.array([1,2,3,4]*2),
            bonds_by_prot={a:[(1,2,True),(3,4,True)] for a in ("a","b")})
        v=np.array([.05,.1,.8,.9,.15,.2,.7,.8])
        raw=audit.b10.pairing_test(sites,v,np.random.default_rng(7))
        shifted=audit.b10.pairing_test(sites,v+.2,np.random.default_rng(7))
        for key in ("observed_mean_within_bond_abs_diff","null_mean","difference_observed_minus_null_mean"):
            self.assertAlmostEqual(raw[key],shifted[key],places=12)

    def test_depth_adaptation_has_known_correct_endpoint(self):
        n=30;prot=np.repeat([0,1],15);group=np.tile(np.repeat(np.arange(5),3),2)
        matrix=np.zeros((n,20));matrix[:,0]=np.nan
        for aa in "DEKR": matrix[:,audit.ALPHABET.index(aa)]=group*.1
        sites=dict(nat=np.zeros(n,dtype=int),prot=prot,grp=group,prot_names=["a","b"])
        result=audit.charge_groups(matrix,sites,{"a":"g1","b":"g2"})
        self.assertAlmostEqual(result["estimate"],.4)

    def test_training_gate_is_hard_failure(self):
        with self.assertRaises(AssertionError): audit.forbidden_fit()


if __name__=="__main__": unittest.main()

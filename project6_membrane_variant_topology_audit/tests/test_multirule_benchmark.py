"""Test dependence-group splits and paired benchmark metrics with synthetic inputs."""
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts/beat_am"))
from m4_benchmark_multirule import conservative_groups,fold_plans,paired_summary,protein_metrics,ABLATIONS


class TestSplits(unittest.TestCase):
    def test_union_retains_links_from_both_old_groupings(self):
        acc = np.array(["a","b","c","d"])
        groups = conservative_groups(acc,["x","x","y","z"],["p","q","q","r"])
        self.assertEqual(len(set(groups[:3])),1)
        self.assertNotEqual(groups[0],groups[3])

    def test_guarded_loocv_holds_relatives_out_together(self):
        groups = np.repeat(np.array(list("abcdefghij")),[2,1,1,1,1,1,1,1,1,1])
        plans = fold_plans(groups,"guarded_loocv")
        self.assertEqual(len(plans),10)
        self.assertEqual(len(plans[0]["test"]),2)
        for p in plans:
            self.assertFalse(set(groups[p["train"]]) & set(groups[p["test"]]))
            for tr,cal in p["inner"]:
                self.assertFalse(set(groups[tr]) & set(groups[np.r_[cal,p["test"]]]))

    def test_fivefold_once_per_variant(self):
        groups = np.repeat(np.array(list("abcdefghijkl")),3)
        plans = fold_plans(groups,"grouped_cv")
        self.assertEqual(len(plans),5)
        np.testing.assert_array_equal(np.sort(np.concatenate([p["test"] for p in plans])),np.arange(len(groups)))


class TestMetrics(unittest.TestCase):
    def records(self):
        rows = []
        for a,g,mem,b,r in [("a","g1",1,.5,.6),("b","g1",1,.4,.5),("c","g2",1,.3,.4),("d","g3",0,.8,.8)]:
            for name,rho in [("raw",b),("refined",r)]:
                rows.append(dict(accession=a,group=g,membrane=mem,model=name,spearman=rho,mse=1-rho))
        return rows

    def test_paired_mean_and_group_count_not_variant_replication(self):
        out = paired_summary(self.records(),"refined","raw",draws=100)
        self.assertAlmostEqual(out["delta_spearman"],.1)
        self.assertEqual((out["proteins"],out["groups"]),(3,2))
        self.assertEqual((out["wins"],out["losses"],out["ties"]),(3,0,0))
        np.testing.assert_allclose(out["ci95_delta_spearman"],[.1,.1])

    def test_same_predictions_give_exact_zero_and_negative_scope(self):
        out = paired_summary(self.records(),"raw","raw",draws=100)
        self.assertEqual(out["delta_spearman"],0)
        self.assertEqual(out["ci95_delta_spearman"],[0,0])
        out = paired_summary(self.records(),"refined","raw",scope="non_membrane",draws=100)
        self.assertEqual(out["delta_spearman"],0)

    def test_assays_receive_equal_weight_within_protein(self):
        y = np.r_[np.arange(20),np.arange(40)].astype(float)
        pred = np.r_[np.arange(20),-np.arange(40)]
        matrix = np.repeat(pred[:,None],len(ABLATIONS),axis=1)
        data = dict(accession=np.repeat("a",60),dms_id=np.repeat(["small","large"],[20,40]),
                    y=y,gene=np.repeat("gene",60),groups=np.repeat("g",60),is_membrane=np.ones(60))
        rows = protein_metrics(matrix,data)
        self.assertTrue(all(abs(row["spearman"])<1e-12 for row in rows))


if __name__ == "__main__":
    unittest.main()

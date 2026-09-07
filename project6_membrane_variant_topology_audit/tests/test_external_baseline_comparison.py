"""Test mutation joins and matched external-comparator aggregation."""
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts/beat_am"))
from m11_compare_external_baselines import (
    assay_metrics, common_mask, map_mutant, paired_summary, protein_metrics,
)


class TestExternalComparison(unittest.TestCase):
    def test_construct_insertion_is_not_a_global_offset(self):
        mapping = np.array([0,1,0,2])
        self.assertEqual(map_mutant("A3V","MGA","MA",mapping),(2,"A","V"))
        self.assertEqual(map_mutant("M1L","MGA","MA",mapping),(1,"M","L"))
        self.assertIsNone(map_mutant("G2A","MGA","MA",mapping))

    def test_native_mismatch_and_complex_mutants_are_rejected(self):
        mapping = np.array([0,1,0,2])
        for token in ("Y1A", "A3V:M1L", "A3A", "A99V", "A2V"):
            self.assertIsNone(map_mutant(token,"MGA","MA",mapping))

    def test_shared_rows_are_intersection_not_individual_coverage(self):
        values = np.array([[1,2,3],[2,np.nan,4],[3,4,np.nan],[4,5,6]])
        mask = common_mask(values,["raw","GEMME","VESPA"],["GEMME","VESPA"])
        np.testing.assert_array_equal(mask,[True,False,False,True])

    def test_macro_average_is_equal_proteins_after_equal_assays(self):
        rows = [dict(accession="A",model="raw",gene="A",membrane=1,group="g1",
                     rows=10000,spearman=.1,rank_mse=.1,coarse_selection_type="Activity"),
                dict(accession="A",model="raw",gene="A",membrane=1,group="g1",
                     rows=20,spearman=.3,rank_mse=.1,coarse_selection_type="Binding"),
                dict(accession="B",model="raw",gene="B",membrane=1,group="g2",
                     rows=20,spearman=.9,rank_mse=.1,coarse_selection_type="Activity")]
        proteins = protein_metrics(rows)
        self.assertAlmostEqual(proteins[0]["spearman"],.2)
        self.assertAlmostEqual(np.mean([r["spearman"] for r in proteins]),.55)
        activity = protein_metrics(rows,"Activity")
        self.assertAlmostEqual(np.mean([r["spearman"] for r in activity]),.5)

    def test_grouped_bootstrap_keeps_related_proteins_together(self):
        rows = []
        for acc,group,delta in (("A","g1",.2),("B","g1",.2),("C","g2",-.1)):
            for model,rho in (("raw",.3),("refined",.3+delta)):
                rows.append(dict(accession=acc,group=group,model=model,spearman=rho,rank_mse=.1))
        result = paired_summary(rows,"refined","raw",draws=1000)
        self.assertEqual(result["groups"],2)
        self.assertEqual(result["proteins"],3)
        self.assertAlmostEqual(result["delta_spearman"],.1)
        np.testing.assert_allclose(result["ci95_delta_spearman"],[-.1,.2],atol=1e-12)
        self.assertEqual((result["wins"],result["losses"]),(2,1))

    def test_paired_comparison_rejects_different_protein_sets(self):
        rows = [dict(accession="A",group="g1",model="raw",spearman=.3,rank_mse=.1),
                dict(accession="B",group="g2",model="refined",spearman=.4,rank_mse=.1)]
        with self.assertRaises(ValueError):
            paired_summary(rows,"refined","raw")

    def test_fitness_to_damage_conversion_has_the_expected_sign(self):
        n = 25
        table = pd.DataFrame(dict(dms_id=["test"]*n,target=np.arange(n,dtype=float),group=["g"]*n))
        metadata = {"test":dict(assay="test",accession="A",gene="A",membrane=1,coarse_selection_type="Activity")}
        fitness = -np.arange(n,dtype=float)
        values = np.column_stack([-fitness,fitness])
        rows, excluded = assay_metrics(table,values,["damage","unconverted_fitness"],np.ones(n,dtype=bool),metadata)
        self.assertFalse(excluded)
        self.assertAlmostEqual(rows[0]["spearman"],1)
        self.assertAlmostEqual(rows[1]["spearman"],-1)
        self.assertAlmostEqual(rows[0]["rank_mse"],0)

    def test_constant_prediction_excludes_assay_for_both_models(self):
        n = 25
        table = pd.DataFrame(dict(dms_id=["test"]*n,target=np.arange(n,dtype=float),group=["g"]*n))
        values = np.column_stack([np.arange(n),np.zeros(n)])
        records, excluded = assay_metrics(table,values,["raw","constant"],np.ones(n,dtype=bool),{})
        self.assertEqual(records,[])
        self.assertIn("constant",excluded[0]["reason"])


if __name__ == "__main__":
    unittest.main()

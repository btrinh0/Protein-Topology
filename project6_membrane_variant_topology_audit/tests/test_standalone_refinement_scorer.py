"""Test standalone feature construction and inference wiring with synthetic heads."""
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts/beat_am"))
sys.path.insert(0,str(ROOT/"tests"))
from test_multirule_refinement import prepared
from m2_standalone_refinement_scorer import (FEATURES,build_tabular,score_positions,
    same_model_probe,refine_position_scores)
from membrane_rule_audit.refinement import ALPHABET,RULES


class FakePCA:
    def transform(self, values):
        return np.zeros((len(values),64),dtype=np.float32)


class FakeHead:
    def predict(self, values):
        self.inputs = values
        return .1 + values[:,FEATURES.index("side")]*.01


def fake_bundle():
    return dict(features=FEATURES.copy(),pca=FakePCA(),model=FakeHead())


class TestStandaloneWiring(unittest.TestCase):
    def test_feature_codes_for_both_sides_and_TM(self):
        x = build_tabular(prepared(),[1,25,40],["A"]*3,["K"]*3)
        np.testing.assert_array_equal(x[:,FEATURES.index("side")],[2,0,1])
        np.testing.assert_array_equal(x[:,FEATURES.index("in_tm")],[0,1,0])
        self.assertEqual(x.shape,(3,18))
        self.assertNotIn("am",FEATURES); self.assertNotIn("esm1b",FEATURES)

    def test_identity_and_fractional_coordinate_rejected(self):
        for p,n,a in ((1,"A","A"),(1.5,"A","C"),(8,"A","C")):
            with self.assertRaises(ValueError): build_tabular(prepared(),[p],[n],[a])

    def test_fake_scoring_has_19_missense_rows_and_native_nan(self):
        bundle = fake_bundle()
        scored = score_positions(bundle,prepared(),[1,25,40],np.arange(1,46),np.zeros((45,1024)))
        self.assertEqual(bundle["model"].inputs.shape,(57,82))
        self.assertTrue(np.isnan(scored["scores"][:,0]).all())
        np.testing.assert_allclose(scored["scores"][:,ALPHABET.index("K")],[.12,.1,.11])

    def test_missing_embeddings_not_zero_filled(self):
        with self.assertRaises(ValueError):
            score_positions(fake_bundle(),prepared(),[1,25],[1],np.zeros((1,1024)))

    def test_scoring_rejects_fractional_positions(self):
        with self.assertRaises(ValueError):
            score_positions(fake_bundle(),prepared(),[1.5],[1],np.zeros((1,1024)))

    def test_probe_uses_supplied_baseline(self):
        raw = np.array([.2,.7]); features = np.zeros((2,4)); features[:,0] = [-.5,.5]
        np.testing.assert_allclose(same_model_probe(raw,features,"positive_inside",.1),[.15,.75])

    def test_complete_four_rule_wiring_and_missing_scores(self):
        rec = prepared(); scored = score_positions(fake_bundle(),rec,[4,8,25,40],np.arange(1,46),np.zeros((45,1024)))
        scored["scores"][0,ALPHABET.index("W")] = np.nan
        out = refine_position_scores(rec,scored,dict(zip(RULES,[.2,.3,.4,.5])))
        self.assertTrue(np.isnan(out["scores"][0,ALPHABET.index("W")]))
        self.assertAlmostEqual(out["scores"][1,ALPHABET.index("S")]-scored["scores"][1,ALPHABET.index("S")],.5)
        self.assertEqual(out["contributions"].shape,(4,20,4))
        with self.assertRaises(ValueError): refine_position_scores(rec,out,{r:0 for r in RULES})

    def test_unmasked_identity_cannot_enter_refinement(self):
        rec = prepared(); scored = score_positions(fake_bundle(),rec,[1],np.arange(1,46),np.zeros((45,1024)))
        scored["scores"][0,0] = 0
        with self.assertRaises(ValueError): refine_position_scores(rec,scored,{r:0 for r in RULES})


if __name__ == "__main__":
    unittest.main()

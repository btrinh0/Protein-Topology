"""Test timing normalization and published-model adapters."""
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np
import torch
from scipy.special import expit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts/beat_am"))
import m13_local_speed_benchmark as speed


class TimingChecks(unittest.TestCase):
    def test_inner_repetitions_are_normalized_before_reporting_latency(self):
        result = speed.statistics([.010,.020,.030],rows=5,inner_loops=10)
        self.assertAlmostEqual(result["median_seconds"],.002)
        self.assertAlmostEqual(result["variants_per_second"],2500)
        self.assertAlmostEqual(result["min_seconds"],.001)

    def test_wrong_output_shape_is_rejected_before_timing(self):
        with self.assertRaises(ValueError):
            speed.time_methods({"bad":(lambda:np.zeros(3),1)},repeats=2,rows=2)

    def test_missing_predictions_are_rejected_before_timing(self):
        with self.assertRaises(ValueError):
            speed.time_methods({"bad":(lambda:np.array([np.nan]),1)},repeats=2,rows=1)


@unittest.skipUnless((speed.ASSETS/"cached_esm2/ACE2_HUMAN_Chan_2020.npz").exists(),"Published cached-input assets required")
class PublishedModelChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        torch.set_num_threads(1)
        cls.student,cls.conservation,cls.light,_ = speed.published_heads()

    def test_vespag_native_forward_matches_explicit_published_weight_formula(self):
        with np.load(speed.ASSETS/"cached_esm2/ACE2_HUMAN_Chan_2020.npz",allow_pickle=False) as archive:
            x = torch.from_numpy(archive["embedding"][:5].astype(np.float32))
        state = self.student.state_dict()
        with torch.inference_mode():
            hidden = torch.nn.functional.linear(x,state["net.0.weight"],state["net.0.bias"])
            hidden = torch.nn.functional.leaky_relu(hidden,negative_slope=.01)
            expected = torch.nn.functional.linear(hidden,state["net.2.weight"],state["net.2.bias"])
            actual = self.student(x)
        torch.testing.assert_close(actual,expected,atol=1e-6,rtol=1e-6)

    def test_vespal_real_conservation_features_and_legacy_lr_are_numerically_valid(self):
        with np.load(speed.EMB/"Q9BYF1.npz",allow_pickle=False) as archive:
            x = torch.from_numpy(archive["embedding"][:20].astype(np.float32)).unsqueeze(0)
        with torch.inference_mode():
            p = self.conservation.extract_probabilities(self.conservation(x))[0].T.numpy()
        np.testing.assert_allclose(p.sum(axis=1),np.ones(20),atol=2e-7)
        self.assertTrue(np.all((p>=0)&(p<=1)))
        inputs = np.column_stack([p,np.full(20,-1.0)])
        for classifier in self.light.values():
            actual = classifier.predict_proba(inputs)[:,1]
            expected = expit(inputs@classifier.coef_[0]+classifier.intercept_[0])
            np.testing.assert_allclose(actual,expected,atol=1e-12,rtol=0)

    def test_timing_preparation_never_reads_experimental_target_arrays(self):
        original = np.load
        reads = []
        class GuardedArchive:
            def __init__(self,archive):
                self.archive = archive
            def __enter__(self):
                return self
            def __exit__(self,*args):
                self.archive.close()
            def __getitem__(self,key):
                reads.append(key)
                if key == "target":
                    raise AssertionError("Timing tried to read experimental target")
                return self.archive[key]
        def guarded(*args,**kwargs):
            return GuardedArchive(original(*args,**kwargs))
        with patch.object(np,"load",side_effect=guarded):
            frame,*_ = speed.prepare()
        self.assertEqual(len(frame),43991)
        self.assertNotIn("target",reads)


if __name__ == "__main__":
    unittest.main()

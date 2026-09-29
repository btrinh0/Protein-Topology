"""Test directional glycosylation constraints with known-answer fixtures."""
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from membrane_rule_audit.glycosylation_constraints import (
    coefficients_with_glycosylation_floor, fit_glycosylation_floor,
)
from membrane_rule_audit.refinement import RULES


def row(accession, contrast, delta, slope, group="train"):
    return dict(accession=accession, group=group, contrast=contrast,
                baseline_delta=delta, feature_delta=slope)


class GlycosylationConstraintsTest(unittest.TestCase):
    def test_smallest_weight_solves_both_directions_and_nonunit_response(self):
        fit = fit_glycosylation_floor(
            [row("a", "A", -.03, .5), row("a", "B", .08, -1.)],
            calibration_groups=["train"], excluded_groups=["test"])
        self.assertAlmostEqual(fit["beta"], .08)
        self.assertAlmostEqual(fit["minimum_signed_calibration_slack"], 0.)
        self.assertEqual(fit["binding_constraints"], [{"accession": "a", "contrast": "B"}])

    def test_compliant_scores_require_no_change(self):
        fit = fit_glycosylation_floor(
            [row("a", "A", .1, 1.), row("a", "B", -.2, -1.)],
            calibration_groups=["train"])
        self.assertEqual(fit["beta"], 0.)

    def test_excluded_group_and_group_overlap_are_rejected(self):
        with self.assertRaises(ValueError):
            fit_glycosylation_floor([row("a", "A", -.1, 1., "test")],
                                   calibration_groups=["train"], excluded_groups=["test"])
        with self.assertRaises(ValueError):
            fit_glycosylation_floor([row("a", "A", -.1, 1.)],
                                   calibration_groups=["train"], excluded_groups=["train"])

    def test_unsupported_constraint_is_not_silently_certified(self):
        for slope in (0., -1.):
            with self.assertRaises(ValueError):
                fit_glycosylation_floor([row("a", "A", -.1, slope)],
                                       calibration_groups=["train"])

    def test_duplicate_protein_contrast_rejected(self):
        with self.assertRaises(ValueError):
            fit_glycosylation_floor([row("a", "A", -.1, 1.)] * 2,
                                   calibration_groups=["train"])

    def test_floor_only_changes_glycosylation_and_does_not_shrink_existing_weight(self):
        fit = fit_glycosylation_floor([row("a", "A", -.1, 1.)],
                                     calibration_groups=["train"])
        coefficients = dict(zip(RULES, [.02, .03, 0., .04]))
        result = coefficients_with_glycosylation_floor(coefficients, fit)
        self.assertEqual(result, dict(zip(RULES, [.02, .03, .1, .04])))
        self.assertEqual(coefficients["glycosylation_loss"], 0.)
        coefficients["glycosylation_loss"] = .2
        self.assertEqual(coefficients_with_glycosylation_floor(coefficients, fit)["glycosylation_loss"], .2)

    def test_matched_site_adapter_and_loss_feature_solve_the_known_contrasts(self):
        root = Path(__file__).resolve().parents[1]
        sys.path.insert(0, str(root / "scripts" / "beat_am"))
        import m33_glycosylation_constraint_repair as repair
        alphabet = repair.audit.ALPHABET
        positions, natives = [2, 3, 5, 6, 8, 9], ["A", "S", "A", "S", "A", "S"]
        scores = np.zeros((6, 20))
        feature = np.zeros_like(scores)
        for i, native in enumerate(natives):
            scores[i, alphabet.index(native)] = np.nan
        scores[0, alphabet.index("P")] = -.1
        feature[0, alphabet.index("P")] = 1.
        scores[1, alphabet.index("T")] = .2
        feature[1] = [float(a not in "ST") for a in alphabet]
        sites = [("a", p, n, "noncytosolic", "A" if n == "A" else "B", int(p < 4))
                 for p, n in zip(positions, natives)]
        index = {("a", p): i for i, p in enumerate(positions)}
        bank, groups = {"a": {"glycosites": [{"position": 1}]}}, {"a": "train"}
        constraints = repair.constraint_rows(scores, feature, sites, index, bank, groups)
        fit = fit_glycosylation_floor(constraints, calibration_groups=["train"])
        self.assertAlmostEqual(fit["beta"], .2)
        after = repair.constraint_rows(scores + fit["beta"] * feature, feature,
                                       sites, index, bank, groups)
        self.assertAlmostEqual(after[0]["baseline_delta"], .1)
        self.assertAlmostEqual(after[1]["baseline_delta"], 0.)
        self.assertEqual(after[0]["n_target"], 1)
        self.assertEqual(after[0]["n_reference"], 2)

    def test_training_floor_does_not_certify_a_more_extreme_held_out_protein(self):
        fit = fit_glycosylation_floor([row("a", "A", -.02, 1.)],
                                     calibration_groups=["train"], excluded_groups=["test"])
        held_out_delta = -.1 + fit["beta"]
        self.assertLess(held_out_delta, 0.)


if __name__ == "__main__":
    unittest.main()

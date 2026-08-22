"""Test placebo rules, disjoint comparison arms and shared-draw contrasts."""
from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from membrane_rule_audit import audit


def old_penalties(table, kind, group):
    """The arm construction as it stood before the group was held out."""
    out = []
    for key, alternates in table["sites"].items():
        reference = key[2]
        if kind == "gain":
            if reference in audit.CHARGED:
                continue
            numerator = [v for a, v in alternates.items() if a in group]
            comparator = [v for a, v in alternates.items() if a not in audit.CHARGED]
        else:
            if reference not in group:
                continue
            numerator = [v for a, v in alternates.items() if a not in audit.CHARGED]
            comparator = [v for a, v in alternates.items() if a in group and a != reference]
        if not numerator or not comparator:
            continue
        info = table["context"][key]
        out.append({**info, "position": key[1], "reference": reference,
                    "penalty": float(np.mean(numerator)) - float(np.mean(comparator))})
    return out


def fixture(seed=7, genes=24, per_cell=3):
    """Two-sided synthetic sites with every alternate present at each position."""
    rng = np.random.default_rng(seed)
    sites, context = {}, {}
    natives = "ANQVSTLMIFGWKRDE"  # Include charged native residues so the fixture exercises loss rules.
    position = 0
    for g in range(genes):
        gene = f"GENE{g}"
        for side in (0, 1):
            for _ in range(per_cell):
                for reference in natives:
                    position += 1
                    key = (f"ACC{g}", position, reference)
                    sites[key] = {a: float(rng.normal(0.0, 1.0))
                                  for a in audit.AMINO_ACIDS if a != reference}
                    context[key] = {"accession": f"ACC{g}", "gene": gene,
                                    "side": side, "distance": 1 + (position % 10)}
    return {"sites": sites, "context": context}


class ArmsStayDisjoint(unittest.TestCase):
    def setUp(self):
        self.table = fixture()

    def test_charge_rules_are_unchanged_by_the_fix(self):
        for rule, (kind, group) in audit.RULES.items():
            if not group <= audit.CHARGED:
                continue
            new = audit.penalties(self.table, kind, group)
            old = old_penalties(self.table, kind, group)
            self.assertEqual(len(new), len(old), rule)
            for a, b in zip(new, old):
                self.assertEqual(a["penalty"], b["penalty"], rule)

    def test_placebo_arms_overlapped_before_the_fix(self):
        group = audit.PLACEBO_VI
        new = audit.penalties(self.table, "gain", group)
        old = old_penalties(self.table, "gain", group)
        self.assertTrue(new, "placebo produced no rows")
        self.assertEqual(len(new), len(old))
        changed = sum(1 for a, b in zip(new, old) if a["penalty"] != b["penalty"])
        self.assertEqual(changed, len(new),
                         "holding the group out of the comparator must move every row")

    def test_placebo_contrast_is_diluted_toward_zero_when_arms_overlap(self):
        group = audit.PLACEBO_VI
        new = np.mean([abs(r["penalty"]) for r in audit.penalties(self.table, "gain", group)])
        old = np.mean([abs(r["penalty"]) for r in old_penalties(self.table, "gain", group)])
        self.assertGreater(new, old)

    def test_every_rule_yields_two_sided_strata(self):
        for rule, (kind, group) in audit.RULES.items():
            block = audit.strata(audit.penalties(self.table, kind, group))
            self.assertGreater(block["diff"].size, 1, rule)


class SharedDrawContrasts(unittest.TestCase):
    def setUp(self):
        self.table = fixture()
        self.blocks = {r: audit.strata(audit.penalties(self.table, k, g))
                       for r, (k, g) in audit.RULES.items()}

    def test_contrast_point_estimate_is_the_difference_of_rule_estimates(self):
        got = audit.paired_contrasts(self.blocks, audit.CONTRASTS, draws=200, seed=11)
        for name, row in got["contrasts"].items():
            left = got["rules"][row["minuend"]]["estimate"]
            right = got["rules"][row["subtrahend"]]["estimate"]
            self.assertAlmostEqual(row["estimate"], left - right, places=12, msg=name)

    def test_shared_draws_beat_independent_draws_for_a_contrast(self):
        """Pairing must not widen an interval relative to treating rules as independent."""
        got = audit.paired_contrasts(self.blocks, audit.CONTRASTS, draws=400, seed=3)
        row = got["contrasts"]["positive_gain_minus_placebo_vi"]
        paired_width = row["ci95"][1] - row["ci95"][0]
        left = audit.point_and_interval(self.blocks["positive_charge_gain"], 400, 3)
        right = audit.point_and_interval(self.blocks["neutral_placebo_gain_vi"], 400, 3)
        independent = float(np.hypot(left["ci95"][1] - left["ci95"][0],
                                     right["ci95"][1] - right["ci95"][0]))
        self.assertLess(paired_width, independent * 1.05)

    def test_rule_estimates_match_the_unpaired_point_estimate(self):
        got = audit.paired_contrasts(self.blocks, audit.CONTRASTS, draws=50, seed=5)
        for rule, row in got["rules"].items():
            solo = audit.point_and_interval(self.blocks[rule], 50, 5)
            self.assertAlmostEqual(row["estimate"], solo["estimate"], places=12, msg=rule)


class MatchedComparator(unittest.TestCase):
    def setUp(self):
        self.table = fixture()
        self.gain = {r: g for r, (k, g) in audit.RULES.items() if k == "gain"}

    def test_explicit_comparator_is_rejected_for_loss_rules(self):
        with self.assertRaises(ValueError):
            audit.penalties(self.table, "loss", audit.POSITIVE, audit.COMMON_COMPARATOR)

    def test_common_comparator_excludes_every_gain_group(self):
        for group in self.gain.values():
            self.assertFalse(audit.COMMON_COMPARATOR & group)

    def test_singleton_group_loses_its_own_native_sites(self):
        """The reason shared_sites is needed at all."""
        q = audit.penalties(self.table, "gain", audit.PLACEBO_Q, audit.COMMON_COMPARATOR)
        vi = audit.penalties(self.table, "gain", audit.PLACEBO_VI, audit.COMMON_COMPARATOR)
        self.assertLess(len(q), len(vi))
        self.assertFalse(any(r["reference"] == "Q" for r in q))

    def test_shared_sites_makes_every_rule_one_population(self):
        rows = audit.shared_sites(
            {r: audit.penalties(self.table, "gain", g, audit.COMMON_COMPARATOR)
             for r, g in self.gain.items()})
        identify = lambda r: (r["accession"], r["position"], r["reference"])
        populations = {frozenset(identify(r) for r in v) for v in rows.values()}
        self.assertEqual(len(populations), 1, "rules still differ in site population")
        self.assertTrue(next(iter(populations)))

    def test_matched_rules_share_comparator_arm_size(self):
        """With one comparator pool, equal-size numerators give equal-size denominators."""
        rows = audit.shared_sites(
            {r: audit.penalties(self.table, "gain", g, audit.COMMON_COMPARATOR)
             for r, g in self.gain.items()})
        counts = {r: len(v) for r, v in rows.items()}
        self.assertEqual(len(set(counts.values())), 1, counts)

    def test_matched_contrast_is_difference_of_matched_rules(self):
        rows = audit.shared_sites(
            {r: audit.penalties(self.table, "gain", g, audit.COMMON_COMPARATOR)
             for r, g in self.gain.items()})
        contrasts = {n: pair for n, pair in audit.CONTRASTS.items()
                     if pair[0] in self.gain and pair[1] in self.gain}
        got = audit.paired_contrasts({r: audit.strata(v) for r, v in rows.items()},
                                     contrasts, draws=150, seed=17)
        self.assertTrue(got["contrasts"])
        for name, row in got["contrasts"].items():
            left = got["rules"][row["minuend"]]["estimate"]
            right = got["rules"][row["subtrahend"]]["estimate"]
            self.assertAlmostEqual(row["estimate"], left - right, places=12, msg=name)


class Verdicts(unittest.TestCase):
    def test_sign_agreement_and_disagreement(self):
        self.assertEqual(audit.verdict(-0.2, [-0.3, -0.1], -1), "consistent")
        self.assertEqual(audit.verdict(+0.2, [+0.1, +0.3], -1), "inconsistent")
        self.assertEqual(audit.verdict(+0.2, [+0.1, +0.3], +1), "consistent")

    def test_interval_covering_zero_is_indeterminate(self):
        self.assertEqual(audit.verdict(0.01, [-0.2, 0.3], -1),
                         "indeterminate: interval includes zero")

    def test_placebo_is_read_against_zero(self):
        self.assertEqual(audit.verdict(0.01, [-0.2, 0.3], 0), "placebo consistent with zero")
        self.assertEqual(audit.verdict(0.3, [0.1, 0.5], 0), "placebo excludes zero")

    def test_floating_point_noise_does_not_decide_a_verdict(self):
        """A planted-zero effect lands near 1e-17 with an interval excluding zero."""
        noise = 2.7755575615628914e-17
        self.assertEqual(
            audit.verdict(noise, [1.850371707708594e-17, 3.700743415417188e-17], 0),
            "placebo consistent with zero")
        self.assertEqual(
            audit.verdict(noise, [1.850371707708594e-17, 3.700743415417188e-17], -1),
            "indeterminate: interval includes zero")

    def test_tolerance_does_not_mask_a_real_effect(self):
        self.assertEqual(audit.verdict(0.012, [0.010, 0.014], 0), "placebo excludes zero")
        self.assertEqual(audit.verdict(-0.012, [-0.014, -0.010], -1), "consistent")

    def test_undeclared_rule_gets_no_verdict(self):
        self.assertEqual(audit.verdict(0.3, [0.1, 0.5], None), "no declared prediction")


if __name__ == "__main__":
    unittest.main()

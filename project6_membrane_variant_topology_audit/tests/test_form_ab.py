"""Test the Form A/B generator, oracle, students, split, guards, and correction.

Run with: .venv_ml/Scripts/python.exe -m unittest discover -s tests -p test_form_ab.py -v
"""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import form_ab_lib as lib  # noqa: E402


def fixture_sites(seed=3, proteins=40, per_protein=120, side_composition=False):
    """Two-sided synthetic geometry; optionally a different native-residue mix on each side."""
    rng = np.random.default_rng(seed)
    protein = np.repeat(np.arange(proteins), per_protein)
    position = np.tile(np.arange(1, per_protein + 1), proteins)
    side = rng.integers(0, 2, len(protein))
    weights = np.ones(20)
    ref = rng.choice(20, len(protein), p=weights / weights.sum())
    if side_composition:
        hydrophobic = np.array([lib.ALPHABET.index(a) for a in "LIVF"])
        flip = (side == 0) & (rng.random(len(protein)) < 0.35)
        ref[flip] = rng.choice(hydrophobic, flip.sum())
    return {"protein": protein, "position": position, "side": side,
            "distance": rng.integers(1, 11, len(protein)),
            "tmd_index": np.zeros(len(protein), dtype=np.int64), "ref_aa": ref}


class GeometryOnlyLoader(unittest.TestCase):
    def test_score_array_is_never_read_and_sites_are_unique(self):
        sites = fixture_sites(proteins=3, per_protein=10)
        rows = {name: np.repeat(values, 19) for name, values in sites.items()}
        rows["alt_aa"] = lib.ALT_TABLE[sites["ref_aa"]].ravel()
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "archive.npz"
            np.savez(path, **rows, am_pathogenicity=np.zeros(len(rows["alt_aa"])),
                     accessions=np.array(["A", "B", "C"]), aa_alphabet=np.array(list(lib.ALPHABET)))
            opened = []
            real_load = np.load

            class Spy:
                def __init__(self, blob):
                    self.blob = blob

                def __enter__(self):
                    return self

                def __exit__(self, *exc):
                    self.blob.close()

                def __getitem__(self, key):
                    opened.append(key)
                    return self.blob[key]

            with mock.patch.object(lib.np, "load", lambda *a, **k: Spy(real_load(*a, **k))):
                got, accessions = lib.load_sites(path)
        self.assertNotIn("am_pathogenicity", opened)
        self.assertEqual(len(got["protein"]), 30)
        self.assertTrue(np.array_equal(got["ref_aa"], sites["ref_aa"]))
        self.assertEqual(list(accessions), ["A", "B", "C"])


class OracleAndEstimator(unittest.TestCase):
    def setUp(self):
        self.sites = fixture_sites()
        self.base = lib.blosum_base()

    def noise_free(self, amp, lam, side_main=0.0, site=True):
        parts = lib.draw_target(self.sites, self.base, amp, lam, side_main, np.random.default_rng(1))
        matrix = parts["clean"] + (parts["site"] if site else 0.0)
        return lib.mask_native(matrix, self.sites["ref_aa"]), parts

    def test_flat_planted_rule_is_recovered_exactly(self):
        scores, parts = self.noise_free(-0.08, None, side_main=0.3)
        result, _, _ = lib.audit_matrix(self.sites, scores, draws=2)
        oracle, _ = lib.oracle_delta(self.sites, parts["h"])
        self.assertAlmostEqual(oracle, -0.08, places=12)
        self.assertAlmostEqual(result["estimate_raw"], oracle, places=12)

    def test_decaying_rule_matches_the_independent_oracle(self):
        scores, parts = self.noise_free(-0.2, 8.0)
        result, _, _ = lib.audit_matrix(self.sites, scores, draws=2)
        oracle, _ = lib.oracle_delta(self.sites, parts["h"])
        self.assertLess(abs(oracle), 0.2)
        self.assertAlmostEqual(result["estimate_raw"], oracle, places=12)

    def test_site_effect_and_side_baseline_cancel(self):
        with_effects, _ = self.noise_free(0.0, None, side_main=0.3, site=True)
        result, _, _ = lib.audit_matrix(self.sites, with_effects, draws=2)
        self.assertLess(abs(result["estimate_raw"]), 1e-12)

    def test_side_blind_predictions_read_zero_even_when_the_rule_is_planted(self):
        parts = lib.draw_target(self.sites, self.base, -0.2, None, 0.3, np.random.default_rng(2))
        target = parts["clean"] + parts["site"] + parts["noise"]
        table = np.array([np.nanmean(lib.mask_native(target, self.sites["ref_aa"])[self.sites["ref_aa"] == r],
                                     axis=0) if (self.sites["ref_aa"] == r).any() else np.full(20, np.nan)
                          for r in range(20)])
        blind = lib.mask_native(table[self.sites["ref_aa"]], self.sites["ref_aa"])
        result, _, _ = lib.audit_matrix(self.sites, blind, draws=2)
        self.assertLess(abs(result["estimate_raw"]), 1e-12)

    def test_analytic_sd_matches_the_empirical_penalty_sd(self):
        sites = fixture_sites(proteins=120, per_protein=120)
        parts = lib.draw_target(sites, self.base, 0.0, None, 0.0, np.random.default_rng(4))
        scores = lib.mask_native(parts["clean"] + parts["site"] + parts["noise"], sites["ref_aa"])
        result, _, _ = lib.audit_matrix(sites, scores, draws=2)
        self.assertAlmostEqual(result["penalty_sd"] / lib.sd_z0(sites, self.base), 1.0, delta=0.03)

    def test_npz_path_matches_the_direct_path(self):
        parts = lib.draw_target(self.sites, self.base, -0.1, None, 0.0, np.random.default_rng(5))
        scores = lib.mask_native(parts["clean"] + parts["site"] + parts["noise"], self.sites["ref_aa"])
        direct, _, _ = lib.audit_matrix(self.sites, scores, draws=99, seed=7)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "s.npz"
            lib.write_scores_npz(path, self.sites, scores, np.array([f"P{i}" for i in range(40)]))
            official = lib.P.run_panel_audit(path, Path(folder) / "out", "t", "higher-is-damaging",
                                             10, 99, 7)["result"]
        self.assertAlmostEqual(direct["estimate_raw"], official["estimate_raw"], places=12)
        np.testing.assert_allclose(direct["ci95_raw"], official["ci95_raw"], atol=1e-12)
        self.assertEqual(direct["n_strata"], official["n_strata"])

    def test_naive_pooled_comparator_leaks_composition_where_the_matched_estimate_does_not(self):
        sites = fixture_sites(proteins=200, side_composition=True)
        parts = lib.draw_target(sites, self.base, 0.0, None, 0.0, np.random.default_rng(6))
        target = lib.mask_native(parts["clean"] + parts["site"] + parts["noise"], sites["ref_aa"])
        table = np.array([np.nanmean(target[sites["ref_aa"] == r], axis=0) for r in range(20)])
        blind = lib.mask_native(table[sites["ref_aa"]], sites["ref_aa"])
        matched, _, _ = lib.audit_matrix(sites, blind, draws=2)
        naive = lib.naive_comparators(sites, blind)
        self.assertLess(abs(matched["estimate_raw"]), 1e-12)
        self.assertGreater(abs(naive["pooled_paired"]), 0.02)


class DesignAlignment(unittest.TestCase):
    def test_rows_targets_and_predictions_line_up(self):
        sites = fixture_sites(proteins=2, per_protein=15)
        matrix = lib.mask_native(np.random.default_rng(9).normal(size=(30, 20)), sites["ref_aa"])
        columns = lib.site_columns(sites, np.zeros((30, 3), dtype=np.int64))
        index = np.array([4, 7, 20])
        X = lib.design_matrix(columns, index, lib.STUDENT_FEATURES["S1"])
        y = lib.row_targets(matrix[index], sites["ref_aa"][index])
        for i, site in enumerate(index):
            for j in range(19):
                row = 19 * i + j
                alt = lib.ALT_TABLE[sites["ref_aa"][site], j]
                self.assertEqual(int(X[row, 0]), sites["ref_aa"][site])
                self.assertEqual(int(X[row, 1]), alt)
                self.assertEqual(y[row], matrix[site, alt])
        back = lib.predictions_to_matrix(y, sites["ref_aa"][index])
        np.testing.assert_array_equal(np.isnan(back), np.isnan(matrix[index]))
        np.testing.assert_allclose(np.nan_to_num(back), np.nan_to_num(matrix[index]))

    def test_side_blind_design_has_no_side_or_distance_column(self):
        sites = fixture_sites(proteins=1, per_protein=5)
        columns = lib.site_columns(sites, np.zeros((5, 3), dtype=np.int64))
        self.assertEqual(lib.design_matrix(columns, np.arange(5), lib.STUDENT_FEATURES["S2"]).shape[1], 2)
        for student, names in lib.STUDENT_FEATURES.items():
            self.assertNotIn("tmd_index", names, student)
        self.assertNotIn("side", lib.STUDENT_FEATURES["S3"])
        self.assertNotIn("distance", lib.STUDENT_FEATURES["S3"])


class ContextFeatures(unittest.TestCase):
    def test_counts_match_a_brute_force_window(self):
        rng = np.random.default_rng(11)
        census = {f"P{i}": {"sequence": "".join(rng.choice(list(lib.ALPHABET), 60))} for i in range(3)}
        positions = np.array([1, 2, 30, 59, 60])
        sites = {"protein": np.repeat(np.arange(3), 5), "position": np.tile(positions, 3)}
        got = lib.context_counts(sites, np.array(list(census)), census)
        for row, (p, pos) in enumerate(zip(sites["protein"], sites["position"])):
            sequence = census[f"P{p}"]["sequence"]
            window = [sequence[j] for j in range(max(0, pos - 6), min(60, pos + 5)) if j != pos - 1]
            expected = [sum(c in letters for c in window) for letters in lib.CONTEXT_CLASSES]
            self.assertEqual(got[row].tolist(), expected, (p, pos))


class JobFiles(unittest.TestCase):
    def test_jobs_are_atomic_resumable_and_rebuild_sorted_tables(self):
        import csv
        import types
        import run_form_a as ra
        with tempfile.TemporaryDirectory() as folder:
            args = types.SimpleNamespace(out=Path(folder))
            row = lambda sc, st, rep: (
                {name: "" for name in ra.FIELDS}
                | {"scenario": sc, "student": st, "replicate": rep})
            for scenario, rep in (("P_neg05", 1), ("S_main", 0), ("N0", 0)):
                ra.save_job(args, scenario, rep, {"estimates": [row(scenario, "S1", rep), row(scenario, "S0", rep)],
                                                  "naive": [], "diag": []})
            self.assertEqual(ra.completed(args), {("P_neg05", 1), ("S_main", 0), ("N0", 0)})
            self.assertFalse(list((Path(folder) / "jobs").glob("*.tmp")))
            ra.write_tables(args)
            with (Path(folder) / "estimates.csv").open(encoding="utf-8") as stream:
                got = [(r["replicate"], r["scenario"], r["student"]) for r in csv.DictReader(stream)]
            self.assertEqual(got, [("0", "N0", "S0"), ("0", "N0", "S1"), ("0", "S_main", "S0"),
                                   ("0", "S_main", "S1"), ("1", "P_neg05", "S0"), ("1", "P_neg05", "S1")])
            self.assertFalse((Path(folder) / "naive_comparators.csv").exists())


class SplitRules(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(13)
        self.census = {}
        for i in range(50):
            sequence = "".join(rng.choice(list(lib.ALPHABET), 220))
            self.census[f"Q{i:03d}"] = {"sequence": sequence, "gene": f"GENE{i // 5}X{i % 5}"}
        self.census["Q900"] = {"sequence": self.census["Q007"]["sequence"], "gene": "GENE1X2"}

    def test_duplicates_share_a_cluster_and_the_split_is_stable(self):
        first = lib.cluster_split(self.census, seed=5)
        second = lib.cluster_split(self.census, seed=5)
        self.assertEqual(first["label"], second["label"])
        self.assertEqual(first["assignment_sha256"], second["assignment_sha256"])
        self.assertEqual(first["label"]["Q007"], first["label"]["Q900"])
        self.assertTrue(set(first["label"].values()) <= {"train", "validation", "test"})

    def test_symbol_families_split_keeps_a_family_together(self):
        got = lib.cluster_split(self.census, kind="symbol")
        families = {}
        for accession, record in self.census.items():
            families.setdefault(record["gene"].rstrip("0123456789"), set()).add(got["label"][accession])
        self.assertTrue(all(len(labels) == 1 for labels in families.values()))

    def test_bucket_boundaries(self):
        self.assertEqual([lib.bucket_label(b) for b in range(10)],
                         ["train"] * 6 + ["validation"] + ["test"] * 3)


class StudentSelection(unittest.TestCase):
    def test_iteration_is_chosen_on_the_validation_rows(self):
        rng = np.random.default_rng(17)
        X = np.column_stack([rng.integers(0, 20, 4000), rng.integers(0, 20, 4000), rng.normal(size=4000)])
        y = (X[:, 0] == X[:, 1]) * 1.0 + 0.3 * X[:, 2] + rng.normal(0, 0.2, 4000)
        Xv, yv = X[:1000], y[:1000]
        model, best, curve = lib.fit_student(X[1000:], y[1000:], Xv, yv, seed=1, max_iter=30)
        self.assertEqual(best, int(np.argmin(curve)) + 1)
        got = lib.predict_at(model, Xv, best)
        self.assertAlmostEqual(float(np.mean((got - yv) ** 2)), float(curve[best - 1]), places=10)


class FreezeGuard(unittest.TestCase):
    """require_frozen_prereg against a throwaway git repository with a fake remote branch."""

    def git(self, *args):
        import subprocess
        return subprocess.run(["git", *args], cwd=self.root, capture_output=True, text=True, check=True)

    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.root = Path(self.folder.name)
        (self.root / "docs").mkdir()
        self.prereg = self.root / "docs" / "FORM_AB_PREREGISTRATION.md"
        self.prereg.write_text("plan\n", encoding="utf-8", newline="\n")
        self.git("init", "-q")
        self.git("config", "user.email", "t@example.org")
        self.git("config", "user.name", "t")
        self.digest = lib.P.sha256(self.prereg)
        patches = [mock.patch.object(lib, "ROOT", self.root), mock.patch.object(lib, "PREREG", self.prereg),
                   mock.patch.object(lib, "frozen_file_drift", lambda text: [])]
        for patch in patches:
            patch.start()
            self.addCleanup(patch.stop)

    def tearDown(self):
        self.folder.cleanup()

    def commit(self):
        self.git("add", "docs/FORM_AB_PREREGISTRATION.md")
        self.git("commit", "-q", "-m", "freeze")
        return self.git("rev-parse", "HEAD").stdout.strip()

    def test_untracked_file_is_refused(self):
        with self.assertRaisesRegex(SystemExit, "not tracked"):
            lib.require_frozen_prereg(self.digest)

    def test_edit_after_commit_is_refused(self):
        self.commit()
        self.prereg.write_text("plan, edited\n", encoding="utf-8", newline="\n")
        with self.assertRaises(SystemExit):
            lib.require_frozen_prereg(lib.P.sha256(self.prereg))

    def test_unpushed_commit_is_refused_unless_waived(self):
        self.commit()
        with self.assertRaisesRegex(SystemExit, "remote branch"):
            lib.require_frozen_prereg(self.digest)
        self.assertFalse(lib.require_frozen_prereg(self.digest, check_push=False)["prereg_pushed"])

    def test_committed_and_pushed_passes_and_records_the_commit(self):
        commit = self.commit()
        self.git("update-ref", "refs/remotes/origin/main", commit)
        got = lib.require_frozen_prereg(self.digest)
        self.assertEqual((got["prereg_commit"], got["prereg_pushed"]), (commit, True))

    def test_hash_other_than_the_recorded_one_is_refused(self):
        commit = self.commit()
        self.git("update-ref", "refs/remotes/origin/main", commit)
        with self.assertRaisesRegex(SystemExit, "does not match"):
            lib.require_frozen_prereg("0" * 64)

    def test_split_assignment_must_be_written_in_the_preregistration(self):
        self.prereg.write_text("splits: " + "a" * 64 + "\n", encoding="utf-8")
        lib.require_split_pinned("a" * 64)
        with self.assertRaises(SystemExit):
            lib.require_split_pinned("b" * 64)


class Guards(unittest.TestCase):
    def test_validated_tool_is_the_release(self):
        info = lib.validated_tool()
        self.assertEqual(len(info["wheel_sha256"]), 64)

    def test_frozen_file_drift_names_exactly_the_unlisted_files(self):
        hashes = lib.frozen_file_hashes()
        self.assertEqual(lib.frozen_file_drift("\n".join(hashes.values())), [])
        victim = "scripts/run_form_b.py"
        text = "\n".join(h for name, h in hashes.items() if name != victim)
        self.assertEqual(lib.frozen_file_drift(text), [victim])

    def test_confirmatory_stage_refuses_a_mismatched_preregistration(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "prereg.md"
            path.write_text("plan", encoding="utf-8")
            with mock.patch.object(lib, "PREREG", path):
                with self.assertRaises(SystemExit):
                    lib.require_frozen_prereg("0" * 64)


def scorer_panel(delta, seed=21, proteins=60):
    """Real-like scorer on fixture geometry: noisy scores with a K/R-gain side contrast of about `delta`."""
    sites = fixture_sites(seed=seed, proteins=proteins)
    parts = lib.draw_target(sites, lib.blosum_base(), delta, None, 0.0, np.random.default_rng(seed))
    scores = lib.mask_native(parts["clean"] + parts["site"] + parts["noise"], sites["ref_aa"])
    return {**sites, "scores": scores, "accessions": np.array([f"P{i}" for i in range(proteins)])}


class Correction(unittest.TestCase):
    def fit(self, panel, declared=-1, margin=0.05, length=8.0):
        result = lib.audit_matrix(panel, panel["scores"], draws=2)[0]
        k = lib.correction_weight(panel, length)
        beta, target = lib.fit_beta(result["estimate_raw"], k, result["penalty_sd"], margin, declared)
        return result, k, beta, target

    def test_refitted_training_set_returns_the_target_to_numerical_precision(self):
        panel = scorer_panel(+0.15)
        result, k, beta, target = self.fit(panel)
        self.assertGreater(beta, 0.0)
        corrected = lib.apply_correction(panel, beta, 8.0)
        again = lib.audit_matrix(panel, corrected, draws=2)[0]
        self.assertAlmostEqual(again["estimate_raw"], target, places=12)

    def test_contrast_is_linear_in_beta_with_slope_K(self):
        panel = scorer_panel(+0.1)
        base, k, _, _ = self.fit(panel)
        for beta in (0.05, 0.4):
            moved = lib.audit_matrix(panel, lib.apply_correction(panel, beta, 8.0), draws=2)[0]
            self.assertAlmostEqual(moved["estimate_raw"], base["estimate_raw"] - beta * k, places=12)

    def test_no_decay_option_uses_a_constant_shift(self):
        panel = scorer_panel(+0.1)
        base, k, _, _ = self.fit(panel, length=None)
        self.assertAlmostEqual(k, 1.0, places=12)
        moved = lib.audit_matrix(panel, lib.apply_correction(panel, 0.2, None), draws=2)[0]
        self.assertAlmostEqual(moved["estimate_raw"], base["estimate_raw"] - 0.2, places=12)

    def test_scorer_already_beyond_the_target_is_left_alone(self):
        panel = scorer_panel(-0.4)
        _, _, beta, _ = self.fit(panel)
        self.assertEqual(beta, 0.0)
        np.testing.assert_array_equal(lib.apply_correction(panel, beta, 8.0), panel["scores"])

    def test_declaring_the_opposite_direction_mirrors_the_rule(self):
        consistent = scorer_panel(-0.4)
        result, k, beta, target = self.fit(consistent, declared=+1)
        self.assertGreater(beta, 0.0)
        self.assertGreater(target, 0.0)
        moved = lib.audit_matrix(consistent, lib.apply_correction(consistent, beta, 8.0, declared=+1), draws=2)[0]
        self.assertAlmostEqual(moved["estimate_raw"], target, places=12)
        already = scorer_panel(+0.4)
        self.assertEqual(self.fit(already, declared=+1)[2], 0.0)

    def test_only_kr_gain_rows_on_uncharged_natives_change(self):
        panel = scorer_panel(+0.1)
        corrected = lib.apply_correction(panel, 0.3, 8.0)
        changed = np.isfinite(panel["scores"]) & (corrected != panel["scores"])
        expected = np.zeros_like(changed)
        expected[np.ix_(np.arange(len(changed)), lib.KR_INDEX)] = True
        expected &= np.isin(panel["ref_aa"], lib.NEUTRAL_INDEX)[:, None] & np.isfinite(panel["scores"])
        np.testing.assert_array_equal(changed, expected)

    def test_shuffling_preserves_group_side_counts_and_defeats_the_correction(self):
        panel = scorer_panel(+0.1, proteins=150)
        shuffled = lib.shuffle_side(panel, np.random.default_rng(3))
        group = panel["protein"] * 20 + panel["ref_aa"]
        for g in np.unique(group)[:50]:
            self.assertEqual(int((panel["side"][group == g] == 0).sum()), int((shuffled[group == g] == 0).sum()))
        base, k, beta, target = self.fit(panel)
        real = lib.audit_matrix(panel, lib.apply_correction(panel, beta, 8.0), draws=2)[0]["estimate_raw"]
        wrong = lib.audit_matrix(panel, lib.apply_correction(panel, beta, 8.0, side=shuffled),
                                 draws=2)[0]["estimate_raw"]
        self.assertLess(abs(wrong - base["estimate_raw"]), 0.25 * abs(real - base["estimate_raw"]))


class PerProteinAndFidelity(unittest.TestCase):
    def test_protein_contrasts_aggregate_to_the_overall_estimate(self):
        panel = scorer_panel(+0.1)
        overall = lib.audit_matrix(panel, panel["scores"], draws=2)[0]["estimate_raw"]
        each = lib.protein_contrasts(panel, min_strata=1)
        self.assertAlmostEqual(float((each["contrast"] * each["weight"]).sum() / each["weight"].sum()),
                               overall, places=12)
        self.assertTrue(np.all(lib.protein_contrasts(panel, min_strata=3)["n_strata"] >= 3))

    def test_correction_moves_every_protein_toward_the_declared_direction_by_construction(self):
        panel = scorer_panel(+0.1)
        before = lib.protein_contrasts(panel, min_strata=1)
        after = lib.protein_contrasts({**panel, "scores": lib.apply_correction(panel, 0.2, 8.0)}, min_strata=1)
        self.assertTrue(np.all(after["contrast"] < before["contrast"]))

    def test_fidelity_of_identical_scores_and_of_a_shift(self):
        x = np.random.default_rng(5).normal(size=(50, 20))
        same = lib.fidelity(x, x.copy())
        self.assertEqual((same["spearman"], same["fraction_changed"], same["mean_abs_change"]), (1.0, 0.0, 0.0))
        shifted = x.copy()
        shifted[:, 3] += 0.001
        got = lib.fidelity(x, shifted)
        self.assertAlmostEqual(got["fraction_changed"], 0.05, places=12)
        self.assertGreater(got["spearman"], 0.999)


def mini_census(proteins=8):
    """Tiny topology + mapped scores for the working-copy report card."""
    rng = np.random.default_rng(31)
    layout = {"A": (7, 25), "N": (9, 27), "Q": (11, 29), "V": (13, 31), "L": (5, 36), "S": (3, 38),
              "K": (8, 32), "R": (10, 33), "D": (12, 34), "E": (14, 35)}
    census, rows = [], []
    for p in range(proteins):
        sequence = list("A" * 50)
        sequence[14:24] = list("L" * 10)
        for ref, positions in layout.items():
            for position in positions:
                sequence[position - 1] = ref
        accession = f"MINI{p:03d}"
        census.append({"uniprot": accession, "gene": accession, "sequence": "".join(sequence),
                       "tmds": [{"start": 15, "end": 24}], "numTMDs": 1,
                       "loops": {"1": {"start": 1, "end": 14, "location": "inside"},
                                 "2": {"start": 25, "end": 50, "location": "outside"}}})
        for ref, positions in layout.items():
            for side, position in enumerate(positions):
                distance = 15 - position if side == 0 else position - 24
                rows.append((p, position, side, distance, lib.ALPHABET.index(ref)))
    array = np.array(rows, dtype=np.int64)
    sites = {"protein": array[:, 0], "position": array[:, 1], "side": array[:, 2], "distance": array[:, 3],
             "tmd_index": np.zeros(len(array), dtype=np.int64), "ref_aa": array[:, 4]}
    scores = lib.mask_native(rng.normal(size=(len(array), 20)), sites["ref_aa"])
    return census, {**sites, "scores": scores, "accessions": np.array([c["uniprot"] for c in census])}


class ReportCard(unittest.TestCase):
    def test_working_copy_carries_the_placebo_rules_the_validated_wheel_lacks(self):
        tool, digest = lib.report_card_tool()
        self.assertEqual(digest, lib.P.sha256(lib.WORKING_AUDIT))
        self.assertIn("neutral_placebo_gain_vi", tool.RULES)
        self.assertEqual(tool.EXPECTED_SIGN["positive_charge_gain"], -1)
        installed = Path(lib.P.__file__).resolve().parent / "audit.py"
        self.assertNotIn("EXPECTED_SIGN", installed.read_text(encoding="utf-8"))

    def test_correction_changes_only_entries_that_involve_the_corrected_rule(self):
        tool, _ = lib.report_card_tool()
        census, panel = mini_census()
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            topology = folder / "census.json"
            topology.write_text(__import__("json").dumps(census), encoding="utf-8")
            with mock.patch.object(lib, "TOPOLOGY", topology), mock.patch.dict(
                    lib.REPORT_CARD, {"draws": 40, "folds": 2, "max_distance": 60, "seed": 3}):
                (folder / "a").mkdir(), (folder / "b").mkdir()
                before = lib.run_report_card(tool, panel, panel["scores"], folder / "a")
                after = lib.run_report_card(tool, panel, lib.apply_correction(panel, 0.5, 8.0), folder / "b")
            diff = lib.card_diff(tool, before, after)
        self.assertGreater(diff["expected_changes"], 0)
        self.assertEqual(diff["n_unexpected"], 0, diff["unexpected_changes"])
        rule = before["rule_audit"]["all"]["rules"]["positive_charge_gain"]["estimate"]
        moved = after["rule_audit"]["all"]["rules"]["positive_charge_gain"]["estimate"]
        self.assertLess(moved, rule)

    def test_card_diff_flags_a_change_to_an_uninvolved_rule(self):
        tool, _ = lib.report_card_tool()
        a = {"rule_audit": {"all": {"rules": {"neutral_placebo_gain_vi": {"estimate": 0.1},
                                              "positive_charge_gain": {"estimate": 0.2}}}}}
        b = {"rule_audit": {"all": {"rules": {"neutral_placebo_gain_vi": {"estimate": 0.2},
                                              "positive_charge_gain": {"estimate": 0.1}}}}}
        diff = lib.card_diff(tool, a, b)
        self.assertEqual(diff["expected_changes"], 1)
        self.assertEqual(diff["n_unexpected"], 1)


class NeuralStudentFunctional(unittest.TestCase):
    """The student's contrasts are linear functionals of its table; they must equal the tool's estimates."""

    def setUp(self):
        import form_b9 as b9
        self.b9 = b9
        self.tool, _ = lib.report_card_tool()
        self.census, self.panel = mini_census()
        self.table = np.random.default_rng(3).normal(size=(20, 20, 2, 60))
        self.masks = b9.rule_masks(self.tool)
        self.window = lib.window_panel(self.panel)
        gene = b9.gene_codes({c["uniprot"]: c for c in self.census}, self.panel, self.window)
        self.coeffs = {r: b9.coefficients(self.window, gene, self.masks[r][2]) for r in self.masks}

    def test_numpy_functional_equals_the_report_card_for_every_rule(self):
        scores = self.b9.student_scores(self.table, self.panel)
        folder = Path(tempfile.mkdtemp())
        topology = folder / "census.json"
        topology.write_text(__import__("json").dumps(self.census), encoding="utf-8")
        with mock.patch.object(lib, "TOPOLOGY", topology), mock.patch.dict(
                lib.REPORT_CARD, {"draws": 20, "folds": 2, "max_distance": 60, "seed": 3}):
            card = lib.run_report_card(self.tool, self.panel, scores, folder)
        rules = card["rule_audit"]["within_10_residues"]["rules"]
        found = self.b9.contrasts_numpy(self.table, self.coeffs, self.masks)
        self.assertGreaterEqual(len(rules), 6)
        for rule, row in rules.items():
            self.assertAlmostEqual(row["estimate"], found[rule], places=11, msg=rule)

    def test_torch_functional_equals_numpy(self):
        import torch
        torch.set_default_dtype(torch.float64)
        problem = self.b9.Problem(self.tool, self.b9.aggregate(self.panel), self.coeffs, self.masks, 1.0)
        got = problem.contrasts(torch.tensor(self.table))
        want = self.b9.contrasts_numpy(self.table, self.coeffs, self.masks)
        for rule in want:
            self.assertAlmostEqual(float(got[rule]), want[rule], places=12, msg=rule)

    def test_archive_dtypes_do_not_overflow_cell_indices(self):
        narrow = {name: self.panel[name].astype(dtype) for name, dtype in
                  (("protein", np.int32), ("position", np.int32), ("side", np.int8), ("distance", np.int8),
                   ("tmd_index", np.int8), ("ref_aa", np.int8))}
        loaded = {**narrow, "scores": self.panel["scores"], "accessions": self.panel["accessions"]}
        widened = lib.subset_panel(loaded, np.ones(len(loaded["protein"]), dtype=bool))
        self.assertEqual(widened["ref_aa"].dtype, np.int64)
        cells = self.b9.aggregate(widened)
        want = self.b9.aggregate(self.panel)
        np.testing.assert_array_equal(cells["n"], want["n"])
        window = lib.window_panel(loaded)
        gene = self.b9.gene_codes({c["uniprot"]: c for c in self.census}, loaded, window)
        for rule in self.masks:
            np.testing.assert_allclose(self.b9.coefficients(window, gene, self.masks[rule][2]),
                                       self.coeffs[rule], atol=1e-15)
        self.assertGreaterEqual(min(float(c.min()) for c in self.coeffs.values()), -1.0)

    def test_aggregated_cells_reproduce_the_rows(self):
        cells = self.b9.aggregate(self.panel)
        finite = np.isfinite(self.panel["scores"])
        self.assertEqual(cells["n"].sum(), finite.sum())
        self.assertAlmostEqual(cells["sum"].sum(), self.panel["scores"][finite].sum(), places=9)
        self.assertAlmostEqual(cells["sumsq"].sum(), (self.panel["scores"][finite] ** 2).sum(), places=9)

    def test_student_table_ordering_matches_cell_ids(self):
        import torch
        torch.set_default_dtype(torch.float64)
        student = self.b9.Student(8, 4)
        flat = student().detach().numpy().reshape(-1)
        r, a, s, d = 3, 8, 1, 17
        table = student().detach().numpy()
        self.assertEqual(flat[self.b9.cell_ids(r, a, s, d + 1)], table[r, a, s, d])


class B9Gate(unittest.TestCase):
    """G11: the target rule's verdict may change; any other rule's verdict may not."""

    def setUp(self):
        import form_b9 as b9
        self.b9 = b9
        self.target = {"estimate": -0.04, "ci95": [-0.05, -0.001]}

    def test_target_rule_change_alone_passes(self):
        gate = self.b9.gate_b10_5(self.target, 1.0, 0.0002,
                                  {"positive_charge_gain": ["inconsistent", "consistent"]}, 1e-16)
        self.assertTrue(gate["no_other_rule_changes_verdict"])
        self.assertTrue(gate["passed"])

    def test_another_rule_changing_verdict_fails(self):
        gate = self.b9.gate_b10_5(self.target, 1.0, 0.0002,
                                  {"positive_charge_gain": ["inconsistent", "consistent"],
                                   "positive_charge_loss": ["indeterminate: interval includes zero", "consistent"]},
                                  1e-16)
        self.assertFalse(gate["no_other_rule_changes_verdict"])
        self.assertFalse(gate["passed"])

    def test_each_numeric_condition_is_enforced(self):
        self.assertFalse(self.b9.gate_b10_5({"estimate": -0.02, "ci95": [-0.03, -0.001]}, 1.0, 0.0,
                                            {}, 0.0)["target_met"])
        self.assertFalse(self.b9.gate_b10_5(self.target, 1.0, 0.02, {}, 0.0)["fidelity_loss_at_most_0.01"])
        self.assertFalse(self.b9.gate_b10_5(self.target, 1.0, 0.0, {}, 1e-6)["functional_matches_tool"])


if __name__ == "__main__":
    unittest.main()

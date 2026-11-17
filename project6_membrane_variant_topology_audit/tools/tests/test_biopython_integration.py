# Copyright (c) 2026, btrinh0.
# SPDX-License-Identifier: LicenseRef-Biopython OR BSD-3-Clause
"""Test the optional SeqRecord API using synthetic scores and annotations.

"""

from copy import deepcopy
import csv
from io import BytesIO
import json
from pathlib import Path
import unittest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition
from Bio.SeqFeature import CompoundLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import SimpleLocation
from Bio.SeqFeature import UncertainPosition
from Bio.SeqRecord import SeqRecord
import numpy as np

from membrane_rule_audit.biopython import audit_seqrecords
from membrane_rule_audit.biopython import topology_from_seqrecord

FIXTURES = Path(__file__).resolve().parents[1] / "examples" / "biopython"


def annotated_record(identifier="ONE", outside="Extracellular"):
    """Make one explicitly annotated synthetic membrane protein."""
    return SeqRecord(
        Seq("AAAAALLLLLLLLLLAAAAA"),
        id=identifier,
        annotations={"molecule_type": "protein"},
        features=[
            SeqFeature(
                SimpleLocation(0, 5),
                type="topological domain",
                qualifiers={"description": "Cytoplasmic"},
            ),
            SeqFeature(SimpleLocation(5, 15), type="transmembrane region"),
            SeqFeature(
                SimpleLocation(15, 20),
                type="topological domain",
                qualifiers={"description": outside},
            ),
        ],
    )


def scored_rows(identifiers=("ONE", "TWO")):
    """Return a hand-computable side difference of -0.30."""
    return [
        (identifier, position, "A", alternate, value if alternate in "KR" else 0.1)
        for identifier in identifiers
        for position, value in ((4, 0.2), (17, 0.5))
        for alternate in "KRN"
    ]


def run(records, rows, **kwargs):
    """Audit a fixture without spending bootstrap draws by default."""
    kwargs.setdefault("draws", 0)
    return audit_seqrecords(
        records, rows, score_direction="higher-is-damaging", **kwargs
    )


class SeqRecordIntegration(unittest.TestCase):
    """Check data contracts and numerical behavior of the installed API."""

    def test_zero_based_features_become_inclusive_coordinates(self):
        topology = topology_from_seqrecord(annotated_record())
        self.assertEqual(topology["tmds"], [{"start": 6, "end": 15}])
        self.assertEqual(
            topology["loops"]["1"], {"start": 1, "end": 5, "location": "inside"}
        )
        self.assertEqual(topology["loops"]["2"]["end"], 20)

    def test_uniprot_xml_scores_recover_known_answer(self):
        with BytesIO((FIXTURES / "synthetic_proteins.xml").read_bytes()) as handle:
            records = list(SeqIO.parse(handle, "uniprot-xml"))
        with (FIXTURES / "toy_scores.tsv").open(encoding="utf-8", newline="") as handle:
            rows = [
                (
                    row["uniprot_id"],
                    int(row["protein_variant"][1:-1]),
                    row["protein_variant"][0],
                    row["protein_variant"][-1],
                    float(row["score"]),
                )
                for row in csv.DictReader(handle, delimiter="\t")
            ]
        report = run(records, rows, draws=99)
        self.assertAlmostEqual(report["result"]["estimate_raw"], -0.30)
        self.assertEqual(report["result"]["n_proteins"], 2)
        self.assertEqual(report["result"]["n_contributing_sites"], 16)
        self.assertEqual(report["result"]["n_strata"], 8)
        json.dumps(report, allow_nan=False)

    def test_score_direction_is_explicit_and_reverses_the_contrast(self):
        records = [annotated_record("ONE"), annotated_record("TWO")]
        report = audit_seqrecords(
            records, scored_rows(), score_direction="higher-is-tolerated", draws=0
        )
        self.assertAlmostEqual(report["result"]["estimate_raw"], 0.30)
        with self.assertRaises(ValueError):
            audit_seqrecords(records, [], score_direction="automatic")
        with self.assertRaises(TypeError):
            audit_seqrecords(records, [])

    def test_harmonic_weighting_has_independent_hand_computed_answer(self):
        records = [annotated_record("ONE"), annotated_record("TWO")]
        # ONE: inside-outside = -0.3, weight 1*1/2 = 1/2.
        # TWO: inside-outside = +0.1, weight 2*1/3 = 2/3.
        rows = scored_rows(("ONE",))
        for position, value in ((3, 0.6), (4, 0.6), (17, 0.5)):
            rows.extend(
                ("TWO", position, "A", alt, value if alt in "KR" else 0.1)
                for alt in "KRN"
            )
        report = run(records, rows)
        expected = (-0.3 * 0.5 + 0.1 * (2 / 3)) / (0.5 + 2 / 3)
        self.assertAlmostEqual(report["result"]["estimate_raw"], expected)
        self.assertEqual(report["result"]["n_contributing_sites"], 5)

    def test_whole_protein_bootstrap_covers_both_distinct_protein_effects(self):
        records = [annotated_record("ONE"), annotated_record("TWO")]
        rows = scored_rows(("ONE",))
        for position, value in ((4, 0.6), (17, 0.5)):
            rows.extend(
                ("TWO", position, "A", alt, value if alt in "KR" else 0.1)
                for alt in "KRN"
            )
        report = run(records, rows, draws=299, seed=7)
        interval = report["result"]["ci95_raw"]
        self.assertAlmostEqual(interval[0], -0.30)
        self.assertAlmostEqual(interval[1], 0.10)
        self.assertEqual(report, run(records, rows, draws=299, seed=7))

    def test_missing_comparator_arm_cannot_generate_a_result(self):
        rows = [row for row in scored_rows(("ONE",)) if row[3] == "K"]
        report = run([annotated_record()], rows)
        self.assertIsNone(report["result"]["estimate_raw"])
        self.assertIsNone(report["result"]["ci95_raw"])

    def test_unknown_side_is_counted_and_excluded(self):
        record = annotated_record(outside="Mitochondrial matrix")
        report = run([record], scored_rows(("ONE",)))
        self.assertEqual(report["ignored_unknown_side_domains"], 1)
        self.assertIsNone(report["result"]["estimate_raw"])

    def test_side_aliases_are_explicit_and_case_insensitive(self):
        for description in ("Extracellular", " lumenal ", "LUMINAL", "Periplasmic"):
            with self.subTest(description=description):
                record = annotated_record(outside=description)
                self.assertEqual(
                    topology_from_seqrecord(record)["loops"]["2"]["location"],
                    "outside",
                )

    def test_reference_mismatch_and_unknown_id_are_errors(self):
        for row in [("ONE", 4, "G", "K", 0.2), ("MISSING", 4, "A", "K", 0.2)]:
            with self.subTest(row=row), self.assertRaises(ValueError):
                run([annotated_record()], [row])

    def test_conflicting_duplicates_are_rejected(self):
        rows = [("ONE", 4, "A", "K", 0.2), ("ONE", 4, "A", "K", 0.3)]
        with self.assertRaises(ValueError):
            run([annotated_record()], rows)

    def test_identical_duplicates_do_not_change_result_or_site_counts(self):
        records = [annotated_record("ONE"), annotated_record("TWO")]
        rows = scored_rows()
        report = run(records, rows + rows)
        self.assertEqual(report["annotated_scored_flank_sites"], 4)
        self.assertEqual(report["result"], run(records, rows)["result"])

    def test_missing_duplicate_keeps_finite_score_in_either_order(self):
        records = [annotated_record("ONE"), annotated_record("TWO")]
        missing = [("ONE", 4, "A", "K", float("nan"))]
        expected = run(records, scored_rows())["result"]
        for rows in (missing + scored_rows(), scored_rows() + missing):
            with self.subTest(missing_first=np.isnan(rows[0][-1])):
                self.assertEqual(run(records, rows)["result"], expected)

    def test_fuzzy_uncertain_compound_remote_locations_are_rejected(self):
        locations = [
            SimpleLocation(BeforePosition(5), 15),
            SimpleLocation(UncertainPosition(5), 15),
            CompoundLocation([SimpleLocation(5, 9), SimpleLocation(11, 15)]),
            SimpleLocation(5, 15, ref="OTHER"),
            SimpleLocation(5, 15, strand=-1),
            None,
        ]
        for location in locations:
            record = annotated_record()
            record.features[1].location = location
            with self.subTest(location=location), self.assertRaises(ValueError):
                topology_from_seqrecord(record)

    def test_overlapping_and_out_of_bounds_annotations_are_rejected(self):
        for location in (SimpleLocation(0, 6), SimpleLocation(0, 21)):
            record = annotated_record()
            record.features[0].location = location
            with self.subTest(location=location), self.assertRaises(ValueError):
                topology_from_seqrecord(record)

    def test_fasta_record_needs_topology_annotations(self):
        with self.assertRaises(ValueError):
            topology_from_seqrecord(SeqRecord(Seq("AAAAAA"), id="ONE"))

    def test_duplicate_record_ids_are_rejected(self):
        record = annotated_record()
        with self.assertRaises(ValueError):
            run([record, deepcopy(record)], [])

    def test_invalid_positions_residues_and_scores_are_rejected(self):
        rows = [
            ("ONE", True, "A", "K", 0.2),
            ("ONE", np.bool_(True), "A", "K", 0.2),
            ("ONE", 4.0, "A", "K", 0.2),
            ("ONE", 0, "A", "K", 0.2),
            ("ONE", 4, "A", "A", 0.2),
            ("ONE", 4, "A", "K", np.inf),
            ("ONE", 4, "A", "K", -np.inf),
            ("ONE", 4, None, "K", 0.2),
            ("ONE", 4, "A", "K", None),
            ("ONE", 21, "A", "K", 0.2),
        ]
        for row in rows:
            with self.subTest(row=row), self.assertRaises(ValueError):
                run([annotated_record()], [row])

    def test_integer_options_do_not_accept_booleans_or_floats(self):
        for option in ("window", "draws", "seed"):
            for value in (True, np.bool_(False), 1.5, -1):
                with (
                    self.subTest(option=option, value=value),
                    self.assertRaises(ValueError),
                ):
                    run([annotated_record()], [], **{option: value})
        with self.assertRaises(ValueError):
            run([annotated_record()], [], window=0)

    def test_singleton_description_list_is_supported(self):
        record = annotated_record()
        record.features[0].qualifiers["description"] = ["Cytoplasmic"]
        self.assertEqual(
            topology_from_seqrecord(record)["loops"]["1"]["location"], "inside"
        )

    def test_ambiguous_description_lists_are_rejected(self):
        record = annotated_record()
        record.features[0].qualifiers["description"] = ["Cytoplasmic", "Lumenal"]
        with self.assertRaises(ValueError):
            topology_from_seqrecord(record)

    def test_empty_inputs_return_unavailable_or_raise_explicitly(self):
        report = run([annotated_record()], [])
        self.assertIsNone(report["result"]["estimate_raw"])
        self.assertEqual(report["annotated_scored_flank_sites"], 0)
        with self.assertRaises(ValueError):
            run([], [])

    def test_single_contributing_protein_has_no_across_protein_interval(self):
        report = run([annotated_record()], scored_rows(("ONE",)), draws=99)
        self.assertAlmostEqual(report["result"]["estimate_raw"], -0.30)
        self.assertIsNone(report["result"]["ci95_raw"])

    def test_native_residues_are_matched_across_sides(self):
        record = annotated_record()
        record.seq = Seq(str(record.seq)[:16] + "N" + str(record.seq)[17:])
        rows = [
            ("ONE", 4, "A", "K", 0.2),
            ("ONE", 4, "A", "N", 0.1),
            ("ONE", 17, "N", "K", 0.5),
            ("ONE", 17, "N", "A", 0.1),
        ]
        report = run([record], rows)
        self.assertIsNone(report["result"]["estimate_raw"])
        self.assertEqual(report["result"]["n_strata"], 0)

    def test_nonadjacent_domains_and_outside_window_sites_are_excluded(self):
        record = annotated_record()
        record.features[0].location = SimpleLocation(0, 4)
        report = run([record], scored_rows(("ONE",)))
        self.assertEqual(report["annotated_scored_flank_sites"], 1)
        self.assertIsNone(report["result"]["estimate_raw"])
        report = run([annotated_record()], scored_rows(("ONE",)), window=1)
        self.assertEqual(report["annotated_scored_flank_sites"], 0)

    def test_input_order_and_generator_inputs_preserve_results(self):
        records = [annotated_record("ONE"), annotated_record("TWO")]
        expected = run(records, scored_rows(), draws=99, seed=7)
        report = run(
            iter(reversed(records)), iter(reversed(scored_rows())), draws=99, seed=7
        )
        self.assertEqual(report, expected)

    def test_records_are_unchanged_and_unrelated_features_ignored(self):
        record = annotated_record()
        record.features.append(SeqFeature(None, type="active site"))
        before = deepcopy(record)
        run([record], scored_rows(("ONE",)))
        self.assertEqual(str(record.seq), str(before.seq))
        self.assertEqual(record.features, before.features)
        self.assertEqual(record.annotations, before.annotations)

    def test_explicit_nonprotein_and_missing_identity_are_rejected(self):
        for molecule_type in ("DNA", "RNA"):
            record = annotated_record()
            record.annotations["molecule_type"] = molecule_type
            with (
                self.subTest(molecule_type=molecule_type),
                self.assertRaises(ValueError),
            ):
                topology_from_seqrecord(record)
        for identifier in ("", "<unknown id>"):
            with self.subTest(identifier=identifier), self.assertRaises(ValueError):
                topology_from_seqrecord(annotated_record(identifier))
        record = SeqRecord(None, id="ONE", features=annotated_record().features)
        with self.assertRaises(ValueError):
            topology_from_seqrecord(record)


if __name__ == "__main__":
    unittest.main()

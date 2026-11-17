# Copyright (c) 2026, btrinh0.
# SPDX-License-Identifier: LicenseRef-Biopython OR BSD-3-Clause
"""Run an offline SeqIO-to-audit example with fully synthetic inputs.

Install the package with its biopython extra before running this example.
"""

import argparse
import csv
import json
from pathlib import Path

from Bio import SeqIO

from membrane_rule_audit.biopython import audit_seqrecords


def main():
    """Parse annotated UniProt XML and a score table, then print the audit."""
    fixtures = Path(__file__).with_name("biopython")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xml", type=Path, default=fixtures / "synthetic_proteins.xml")
    parser.add_argument("--scores", type=Path, default=fixtures / "toy_scores.tsv")
    parser.add_argument("--draws", type=int, default=1999)
    args = parser.parse_args()
    with args.xml.open("rb") as handle:
        records = list(SeqIO.parse(handle, "uniprot-xml"))
    with args.scores.open(encoding="utf-8", newline="") as handle:
        scores = [
            (
                row["uniprot_id"],
                int(row["protein_variant"][1:-1]),
                row["protein_variant"][0],
                row["protein_variant"][-1],
                float(row["score"]),
            )
            for row in csv.DictReader(handle, delimiter="\t")
        ]
    report = audit_seqrecords(
        records, scores, score_direction="higher-is-damaging", draws=args.draws
    )
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

"""Extract ESM1b scores at TMD flanks, in the same shape as the AlphaMissense extract.

Source: ALL_hum_isoforms_ESM1b_LLR.zip from the esm_variants portal (Brandes et al.,
Nat Genet 2023), one CSV per UniProt accession. Each file is a matrix: columns are
"<wild-type residue> <position>", rows are the 20 possible residues, values are
log-likelihood ratios.

Two conventions differ from AlphaMissense and are normalised here:

* **Sign.** ESM1b's LLR is negative for damaging substitutions; AlphaMissense's
  pathogenicity is positive for them. This writes `score = -LLR`, so in the output
  higher always means more damaging and every downstream estimator is unchanged.
* **Alphabet coverage.** ESM1b scores all 19 substitutions per residue, while the
  AlphaMissense genomic table only covers those reachable by one nucleotide change.
  Everything is kept here; restricting to a common set is done at comparison time.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import zipfile
from array import array
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from extract_am_flank_scores import AA_ALPHABET, AA_INDEX, build_flank_map
from run_gate1 import read_topology, sha256

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TOPOLOGY = ROOT / "data" / "raw" / "topology_viewer.html"
DEFAULT_ZIP = ROOT / "data" / "raw" / "ESM1b_ALL_hum_isoforms_LLR.zip"
DEFAULT_OUT = ROOT / "data" / "processed" / "gate0"

MEMBER_TEMPLATE = "content/ALL_hum_isoforms_ESM1b_LLR/{accession}_LLR.csv"


def parse_matrix(handle: io.TextIOBase) -> tuple[dict[int, str], dict[str, list[str]]]:
    """Return {position: wild-type residue} and {alt residue: row of values}."""
    reader = csv.reader(handle)
    header = next(reader)
    positions: dict[int, str] = {}
    for column, label in enumerate(header[1:], start=1):
        parts = label.strip().split()
        if len(parts) == 2 and parts[1].isdigit():
            positions[column] = parts[0]
    rows = {row[0].strip(): row for row in reader if row and row[0].strip()}
    return positions, rows


def run(args: argparse.Namespace) -> dict:
    topology_path, zip_path = Path(args.topology), Path(args.zip)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    proteins = read_topology(topology_path)
    accessions = sorted(proteins)
    accession_index = {accession: i for i, accession in enumerate(accessions)}
    flank_maps = {a: build_flank_map(r, args.max_distance) for a, r in proteins.items()}

    protein_col, position_col = array("i"), array("i")
    side_col, distance_col, tmd_col = array("b"), array("b"), array("b")
    ref_col, alt_col, score_col = array("b"), array("b"), array("f")
    qc: dict[str, int] = defaultdict(int)

    with zipfile.ZipFile(zip_path) as archive:
        available = set(archive.namelist())
        for accession in accessions:
            flanks = flank_maps[accession]
            if not flanks:
                continue
            member = MEMBER_TEMPLATE.format(accession=accession)
            if member not in available:
                qc["census_proteins_absent_from_esm1b"] += 1
                continue
            qc["census_proteins_found"] += 1
            sequence = proteins[accession]["sequence"]
            with archive.open(member) as raw:
                with io.TextIOWrapper(raw, encoding="utf-8", errors="replace") as handle:
                    header_positions, rows = parse_matrix(handle)

            # The header carries a position index; column labels are the source of
            # truth for which residue ESM1b thinks is wild type there.
            label_by_column = {}
            with archive.open(member) as raw:
                with io.TextIOWrapper(raw, encoding="utf-8", errors="replace") as handle:
                    header = next(csv.reader(handle))
            for column, label in enumerate(header[1:], start=1):
                parts = label.strip().split()
                if len(parts) == 2 and parts[1].isdigit():
                    label_by_column[column] = (parts[0], int(parts[1]))

            for column, (ref_aa, position) in label_by_column.items():
                flank = flanks.get(position)
                if flank is None:
                    continue
                if position > len(sequence) or sequence[position - 1] != ref_aa:
                    qc["reference_sequence_mismatch"] += 1
                    continue
                if ref_aa not in AA_INDEX:
                    continue
                side_idx, distance, tmd_index = flank
                for alt_aa, row in rows.items():
                    if alt_aa == ref_aa or alt_aa not in AA_INDEX:
                        continue
                    if column >= len(row):
                        continue
                    raw_value = row[column].strip()
                    if not raw_value:
                        continue
                    try:
                        llr = float(raw_value)
                    except ValueError:
                        continue
                    protein_col.append(accession_index[accession])
                    position_col.append(position)
                    side_col.append(side_idx)
                    distance_col.append(distance)
                    tmd_col.append(min(tmd_index, 127))
                    ref_col.append(AA_INDEX[ref_aa])
                    alt_col.append(AA_INDEX[alt_aa])
                    score_col.append(-llr)  # higher = more damaging
                    qc["flank_rows_kept"] += 1

    out_path = out_dir / "esm1b_flank_scores.npz"
    np.savez_compressed(
        out_path,
        protein=np.frombuffer(protein_col, dtype=np.int32),
        position=np.frombuffer(position_col, dtype=np.int32),
        side=np.frombuffer(side_col, dtype=np.int8),
        distance=np.frombuffer(distance_col, dtype=np.int8),
        tmd_index=np.frombuffer(tmd_col, dtype=np.int8),
        ref_aa=np.frombuffer(ref_col, dtype=np.int8),
        alt_aa=np.frombuffer(alt_col, dtype=np.int8),
        am_pathogenicity=np.frombuffer(score_col, dtype=np.float32),
        accessions=np.array(accessions),
        aa_alphabet=np.array(list(AA_ALPHABET)),
        sides=np.array(["cytosolic", "outer"]),
    )
    summary = {
        "analysis": "ESM1b scores at TMD flanks (no clinical labels)",
        "score_convention": "negated log-likelihood ratio; higher = more damaging",
        "alphabet_coverage": "all 19 substitutions per residue",
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "max_distance_residues": args.max_distance,
        "quality_control": dict(qc),
        "output": {"path": str(out_path), "bytes": out_path.stat().st_size},
        "data_files": {"esm1b_zip": {"sha256": sha256(zip_path)}},
    }
    (out_dir / "esm1b_extract_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k != "data_files"}, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", default=str(DEFAULT_TOPOLOGY))
    parser.add_argument("--zip", default=str(DEFAULT_ZIP))
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    parser.add_argument("--max-distance", type=int, default=60)
    run(parser.parse_args())


if __name__ == "__main__":
    main()

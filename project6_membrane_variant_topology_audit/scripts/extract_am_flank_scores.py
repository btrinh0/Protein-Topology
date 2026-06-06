"""Extract AlphaMissense scores at TMD-flank positions of the human ER membrane census.

Label-free: this script never reads ClinVar. It produces the substrate for Gate 0,
which asks whether AlphaMissense itself encodes the positive-inside rule.

Output is a compact .npz with one row per (canonical protein position, alternate
amino acid) inside `--max-distance` residues of an annotated TMD edge.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from array import array
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from run_gate1 import read_topology, sha256

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TOPOLOGY = ROOT / "data" / "raw" / "topology_viewer.html"
DEFAULT_AM = ROOT / "data" / "raw" / "AlphaMissense_hg38.tsv.gz"
DEFAULT_OUT = ROOT / "data" / "processed" / "gate0"

AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {aa: i for i, aa in enumerate(AA_ALPHABET)}
SIDE_INDEX = {"cytosolic": 0, "outer": 1}


def build_flank_map(protein: dict, max_distance: int) -> dict[int, tuple[int, int, int]]:
    """Map loop position -> (side index, distance to nearest TMD edge, 1-based TMD index).

    Boundaries are taken from TMD coordinates rather than loop numbering, because 369
    census proteins have gaps in their loop keys where two TMDs abut directly.
    """
    ends = {tmd["end"]: i for i, tmd in enumerate(protein["tmds"], start=1)}
    starts = {tmd["start"]: i for i, tmd in enumerate(protein["tmds"], start=1)}
    flanks: dict[int, tuple[int, int, int]] = {}
    for loop in protein["loops"].values():
        side = loop.get("location", "").strip().lower()
        if side not in {"inside", "outside"}:
            continue
        side_idx = SIDE_INDEX["cytosolic" if side == "inside" else "outer"]
        start, end = int(loop["start"]), int(loop["end"])
        borders = []
        if start - 1 in ends:
            borders.append((start - 1, ends[start - 1], -1))
        if end + 1 in starts:
            borders.append((end + 1, starts[end + 1], 1))
        if not borders:
            continue
        for position in range(start, end + 1):
            best = None
            for boundary, tmd_index, direction in borders:
                distance = (position - boundary) if direction == -1 else (boundary - position)
                if 0 < distance <= max_distance and (best is None or distance < best[0]):
                    best = (distance, tmd_index)
            if best is not None:
                flanks[position] = (side_idx, best[0], best[1])
    return flanks


def run(args: argparse.Namespace) -> dict:
    topology_path = Path(args.topology)
    am_path = Path(args.alphamissense)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    proteins = read_topology(topology_path)
    accessions = sorted(proteins)
    accession_index = {accession: i for i, accession in enumerate(accessions)}

    flank_maps = {
        accession: build_flank_map(record, args.max_distance)
        for accession, record in proteins.items()
    }
    total_flank_positions = sum(len(flanks) for flanks in flank_maps.values())

    protein_col = array("i")
    position_col = array("i")
    side_col = array("b")
    distance_col = array("b")
    tmd_col = array("b")
    ref_col = array("b")
    alt_col = array("b")
    score_col = array("f")

    qc: dict[str, int] = defaultdict(int)
    with gzip.open(am_path, "rt", encoding="utf-8-sig", newline="") as handle:
        for line in handle:
            if line.startswith("#CHROM\t"):
                break
        else:
            raise ValueError("AlphaMissense file is missing the #CHROM header")

        for line in handle:
            qc["alphamissense_rows_scanned"] += 1
            if args.limit and qc["alphamissense_rows_scanned"] > args.limit:
                break
            fields = line.rstrip("\n").split("\t")
            accession = fields[5]
            flanks = flank_maps.get(accession)
            if not flanks:
                continue
            qc["rows_in_census_proteins"] += 1
            variant = fields[7]
            try:
                position = int(variant[1:-1])
            except ValueError:
                qc["non_simple_protein_change"] += 1
                continue
            flank = flanks.get(position)
            if flank is None:
                continue
            ref_aa, alt_aa = variant[0], variant[-1]
            sequence = proteins[accession]["sequence"]
            if position > len(sequence) or sequence[position - 1] != ref_aa:
                qc["reference_sequence_mismatch"] += 1
                continue
            if ref_aa not in AA_INDEX or alt_aa not in AA_INDEX:
                qc["non_standard_amino_acid"] += 1
                continue
            side_idx, distance, tmd_index = flank
            protein_col.append(accession_index[accession])
            position_col.append(position)
            side_col.append(side_idx)
            distance_col.append(distance)
            tmd_col.append(min(tmd_index, 127))
            ref_col.append(AA_INDEX[ref_aa])
            alt_col.append(AA_INDEX[alt_aa])
            score_col.append(float(fields[8]))
            qc["flank_rows_kept"] += 1

    out_path = out_dir / "am_flank_scores.npz"
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
        "analysis": "Gate 0 substrate: AlphaMissense scores at TMD flanks (no clinical labels)",
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "max_distance_residues": args.max_distance,
        "topology_proteins": len(proteins),
        "flank_positions_in_census": total_flank_positions,
        "quality_control": dict(qc),
        "output": {"path": str(out_path), "bytes": out_path.stat().st_size},
        "data_files": {
            "topology_viewer_html": {
                "path": str(topology_path), "sha256": sha256(topology_path),
            },
            "alphamissense_hg38": {
                "path": str(am_path), "sha256": sha256(am_path),
            },
        },
    }
    summary_path = out_dir / "extract_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k != "data_files"}, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", default=str(DEFAULT_TOPOLOGY))
    parser.add_argument("--alphamissense", default=str(DEFAULT_AM))
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    parser.add_argument("--max-distance", type=int, default=60)
    parser.add_argument("--limit", type=int, default=0, help="stop after N AlphaMissense rows (smoke test)")
    run(parser.parse_args())


if __name__ == "__main__":
    main()

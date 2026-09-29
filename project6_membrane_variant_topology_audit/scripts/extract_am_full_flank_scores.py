"""Extract AlphaMissense scores at TMD flanks from the ALL-substitutions table.

The genomic hg38 table only contains substitutions reachable by a single nucleotide
change, about 6 of the 19 per residue, and which 6 depends on the codon. That is a
real limitation of the earlier audit: the charge-preserving comparator K<->R needs
an AGA/AGG arginine codon, so CGN arginines were dropped entirely.

AlphaMissense_aa_substitutions.tsv.gz carries all 19 substitutions per residue, so
this removes the codon selection and puts AlphaMissense on the same alphabet as
ESM1b and a saturating deep mutational scan.

Output matches the schema of extract_am_flank_scores.py so every downstream
estimator runs unchanged.
"""

from __future__ import annotations

import argparse
import gzip
import json
from array import array
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from extract_am_flank_scores import AA_ALPHABET, AA_INDEX, build_flank_map
from run_gate1 import read_topology, sha256

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TOPOLOGY = ROOT / "data" / "raw" / "topology_viewer.html"
DEFAULT_TABLE = ROOT / "data" / "raw" / "AlphaMissense_aa_substitutions.tsv.gz"
DEFAULT_OUT = ROOT / "data" / "processed" / "gate0"


def run(args: argparse.Namespace) -> dict:
    topology_path, table_path = Path(args.topology), Path(args.table)
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

    with gzip.open(table_path, "rt", encoding="utf-8-sig", newline="") as handle:
        for line in handle:
            if line.startswith("uniprot_id\t"):
                break
            if not line.startswith("#"):
                raise ValueError("Unexpected preamble before the header row")
        for line in handle:
            qc["rows_scanned"] += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            accession = fields[0]
            flanks = flank_maps.get(accession)
            if not flanks:
                continue
            qc["rows_in_census_proteins"] += 1
            variant = fields[1]
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
            score_col.append(float(fields[2]))
            qc["flank_rows_kept"] += 1

    out_path = out_dir / "am_full_flank_scores.npz"
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
        "analysis": "AlphaMissense all-substitutions scores at TMD flanks",
        "alphabet_coverage": "all 19 substitutions per residue",
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "max_distance_residues": args.max_distance,
        "quality_control": dict(qc),
        "output": {"path": str(out_path), "bytes": out_path.stat().st_size},
        "data_files": {"alphamissense_aa_substitutions": {"sha256": sha256(table_path)}},
    }
    (out_dir / "am_full_extract_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k != "data_files"}, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", default=str(DEFAULT_TOPOLOGY))
    parser.add_argument("--table", default=str(DEFAULT_TABLE))
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    parser.add_argument("--max-distance", type=int, default=60)
    run(parser.parse_args())


if __name__ == "__main__":
    main()

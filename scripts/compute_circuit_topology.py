from __future__ import annotations

import argparse
import csv
import gzip
import math
from pathlib import Path
from typing import Iterable, TextIO

from Bio.PDB import MMCIFParser, PDBParser, Polypeptide


REPRESENTATIVE_ATOMS = ("CB", "CA")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute first-pass circuit-topology summaries for protein chains."
    )
    parser.add_argument(
        "--chain-manifest",
        type=Path,
        default=Path("data/processed/chain_manifest.csv"),
        help="Chain manifest produced by ingest_structures.py.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/processed/chain_topology_manifest.csv"),
        help="CSV path for chain topology summaries.",
    )
    parser.add_argument(
        "--distance-threshold",
        type=float,
        default=8.0,
        help="Representative-atom distance threshold in angstroms.",
    )
    parser.add_argument(
        "--min-sequence-separation",
        type=int,
        default=4,
        help="Minimum absolute sequence separation for residue contacts.",
    )
    return parser.parse_args()


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def detect_format(path: Path) -> str:
    name = path.name.lower()
    if name.endswith((".cif", ".mmcif", ".cif.gz", ".mmcif.gz")):
        return "mmcif"
    if name.endswith((".pdb", ".ent", ".pdb.gz", ".ent.gz")):
        return "pdb"
    raise ValueError(f"Unsupported structure format: {path}")


def load_structure(path: Path, file_format: str, structure_id: str):
    parser = MMCIFParser(QUIET=True) if file_format == "mmcif" else PDBParser(QUIET=True)
    with open_text(path) as handle:
        return parser.get_structure(structure_id, handle)


def read_chain_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def representative_atom(residue):
    for atom_name in REPRESENTATIVE_ATOMS:
        if atom_name == "CA" and residue.get_resname() != "GLY" and "CB" in residue:
            continue
        if atom_name in residue:
            return residue[atom_name]
    return None


def chain_residues(chain) -> list[tuple[int, object]]:
    residues: list[tuple[int, object]] = []
    for residue in chain:
        if residue.id[0] != " ":
            continue
        if not Polypeptide.is_aa(residue, standard=False):
            continue
        atom = representative_atom(residue)
        if atom is None:
            continue
        residues.append((residue.id[1], atom))
    return residues


def contact_pairs(
    residues: list[tuple[int, object]],
    distance_threshold: float,
    min_sequence_separation: int,
) -> list[tuple[int, int]]:
    contacts: list[tuple[int, int]] = []
    threshold_sq = distance_threshold * distance_threshold

    for i in range(len(residues)):
        seq_i, atom_i = residues[i]
        coord_i = atom_i.coord
        for j in range(i + 1, len(residues)):
            seq_j, atom_j = residues[j]
            if abs(seq_j - seq_i) < min_sequence_separation:
                continue
            diff = coord_i - atom_j.coord
            dist_sq = float(diff[0] ** 2 + diff[1] ** 2 + diff[2] ** 2)
            if dist_sq <= threshold_sq:
                left, right = sorted((seq_i, seq_j))
                contacts.append((left, right))

    return contacts


def relation(contact_a: tuple[int, int], contact_b: tuple[int, int]) -> str | None:
    i, j = contact_a
    k, l = contact_b

    if len({i, j, k, l}) < 4:
        return None
    if j < k or l < i:
        return "series"
    if (i < k < l < j) or (k < i < j < l):
        return "parallel"
    if (i < k < j < l) or (k < i < l < j):
        return "cross"
    return None


def psx_counts(contacts: list[tuple[int, int]]) -> tuple[int, int, int, int]:
    series = 0
    parallel = 0
    cross = 0
    pair_count = 0

    for idx, contact_a in enumerate(contacts):
        for contact_b in contacts[idx + 1 :]:
            rel = relation(contact_a, contact_b)
            if rel is None:
                continue
            pair_count += 1
            if rel == "series":
                series += 1
            elif rel == "parallel":
                parallel += 1
            elif rel == "cross":
                cross += 1

    return series, parallel, cross, pair_count


def safe_fraction(count: int, total: int) -> str:
    if total == 0:
        return ""
    return f"{count / total:.6f}"


def main() -> int:
    args = parse_args()
    rows = read_chain_manifest(args.chain_manifest)
    by_structure: dict[tuple[str, str], dict[str, str]] = {
        (row["structure_id"], row["chain_id"]): row for row in rows
    }
    structures = sorted({(row["structure_id"], row["source_path"], row["file_format"]) for row in rows})

    topology_rows: list[dict[str, str]] = []

    for structure_id, source_path, file_format in structures:
        structure = load_structure(Path(source_path), file_format or detect_format(Path(source_path)), structure_id)
        first_model = next(structure.get_models())

        for chain in first_model:
            key = (structure_id, chain.id)
            if key not in by_structure:
                continue

            residues = chain_residues(chain)
            contacts = contact_pairs(
                residues,
                distance_threshold=args.distance_threshold,
                min_sequence_separation=args.min_sequence_separation,
            )
            series, parallel, cross, pair_count = psx_counts(contacts)

            topology_rows.append(
                {
                    "structure_id": structure_id,
                    "chain_id": chain.id,
                    "source_path": source_path,
                    "file_format": file_format,
                    "eligible_residue_count": str(len(residues)),
                    "contact_count": str(len(contacts)),
                    "contact_pair_count": str(pair_count),
                    "series_pair_count": str(series),
                    "parallel_pair_count": str(parallel),
                    "cross_pair_count": str(cross),
                    "s_fraction": safe_fraction(series, pair_count),
                    "p_fraction": safe_fraction(parallel, pair_count),
                    "x_fraction": safe_fraction(cross, pair_count),
                }
            )

    write_csv(
        args.output,
        [
            "structure_id",
            "chain_id",
            "source_path",
            "file_format",
            "eligible_residue_count",
            "contact_count",
            "contact_pair_count",
            "series_pair_count",
            "parallel_pair_count",
            "cross_pair_count",
            "s_fraction",
            "p_fraction",
            "x_fraction",
        ],
        topology_rows,
    )

    print(f"Wrote topology summaries for {len(topology_rows)} chains.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

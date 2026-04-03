from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path
from typing import Iterable, TextIO

from Bio.PDB import MMCIFParser, PDBParser, Polypeptide
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.parse_pdb_header import parse_pdb_header


SUPPORTED_SUFFIXES = (
    ".cif",
    ".mmcif",
    ".pdb",
    ".ent",
    ".cif.gz",
    ".mmcif.gz",
    ".pdb.gz",
    ".ent.gz",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build structure and chain manifests from local PDB/mmCIF files."
    )
    parser.add_argument("input_dir", type=Path, help="Directory containing structure files.")
    parser.add_argument(
        "--structure-output",
        type=Path,
        default=Path("data/processed/structure_manifest.csv"),
        help="CSV path for structure-level metadata.",
    )
    parser.add_argument(
        "--chain-output",
        type=Path,
        default=Path("data/processed/chain_manifest.csv"),
        help="CSV path for chain-level metadata.",
    )
    parser.add_argument(
        "--failures-output",
        type=Path,
        default=Path("data/processed/ingest_failures.csv"),
        help="CSV path for files that failed to parse.",
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


def supported_files(root: Path) -> list[Path]:
    return sorted(
        path
        for path in root.rglob("*")
        if path.is_file() and path.name.lower().endswith(SUPPORTED_SUFFIXES)
    )


def stem_without_compression(path: Path) -> str:
    if path.suffix == ".gz":
        return Path(path.stem).stem
    return path.stem


def first_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, list):
        for item in value:
            if item not in (None, "?", "."):
                return str(item)
        return ""
    if value in ("?", "."):
        return ""
    return str(value)


def list_value(value: object) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item) for item in value if item not in (None, "?", ".")]
    if value in ("?", "."):
        return []
    return [str(value)]


def load_mmcif_metadata(path: Path) -> dict[str, str]:
    with open_text(path) as handle:
        mmcif = MMCIF2Dict(handle)

    methods = list_value(mmcif.get("_exptl.method"))
    resolution = first_value(mmcif.get("_refine.ls_d_res_high"))
    if not resolution:
        resolution = first_value(mmcif.get("_em_3d_reconstruction.resolution"))

    revision_dates = list_value(mmcif.get("_pdbx_audit_revision_history.revision_date"))
    release_date = revision_dates[0] if revision_dates else ""

    return {
        "structure_id": first_value(mmcif.get("_entry.id")).lower() or stem_without_compression(path).lower(),
        "experimental_method": ";".join(methods),
        "resolution": resolution,
        "deposition_date": first_value(mmcif.get("_pdbx_database_status.recvd_initial_deposition_date")),
        "release_date": release_date,
        "title": first_value(mmcif.get("_struct.title")),
    }


def load_pdb_metadata(path: Path) -> dict[str, str]:
    with open_text(path) as handle:
        header = parse_pdb_header(handle)

    return {
        "structure_id": str(header.get("idcode") or stem_without_compression(path)).lower(),
        "experimental_method": str(header.get("structure_method") or ""),
        "resolution": "" if header.get("resolution") is None else str(header.get("resolution")),
        "deposition_date": str(header.get("deposition_date") or ""),
        "release_date": str(header.get("release_date") or ""),
        "title": str(header.get("name") or ""),
    }


def load_structure(path: Path, file_format: str, structure_id: str):
    parser = MMCIFParser(QUIET=True) if file_format == "mmcif" else PDBParser(QUIET=True)
    with open_text(path) as handle:
        return parser.get_structure(structure_id, handle)


def chain_stats(structure) -> list[dict[str, str]]:
    first_model = next(structure.get_models())
    rows: list[dict[str, str]] = []

    for chain in first_model:
        residues = [
            residue
            for residue in chain
            if residue.id[0] == " " and Polypeptide.is_aa(residue, standard=False)
        ]
        if not residues:
            continue

        seq_numbers = [residue.id[1] for residue in residues]
        rows.append(
            {
                "chain_id": chain.id,
                "residue_count": str(len(residues)),
                "resolved_residue_count": str(len(residues)),
                "min_resseq": str(min(seq_numbers)),
                "max_resseq": str(max(seq_numbers)),
            }
        )

    return rows


def structure_row(path: Path, file_format: str, metadata: dict[str, str], chains: list[dict[str, str]], structure) -> dict[str, str]:
    protein_residue_count = sum(int(row["residue_count"]) for row in chains)
    chain_ids = ",".join(row["chain_id"] for row in chains)

    return {
        "structure_id": metadata["structure_id"],
        "source_path": str(path.resolve()),
        "file_format": file_format,
        "experimental_method": metadata["experimental_method"],
        "resolution": metadata["resolution"],
        "deposition_date": metadata["deposition_date"],
        "release_date": metadata["release_date"],
        "title": metadata["title"],
        "model_count": str(sum(1 for _ in structure.get_models())),
        "protein_chain_count": str(len(chains)),
        "protein_residue_count": str(protein_residue_count),
        "protein_chain_ids": chain_ids,
    }


def write_csv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    files = supported_files(args.input_dir)
    if not files:
        print(f"No supported structure files found in {args.input_dir}", file=sys.stderr)
        return 1

    structure_rows: list[dict[str, str]] = []
    chain_rows: list[dict[str, str]] = []
    failure_rows: list[dict[str, str]] = []

    for path in files:
        try:
            file_format = detect_format(path)
            metadata = load_mmcif_metadata(path) if file_format == "mmcif" else load_pdb_metadata(path)
            structure = load_structure(path, file_format, metadata["structure_id"])
            chains = chain_stats(structure)
            if not chains:
                raise ValueError("No protein chains detected in first model.")

            structure_rows.append(structure_row(path, file_format, metadata, chains, structure))
            for row in chains:
                chain_rows.append(
                    {
                        "structure_id": metadata["structure_id"],
                        "chain_id": row["chain_id"],
                        "source_path": str(path.resolve()),
                        "file_format": file_format,
                        "residue_count": row["residue_count"],
                        "resolved_residue_count": row["resolved_residue_count"],
                        "min_resseq": row["min_resseq"],
                        "max_resseq": row["max_resseq"],
                    }
                )
        except Exception as exc:
            failure_rows.append(
                {
                    "source_path": str(path.resolve()),
                    "error": str(exc),
                }
            )

    write_csv(
        args.structure_output,
        [
            "structure_id",
            "source_path",
            "file_format",
            "experimental_method",
            "resolution",
            "deposition_date",
            "release_date",
            "title",
            "model_count",
            "protein_chain_count",
            "protein_residue_count",
            "protein_chain_ids",
        ],
        structure_rows,
    )
    write_csv(
        args.chain_output,
        [
            "structure_id",
            "chain_id",
            "source_path",
            "file_format",
            "residue_count",
            "resolved_residue_count",
            "min_resseq",
            "max_resseq",
        ],
        chain_rows,
    )
    write_csv(args.failures_output, ["source_path", "error"], failure_rows)

    print(
        f"Ingested {len(structure_rows)} structures, "
        f"{len(chain_rows)} protein chains, "
        f"{len(failure_rows)} failures."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

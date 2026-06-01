from __future__ import annotations

"""Evaluate TM-align similarity on the locked V2 matched pairs."""

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

DEPS = Path(__file__).parent.parent / ".deps" / "tmtools312"
sys.path.insert(0, str(DEPS))
sys.path.insert(1, str(Path(__file__).parent))

from Bio.PDB import Polypeptide
from tmtools import tm_align

from compute_domain_topology_v2 import amino_acid_residues, detect_format, load_structure, select_domain_residues


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate TM-align scores on V2 matched pairs.")
    parser.add_argument("--pairs", type=Path, default=Path("data/processed/v2/structural_coherence_pdbredo_primary_pairs_v2.csv"))
    parser.add_argument("--topology", type=Path, default=Path("data/processed/v2/domain_topology_pdbredo_primary_alpha_beta_v2.csv"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/tmalign_structure_pdbredo_primary_pairs_v2.csv"))
    parser.add_argument("--summary-output", type=Path, default=Path("data/processed/v2/tmalign_structure_pdbredo_primary_summary_v2.json"))
    parser.add_argument("--errors-output", type=Path, default=Path("data/processed/v2/tmalign_structure_pdbredo_primary_errors_v2.csv"))
    parser.add_argument("--max-pairs", type=int, default=0)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def domain_coordinates(row: dict[str, str]) -> tuple[str, np.ndarray]:
    structure_path = Path(row["source_path"])
    structure = load_structure(structure_path, detect_format(structure_path), row["structure_id"])
    model = next(structure.get_models())
    residues = amino_acid_residues(model[row["chain_id"]])
    selected = select_domain_residues(residues, json.loads(row["domain_segments_json"]), row["chain_id"])
    sequence: list[str] = []
    coordinates: list[list[float]] = []
    for residue, _ in selected:
        if "CA" not in residue:
            continue
        sequence.append(Polypeptide.protein_letters_3to1.get(residue.get_resname().upper(), "X"))
        coordinates.append([float(value) for value in residue["CA"].coord])
    if len(sequence) < 3:
        raise ValueError("fewer than three C-alpha coordinates")
    return "".join(sequence), np.asarray(coordinates, dtype=np.float64)


def aligned_length(seq_x: str, seq_y: str) -> int:
    return sum(left != "-" and right != "-" for left, right in zip(seq_x, seq_y))


def mean(values: list[float]) -> float | None:
    return sum(values) / len(values) if values else None


def median(values: list[float]) -> float | None:
    if not values:
        return None
    values = sorted(values)
    middle = len(values) // 2
    return values[middle] if len(values) % 2 else (values[middle - 1] + values[middle]) / 2.0


def main() -> int:
    args = parse_args()
    pairs = read_csv(args.pairs)
    if args.max_pairs > 0:
        pairs = pairs[: args.max_pairs]
    topology = {row["domain_uid"]: row for row in read_csv(args.topology)}
    bundles: dict[str, tuple[str, np.ndarray]] = {}
    output_rows: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
    grouped: dict[tuple[str, str, str], list[dict[str, float]]] = defaultdict(list)
    for pair in pairs:
        left_uid, right_uid = pair["left_domain_uid"], pair["right_domain_uid"]
        try:
            if left_uid not in bundles:
                bundles[left_uid] = domain_coordinates(topology[left_uid])
            if right_uid not in bundles:
                bundles[right_uid] = domain_coordinates(topology[right_uid])
            left_sequence, left_coordinates = bundles[left_uid]
            right_sequence, right_coordinates = bundles[right_uid]
            result = tm_align(left_coordinates, right_coordinates, left_sequence, right_sequence)
            aligned_count = aligned_length(result.seqxA, result.seqyA)
            tm_mean = (float(result.tm_norm_chain1) + float(result.tm_norm_chain2)) / 2.0
            output_rows.append({
                **pair,
                "tm_norm_chain1": f"{float(result.tm_norm_chain1):.6f}",
                "tm_norm_chain2": f"{float(result.tm_norm_chain2):.6f}",
                "tm_score_mean": f"{tm_mean:.6f}",
                "tm_align_rmsd_angstrom": f"{float(result.rmsd):.6f}",
                "aligned_residue_count": str(aligned_count),
                "aligned_length_fraction": f"{aligned_count / max(max(len(left_sequence), len(right_sequence)), 1):.6f}",
            })
            grouped[(pair["threshold"], pair["split"], pair["pair_type"])].append({
                "tm1": float(result.tm_norm_chain1),
                "tm2": float(result.tm_norm_chain2),
                "tm_mean": tm_mean,
                "rmsd": float(result.rmsd),
            })
        except Exception as exc:
            errors.append({
                "threshold": pair.get("threshold", ""),
                "split": pair.get("split", ""),
                "pair_type": pair.get("pair_type", ""),
                "left_domain_uid": left_uid,
                "right_domain_uid": right_uid,
                "error": str(exc),
            })
    summary_rows: list[dict[str, object]] = []
    for (threshold, split, pair_type), values in sorted(grouped.items()):
        summary_rows.append({
            "threshold": float(threshold),
            "split": split,
            "pair_type": pair_type,
            "pair_count": len(values),
            "mean_tm_norm_chain1": mean([row["tm1"] for row in values]),
            "mean_tm_norm_chain2": mean([row["tm2"] for row in values]),
            "mean_tm_score": mean([row["tm_mean"] for row in values]),
            "median_tm_score": median([row["tm_mean"] for row in values]),
            "mean_tm_align_rmsd_angstrom": mean([row["rmsd"] for row in values]),
        })
    for threshold in sorted({row["threshold"] for row in pairs}):
        for split in sorted({row["split"] for row in pairs}):
            within = next((row for row in summary_rows if row["threshold"] == float(threshold) and row["split"] == split and row["pair_type"] == "within"), None)
            between = next((row for row in summary_rows if row["threshold"] == float(threshold) and row["split"] == split and row["pair_type"] == "between"), None)
            if within and between:
                summary_rows.append({
                    "threshold": float(threshold),
                    "split": split,
                    "pair_type": "effect_between_minus_within",
                    "pair_count": min(within["pair_count"], between["pair_count"]),
                    "mean_tm_score": between["mean_tm_score"] - within["mean_tm_score"],
                    "median_tm_score": between["median_tm_score"] - within["median_tm_score"],
                    "mean_tm_align_rmsd_angstrom": between["mean_tm_align_rmsd_angstrom"] - within["mean_tm_align_rmsd_angstrom"],
                })
    fields = list(output_rows[0].keys()) if output_rows else list(pairs[0].keys()) + ["tm_norm_chain1", "tm_norm_chain2", "tm_score_mean", "tm_align_rmsd_angstrom", "aligned_residue_count", "aligned_length_fraction"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output_rows)
    args.errors_output.parent.mkdir(parents=True, exist_ok=True)
    with args.errors_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["threshold", "split", "pair_type", "left_domain_uid", "right_domain_uid", "error"])
        writer.writeheader()
        writer.writerows(errors)
    summary = {
        "pairs_input": str(args.pairs),
        "topology_input": str(args.topology),
        "pairs_requested": len(pairs),
        "pairs_written": len(output_rows),
        "errors": len(errors),
        "method": "tmtools_TM_align_20210224",
        "status": "pass" if not errors and output_rows else "fail",
        "results": summary_rows,
        "notes": [
            "TM-score is reported normalized by each input chain and as their arithmetic mean.",
            "The pair list is inherited from the endpoint-matched representation check and was not selected using TM-score.",
        ],
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Evaluated {len(output_rows)} TM-align pairs with {len(errors)} errors.")
    return 0 if summary["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())

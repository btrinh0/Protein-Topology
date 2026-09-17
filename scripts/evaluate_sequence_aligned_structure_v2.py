from __future__ import annotations

"""Evaluate sequence-aligned C-alpha similarity for V2 comparison pairs."""

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio.Align import PairwiseAligner
from Bio.PDB import Polypeptide

from compute_domain_topology_v2 import amino_acid_residues, detect_format, load_structure, select_domain_residues


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate sequence-aligned C-alpha RMSD on V2 matched pairs.")
    parser.add_argument("--pairs", type=Path, default=Path("data/processed/v2/structural_coherence_pydssp_pairs_v2.csv"))
    parser.add_argument("--topology", type=Path, default=Path("data/processed/v2/domain_topology_pydssp_calibration_pool_v2.csv"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/sequence_aligned_structure_pydssp_pairs_v2.csv"))
    parser.add_argument("--summary-output", type=Path, default=Path("data/processed/v2/sequence_aligned_structure_pydssp_summary_v2.json"))
    parser.add_argument("--errors-output", type=Path, default=Path("data/processed/v2/sequence_aligned_structure_pydssp_errors_v2.csv"))
    parser.add_argument("--max-pairs", type=int, default=0)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def domain_coordinates(row: dict[str, str]):
    structure_path = Path(row["source_path"])
    structure = load_structure(structure_path, detect_format(structure_path), row["structure_id"])
    model = next(structure.get_models())
    residues = amino_acid_residues(model[row["chain_id"]])
    selected = select_domain_residues(residues, json.loads(row["domain_segments_json"]), row["chain_id"])
    sequence = []
    coordinates = []
    for residue, _ in selected:
        if "CA" not in residue:
            continue
        sequence.append(Polypeptide.protein_letters_3to1.get(residue.get_resname().upper(), "X"))
        coordinates.append([float(value) for value in residue["CA"].coord])
    if len(sequence) < 3:
        raise ValueError("fewer than three C-alpha coordinates")
    return "".join(sequence), np.asarray(coordinates, dtype=np.float64)


def aligned_coordinates(sequence_left: str, coords_left: np.ndarray, sequence_right: str, coords_right: np.ndarray):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(sequence_left, sequence_right)[0]
    pairs = []
    for (left_start, left_end), (right_start, right_end) in zip(*alignment.aligned):
        width = min(int(left_end - left_start), int(right_end - right_start))
        pairs.extend((left_start + offset, right_start + offset) for offset in range(width))
    if len(pairs) < 3:
        raise ValueError("fewer than three aligned residues")
    left_indices, right_indices = zip(*pairs)
    return coords_left[list(left_indices)], coords_right[list(right_indices)], len(pairs)


def kabsch_rmsd(left: np.ndarray, right: np.ndarray) -> float:
    left_centered = left - left.mean(axis=0)
    right_centered = right - right.mean(axis=0)
    covariance = left_centered.T @ right_centered
    u, _, vt = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(vt.T @ u.T)
    rotation = vt.T @ correction @ u.T
    aligned = left_centered @ rotation
    return math.sqrt(float(np.mean(np.sum((aligned - right_centered) ** 2, axis=1))))


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
    pair_rows = read_csv(args.pairs)
    if args.max_pairs > 0:
        pair_rows = pair_rows[: args.max_pairs]
    topology_rows = {row["domain_uid"]: row for row in read_csv(args.topology)}
    bundles: dict[str, tuple[str, np.ndarray]] = {}
    errors: list[dict[str, str]] = []
    output_rows: list[dict[str, str]] = []
    grouped: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    for pair in pair_rows:
        left_uid = pair["left_domain_uid"]
        right_uid = pair["right_domain_uid"]
        try:
            if left_uid not in bundles:
                bundles[left_uid] = domain_coordinates(topology_rows[left_uid])
            if right_uid not in bundles:
                bundles[right_uid] = domain_coordinates(topology_rows[right_uid])
            left_sequence, left_coordinates = bundles[left_uid]
            right_sequence, right_coordinates = bundles[right_uid]
            aligned_left, aligned_right, aligned_count = aligned_coordinates(left_sequence, left_coordinates, right_sequence, right_coordinates)
            rmsd = kabsch_rmsd(aligned_left, aligned_right)
            identity = sum(a == b for a, b in zip(left_sequence, right_sequence)) / max(min(len(left_sequence), len(right_sequence)), 1)
            coverage = aligned_count / max(max(len(left_sequence), len(right_sequence)), 1)
            output_rows.append({
                **pair,
                "aligned_residue_count": str(aligned_count),
                "aligned_length_fraction": f"{coverage:.6f}",
                "sequence_identity_short_alignment": f"{identity:.6f}",
                "kabsch_ca_rmsd_angstrom": f"{rmsd:.6f}",
            })
            grouped[(pair["threshold"], pair["split"], pair["pair_type"])].append(rmsd)
        except Exception as exc:
            errors.append({"threshold": pair.get("threshold", ""), "split": pair.get("split", ""), "pair_type": pair.get("pair_type", ""), "left_domain_uid": left_uid, "right_domain_uid": right_uid, "error": str(exc)})
    summary = []
    for (threshold, split, pair_type), values in sorted(grouped.items()):
        summary.append({
            "threshold": float(threshold),
            "split": split,
            "pair_type": pair_type,
            "pair_count": len(values),
            "mean_kabsch_ca_rmsd_angstrom": mean(values),
            "median_kabsch_ca_rmsd_angstrom": median(values),
        })
    for threshold in sorted({row["threshold"] for row in pair_rows}):
        for split in sorted({row["split"] for row in pair_rows}):
            within = next((row for row in summary if row["threshold"] == float(threshold) and row["split"] == split and row["pair_type"] == "within"), None)
            between = next((row for row in summary if row["threshold"] == float(threshold) and row["split"] == split and row["pair_type"] == "between"), None)
            if within and between:
                summary.append({
                    "threshold": float(threshold),
                    "split": split,
                    "pair_type": "effect_between_minus_within",
                    "pair_count": min(within["pair_count"], between["pair_count"]),
                    "mean_kabsch_ca_rmsd_angstrom": between["mean_kabsch_ca_rmsd_angstrom"] - within["mean_kabsch_ca_rmsd_angstrom"],
                    "median_kabsch_ca_rmsd_angstrom": between["median_kabsch_ca_rmsd_angstrom"] - within["median_kabsch_ca_rmsd_angstrom"],
                })
    fields = list(output_rows[0].keys()) if output_rows else list(pair_rows[0].keys()) + ["aligned_residue_count", "aligned_length_fraction", "sequence_identity_short_alignment", "kabsch_ca_rmsd_angstrom"]
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
    result = {
        "pairs_input": str(args.pairs),
        "topology_input": str(args.topology),
        "pairs_requested": len(pair_rows),
        "pairs_written": len(output_rows),
        "errors": len(errors),
        "method": "global_sequence_alignment_plus_Kabsch_C_alpha_RMSD",
        "status": "sequence_aligned_geometry_proxy_only; TM-align_or_Foldseek_not_run",
        "results": summary,
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Evaluated {len(output_rows)} sequence-aligned structural pairs with {len(errors)} errors.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

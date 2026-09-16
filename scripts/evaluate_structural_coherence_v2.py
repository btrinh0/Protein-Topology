from __future__ import annotations

"""Evaluate held-out topology coherence with normalized C-alpha geometry."""

import argparse
import csv
import json
import math
import random
from collections import defaultdict
from itertools import combinations
from pathlib import Path

from calibrate_topology_signature_v2 import macro_labels, parse_signature
from compute_domain_topology_v2 import amino_acid_residues, detect_format, load_structure, select_domain_residues


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate V2 topology coherence with a coordinate geometry proxy.")
    parser.add_argument("--split", type=Path, default=Path("data/processed/v2/calibration_split_v2.csv"))
    parser.add_argument("--topology", type=Path, default=Path("data/processed/v2/domain_topology_alpha_beta_70_160_calibration_pool.csv"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/structural_coherence_pairs_v2.csv"))
    parser.add_argument("--summary-output", type=Path, default=Path("data/processed/v2/structural_coherence_summary_v2.json"))
    parser.add_argument("--errors-output", type=Path, default=Path("data/processed/v2/structural_coherence_errors_v2.csv"))
    parser.add_argument("--samples", type=int, default=32)
    parser.add_argument("--max-pairs", type=int, default=2500)
    parser.add_argument("--seed", default="tct-v2-structural-coherence-2026-08-26")
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def resampled_geometry(coords: list[tuple[float, float, float]], samples: int) -> tuple[float, ...] | None:
    if len(coords) < 4 or samples < 4:
        return None
    centre = tuple(sum(point[index] for point in coords) / len(coords) for index in range(3))
    centered = [tuple(point[index] - centre[index] for index in range(3)) for point in coords]
    radius = math.sqrt(sum(sum(value * value for value in point) for point in centered) / len(centered))
    if radius <= 1e-8:
        return None
    resampled: list[tuple[float, float, float]] = []
    for sample_index in range(samples):
        position = sample_index * (len(centered) - 1) / (samples - 1)
        left = int(position)
        right = min(left + 1, len(centered) - 1)
        fraction = position - left
        resampled.append(tuple(
            (1.0 - fraction) * centered[left][axis] + fraction * centered[right][axis]
            for axis in range(3)
        ))
    values: list[float] = []
    for left_index in range(samples):
        for right_index in range(left_index + 1, samples):
            delta = tuple(resampled[left_index][axis] - resampled[right_index][axis] for axis in range(3))
            values.append(math.sqrt(sum(value * value for value in delta)) / radius)
    return tuple(values)


def load_descriptor(row: dict[str, str], samples: int):
    structure_path = Path(row["source_path"])
    structure = load_structure(structure_path, detect_format(structure_path), row["structure_id"])
    model = next(structure.get_models())
    chain = model[row["chain_id"]]
    residues = amino_acid_residues(chain)
    segments = json.loads(row["domain_segments_json"])
    selected = select_domain_residues(residues, segments, row["chain_id"])
    coords = [tuple(float(value) for value in residue["CA"].coord) for residue, _ in selected if "CA" in residue]
    return resampled_geometry(coords, samples), len(coords)


def numeric(row: dict[str, str], field: str, default: float = 0.0) -> float:
    try:
        return float(row.get(field, ""))
    except (TypeError, ValueError):
        return default


def covariate_vector(row: dict[str, str]) -> tuple[float, ...]:
    return (
        numeric(row, "resolved_residue_count") / 100.0,
        numeric(row, "sse_count") / 10.0,
        numeric(row, "contact_density"),
        numeric(row, "contact_order"),
    )


def covariate_distance(left: dict[str, str], right: dict[str, str]) -> float:
    values = [a - b for a, b in zip(covariate_vector(left), covariate_vector(right))]
    return math.sqrt(sum(value * value for value in values))


def geometry_distance(left: tuple[float, ...], right: tuple[float, ...]) -> float:
    if len(left) != len(right):
        return 1.0
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right)) / len(left))


def pair_rows(rows: list[dict[str, str]], labels: dict[str, str], max_pairs: int, seed: str) -> tuple[list[tuple[int, int, str, float]], list[tuple[int, int, str, float]]]:
    grouped: dict[str, list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        label = labels.get(row.get("topology_signature_exact", ""), "")
        if label:
            grouped[label].append(index)
    within = [
        (left, right, label, covariate_distance(rows[left], rows[right]))
        for label, members in sorted(grouped.items())
        for left, right in combinations(members, 2)
    ]
    within.sort(key=lambda pair: (pair[3], rows[pair[0]]["domain_uid"], rows[pair[1]]["domain_uid"]))
    within = within[:max_pairs]
    within_labels = {index for pair in within for index in pair[:2]}
    rng = random.Random(seed)
    indices = list(range(len(rows)))
    rng.shuffle(indices)
    between: list[tuple[int, int, str, float]] = []
    for position, left in enumerate(indices):
        left_label = labels.get(rows[left].get("topology_signature_exact", ""), "")
        if not left_label:
            continue
        for right in indices[position + 1 :]:
            right_label = labels.get(rows[right].get("topology_signature_exact", ""), "")
            if right_label and right_label != left_label:
                between.append((left, right, "", covariate_distance(rows[left], rows[right])))
    between.sort(key=lambda pair: (pair[3], rows[pair[0]]["domain_uid"], rows[pair[1]]["domain_uid"]))
    if within_labels:
        target = len(within)
    else:
        target = max_pairs
    return within, between[:target]


def mean(values: list[float]) -> float | None:
    return sum(values) / len(values) if values else None


def median(values: list[float]) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    middle = len(ordered) // 2
    return ordered[middle] if len(ordered) % 2 else (ordered[middle - 1] + ordered[middle]) / 2.0


def main() -> int:
    args = parse_args()
    split_rows = read_csv(args.split)
    topology_rows = {row["domain_uid"]: row for row in read_csv(args.topology)}
    rows = [
        {**row, **topology_rows.get(row["domain_uid"], {})}
        for row in split_rows
        if row.get("eligible_precalibration") == "true"
    ]
    descriptors: dict[str, tuple[float, ...]] = {}
    descriptor_lengths: dict[str, int] = {}
    errors: list[dict[str, str]] = []
    for index, row in enumerate(rows, start=1):
        try:
            descriptor, length = load_descriptor(row, args.samples)
            if descriptor is None:
                raise ValueError("insufficient C-alpha coordinates for geometry descriptor")
            descriptors[row["domain_uid"]] = descriptor
            descriptor_lengths[row["domain_uid"]] = length
        except Exception as exc:
            errors.append({"domain_uid": row["domain_uid"], "structure_id": row["structure_id"], "error": str(exc)})
        if index % 100 == 0 or index == len(rows):
            print(f"Loaded geometry descriptors for {index}/{len(rows)} domains.")

    usable_rows = [row for row in rows if row["domain_uid"] in descriptors]
    signatures = {
        row["topology_signature_exact"]: parse_signature(row["topology_signature_exact"])
        for row in usable_rows
        if row.get("topology_signature_exact")
    }
    pair_output: list[dict[str, str]] = []
    summary: list[dict[str, object]] = []
    for threshold in (0.10, 0.20):
        labels = macro_labels(signatures, threshold)
        for split in ("calibration", "held_out"):
            subset = [row for row in usable_rows if row["split"] == split]
            within, between = pair_rows(subset, labels, args.max_pairs, f"{args.seed}:{threshold}:{split}")
            distances: dict[str, list[float]] = {"within": [], "between": []}
            covariates: dict[str, list[float]] = {"within": [], "between": []}
            for pair_type, pairs in (("within", within), ("between", between)):
                for left_index, right_index, label, cov_distance in pairs:
                    left = subset[left_index]
                    right = subset[right_index]
                    left_descriptor = descriptors[left["domain_uid"]]
                    right_descriptor = descriptors[right["domain_uid"]]
                    structural_distance = geometry_distance(left_descriptor, right_descriptor)
                    distances[pair_type].append(structural_distance)
                    covariates[pair_type].append(cov_distance)
                    pair_output.append({
                        "threshold": f"{threshold:.2f}",
                        "split": split,
                        "pair_type": pair_type,
                        "topology_macroclass": label,
                        "left_domain_uid": left["domain_uid"],
                        "right_domain_uid": right["domain_uid"],
                        "left_structure_id": left["structure_id"],
                        "right_structure_id": right["structure_id"],
                        "left_ca_count": str(descriptor_lengths[left["domain_uid"]]),
                        "right_ca_count": str(descriptor_lengths[right["domain_uid"]]),
                        "covariate_distance": f"{cov_distance:.6f}",
                        "normalized_ca_geometry_distance": f"{structural_distance:.6f}",
                    })
            summary.append({
                "threshold": threshold,
                "split": split,
                "within_pair_count": len(distances["within"]),
                "between_pair_count": len(distances["between"]),
                "within_mean_normalized_ca_geometry_distance": mean(distances["within"]),
                "between_mean_normalized_ca_geometry_distance": mean(distances["between"]),
                "within_median_normalized_ca_geometry_distance": median(distances["within"]),
                "between_median_normalized_ca_geometry_distance": median(distances["between"]),
                "within_mean_covariate_distance": mean(covariates["within"]),
                "between_mean_covariate_distance": mean(covariates["between"]),
                "coherence_effect_between_minus_within": (
                    mean(distances["between"]) - mean(distances["within"])
                    if distances["within"] and distances["between"] else None
                ),
            })

    fields = [
        "threshold", "split", "pair_type", "topology_macroclass", "left_domain_uid", "right_domain_uid",
        "left_structure_id", "right_structure_id", "left_ca_count", "right_ca_count", "covariate_distance",
        "normalized_ca_geometry_distance",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(pair_output)
    args.errors_output.parent.mkdir(parents=True, exist_ok=True)
    with args.errors_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["domain_uid", "structure_id", "error"])
        writer.writeheader()
        writer.writerows(errors)
    output_summary = {
        "split_input": str(args.split),
        "topology_input": str(args.topology),
        "descriptor": "centered_resampled_C_alpha_pairwise_distance_matrix_normalized_by_radius_of_gyration",
        "samples": args.samples,
        "eligible_rows": len(rows),
        "usable_rows": len(usable_rows),
        "descriptor_errors": len(errors),
        "candidate_thresholds": [0.10, 0.20],
        "results": summary,
        "status": "coordinate_geometry_proxy_only; TM-align_or_Foldseek_not_run",
        "interpretation_note": "This descriptor is an independent geometry check, not a replacement for sequence-aligned structural similarity.",
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(output_summary, indent=2) + "\n", encoding="utf-8")
    print(f"Evaluated {len(summary)} structural-coherence conditions from {len(usable_rows)} usable domains.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

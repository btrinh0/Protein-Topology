from __future__ import annotations

"""Test provisional topology signatures under contact-parameter perturbations."""

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

from calibrate_topology_signature_v2 import macro_labels, normalized_mutual_information, parse_signature, signature_distance
from compute_domain_topology_v2 import (
    SSESegment,
    amino_acid_residues,
    contact_pairs,
    detect_format,
    load_structure,
    select_domain_residues,
    topology_signature,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Assess V2 topology signature perturbation stability.")
    parser.add_argument("--split", type=Path, default=Path("data/processed/v2/calibration_split_v2.csv"))
    parser.add_argument("--topology", type=Path, default=Path("data/processed/v2/domain_topology_alpha_beta_70_160_calibration_pool.csv"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/topology_perturbation_v2.csv"))
    parser.add_argument("--summary-output", type=Path, default=Path("data/processed/v2/topology_perturbation_summary_v2.json"))
    parser.add_argument("--errors-output", type=Path, default=Path("data/processed/v2/topology_perturbation_errors_v2.csv"))
    parser.add_argument("--min-sequence-separation", type=int, default=4)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sse_segments(row: dict[str, str]) -> list[SSESegment]:
    return [SSESegment(str(item["state"]), int(item["start_position"]), int(item["end_position"])) for item in json.loads(row["sse_segments_json"])]


def domain_signatures(rows: list[dict[str, str]], distances: tuple[float, ...], supports: tuple[int, ...], min_sequence_separation: int):
    signatures: dict[tuple[float, int], dict[str, str]] = {
        (distance, support): {} for distance in distances for support in supports
    }
    errors: list[dict[str, str]] = []
    rows_by_structure: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        rows_by_structure[row["structure_id"]].append(row)
    for structure_id in sorted(rows_by_structure):
        structure_rows = rows_by_structure[structure_id]
        structure_path = Path(structure_rows[0]["source_path"])
        try:
            structure = load_structure(structure_path, detect_format(structure_path), structure_id)
            model = next(structure.get_models())
        except Exception as exc:
            errors.extend({"domain_uid": row["domain_uid"], "structure_id": structure_id, "error": str(exc)} for row in structure_rows)
            continue
        chain_cache = {chain.id: amino_acid_residues(chain) for chain in model}
        for row in structure_rows:
            try:
                segments = json.loads(row["domain_segments_json"])
                selected = select_domain_residues(chain_cache.get(row["chain_id"], []), segments, row["chain_id"])
                sse = sse_segments(row)
                for distance in distances:
                    contacts, eligible = contact_pairs(selected, distance, min_sequence_separation)
                    eligible_indices = [selected.index(item) for item in eligible]
                    for support in supports:
                        signature, _, status = topology_signature(sse, eligible_indices, contacts, support)
                        signatures[(distance, support)][row["domain_uid"]] = signature if status == "computed" else ""
            except Exception as exc:
                errors.append({"domain_uid": row["domain_uid"], "structure_id": structure_id, "error": str(exc)})
    return signatures, errors


def main() -> int:
    args = parse_args()
    split_rows = read_csv(args.split)
    topology_rows = {row["domain_uid"]: row for row in read_csv(args.topology)}
    rows = [
        {**row, **topology_rows.get(row["domain_uid"], {})}
        for row in split_rows
        if row.get("eligible_precalibration") == "true"
    ]
    distances = (7.0, 8.0, 9.0)
    supports = (1, 2, 3)
    signatures, errors = domain_signatures(rows, distances, supports, args.min_sequence_separation)
    baseline_condition = (8.0, 2)
    baseline = signatures[baseline_condition]
    baseline_parsed = {uid: parse_signature(value) for uid, value in baseline.items() if value}
    output_rows: list[dict[str, str]] = []
    summary: list[dict[str, object]] = []
    for condition in [(distance, support) for distance in distances for support in supports]:
        distance, support = condition
        current = signatures[condition]
        current_parsed = {uid: parse_signature(value) for uid, value in current.items() if value}
        exact_matches = []
        signature_distances = []
        for row in rows:
            uid = row["domain_uid"]
            base_value = baseline.get(uid, "")
            current_value = current.get(uid, "")
            exact_match = bool(base_value and current_value and base_value == current_value)
            exact_matches.append(exact_match)
            distance_value = signature_distance(baseline_parsed[uid], current_parsed[uid]) if uid in baseline_parsed and uid in current_parsed else 1.0
            signature_distances.append(distance_value)
            output_rows.append({
                "domain_uid": uid,
                "split": row["split"],
                "distance_threshold_angstrom": f"{distance:.1f}",
                "min_sse_contact_support": str(support),
                "baseline_signature": base_value,
                "perturbed_signature": current_value,
                "exact_match_baseline": "true" if exact_match else "false",
                "signature_distance_baseline": f"{distance_value:.6f}",
            })
        current_labels = macro_labels({value: parse_signature(value) for value in current.values() if value}, 0.10)
        baseline_labels = macro_labels({value: parse_signature(value) for value in baseline.values() if value}, 0.10)
        for split in ("calibration", "held_out"):
            subset = [row for row in rows if row["split"] == split]
            subset_uids = [row["domain_uid"] for row in subset]
            base_labels = [baseline_labels.get(baseline.get(uid, ""), "") for uid in subset_uids]
            perturbed_labels = [current_labels.get(current.get(uid, ""), "") for uid in subset_uids]
            match_values = [value for row, value in zip(rows, exact_matches) if row["split"] == split]
            distance_values = [value for row, value in zip(rows, signature_distances) if row["split"] == split]
            summary.append({
                "distance_threshold_angstrom": distance,
                "min_sse_contact_support": support,
                "split": split,
                "row_count": len(subset),
                "exact_signature_match_fraction": sum(match_values) / len(match_values) if match_values else 0.0,
                "mean_signature_distance_from_baseline": sum(distance_values) / len(distance_values) if distance_values else None,
                "macroclass_nmi_to_baseline": normalized_mutual_information(base_labels, perturbed_labels),
                "baseline_nonempty_signature_count": sum(bool(baseline.get(uid, "")) for uid in subset_uids),
                "perturbed_nonempty_signature_count": sum(bool(current.get(uid, "")) for uid in subset_uids),
            })
    fields = [
        "domain_uid", "split", "distance_threshold_angstrom", "min_sse_contact_support", "baseline_signature",
        "perturbed_signature", "exact_match_baseline", "signature_distance_baseline",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output_rows)
    args.errors_output.parent.mkdir(parents=True, exist_ok=True)
    with args.errors_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["domain_uid", "structure_id", "error"])
        writer.writeheader()
        writer.writerows(errors)
    output_summary = {
        "split_input": str(args.split),
        "topology_input": str(args.topology),
        "baseline_condition": {"distance_threshold_angstrom": 8.0, "min_sse_contact_support": 2},
        "conditions": [{"distance_threshold_angstrom": distance, "min_sse_contact_support": support} for distance in distances for support in supports],
        "rows": len(rows),
        "errors": len(errors),
        "results": summary,
        "sse_assignment_status": "provisional_author_annotation; secondary-assignment sensitivity remains pending",
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(output_summary, indent=2) + "\n", encoding="utf-8")
    print(f"Evaluated {len(distances) * len(supports)} perturbation conditions for {len(rows)} domains.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

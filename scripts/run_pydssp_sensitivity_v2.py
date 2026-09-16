from __future__ import annotations

"""Run a simplified DSSP-like secondary-structure sensitivity pass for V2."""

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(1, str(Path(__file__).parent.parent / ".deps" / "pydssp-source"))

import pydssp
from calibrate_topology_signature_v2 import parse_signature, signature_distance
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
    parser = argparse.ArgumentParser(description="Compute V2 topology with simplified PyDSSP assignments.")
    parser.add_argument("--split", type=Path, default=Path("data/processed/v2/calibration_split_v2.csv"))
    parser.add_argument("--topology", type=Path, default=Path("data/processed/v2/domain_topology_alpha_beta_70_160_calibration_pool.csv"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/domain_topology_pydssp_sensitivity_v2.csv"))
    parser.add_argument("--summary-output", type=Path, default=Path("data/processed/v2/domain_topology_pydssp_sensitivity_summary_v2.json"))
    parser.add_argument("--errors-output", type=Path, default=Path("data/processed/v2/domain_topology_pydssp_sensitivity_errors_v2.csv"))
    parser.add_argument("--distance-threshold", type=float, default=8.0)
    parser.add_argument("--min-sequence-separation", type=int, default=4)
    parser.add_argument("--min-sse-length", type=int, default=3)
    parser.add_argument("--min-sse-contact-support", type=int, default=2)
    parser.add_argument("--max-domains", type=int, default=0)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def author_states(row: dict[str, str], length: int) -> list[str]:
    states = ["L"] * length
    for item in json.loads(row["sse_segments_json"]):
        for position in range(int(item["start_position"]), int(item["end_position"]) + 1):
            if 0 <= position < length:
                states[position] = str(item["state"])
    return states


def backbone_coordinates(selected):
    positions = []
    coordinates = []
    for position, (residue, _) in enumerate(selected):
        if all(name in residue for name in ("N", "CA", "C", "O")):
            positions.append(position)
            coordinates.append([[float(value) for value in residue[name].coord] for name in ("N", "CA", "C", "O")])
    return positions, np.asarray(coordinates, dtype=np.float64)


def pydssp_states(selected: list[tuple[object, object]]) -> tuple[list[str], int]:
    positions, coordinates = backbone_coordinates(selected)
    if len(positions) < 8:
        raise ValueError("fewer than eight complete backbone residues")
    assigned = pydssp.assign(coordinates, out_type="c3")
    states = ["L"] * len(selected)
    for position, state in zip(positions, assigned):
        value = str(state)
        states[position] = value if value in {"H", "E"} else "L"
    return states, len(positions)


def states_to_segments(states: list[str], min_length: int) -> list[SSESegment]:
    segments: list[SSESegment] = []
    index = 0
    while index < len(states):
        if states[index] not in {"H", "E"}:
            index += 1
            continue
        end = index
        while end + 1 < len(states) and states[end + 1] == states[index]:
            end += 1
        if end - index + 1 >= min_length:
            segments.append(SSESegment(states[index], index, end))
        index = end + 1
    return segments


def main() -> int:
    args = parse_args()
    split_rows = read_csv(args.split)
    topology_rows = {row["domain_uid"]: row for row in read_csv(args.topology)}
    rows = [
        {**row, **topology_rows.get(row["domain_uid"], {})}
        for row in split_rows
        if row.get("eligible_precalibration") == "true"
    ]
    if args.max_domains > 0:
        rows = rows[: args.max_domains]
    rows_by_structure: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        rows_by_structure[row["structure_id"]].append(row)
    output_rows: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
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
                states, complete_backbone_count = pydssp_states(selected)
                sse = states_to_segments(states, args.min_sse_length)
                contacts, eligible = contact_pairs(selected, args.distance_threshold, args.min_sequence_separation)
                eligible_indices = [selected.index(item) for item in eligible]
                signature, edge_count, signature_status = topology_signature(
                    sse, eligible_indices, contacts, args.min_sse_contact_support
                )
                author_signature = row.get("topology_signature_exact", "")
                author_distance = signature_distance(parse_signature(author_signature), parse_signature(signature)) if author_signature and signature else 1.0
                author_full_states = author_states(row, len(selected))
                state_agreement = sum(left == right for left, right in zip(author_full_states, states)) / max(len(states), 1)
                candidate_status = "pydssp_alpha_beta_candidate" if any(segment.state == "H" for segment in sse) and any(segment.state == "E" for segment in sse) else "not_alpha_beta_by_pydssp"
                output_rows.append({
                    "domain_uid": row["domain_uid"],
                    "source_domain_id": row["source_domain_id"],
                    "structure_id": structure_id,
                    "chain_id": row["chain_id"],
                    "split": row["split"],
                    "sse_assignment_method": "pydssp_simplified_c3",
                    "sse_assignment_status": "sensitivity_only",
                    "analysis_candidate_status": candidate_status,
                    "complete_backbone_residue_count": str(complete_backbone_count),
                    "backbone_coverage_fraction": f"{complete_backbone_count / max(len(selected), 1):.6f}",
                    "sse_count": str(len(sse)),
                    "sse_order": "-".join(segment.state for segment in sse),
                    "sse_segments_json": json.dumps([segment.__dict__ for segment in sse], separators=(",", ":")),
                    "topology_signature_exact": signature if signature_status == "computed" else "",
                    "sse_contact_edge_count": str(edge_count),
                    "author_signature_distance": f"{author_distance:.6f}",
                    "author_state_agreement_fraction": f"{state_agreement:.6f}",
                    "contact_definition_id": row.get("contact_definition_id", ""),
                    "quality_flags": "simplified_pydssp_sensitivity;not_primary_assignment",
                })
            except Exception as exc:
                errors.append({"domain_uid": row["domain_uid"], "structure_id": structure_id, "error": str(exc)})
    fields = list(output_rows[0].keys()) if output_rows else ["domain_uid", "structure_id"]
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
    agreement_values = [float(row["author_state_agreement_fraction"]) for row in output_rows]
    signature_distances = [float(row["author_signature_distance"]) for row in output_rows]
    summary = {
        "split_input": str(args.split),
        "topology_input": str(args.topology),
        "rows_requested": len(rows),
        "rows_written": len(output_rows),
        "errors": len(errors),
        "assignment_method": "pydssp_simplified_c3",
        "assignment_role": "sensitivity_only;not_primary_assignment",
        "state_agreement_mean": sum(agreement_values) / len(agreement_values) if agreement_values else None,
        "state_agreement_median": sorted(agreement_values)[len(agreement_values) // 2] if agreement_values else None,
        "signature_distance_mean": sum(signature_distances) / len(signature_distances) if signature_distances else None,
        "pydssp_alpha_beta_candidates": sum(row["analysis_candidate_status"] == "pydssp_alpha_beta_candidate" for row in output_rows),
        "notes": [
            "PyDSSP is a simplified C3 implementation and is used only for sensitivity analysis.",
            "Missing backbone atoms are excluded from the assignment and recorded as coverage.",
            "A final primary assignment still requires a declared DSSP or STRIDE executable and reference validation.",
        ],
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(output_rows)} PyDSSP sensitivity rows with {len(errors)} errors.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

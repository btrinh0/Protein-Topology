from __future__ import annotations

"""Build a V2 topology table from cached PDB-REDO DSSP annotations."""

import argparse
import csv
import hashlib
import json
import sys
import urllib.error
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from Bio.PDB.DSSP import make_dssp_dict

sys.path.insert(0, str(Path(__file__).parent))

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


PDB_REDO_URL = "https://pdb-redo.eu/dssp/db/{structure_id}/legacy"
STATE_MAP = {"H": "H", "G": "H", "I": "H", "E": "E", "B": "E"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute V2 topology with precomputed PDB-REDO DSSP states.")
    parser.add_argument("--split", type=Path, default=Path("data/processed/v2/calibration_split_v2.csv"))
    parser.add_argument("--topology", type=Path, default=Path("data/processed/v2/domain_topology_alpha_beta_70_160_calibration_pool.csv"))
    parser.add_argument("--cache-dir", type=Path, default=Path("data/raw/v2/pdbredo_dssp"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/domain_topology_pdbredo_dssp_v2.csv"))
    parser.add_argument("--summary-output", type=Path, default=Path("data/processed/v2/domain_topology_pdbredo_dssp_summary_v2.json"))
    parser.add_argument("--errors-output", type=Path, default=Path("data/processed/v2/domain_topology_pdbredo_dssp_errors_v2.csv"))
    parser.add_argument("--distance-threshold", type=float, default=8.0)
    parser.add_argument("--min-sequence-separation", type=int, default=4)
    parser.add_argument("--min-sse-length", type=int, default=3)
    parser.add_argument("--min-sse-contact-support", type=int, default=2)
    parser.add_argument("--max-structures", type=int, default=0)
    parser.add_argument("--download", action="store_true")
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


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


def cache_path(cache_dir: Path, structure_id: str) -> Path:
    return cache_dir / f"{structure_id.lower()}.dssp"


def download_annotation(structure_id: str, destination: Path) -> dict[str, object]:
    url = PDB_REDO_URL.format(structure_id=structure_id.lower())
    request = urllib.request.Request(url, headers={"User-Agent": "TCT-v2/1.0"})
    with urllib.request.urlopen(request, timeout=90) as response:
        payload = response.read()
        content_type = response.headers.get("Content-Type", "")
    if not payload or b"TOTAL NUMBER OF RESIDUES" not in payload:
        raise ValueError("PDB-REDO response is not a legacy DSSP file")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(payload)
    metadata = {
        "structure_id": structure_id,
        "url": url,
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "bytes": len(payload),
        "content_type": content_type,
    }
    destination.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return metadata


def load_annotation(structure_id: str, cache_dir: Path, download: bool) -> tuple[dict[tuple[str, tuple[str, int, str]], tuple], dict[str, object]]:
    path = cache_path(cache_dir, structure_id)
    metadata_path = path.with_suffix(".json")
    if not path.exists():
        if not download:
            raise FileNotFoundError(f"missing cached PDB-REDO annotation: {path}")
        metadata = download_annotation(structure_id, path)
    else:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8")) if metadata_path.exists() else {
            "structure_id": structure_id,
            "url": PDB_REDO_URL.format(structure_id=structure_id.lower()),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "bytes": path.stat().st_size,
        }
    dssp_values, _ = make_dssp_dict(str(path))
    return dssp_values, metadata


def dssp_states(selected, dssp_values) -> tuple[list[str], int, int]:
    states: list[str] = []
    matched = 0
    for residue, _ in selected:
        key = (residue.get_parent().id, residue.id)
        value = dssp_values.get(key)
        if value is None:
            states.append("L")
            continue
        matched += 1
        states.append(STATE_MAP.get(str(value[1]).upper(), "L"))
    return states, matched, len(selected)


def main() -> int:
    args = parse_args()
    split_rows = read_csv(args.split)
    topology_rows = {row["domain_uid"]: row for row in read_csv(args.topology)}
    rows = [
        {**topology_rows.get(row["domain_uid"], {}), **row}
        for row in split_rows
        if row.get("eligible_precalibration") == "true" and row["domain_uid"] in topology_rows
    ]
    rows_by_structure: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        rows_by_structure[row["structure_id"]].append(row)
    structure_ids = sorted(rows_by_structure)
    if args.max_structures > 0:
        structure_ids = structure_ids[: args.max_structures]
    output_rows: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
    annotation_metadata: list[dict[str, object]] = []
    for structure_id in structure_ids:
        structure_rows = rows_by_structure[structure_id]
        try:
            dssp_values, metadata = load_annotation(structure_id, args.cache_dir, args.download)
            structure_path = Path(structure_rows[0]["source_path"])
            structure = load_structure(structure_path, detect_format(structure_path), structure_id)
            model = next(structure.get_models())
        except (OSError, ValueError, urllib.error.URLError, urllib.error.HTTPError) as exc:
            errors.extend({"domain_uid": row["domain_uid"], "structure_id": structure_id, "error": str(exc)} for row in structure_rows)
            continue
        annotation_metadata.append(metadata)
        chain_cache = {chain.id: amino_acid_residues(chain) for chain in model}
        for row in structure_rows:
            try:
                segments = json.loads(row["domain_segments_json"])
                selected = select_domain_residues(chain_cache.get(row["chain_id"], []), segments, row["chain_id"])
                states, matched_count, selected_count = dssp_states(selected, dssp_values)
                sse = states_to_segments(states, args.min_sse_length)
                contacts, eligible = contact_pairs(selected, args.distance_threshold, args.min_sequence_separation)
                eligible_indices = [selected.index(item) for item in eligible]
                signature, edge_count, signature_status = topology_signature(
                    sse, eligible_indices, contacts, args.min_sse_contact_support
                )
                author_signature = row.get("topology_signature_exact", "")
                author_distance = signature_distance(parse_signature(author_signature), parse_signature(signature)) if author_signature and signature else 1.0
                coverage = matched_count / max(selected_count, 1)
                candidate_status = "pdbredo_alpha_beta_candidate" if coverage == 1.0 and any(segment.state == "H" for segment in sse) and any(segment.state == "E" for segment in sse) else "not_alpha_beta_by_pdbredo"
                merged = dict(row)
                merged.update({
                    "analysis_candidate_status": candidate_status,
                    "sse_assignment_method": "pdbredo_dssp_4.6.1",
                    "sse_assignment_status": "precomputed_primary_candidate" if coverage == 1.0 else "partial_assignment",
                    "sse_count": str(len(sse)),
                    "sse_order": "-".join(segment.state for segment in sse),
                    "sse_segments_json": json.dumps([segment.__dict__ for segment in sse], separators=(",", ":")),
                    "topology_signature_exact": signature if signature_status == "computed" else "",
                    "sse_contact_edge_count": str(edge_count),
                    "dssp_residue_count": str(matched_count),
                    "dssp_coverage_fraction": f"{coverage:.6f}",
                    "author_signature_distance": f"{author_distance:.6f}",
                    "quality_flags": f"{row.get('quality_flags', '')};pdbredo_dssp_precomputed;source={PDB_REDO_URL.format(structure_id=structure_id.lower())}",
                })
                output_rows.append(merged)
            except Exception as exc:
                errors.append({"domain_uid": row["domain_uid"], "structure_id": structure_id, "error": str(exc)})
    fields = list(output_rows[0].keys()) if output_rows else list(topology_rows[next(iter(topology_rows))].keys())
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
    coverage_values = [float(row["dssp_coverage_fraction"]) for row in output_rows if row.get("dssp_coverage_fraction")]
    summary = {
        "split_input": str(args.split),
        "topology_input": str(args.topology),
        "cache_dir": str(args.cache_dir),
        "rows_requested": len(rows),
        "structures_requested": len(structure_ids),
        "rows_written": len(output_rows),
        "errors": len(errors),
        "assignment_method": "pdbredo_dssp_4.6.1",
        "assignment_role": "precomputed_primary_candidate",
        "structures_with_annotations": len(annotation_metadata),
        "complete_domain_coverage_fraction": sum(value == 1.0 for value in coverage_values) / max(len(coverage_values), 1),
        "domain_coverage_mean": sum(coverage_values) / len(coverage_values) if coverage_values else None,
        "pdbredo_alpha_beta_candidates": sum(row["analysis_candidate_status"] == "pdbredo_alpha_beta_candidate" for row in output_rows),
        "annotation_metadata": annotation_metadata,
        "notes": [
            "PDB-REDO DSSP legacy annotations are cached byte-for-byte with URL, retrieval time, and SHA-256 metadata.",
            "DSSP H/G/I states are reduced to H and E/B states to E; all other states are loop.",
            "Domains with incomplete residue coverage remain visible but are not primary alpha/beta candidates.",
            "Structural similarity validation remains a separate gate; this source does not replace TM-align or Foldseek.",
        ],
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(output_rows)} PDB-REDO DSSP rows with {len(errors)} errors.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

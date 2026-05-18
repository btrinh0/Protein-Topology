from __future__ import annotations

"""Compare the V2 audit sample against the frozen raw classification sources."""

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

from build_domain_mapping_v2 import parse_cath_boundaries, parse_scope_cla, parse_segment_string


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare raw-source evidence for the V2 mapping audit.")
    parser.add_argument("--audit", type=Path, default=Path("data/processed/v2/domain_mapping_audit_sample_v2.csv"))
    parser.add_argument("--classifications-dir", type=Path, default=Path("data/raw/classifications"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/domain_mapping_audit_raw_evidence_v2.csv"))
    parser.add_argument("--summary-output", type=Path, default=Path("data/processed/v2/domain_mapping_audit_raw_evidence_summary_v2.json"))
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def segment_key(segment) -> tuple:
    return (
        segment.chain_id,
        segment.start_resseq,
        segment.start_icode,
        segment.end_resseq,
        segment.end_icode,
        segment.is_full_chain,
    )


def json_segment_key(item: dict[str, object]) -> tuple:
    return (
        str(item.get("chain_id", "")),
        item.get("start_resseq"),
        str(item.get("start_icode", "")),
        item.get("end_resseq"),
        str(item.get("end_icode", "")),
        bool(item.get("is_full_chain", False)),
    )


def parse_cath_list(path: Path, needed: set[str]) -> dict[str, tuple[str, str, str, str]]:
    result: dict[str, tuple[str, str, str, str]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            if len(fields) >= 5 and fields[0] in needed:
                result[fields[0]] = tuple(fields[1:5])
    return result


def parse_ecod_target_candidates(path: Path, targets: set[tuple[str, str]]) -> dict[tuple[str, str], set[str]]:
    result: dict[tuple[str, str], set[str]] = defaultdict(set)
    header: list[str] | None = None
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.lstrip("#").strip()
            if header is None:
                if stripped.startswith("uid\t"):
                    header = stripped.split("\t")
                continue
            if line.startswith("#") or not line.strip():
                continue
            values = line.rstrip("\n").split("\t")
            if len(values) < len(header):
                continue
            row = dict(zip(header, values))
            pdb_id = row.get("pdb", "").lower()
            segments, _ = parse_segment_string(row.get("pdb_range", ""))
            chain_ids = {segment.chain_id for segment in segments}
            fallback_chain = row.get("chain", "").strip()
            if not chain_ids and fallback_chain:
                chain_ids = {fallback_chain}
            for chain_id in chain_ids:
                if (pdb_id, chain_id) in targets:
                    result[(pdb_id, chain_id)].add(row.get("ecod_domain_id", ""))
    return result


def main() -> int:
    args = parse_args()
    audit_rows = read_csv(args.audit)
    classifications = args.classifications_dir
    scope_rows = {row["source_domain_id"]: row for row in parse_scope_cla(classifications / "scope_cla.txt")}
    cath_boundaries = parse_cath_boundaries(classifications / "cath_domain_boundaries_v4_3_0.txt")
    cath_needed = {row["cath_domain_id"] for row in audit_rows if row.get("cath_domain_id")}
    cath_hierarchy = parse_cath_list(classifications / "cath_domain_list_v4_3_0.txt", cath_needed)
    ecod_targets = {(row["structure_id"].lower(), row["chain_id"]) for row in audit_rows}
    ecod_candidates = parse_ecod_target_candidates(classifications / "ecod_domains.txt", ecod_targets)
    output_rows: list[dict[str, str]] = []
    for row in audit_rows:
        source = scope_rows.get(row["source_domain_id"])
        requested_segments = [json_segment_key(item) for item in json.loads(row["domain_segments_json"])]
        raw_scope_segments = [segment_key(segment) for segment in source["segments"]] if source else []
        scope_match = bool(source and source["structure_id"] == row["structure_id"] and raw_scope_segments == requested_segments)
        raw_cath = cath_boundaries.get((row["structure_id"].lower(), row["chain_id"]), [])
        cath_ids = {
            candidate["candidate_id"][:-2] + "00" if candidate.get("is_whole_chain") else candidate["candidate_id"]
            for candidate in raw_cath
        }
        cath_id = row.get("cath_domain_id", "")
        cath_present = not cath_id or cath_id in cath_ids
        ecod_id = row.get("ecod_domain_id", "")
        ecod_present = not ecod_id or ecod_id in ecod_candidates.get((row["structure_id"].lower(), row["chain_id"]), set())
        hierarchy = cath_hierarchy.get(cath_id)
        hierarchy_match = not cath_id or (hierarchy is not None and hierarchy == (
            row.get("cath_class_id", ""),
            row.get("cath_architecture_id", "").split(".")[-1] if row.get("cath_architecture_id", "") else "",
            row.get("cath_topology_id", "").split(".")[-1] if row.get("cath_topology_id", "") else "",
            row.get("cath_homology_id", "").split(".")[-1] if row.get("cath_homology_id", "") else "",
        ))
        output_rows.append({
            "audit_id": row["audit_id"],
            "domain_uid": row["domain_uid"],
            "scope_raw_row_found": "true" if source else "false",
            "scope_raw_boundary_matches_v2": "true" if scope_match else "false",
            "cath_raw_candidate_count": str(len(raw_cath)),
            "cath_selected_id_present_in_raw_boundaries": "true" if cath_present else "false",
            "ecod_raw_candidate_count_for_chain": str(len(ecod_candidates.get((row["structure_id"].lower(), row["chain_id"]), set()))),
            "ecod_selected_id_present_in_raw_domains": "true" if ecod_present else "false",
            "cath_hierarchy_present_in_static_list": "true" if hierarchy is not None or not cath_id else "false",
            "cath_hierarchy_matches_v2_fields": "true" if hierarchy_match else "false",
            "manual_audit_status": "pending",
        })
    fields = list(output_rows[0].keys()) if output_rows else ["audit_id", "domain_uid"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output_rows)
    summary = {
        "audit_input": str(args.audit),
        "rows": len(output_rows),
        "scope_boundary_match": sum(row["scope_raw_boundary_matches_v2"] == "true" for row in output_rows),
        "cath_selected_id_presence": sum(row["cath_selected_id_present_in_raw_boundaries"] == "true" for row in output_rows),
        "ecod_selected_id_presence": sum(row["ecod_selected_id_present_in_raw_domains"] == "true" for row in output_rows),
        "cath_hierarchy_match": sum(row["cath_hierarchy_matches_v2_fields"] == "true" for row in output_rows),
        "status": "raw_source_precheck_only;manual_review_required",
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Compared {len(output_rows)} audit rows against raw sources; manual status remains pending.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

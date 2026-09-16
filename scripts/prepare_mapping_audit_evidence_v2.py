from __future__ import annotations

"""Prepare objective coordinate and candidate evidence for the V2 audit sample."""

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

from compute_domain_topology_v2 import amino_acid_residues, detect_format, load_structure, select_domain_residues


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare non-authoritative evidence for the V2 mapping audit sample.")
    parser.add_argument("--audit", type=Path, default=Path("data/processed/v2/domain_mapping_audit_sample_v2.csv"))
    parser.add_argument("--chain-manifest", type=Path, default=Path("data/processed/chain_manifest.csv"))
    parser.add_argument("--candidates", type=Path, default=Path("data/processed/v2/domain_mapping_candidates_v2.csv"))
    parser.add_argument("--mapping", type=Path, default=Path("data/processed/v2/domain_mapping_v2.csv"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/domain_mapping_audit_coordinate_evidence_v2.csv"))
    parser.add_argument("--summary-output", type=Path, default=Path("data/processed/v2/domain_mapping_audit_evidence_summary_v2.json"))
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def span_count(segments: list[dict[str, object]], chain_residue_count: int) -> int:
    if any(segment.get("is_full_chain") for segment in segments):
        return chain_residue_count
    return sum(abs(int(segment["end_resseq"]) - int(segment["start_resseq"])) + 1 for segment in segments)


def main() -> int:
    args = parse_args()
    audit_rows = read_csv(args.audit)
    chain_manifest = {(row["structure_id"], row["chain_id"]): row for row in read_csv(args.chain_manifest)}
    mapping_rows = {row["domain_uid"]: row for row in read_csv(args.mapping)}
    candidate_rows: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in read_csv(args.candidates):
        candidate_rows[(row["domain_uid"], row["candidate_source"])].append(row)
    evidence_rows: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
    for row in audit_rows:
        manifest = chain_manifest.get((row["structure_id"], row["chain_id"]))
        selected_count = 0
        chain_count = 0
        coordinate_status = "unavailable"
        coverage = ""
        try:
            if not manifest:
                raise ValueError("chain_not_in_manifest")
            structure_path = Path(manifest["source_path"])
            structure = load_structure(structure_path, detect_format(structure_path), row["structure_id"])
            model = next(structure.get_models())
            chain = model[row["chain_id"]]
            residues = amino_acid_residues(chain)
            segments = json.loads(row["domain_segments_json"])
            selected_count = len(select_domain_residues(residues, segments, row["chain_id"]))
            chain_count = len(residues)
            span = span_count(segments, chain_count)
            coverage = f"{selected_count / span:.6f}" if span else ""
            coordinate_status = "pass" if selected_count > 0 and span > 0 and 0.0 < selected_count / span <= 1.05 else "review"
        except Exception as exc:
            errors.append({"audit_id": row["audit_id"], "domain_uid": row["domain_uid"], "error": str(exc)})
        mapping = mapping_rows.get(row["domain_uid"], {})
        cath = candidate_rows.get((row["domain_uid"], "CATH"), [])
        ecod = candidate_rows.get((row["domain_uid"], "ECOD"), [])
        cath_selected = [candidate["candidate_id"] for candidate in cath if candidate.get("candidate_selected") == "true"]
        ecod_selected = [candidate["candidate_id"] for candidate in ecod if candidate.get("candidate_selected") == "true"]
        evidence_rows.append({
            "audit_id": row["audit_id"],
            "domain_uid": row["domain_uid"],
            "structure_id": row["structure_id"],
            "chain_id": row["chain_id"],
            "source_boundary_status": row["source_boundary_status"],
            "coordinate_chain_residue_count": str(chain_count),
            "coordinate_selected_residue_count": str(selected_count),
            "coordinate_coverage_fraction": coverage,
            "coordinate_slice_precheck": coordinate_status,
            "cath_candidate_count": str(len(cath)),
            "cath_selected_candidate_ids": ";".join(cath_selected),
            "cath_selected_id_present_in_candidates": "true" if not row.get("cath_domain_id") or row["cath_domain_id"] in [candidate["candidate_id"] for candidate in cath] else "false",
            "ecod_candidate_count": str(len(ecod)),
            "ecod_selected_candidate_ids": ";".join(ecod_selected),
            "ecod_selected_id_present_in_candidates": "true" if not row.get("ecod_domain_id") or row["ecod_domain_id"] in [candidate["candidate_id"] for candidate in ecod] else "false",
            "mapping_status_from_v2": mapping.get("mapping_status", ""),
            "manual_audit_status": "pending",
        })
    fields = list(evidence_rows[0].keys()) if evidence_rows else ["audit_id", "domain_uid"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(evidence_rows)
    summary = {
        "audit_input": str(args.audit),
        "rows": len(evidence_rows),
        "coordinate_precheck_pass": sum(row["coordinate_slice_precheck"] == "pass" for row in evidence_rows),
        "coordinate_precheck_review_or_unavailable": sum(row["coordinate_slice_precheck"] != "pass" for row in evidence_rows),
        "cath_selected_id_presence_failures": sum(row["cath_selected_id_present_in_candidates"] == "false" for row in evidence_rows),
        "ecod_selected_id_presence_failures": sum(row["ecod_selected_id_present_in_candidates"] == "false" for row in evidence_rows),
        "errors": len(errors),
        "status": "objective_precheck_only; manual_review_required",
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Prepared objective evidence for {len(evidence_rows)} audit rows; manual status remains pending.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

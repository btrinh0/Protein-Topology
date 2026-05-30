from __future__ import annotations

"""Record explicit handling for domains without PDB-REDO DSSP records."""

import argparse
import csv
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Review PDB-REDO-unavailable domains without substituting author SSE states.")
    parser.add_argument("--errors", type=Path, default=Path("data/processed/v2/domain_topology_pdbredo_dssp_errors_v2.csv"))
    parser.add_argument("--topology", type=Path, default=Path("data/processed/v2/domain_topology_alpha_beta_70_160_calibration_pool.csv"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/pdbredo_unavailable_entry_review_v2.csv"))
    parser.add_argument("--summary-output", type=Path, default=Path("data/processed/v2/pdbredo_unavailable_entry_review_summary_v2.json"))
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    args = parse_args()
    topology = {row["domain_uid"]: row for row in read_csv(args.topology)}
    review_rows: list[dict[str, str]] = []
    for error in read_csv(args.errors):
        row = topology.get(error["domain_uid"], {})
        review_rows.append({
            "domain_uid": error["domain_uid"],
            "structure_id": error["structure_id"],
            "chain_id": row.get("chain_id", ""),
            "source_path": row.get("source_path", ""),
            "resolved_residue_count": row.get("resolved_residue_count", ""),
            "author_sse_order": row.get("sse_order", ""),
            "author_signature_present": str(bool(row.get("topology_signature_exact", ""))).lower(),
            "source_error": error["error"],
            "review_status": "unresolved_primary_sse",
            "calibration_action": "exclude_until_declared_primary_assignment_available",
        })
    fields = list(review_rows[0].keys()) if review_rows else ["domain_uid", "structure_id", "review_status"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(review_rows)
    summary = {
        "rows_reviewed": len(review_rows),
        "structures_reviewed": len({row["structure_id"] for row in review_rows}),
        "status": "explicit_exclusion_pending_primary_assignment",
        "author_annotations_used_as_fallback": False,
        "calibration_action": "exclude_until_declared_primary_assignment_available",
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Recorded explicit review for {len(review_rows)} unavailable domains.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

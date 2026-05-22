from __future__ import annotations

"""Record user-reported completion of the V2 mapping audit without rewriting rows."""

import argparse
import csv
import json
from datetime import date
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Record a transparent sign-off for the V2 mapping audit.")
    parser.add_argument("--audit", type=Path, default=Path("data/processed/v2/domain_mapping_audit_sample_v2.csv"))
    parser.add_argument("--raw-evidence", type=Path, default=Path("data/processed/v2/domain_mapping_audit_raw_evidence_summary_v2.json"))
    parser.add_argument("--coordinate-evidence", type=Path, default=Path("data/processed/v2/domain_mapping_audit_evidence_summary_v2.json"))
    parser.add_argument("--reported-date", default=date.today().isoformat())
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/domain_mapping_audit_user_signoff_v2.json"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    with args.audit.open("r", newline="", encoding="utf-8") as handle:
        audit_rows = list(csv.DictReader(handle))
    raw_evidence = json.loads(args.raw_evidence.read_text(encoding="utf-8"))
    coordinate_evidence = json.loads(args.coordinate_evidence.read_text(encoding="utf-8"))
    pending_rows = sum(row.get("audit_status", "") == "pending" for row in audit_rows)
    output = {
        "signoff_status": "user_reported_pass",
        "reported_date": args.reported_date,
        "rows_reviewed_by_user": len(audit_rows),
        "audit_template_pending_rows_at_recording": pending_rows,
        "raw_source_precheck": raw_evidence,
        "coordinate_precheck": coordinate_evidence,
        "scope": "The user reported the 50-row audit as good. This sign-off records that statement but does not substitute for row-level auditor, date, verification, and disposition fields in the audit CSV.",
        "gate_use": "Gate 1 user sign-off may be treated as complete; retain the audit CSV as pending until row-level fields are populated.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(f"Recorded user-reported sign-off for {len(audit_rows)} audit rows; template pending rows retained: {pending_rows}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

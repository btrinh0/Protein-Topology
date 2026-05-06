from __future__ import annotations

"""Build the quality-filtered V2 calibration pool without changing raw outputs."""

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Filter provisional topology rows for calibration.")
    parser.add_argument("--input", type=Path, default=Path("data/processed/v2/domain_topology_alpha_beta_70_160_provisional.csv"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/domain_topology_alpha_beta_70_160_calibration_pool.csv"))
    parser.add_argument("--summary-output", type=Path, default=Path("data/processed/v2/domain_topology_calibration_pool_summary_v2.json"))
    parser.add_argument("--min-coverage", type=float, default=0.80)
    parser.add_argument("--max-coverage", type=float, default=1.00)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    args = parse_args()
    rows = read_csv(args.input)
    selected: list[dict[str, str]] = []
    exclusions: Counter[str] = Counter()
    for row in rows:
        reasons: list[str] = []
        if row.get("topology_status") != "computed":
            reasons.append("topology_not_computed")
        if not row.get("topology_signature_exact"):
            reasons.append("missing_exact_signature")
        try:
            coverage = float(row["coordinate_coverage_fraction"])
            if coverage < args.min_coverage:
                reasons.append("low_coordinate_coverage")
            if coverage > args.max_coverage:
                reasons.append("coordinate_coverage_above_one")
        except (KeyError, TypeError, ValueError):
            reasons.append("invalid_coordinate_coverage")
        try:
            resolved_count = int(row["resolved_residue_count"])
            if not 70 <= resolved_count <= 160:
                reasons.append("resolved_length_outside_70_160")
        except (KeyError, TypeError, ValueError):
            reasons.append("invalid_resolved_length")
        if "H" not in row.get("sse_order", "") or "E" not in row.get("sse_order", ""):
            reasons.append("not_alpha_beta_by_provisional_sse")
        if reasons:
            exclusions.update(reasons)
            continue
        selected.append(row)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(selected)
    summary = {
        "input": str(args.input),
        "output": str(args.output),
        "input_rows": len(rows),
        "selected_rows": len(selected),
        "selected_structures": len({row["structure_id"] for row in selected}),
        "excluded_rows": len(rows) - len(selected),
        "exclusion_reason_counts": dict(sorted(exclusions.items())),
        "coverage_window": [args.min_coverage, args.max_coverage],
        "length_window": [70, 160],
        "status": "provisional_calibration_pool; author_SSE_assignment_and_macroclass_pending",
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(selected)} calibration-pool rows; excluded {len(rows) - len(selected)} rows.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations

"""Validate V2 topology outputs before representation calibration."""

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate domain-sliced V2 topology output.")
    parser.add_argument("--topology", type=Path, default=Path("data/processed/v2/domain_topology_alpha_beta_70_160_provisional.csv"))
    parser.add_argument("--errors", type=Path, default=Path("data/processed/v2/domain_topology_alpha_beta_70_160_errors.csv"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/domain_topology_validation_v2.json"))
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    args = parse_args()
    rows = read_csv(args.topology)
    errors = read_csv(args.errors) if args.errors.exists() else []
    duplicate_domains = len(rows) - len({row["domain_uid"] for row in rows})
    missing_signatures = sum(not row.get("topology_signature_exact") for row in rows if row.get("topology_status") == "computed")
    coverage_values = []
    invalid_coverage = 0
    low_coverage = 0
    for row in rows:
        try:
            value = float(row["coordinate_coverage_fraction"])
            coverage_values.append(value)
            invalid_coverage += not 0.0 <= value <= 1.0
            low_coverage += value < 0.8
        except (KeyError, TypeError, ValueError):
            invalid_coverage += 1
    result = {
        "topology_input": str(args.topology),
        "rows": len(rows),
        "error_rows": len(errors),
        "duplicate_domain_rows": duplicate_domains,
        "topology_status": dict(sorted(Counter(row.get("topology_status", "") for row in rows).items())),
        "analysis_candidate_status": dict(sorted(Counter(row.get("analysis_candidate_status", "") for row in rows).items())),
        "sse_assignment_status": dict(sorted(Counter(row.get("sse_assignment_status", "") for row in rows).items())),
        "macroclass_status": dict(sorted(Counter(row.get("macroclass_status", "") for row in rows).items())),
        "missing_exact_signatures": missing_signatures,
        "invalid_coordinate_coverage_rows": invalid_coverage,
        "low_coordinate_coverage_rows_below_0_8": low_coverage,
        "coverage_range": [min(coverage_values), max(coverage_values)] if coverage_values else [],
        "gate_status": "not_ready_for_occupancy" if (errors or duplicate_domains or missing_signatures or invalid_coverage) else "ready_for_calibration_review",
        "notes": [
            "Author secondary-structure annotations remain provisional.",
            "Coordinate coverage outside [0, 1] is retained as a diagnostic and must be resolved before final analysis.",
            "A blank topology_macroclass is required until calibration is complete.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations

"""Merge simplified PyDSSP sensitivity assignments into a V2 topology table."""

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a non-primary PyDSSP sensitivity topology table.")
    parser.add_argument("--topology", type=Path, default=Path("data/processed/v2/domain_topology_alpha_beta_70_160_calibration_pool.csv"))
    parser.add_argument("--sensitivity", type=Path, default=Path("data/processed/v2/domain_topology_pydssp_sensitivity_v2.csv"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/domain_topology_pydssp_calibration_pool_v2.csv"))
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    args = parse_args()
    topology_rows = read_csv(args.topology)
    sensitivity_rows = {row["domain_uid"]: row for row in read_csv(args.sensitivity)}
    merged: list[dict[str, str]] = []
    for row in topology_rows:
        sensitivity = sensitivity_rows.get(row["domain_uid"])
        if sensitivity is None:
            continue
        merged_row = dict(row)
        for field in (
            "analysis_candidate_status", "sse_assignment_method", "sse_assignment_status", "sse_count",
            "sse_order", "sse_segments_json", "topology_signature_exact", "sse_contact_edge_count",
        ):
            merged_row[field] = sensitivity.get(field, merged_row.get(field, ""))
        merged_row["topology_macroclass"] = ""
        merged_row["macroclass_status"] = "awaiting_calibration"
        merged_row["quality_flags"] = f"{row.get('quality_flags', '')};simplified_pydssp_sensitivity"
        merged.append(merged_row)
    fields = list(merged[0].keys()) if merged else list(topology_rows[0].keys())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(merged)
    print(f"Wrote {len(merged)} PyDSSP sensitivity topology rows.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations

"""Assign calibrated V2 macroclasses only after the calibration gate is open."""

import argparse
import csv
import json
from pathlib import Path

from calibrate_topology_signature_v2 import macro_labels, parse_signature


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the calibrated V2 topology-signature table.")
    parser.add_argument("--topology", type=Path, default=Path("data/processed/v2/domain_topology_alpha_beta_70_160_calibration_pool.csv"))
    parser.add_argument("--calibration-summary", type=Path, default=Path("data/processed/v2/topology_signature_calibration_summary_v2.json"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/topology_signature_v2.csv"))
    parser.add_argument("--threshold", type=float, default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    summary = json.loads(args.calibration_summary.read_text(encoding="utf-8"))
    if not str(summary.get("decision", "")).startswith("assign_macroclass"):
        raise SystemExit("Calibration gate is closed; no macroclass was assigned.")
    if args.threshold is None:
        raise SystemExit("--threshold is required after the calibration gate is opened.")
    with args.topology.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    signatures = {
        row["topology_signature_exact"]: parse_signature(row["topology_signature_exact"])
        for row in rows
        if row.get("topology_signature_exact")
    }
    labels = macro_labels(signatures, args.threshold)
    for row in rows:
        row["topology_macroclass"] = labels.get(row.get("topology_signature_exact", ""), "")
        row["macroclass_status"] = f"calibrated_threshold_{args.threshold:.2f}"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote calibrated signatures for {len(rows)} domains.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

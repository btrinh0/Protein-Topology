from __future__ import annotations

"""Assign calibration-selected candidate macroclasses without opening occupancy analysis."""

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

from calibrate_topology_signature_v2 import macro_labels


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Assign deterministic candidate macroclasses from calibration-only metrics.")
    parser.add_argument("--topology", type=Path, default=Path("data/processed/v2/domain_topology_pdbredo_primary_alpha_beta_v2.csv"))
    parser.add_argument("--calibration", type=Path, default=Path("data/processed/v2/topology_signature_calibration_pdbredo_primary_complete_summary_v2.json"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/domain_topology_pdbredo_calibrated_candidate_v2.csv"))
    parser.add_argument("--summary-output", type=Path, default=Path("data/processed/v2/domain_topology_pdbredo_calibrated_candidate_summary_v2.json"))
    parser.add_argument("--min-recurring-classes", type=int, default=10)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    args = parse_args()
    rows = read_csv(args.topology)
    calibration = json.loads(args.calibration.read_text(encoding="utf-8"))
    candidates = [
        row for row in calibration["candidate_metrics"]
        if row["split"] == "calibration"
        and row["recurring_class_count"] >= args.min_recurring_classes
        and row["within_mean_feature_distance"] < row["between_mean_feature_distance"]
    ]
    if not candidates:
        raise SystemExit("No calibration threshold met the recurring-class and separation requirements.")
    selected = max(candidates, key=lambda row: (row["nmi_scope_fold"] + row["nmi_cath_topology"], row["nmi_scope_fold"], row["nmi_cath_topology"], -row["threshold"]))
    threshold = float(selected["threshold"])
    signatures = {
        row["topology_signature_exact"]: __import__("calibrate_topology_signature_v2").parse_signature(row["topology_signature_exact"])
        for row in rows if row.get("topology_signature_exact")
    }
    labels = macro_labels(signatures, threshold)
    output_rows: list[dict[str, str]] = []
    for row in rows:
        merged = dict(row)
        merged["topology_macroclass"] = labels.get(row.get("topology_signature_exact", ""), "")
        merged["macroclass_status"] = "calibrated_candidate"
        merged["macroclass_selection_threshold"] = f"{threshold:.2f}"
        merged["quality_flags"] = f"{row.get('quality_flags', '')};macroclass_calibration_candidate"
        output_rows.append(merged)
    fields = list(output_rows[0].keys()) if output_rows else list(rows[0].keys())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output_rows)
    counts = Counter(row["topology_macroclass"] for row in output_rows if row.get("topology_macroclass"))
    summary = {
        "topology_input": str(args.topology),
        "calibration_input": str(args.calibration),
        "rows_written": len(output_rows),
        "selected_threshold": threshold,
        "selection_rule": "calibration split only: maximize nmi_scope_fold+nmi_cath_topology among thresholds with recurring_class_count>=10 and within_feature_distance<between_feature_distance",
        "calibration_candidates_considered": candidates,
        "macroclass_count": len(counts),
        "recurring_macroclass_count": sum(count >= 2 for count in counts.values()),
        "macroclass_status": "calibrated_candidate; occupancy remains gate-blocked",
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Assigned {len(counts)} candidate macroclasses at threshold {threshold:.2f} across {len(output_rows)} rows.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations

"""Promote calibration-selected macroclass candidates into the guarded final table."""

import argparse
import csv
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Promote macroclasses only after the representation gate passes.")
    parser.add_argument("--topology", type=Path, default=Path("data/processed/v2/domain_topology_pdbredo_calibrated_candidate_v2.csv"))
    parser.add_argument("--gate", type=Path, default=Path("data/processed/v2/representation_gate_pdbredo_primary_v2.json"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/domain_topology_v2_final.csv"))
    parser.add_argument("--summary-output", type=Path, default=Path("data/processed/v2/domain_topology_v2_final_summary.json"))
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    args = parse_args()
    gate = json.loads(args.gate.read_text(encoding="utf-8"))
    if not str(gate.get("gate_status", "")).startswith("passed"):
        raise SystemExit(f"Representation gate is {gate.get('gate_status', '')!r}; final topology remains blocked.")
    rows = read_csv(args.topology)
    if not rows or any(not row.get("topology_macroclass") for row in rows):
        raise SystemExit("Candidate topology table contains empty macroclass labels.")
    promoted: list[dict[str, str]] = []
    for row in rows:
        promoted_row = dict(row)
        promoted_row["macroclass_status"] = "calibrated"
        promoted_row["quality_flags"] = f"{row.get('quality_flags', '')};macroclass_gate_passed"
        promoted.append(promoted_row)
    fields = list(promoted[0].keys())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(promoted)
    summary = {
        "topology_input": str(args.topology),
        "gate_input": str(args.gate),
        "gate_status": gate["gate_status"],
        "rows_written": len(promoted),
        "macroclass_count": len({row["topology_macroclass"] for row in promoted}),
        "macroclass_status": "calibrated",
        "sparsity_status": "not_assigned",
        "status": "final_topology_table_ready_for_guarded_occupancy",
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Promoted {len(promoted)} rows into the final topology table.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

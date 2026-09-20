from __future__ import annotations

"""Filter the PDB-REDO table to complete primary-source alpha/beta candidates."""

import argparse
import csv
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Keep only complete PDB-REDO DSSP alpha/beta candidates for calibration.")
    parser.add_argument("--input", type=Path, default=Path("data/processed/v2/domain_topology_pdbredo_dssp_v2.csv"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/domain_topology_pdbredo_primary_alpha_beta_v2.csv"))
    parser.add_argument("--summary-output", type=Path, default=Path("data/processed/v2/domain_topology_pdbredo_primary_alpha_beta_summary_v2.json"))
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    args = parse_args()
    rows = read_csv(args.input)
    candidates = [
        row for row in rows
        if row.get("sse_assignment_status") == "precomputed_primary_candidate"
        and row.get("analysis_candidate_status") == "pdbredo_alpha_beta_candidate"
        and row.get("dssp_coverage_fraction") == "1.000000"
        and row.get("topology_signature_exact")
    ]
    fields = list(rows[0].keys()) if rows else []
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(candidates)
    summary = {
        "input": str(args.input),
        "input_rows": len(rows),
        "candidate_rows": len(candidates),
        "excluded_rows": len(rows) - len(candidates),
        "selection_rule": "precomputed_primary_candidate AND pdbredo_alpha_beta_candidate AND dssp_coverage_fraction=1.000000 AND nonempty_signature",
        "status": "complete_primary_alpha_beta_pool",
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(candidates)} complete primary alpha/beta rows from {len(rows)} source rows.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations

"""Assess CATH topology as a database-level fallback class definition."""

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate CATH topology fallback occupancy and breadth.")
    parser.add_argument("--topology", type=Path, default=Path("data/processed/v2/domain_topology_v2_final.csv"))
    parser.add_argument("--mapping", type=Path, default=Path("data/processed/v2/domain_mapping_breadth_repaired_v2.csv"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/cath_topology_fallback_v2.csv"))
    parser.add_argument("--summary-output", type=Path, default=Path("data/processed/v2/cath_topology_fallback_summary_v2.json"))
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def unique_count(rows: list[dict[str, str]], field: str) -> int:
    return len({row.get(field, "") for row in rows if row.get(field, "")})


def main() -> int:
    args = parse_args()
    topology = read_csv(args.topology)
    mapping = {row["domain_uid"]: row for row in read_csv(args.mapping)}
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for source in topology:
        row = {**source, **mapping.get(source["domain_uid"], {})}
        cath_id = row.get("cath_topology_id", "")
        if cath_id:
            grouped[cath_id].append(row)
    output_rows = []
    for cath_id, rows in sorted(grouped.items()):
        output_rows.append({
            "cath_topology_id": cath_id,
            "n_domain_rows": str(len(rows)),
            "n_structures": str(unique_count(rows, "structure_id")),
            "n_scope_superfamilies": str(unique_count(rows, "scope_superfamily_id")),
            "n_cath_homologies": str(unique_count(rows, "cath_homology_id")),
            "n_ecod_h_groups": str(unique_count(rows, "ecod_h_id")),
            "n_explicit_source_boundaries": str(sum(row.get("source_boundary_status") == "explicit" for row in rows)),
            "fallback_breadth_status": "complete_three_source" if all(unique_count(rows, field) > 0 for field in ("scope_superfamily_id", "cath_homology_id", "ecod_h_id")) else "incomplete_three_source",
        })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(output_rows[0].keys()) if output_rows else ["cath_topology_id"]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output_rows)
    summary = {
        "topology_input": str(args.topology),
        "mapping_input": str(args.mapping),
        "output": str(args.output),
        "cath_topology_class_count": len(output_rows),
        "recurring_class_count": sum(int(row["n_domain_rows"]) >= 2 for row in output_rows),
        "recurring_class_count_ge_three": sum(int(row["n_domain_rows"]) >= 3 for row in output_rows),
        "complete_three_source_class_count": sum(row["fallback_breadth_status"] == "complete_three_source" for row in output_rows),
        "status": "fallback_diagnostic_only",
        "primary_topology_changed": False,
        "interpretation": "CATH topology is reported as a permitted database-level fallback, not substituted for the calibrated contact topology without explicit preregistration.",
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Assessed {len(output_rows)} CATH topology fallback classes.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

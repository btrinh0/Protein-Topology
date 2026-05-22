from __future__ import annotations

"""Build the V2 topology occupancy atlas only after representation validation."""

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a guarded V2 occupancy atlas.")
    parser.add_argument("--topology", type=Path, default=Path("data/processed/v2/domain_topology_v2_final.csv"))
    parser.add_argument("--mapping", type=Path, default=Path("data/processed/v2/domain_mapping_v2.csv"))
    parser.add_argument("--gate", type=Path, default=Path("data/processed/v2/representation_gate_v2.json"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/occupancy_atlas_v2.csv"))
    parser.add_argument("--summary-output", type=Path, default=Path("data/processed/v2/occupancy_atlas_summary_v2.json"))
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def unique_count(rows: list[dict[str, str]], field: str) -> int:
    return len({row.get(field, "") for row in rows if row.get(field, "")})


def main() -> int:
    args = parse_args()
    if not args.gate.exists():
        raise SystemExit(f"Representation gate file not found: {args.gate}")
    gate = json.loads(args.gate.read_text(encoding="utf-8"))
    gate_status = str(gate.get("gate_status", ""))
    if not gate_status.startswith("passed"):
        raise SystemExit(f"Representation gate is {gate_status!r}; occupancy atlas remains blocked.")
    if not args.topology.exists():
        raise SystemExit(f"Final topology table not found: {args.topology}")
    topology_rows = read_csv(args.topology)
    mapping_rows = {row["domain_uid"]: row for row in read_csv(args.mapping)}
    eligible = []
    for row in topology_rows:
        macroclass = row.get("topology_macroclass", "")
        if not macroclass or row.get("macroclass_status") not in {"calibrated", "final"}:
            continue
        eligible.append({**row, **mapping_rows.get(row["domain_uid"], {})})
    if not eligible:
        raise SystemExit("No final calibrated topology macroclasses are available.")
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in eligible:
        grouped[row["topology_macroclass"]].append(row)
    atlas_rows = []
    for macroclass, rows in sorted(grouped.items()):
        atlas_rows.append({
            "topology_macroclass": macroclass,
            "n_domain_rows": str(len(rows)),
            "n_structures": str(unique_count(rows, "structure_id")),
            "n_astral40_clusters": str(unique_count(rows, "astral_cluster_id")),
            "n_scope_superfamilies": str(unique_count(rows, "scope_superfamily_id")),
            "n_cath_topologies": str(unique_count(rows, "cath_topology_id")),
            "n_cath_homologies": str(unique_count(rows, "cath_homology_id")),
            "n_ecod_h_groups": str(unique_count(rows, "ecod_h_id")),
            "n_ecod_t_groups": str(unique_count(rows, "ecod_t_id")),
            "n_explicit_source_boundaries": str(sum(row.get("source_boundary_status") == "explicit" for row in rows)),
            "n_whole_chain_source_boundaries": str(sum(row.get("source_boundary_status") == "whole_chain" for row in rows)),
            "n_rows_with_scope_breadth": str(sum(bool(row.get("scope_superfamily_id")) for row in rows)),
            "n_rows_with_cath_breadth": str(sum(bool(row.get("cath_homology_id")) for row in rows)),
            "n_rows_with_ecod_breadth": str(sum(bool(row.get("ecod_h_id")) for row in rows)),
            "sparsity_label": "",
            "sparsity_rule_status": "awaiting_preregistered_breadth_threshold",
        })
    fields = list(atlas_rows[0].keys())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(atlas_rows)
    summary = {
        "topology_input": str(args.topology),
        "mapping_input": str(args.mapping),
        "gate_input": str(args.gate),
        "gate_status": gate_status,
        "eligible_domain_rows": len(eligible),
        "macroclass_count": len(atlas_rows),
        "sparsity_labels_assigned": False,
        "status": "occupancy_counts_ready; sparse_common_threshold_not_assigned",
        "rule": "Independent evolutionary breadth thresholds must be preregistered before target selection.",
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote occupancy counts for {len(atlas_rows)} calibrated macroclasses.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

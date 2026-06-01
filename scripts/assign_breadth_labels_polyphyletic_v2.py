from __future__ import annotations

"""Build a strict polyphyletic breadth sensitivity atlas."""

import argparse
import csv
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Label strict polyphyletic common classes and singleton sparse classes.")
    parser.add_argument("--atlas", type=Path, default=Path("data/processed/v2/occupancy_atlas_breadth_repaired_v2.csv"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/occupancy_atlas_polyphyletic_breadth_v2.csv"))
    parser.add_argument("--summary-output", type=Path, default=Path("data/processed/v2/occupancy_breadth_polyphyletic_summary_v2.json"))
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def integer(row: dict[str, str], field: str) -> int:
    try:
        return int(row.get(field, "0"))
    except (TypeError, ValueError):
        return 0


def main() -> int:
    args = parse_args()
    rows = read_csv(args.atlas)
    common_ids = {
        row["topology_macroclass"] for row in rows
        if integer(row, "n_scope_superfamilies") > 1
        and integer(row, "n_cath_homologies") > 1
        and integer(row, "n_ecod_h_groups") > 0
    }
    sparse_ids = {
        row["topology_macroclass"] for row in rows
        if integer(row, "n_astral40_clusters") == 1
        and integer(row, "n_scope_superfamilies") == 1
        and integer(row, "n_cath_homologies") == 1
        and integer(row, "n_ecod_h_groups") == 1
        and integer(row, "n_explicit_source_boundaries") == integer(row, "n_domain_rows")
    }
    if len(common_ids) < 6 or len(sparse_ids) < 6:
        raise SystemExit("Strict polyphyletic rule did not produce six classes in both arms.")
    output_rows = []
    for source in rows:
        row = dict(source)
        class_id = row["topology_macroclass"]
        row["sparsity_label"] = "common" if class_id in common_ids else "sparse" if class_id in sparse_ids else ""
        row["sparsity_rule_status"] = "strict_polyphyletic_candidate"
        output_rows.append(row)
    fields = list(output_rows[0].keys())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output_rows)
    summary = {
        "atlas_input": str(args.atlas),
        "output": str(args.output),
        "strict_common_class_count": len(common_ids),
        "strict_sparse_singleton_class_count": len(sparse_ids),
        "common_rule": "n_scope_superfamilies>1, n_cath_homologies>1, and nonzero ECOD breadth",
        "sparse_rule": "all four breadth counts equal one and all rows have explicit source boundaries",
        "sparse_common_labels_assigned": True,
        "status": "candidate_ready_for_six_pair_target_matching;secondary_breadth_lane",
        "primary_atlas_changed": False,
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Assigned {len(sparse_ids)} sparse and {len(common_ids)} strict common candidate classes.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

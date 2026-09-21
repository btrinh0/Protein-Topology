from __future__ import annotations

"""Assess breadth coverage for the separately versioned repaired mapping lane."""

import argparse
import csv
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Check repaired evolutionary-breadth coverage.")
    parser.add_argument("--atlas", type=Path, default=Path("data/processed/v2/occupancy_atlas_breadth_repaired_v2.csv"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/occupancy_breadth_resolution_breadth_repaired_v2.json"))
    parser.add_argument("--minimum-classes-per-arm", type=int, default=8)
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
    complete = [row for row in rows if integer(row, "n_scope_superfamilies") > 0 and integer(row, "n_cath_homologies") > 0]
    scope_broad = [row for row in complete if integer(row, "n_scope_superfamilies") > 1]
    cath_broad = [row for row in complete if integer(row, "n_cath_homologies") > 1]
    ecod_covered = [row for row in rows if integer(row, "n_ecod_h_groups") > 0]
    recurring = [row for row in rows if integer(row, "n_domain_rows") >= 2]
    blocking_reasons = ["Most macroclasses have one observed SCOPe superfamily and one observed CATH homology."]
    if ecod_covered:
        blocking_reasons.append("ECOD breadth is covered by local x_name:h_name tuples; these are a sensitivity field, not official numeric ECOD identifiers.")
    else:
        blocking_reasons.append("No macroclass has a nonzero ECOD H-group breadth in the supplied mapping.")
    blocking_reasons.append("Raw domain-row occupancy is retained descriptively but is not substituted for independent evolutionary breadth.")
    summary = {
        "atlas_input": str(args.atlas),
        "macroclass_count": len(rows),
        "complete_scope_cath_coverage_classes": len(complete),
        "scope_superfamily_breadth_gt_one_classes": len(scope_broad),
        "cath_homology_breadth_gt_one_classes": len(cath_broad),
        "ecod_h_group_covered_classes": len(ecod_covered),
        "recurring_macroclasses": len(recurring),
        "minimum_classes_per_arm": args.minimum_classes_per_arm,
        "sparse_common_label_rule": "requires at least the configured number of classes in each arm with independently resolved SCOPe/CATH breadth; ties cannot be broken by raw PDB count alone",
        "sparse_common_labels_assigned": False,
        "status": "breadth_resolution_insufficient_for_sparse_common_labels" if min(len(scope_broad), len(cath_broad)) < args.minimum_classes_per_arm else "breadth_resolution_candidate",
        "blocking_reasons": blocking_reasons,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Repaired-lane breadth status: {summary['status']}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations

"""Select preregistered sparse/common topology pairs after the V2 atlas gate."""

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a guarded sparse/common target manifest.")
    parser.add_argument("--topology", type=Path, default=Path("data/processed/v2/domain_topology_v2_final.csv"))
    parser.add_argument("--atlas", type=Path, default=Path("data/processed/v2/occupancy_atlas_v2.csv"))
    parser.add_argument("--gate", type=Path, default=Path("data/processed/v2/representation_gate_v2.json"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/matched_target_pairs_v1.csv"))
    parser.add_argument("--summary-output", type=Path, default=Path("data/processed/v2/matched_target_pairs_summary_v1.json"))
    parser.add_argument("--pair-count", type=int, default=8)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def value(row: dict[str, str], field: str) -> float:
    try:
        return float(row.get(field, ""))
    except (TypeError, ValueError):
        return 0.0


def class_features(rows: list[dict[str, str]]) -> tuple[float, ...]:
    fields = ("resolved_residue_count", "sse_count", "contact_density", "contact_order", "sse_contact_edge_count")
    return tuple(sum(value(row, field) for row in rows) / max(len(rows), 1) for field in fields)


def distance(left: tuple[float, ...], right: tuple[float, ...], scales: tuple[float, ...]) -> float:
    return math.sqrt(sum(((a - b) / scale) ** 2 for a, b, scale in zip(left, right, scales)))


def main() -> int:
    args = parse_args()
    if args.pair_count <= 0:
        raise SystemExit("--pair-count must be positive.")
    gate = json.loads(args.gate.read_text(encoding="utf-8")) if args.gate.exists() else {}
    gate_status = str(gate.get("gate_status", ""))
    if not gate_status.startswith("passed"):
        raise SystemExit(f"Representation gate is {gate_status!r}; target selection remains blocked.")
    if not args.atlas.exists():
        raise SystemExit(f"Occupancy atlas not found: {args.atlas}")
    atlas = read_csv(args.atlas)
    if not any(row.get("sparsity_label") in {"sparse", "common"} for row in atlas):
        raise SystemExit("Occupancy atlas has no preregistered sparse/common labels.")
    topology = read_csv(args.topology)
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in topology:
        if row.get("topology_macroclass"):
            grouped[row["topology_macroclass"]].append(row)
    atlas_by_class = {row["topology_macroclass"]: row for row in atlas}
    sparse = [row for row in atlas if row.get("sparsity_label") == "sparse" and row["topology_macroclass"] in grouped]
    common = [row for row in atlas if row.get("sparsity_label") == "common" and row["topology_macroclass"] in grouped]
    scales = (20.0, 2.0, 0.25, 0.25, 5.0)
    candidates = []
    for sparse_row in sparse:
        sparse_features = class_features(grouped[sparse_row["topology_macroclass"]])
        for common_row in common:
            common_features = class_features(grouped[common_row["topology_macroclass"]])
            candidates.append((
                distance(sparse_features, common_features, scales),
                sparse_row["topology_macroclass"],
                common_row["topology_macroclass"],
            ))
    candidates.sort()
    selected = []
    used_sparse: set[str] = set()
    used_common: set[str] = set()
    for covariate_distance, sparse_class, common_class in candidates:
        if sparse_class in used_sparse or common_class in used_common:
            continue
        selected.append({
            "pair_id": f"TCT-V2-PAIR-{len(selected) + 1:02d}",
            "sparse_topology_macroclass": sparse_class,
            "common_topology_macroclass": common_class,
            "covariate_distance": f"{covariate_distance:.6f}",
            "selection_status": "frozen",
        })
        used_sparse.add(sparse_class)
        used_common.add(common_class)
        if len(selected) >= args.pair_count:
            break
    if len(selected) < args.pair_count:
        raise SystemExit(f"Only {len(selected)} balanced pairs were available; {args.pair_count} required.")
    fields = list(selected[0].keys())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(selected)
    summary = {
        "topology_input": str(args.topology),
        "atlas_input": str(args.atlas),
        "gate_input": str(args.gate),
        "pair_count_requested": args.pair_count,
        "pair_count_written": len(selected),
        "selection_status": "frozen",
        "matching_features": ["resolved_residue_count", "sse_count", "contact_density", "contact_order", "sse_contact_edge_count"],
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(selected)} frozen sparse/common topology pairs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

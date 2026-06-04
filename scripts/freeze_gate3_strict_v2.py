from __future__ import annotations

"""Freeze the strict six-pair Gate 3 manifest with provenance hashes."""

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate and freeze the strict polyphyletic Gate 3 lane.")
    parser.add_argument("--pairs", type=Path, default=Path("data/processed/v2/matched_target_pairs_polyphyletic_v2.csv"))
    parser.add_argument("--atlas", type=Path, default=Path("data/processed/v2/occupancy_atlas_polyphyletic_breadth_v2.csv"))
    parser.add_argument("--topology", type=Path, default=Path("data/processed/v2/domain_topology_v2_final.csv"))
    parser.add_argument("--gate", type=Path, default=Path("data/processed/v2/representation_gate_pdbredo_primary_v2.json"))
    parser.add_argument("--mapping", type=Path, default=Path("data/processed/v2/domain_mapping_breadth_repaired_v2.csv"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/matched_target_pairs_gate3_strict_v2.csv"))
    parser.add_argument("--summary-output", type=Path, default=Path("data/processed/v2/gate3_strict_freeze_summary_v2.json"))
    parser.add_argument("--decision-output", type=Path, default=Path("data/processed/v2/gate3_lane_decision_strict_v2.json"))
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def integer(row: dict[str, str], field: str) -> int:
    try:
        return int(row.get(field, "0"))
    except (TypeError, ValueError):
        return 0


def main() -> int:
    args = parse_args()
    for path in (args.pairs, args.atlas, args.topology, args.gate, args.mapping):
        if not path.exists():
            raise SystemExit(f"Required input not found: {path}")
    gate = json.loads(args.gate.read_text(encoding="utf-8"))
    if not str(gate.get("gate_status", "")).startswith("passed"):
        raise SystemExit(f"Representation gate is {gate.get('gate_status')!r}.")
    pairs = read_csv(args.pairs)
    atlas = read_csv(args.atlas)
    topology = read_csv(args.topology)
    mapping = read_csv(args.mapping)
    if len(pairs) != 6:
        raise SystemExit(f"Strict primary lane requires exactly six pairs; found {len(pairs)}.")
    expected_ids = [f"TCT-V2-PAIR-{index:02d}" for index in range(1, 7)]
    if [row.get("pair_id") for row in pairs] != expected_ids:
        raise SystemExit("Pair IDs are not the expected sequential six-pair manifest.")
    atlas_by_class = {row.get("topology_macroclass", ""): row for row in atlas}
    topology_classes = {row.get("topology_macroclass", "") for row in topology if row.get("topology_macroclass")}
    sparse_classes = set()
    common_classes = set()
    for row in atlas:
        label = row.get("sparsity_label", "")
        class_id = row.get("topology_macroclass", "")
        if label == "sparse":
            sparse_classes.add(class_id)
            strict_sparse = all(integer(row, field) == 1 for field in ("n_astral40_clusters", "n_scope_superfamilies", "n_cath_homologies", "n_ecod_h_groups"))
            strict_sparse = strict_sparse and integer(row, "n_explicit_source_boundaries") == integer(row, "n_domain_rows")
            if not strict_sparse:
                raise SystemExit(f"Sparse class {class_id} fails the strict singleton rule.")
        elif label == "common":
            common_classes.add(class_id)
            strict_common = integer(row, "n_scope_superfamilies") > 1 and integer(row, "n_cath_homologies") > 1 and integer(row, "n_ecod_h_groups") > 0
            if not strict_common:
                raise SystemExit(f"Common class {class_id} fails the strict polyphyletic rule.")
    if len(common_classes) != 6 or len(sparse_classes) < 6:
        raise SystemExit("The strict atlas does not contain six common and six sparse classes.")
    if sparse_classes & common_classes:
        raise SystemExit("Sparse and common class sets overlap.")
    for pair in pairs:
        sparse = pair.get("sparse_topology_macroclass", "")
        common = pair.get("common_topology_macroclass", "")
        if sparse not in sparse_classes or common not in common_classes:
            raise SystemExit(f"Pair {pair.get('pair_id')} references an invalid strict label.")
        if sparse not in topology_classes or common not in topology_classes:
            raise SystemExit(f"Pair {pair.get('pair_id')} references a class absent from topology rows.")
    if len({row.get("sparse_topology_macroclass") for row in pairs}) != 6 or len({row.get("common_topology_macroclass") for row in pairs}) != 6:
        raise SystemExit("Strict manifest reuses a sparse or common class.")
    fields = list(pairs[0].keys())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(pairs)
    inputs = {name: {"path": str(path), "sha256": sha256(path)} for name, path in {
        "pairs_source": args.pairs,
        "atlas": args.atlas,
        "topology": args.topology,
        "representation_gate": args.gate,
        "breadth_mapping": args.mapping,
    }.items()}
    decision = {
        "decision_status": "user_directed_primary_lane",
        "decision_date": datetime.now(timezone.utc).date().isoformat(),
        "primary_lane": "strict_polyphyletic_six_pair",
        "primary_pair_count": 6,
        "primary_common_rule": "n_scope_superfamilies>1, n_cath_homologies>1, and nonzero ECOD breadth",
        "primary_sparse_rule": "all four breadth counts equal one and all rows have explicit source boundaries",
        "sensitivity_lane": "quantile_breadth_eight_pair",
        "fallback_lane": "CATH topology database-level sensitivity only",
        "design_campaign_status": "not_started_until_manifest_freeze",
    }
    args.decision_output.parent.mkdir(parents=True, exist_ok=True)
    args.decision_output.write_text(json.dumps(decision, indent=2) + "\n", encoding="utf-8")
    summary = {
        "freeze_status": "passed",
        "manifest": str(args.output),
        "pair_count": len(pairs),
        "common_class_count": len(common_classes),
        "sparse_class_count": len(sparse_classes),
        "selection_status": "frozen_primary_gate3",
        "representation_gate_status": gate.get("gate_status", ""),
        "inputs": inputs,
        "decision_file": str(args.decision_output),
        "sensitivity_manifest": "data/processed/v2/matched_target_pairs_quantile_breadth_v2.csv",
        "validation_checks": {
            "sequential_pair_ids": True,
            "strict_common_rule": True,
            "strict_sparse_rule": True,
            "unique_classes_per_arm": True,
            "classes_present_in_topology": True,
            "representation_gate_passed": True,
        },
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Frozen strict Gate 3 manifest with {len(pairs)} pairs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

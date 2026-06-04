from __future__ import annotations

"""Prepare blinded, equal-budget Gate 4 pilot targets from the frozen pairs."""

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


FEATURES = ("resolved_residue_count", "sse_count", "contact_density", "contact_order", "sse_contact_edge_count")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare neutral design-pilot targets and a separate blinding key.")
    parser.add_argument("--pairs", type=Path, default=Path("data/processed/v2/matched_target_pairs_gate3_strict_v2.csv"))
    parser.add_argument("--topology", type=Path, default=Path("data/processed/v2/domain_topology_v2_final.csv"))
    parser.add_argument("--decision", type=Path, default=Path("data/processed/v2/gate3_lane_decision_strict_v2.json"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/design_pilot_targets_blinded_v2.csv"))
    parser.add_argument("--key-output", type=Path, default=Path("data/processed/v2/design_pilot_blinding_key_v2.json"))
    parser.add_argument("--summary-output", type=Path, default=Path("data/processed/v2/design_pilot_preparation_summary_v2.json"))
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def number(row: dict[str, str], field: str) -> float:
    try:
        return float(row.get(field, ""))
    except (TypeError, ValueError):
        return 0.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def choose_representative(rows: list[dict[str, str]]) -> dict[str, str]:
    medians = {field: sorted(number(row, field) for row in rows)[len(rows) // 2] for field in FEATURES}
    scales = {field: max(abs(number(row, field) - medians[field]) for row in rows) or 1.0 for field in FEATURES}
    return min(
        rows,
        key=lambda row: (
            sum(((number(row, field) - medians[field]) / scales[field]) ** 2 for field in FEATURES),
            row.get("domain_uid", ""),
        ),
    )


def main() -> int:
    args = parse_args()
    for path in (args.pairs, args.topology, args.decision):
        if not path.exists():
            raise SystemExit(f"Required input not found: {path}")
    decision = json.loads(args.decision.read_text(encoding="utf-8"))
    if decision.get("primary_lane") != "strict_polyphyletic_six_pair":
        raise SystemExit("The Gate 3 decision does not select the strict six-pair lane.")
    pairs = read_csv(args.pairs)
    topology = read_csv(args.topology)
    by_class: dict[str, list[dict[str, str]]] = {}
    for row in topology:
        class_id = row.get("topology_macroclass", "")
        if class_id:
            by_class.setdefault(class_id, []).append(row)
    targets: list[dict[str, str]] = []
    key: dict[str, dict[str, str]] = {}
    for pair in pairs:
        for arm_code, label_field in (("A", "sparse_topology_macroclass"), ("B", "common_topology_macroclass")):
            class_id = pair.get(label_field, "")
            candidates = by_class.get(class_id, [])
            if not candidates:
                raise SystemExit(f"No topology row exists for {class_id}.")
            representative = choose_representative(candidates)
            target_id = f"TCT-V2-TGT-{len(targets) + 1:02d}"
            targets.append({
                "target_id": target_id,
                "pair_id": pair["pair_id"],
                "arm_code": arm_code,
                "domain_uid": representative["domain_uid"],
                "structure_id": representative["structure_id"],
                "chain_id": representative["chain_id"],
                "source_path": representative["source_path"],
                "domain_segments_json": representative["domain_segments_json"],
                "resolved_residue_count": representative["resolved_residue_count"],
                "sse_count": representative["sse_count"],
                "sse_order": representative["sse_order"],
                "topology_signature_exact": representative["topology_signature_exact"],
                "template_selection_status": "deterministic_class_median_covariate_representative",
                "design_budget_status": "equal_budget_required",
            })
            key[target_id] = {
                "pair_id": pair["pair_id"],
                "arm_code": arm_code,
                "sparsity_label": "sparse" if arm_code == "A" else "common",
                "topology_macroclass": class_id,
                "domain_uid": representative["domain_uid"],
                "structure_id": representative["structure_id"],
            }
    fields = list(targets[0].keys())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(targets)
    args.key_output.parent.mkdir(parents=True, exist_ok=True)
    args.key_output.write_text(json.dumps(key, indent=2) + "\n", encoding="utf-8")
    summary = {
        "preparation_status": "pass",
        "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
        "target_count": len(targets),
        "pair_count": len(pairs),
        "targets_per_arm": len(targets) // 2,
        "selection_rule": "one deterministic representative per arm per pair, minimizing scaled distance to the class median across the five frozen matching features",
        "blinding": "neutral target IDs and arm codes are in the pilot manifest; sparsity labels are stored only in the separate key",
        "design_budget": {
            "sequence_samples_per_target": 16,
            "backbone_generation_attempts_per_target": 8,
            "identical_filtering_and_scoring": True,
            "positive_controls": 2,
            "negative_controls": 2,
        },
        "inputs": {
            "pairs": {"path": str(args.pairs), "sha256": sha256(args.pairs)},
            "topology": {"path": str(args.topology), "sha256": sha256(args.topology)},
            "decision": {"path": str(args.decision), "sha256": sha256(args.decision)},
        },
        "pilot_manifest": str(args.output),
        "blinding_key": str(args.key_output),
        "generation_status": "pending_toolchain_preflight",
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Prepared {len(targets)} blinded pilot targets across {len(pairs)} pairs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

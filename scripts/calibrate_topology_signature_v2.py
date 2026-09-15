from __future__ import annotations

"""Evaluate provisional SSE signature merge thresholds without assigning a final class."""

import argparse
import csv
import hashlib
import json
import math
import random
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Calibrate V2 topology-signature merge candidates.")
    parser.add_argument(
        "--split",
        type=Path,
        default=Path("data/processed/v2/calibration_split_v2.csv"),
    )
    parser.add_argument(
        "--topology",
        type=Path,
        default=Path("data/processed/v2/domain_topology_alpha_beta_70_160_provisional.csv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/processed/v2/topology_signature_calibration_candidates_v2.csv"),
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=Path("data/processed/v2/topology_signature_calibration_summary_v2.json"),
    )
    parser.add_argument("--seed", default="tct-v2-signature-calibration-2026-08-23")
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def parse_signature(value: str) -> tuple[tuple[str, ...], frozenset[tuple[int, int, str]]]:
    if "|" not in value:
        return tuple(), frozenset()
    nodes_text, edges_text = value.split("|", 1)
    nodes = tuple(node for node in nodes_text.strip().split("-") if node)
    edges: set[tuple[int, int, str]] = set()
    for edge in edges_text.strip().split(";"):
        if not edge or ":" not in edge:
            continue
        endpoints, label = edge.split(":", 1)
        try:
            left, right = endpoints.split("-", 1)
            edges.add((int(left), int(right), label))
        except ValueError:
            continue
    return nodes, frozenset(edges)


def signature_distance(left: tuple[tuple[str, ...], frozenset[tuple[int, int, str]]], right: tuple[tuple[str, ...], frozenset[tuple[int, int, str]]]) -> float:
    left_nodes, left_edges = left
    right_nodes, right_edges = right
    if not left_nodes or not right_nodes or len(left_nodes) != len(right_nodes):
        return 1.0
    node_distance = sum(a != b for a, b in zip(left_nodes, right_nodes)) / len(left_nodes)
    edge_union = left_edges | right_edges
    edge_distance = len(left_edges ^ right_edges) / max(len(edge_union), 1)
    return 0.5 * node_distance + 0.5 * edge_distance


class UnionFind:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))

    def find(self, value: int) -> int:
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left: int, right: int) -> None:
        left_root, right_root = self.find(left), self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root


def macro_labels(signatures: dict[str, tuple[tuple[str, ...], frozenset[tuple[int, int, str]]]], threshold: float) -> dict[str, str]:
    keys = sorted(signatures)
    union_find = UnionFind(len(keys))
    for left_index, left_key in enumerate(keys):
        for right_index in range(left_index + 1, len(keys)):
            right_key = keys[right_index]
            if signature_distance(signatures[left_key], signatures[right_key]) <= threshold:
                union_find.union(left_index, right_index)
    roots: dict[int, list[str]] = {}
    for index, key in enumerate(keys):
        roots.setdefault(union_find.find(index), []).append(key)
    labels: dict[str, str] = {}
    for members in roots.values():
        token = hashlib.sha256("|".join(members).encode("utf-8")).hexdigest()[:16]
        for member in members:
            labels[member] = f"macro_{token}"
    return labels


def entropy(values: list[str]) -> float:
    total = len(values)
    if total == 0:
        return 0.0
    counts = Counter(values)
    return -sum((count / total) * math.log(count / total) for count in counts.values())


def normalized_mutual_information(left: list[str], right: list[str]) -> float:
    pairs = [(a, b) for a, b in zip(left, right) if a and b]
    if not pairs:
        return 0.0
    left_values, right_values = zip(*pairs)
    total = len(pairs)
    joint = Counter(pairs)
    left_counts = Counter(left_values)
    right_counts = Counter(right_values)
    mutual_information = sum(
        (count / total) * math.log((count * total) / (left_counts[a] * right_counts[b]))
        for (a, b), count in joint.items()
    )
    denominator = math.sqrt(entropy(list(left_values)) * entropy(list(right_values)))
    return mutual_information / denominator if denominator else 0.0


def feature_distance(left: dict[str, str], right: dict[str, str]) -> float | None:
    values: list[float] = []
    for field, scale in (("sse_count", 10.0), ("contact_density", 1.0), ("contact_order", 1.0)):
        try:
            values.append(abs(float(left[field]) - float(right[field])) / scale)
        except (KeyError, TypeError, ValueError):
            return None
    return math.sqrt(sum(value * value for value in values))


def coherence_proxy(rows: list[dict[str, str]], labels: dict[str, str], seed: str) -> dict[str, float | int | None]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        label = labels.get(row.get("topology_signature_exact", ""), "")
        if label:
            grouped.setdefault(label, []).append(row)
    within: list[float] = []
    for members in grouped.values():
        for left_index, left in enumerate(members):
            for right in members[left_index + 1 :]:
                distance = feature_distance(left, right)
                if distance is not None:
                    within.append(distance)
    rng = random.Random(seed)
    shuffled = list(rows)
    rng.shuffle(shuffled)
    between: list[float] = []
    for left_index, left in enumerate(shuffled):
        left_label = labels.get(left.get("topology_signature_exact", ""), "")
        for right in shuffled[left_index + 1 :]:
            if len(between) >= max(len(within), 1):
                break
            right_label = labels.get(right.get("topology_signature_exact", ""), "")
            if left_label and right_label and left_label != right_label:
                distance = feature_distance(left, right)
                if distance is not None:
                    between.append(distance)
        if len(between) >= max(len(within), 1):
            break
    return {
        "within_pair_count": len(within),
        "between_pair_count": len(between),
        "within_mean_feature_distance": sum(within) / len(within) if within else None,
        "between_mean_feature_distance": sum(between) / len(between) if between else None,
    }


def main() -> int:
    args = parse_args()
    split_rows = read_csv(args.split)
    topology_rows = {row["domain_uid"]: row for row in read_csv(args.topology)}
    rows = [
        {**row, **topology_rows.get(row["domain_uid"], {})}
        for row in split_rows
        if row.get("eligible_precalibration") == "true"
    ]
    signatures = {
        row["topology_signature_exact"]: parse_signature(row["topology_signature_exact"])
        for row in rows
        if row.get("topology_signature_exact")
    }
    thresholds = (0.0, 0.10, 0.20, 0.30)
    candidate_rows: list[dict[str, str]] = []
    candidate_summary: list[dict[str, object]] = []
    for threshold in thresholds:
        labels = macro_labels(signatures, threshold)
        for split in ("calibration", "held_out"):
            subset = [row for row in rows if row["split"] == split]
            subset_labels = [labels.get(row.get("topology_signature_exact", ""), "") for row in subset]
            class_counts = Counter(label for label in subset_labels if label)
            scope_labels = [row.get("scope_fold_id", "") for row in subset]
            cath_labels = [row.get("cath_topology_id", "") for row in subset]
            baseline = [row.get("sse_order", "") for row in subset]
            coherence = coherence_proxy(subset, labels, f"{args.seed}:{threshold}:{split}")
            candidate_rows.append(
                {
                    "threshold": f"{threshold:.2f}",
                    "split": split,
                    "row_count": str(len(subset)),
                    "macroclass_count": str(len(class_counts)),
                    "macroclass_singleton_fraction": f"{sum(count == 1 for count in class_counts.values()) / max(len(class_counts), 1):.6f}",
                    "macroclass_recurrence_fraction": f"{sum(count >= 2 for count in class_counts.values()) / max(len(class_counts), 1):.6f}",
                    "nmi_scope_fold": f"{normalized_mutual_information(subset_labels, scope_labels):.6f}",
                    "nmi_cath_topology": f"{normalized_mutual_information(subset_labels, cath_labels):.6f}",
                    "baseline_nmi_scope_fold": f"{normalized_mutual_information(baseline, scope_labels):.6f}",
                    "baseline_nmi_cath_topology": f"{normalized_mutual_information(baseline, cath_labels):.6f}",
                    "within_mean_feature_distance": "" if coherence["within_mean_feature_distance"] is None else f"{coherence['within_mean_feature_distance']:.6f}",
                    "between_mean_feature_distance": "" if coherence["between_mean_feature_distance"] is None else f"{coherence['between_mean_feature_distance']:.6f}",
                    "within_pair_count": str(coherence["within_pair_count"]),
                    "between_pair_count": str(coherence["between_pair_count"]),
                }
            )
            candidate_summary.append(
                {
                    "threshold": threshold,
                    "split": split,
                    "macroclass_count": len(class_counts),
                    "recurring_class_count": sum(count >= 2 for count in class_counts.values()),
                    "nmi_scope_fold": normalized_mutual_information(subset_labels, scope_labels),
                    "nmi_cath_topology": normalized_mutual_information(subset_labels, cath_labels),
                    "baseline_nmi_scope_fold": normalized_mutual_information(baseline, scope_labels),
                    "baseline_nmi_cath_topology": normalized_mutual_information(baseline, cath_labels),
                    **coherence,
                }
            )

    fields = list(candidate_rows[0].keys())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(candidate_rows)
    summary = {
        "split_input": str(args.split),
        "topology_input": str(args.topology),
        "thresholds_evaluated": list(thresholds),
        "eligible_rows": len(rows),
        "candidate_metrics": candidate_summary,
        "structural_coherence_status": "proxy_only; TM-align or Foldseek similarity was not run",
        "sse_assignment_status": "provisional_author_annotation; DSSP/STRIDE sensitivity is required",
        "decision": "do_not_assign_macroclass_before_secondary_sse_and_structural_similarity_checks",
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Evaluated {len(thresholds)} merge thresholds across calibration and held-out rows.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations

"""Create the pending 50-domain manual audit for V2 domain mapping.

The script never marks a mapping as valid.  It selects a deterministic,
status-stratified sample and produces a template that must be completed from
the raw source annotations and coordinate records.
"""

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


AUDIT_FIELDS = [
    "audit_id",
    "audit_stratum",
    "audit_rank_within_stratum",
    "domain_uid",
    "source_domain_id",
    "structure_id",
    "chain_id",
    "domain_segments_json",
    "source_boundary_status",
    "cath_mapping_status",
    "cath_domain_id",
    "cath_mapping_overlap_score",
    "cath_mapping_basis",
    "ecod_mapping_status",
    "ecod_domain_id",
    "ecod_mapping_overlap_score",
    "ecod_mapping_basis",
    "quality_flags",
    "audit_status",
    "auditor",
    "audit_date",
    "scope_boundary_verified",
    "coordinate_slice_verified",
    "cath_candidate_set_verified",
    "cath_selection_verified",
    "ecod_candidate_set_verified",
    "ecod_selection_verified",
    "hierarchy_fields_verified",
    "audit_disposition",
    "audit_notes",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a status-stratified pending manual audit for V2 domain mappings."
    )
    parser.add_argument(
        "--mapping",
        type=Path,
        default=Path("data/processed/v2/domain_mapping_v2.csv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/processed/v2/domain_mapping_audit_sample_v2.csv"),
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=Path("data/processed/v2/domain_mapping_audit_summary_v2.json"),
    )
    parser.add_argument("--sample-size", type=int, default=50)
    parser.add_argument("--seed", default="tct-v2-domain-mapping-audit-2026-08-18")
    return parser.parse_args()


def stable_order_key(seed: str, domain_uid: str) -> str:
    return hashlib.sha256(f"{seed}:{domain_uid}".encode("utf-8")).hexdigest()


def audit_stratum(row: dict[str, str]) -> str:
    return " | ".join(
        [
            row.get("source_boundary_status", ""),
            row.get("cath_mapping_status", ""),
            row.get("ecod_mapping_status", ""),
        ]
    )


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=AUDIT_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    if args.sample_size <= 0:
        raise SystemExit("--sample-size must be positive.")
    with args.mapping.open("r", newline="", encoding="utf-8") as handle:
        mapping_rows = list(csv.DictReader(handle))
    if not mapping_rows:
        raise SystemExit(f"No rows found in {args.mapping}")

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in mapping_rows:
        grouped[audit_stratum(row)].append(row)
    for stratum, rows in grouped.items():
        rows.sort(key=lambda row: stable_order_key(args.seed, row["domain_uid"]))

    selected: list[dict[str, str]] = []
    positions = {stratum: 0 for stratum in grouped}
    ordered_strata = sorted(grouped, key=lambda stratum: (len(grouped[stratum]), stratum))
    while len(selected) < min(args.sample_size, len(mapping_rows)):
        added = False
        for stratum in ordered_strata:
            position = positions[stratum]
            if position >= len(grouped[stratum]):
                continue
            selected.append(grouped[stratum][position])
            positions[stratum] += 1
            added = True
            if len(selected) == min(args.sample_size, len(mapping_rows)):
                break
        if not added:
            break

    ranks = Counter()
    audit_rows: list[dict[str, str]] = []
    for index, row in enumerate(selected, start=1):
        stratum = audit_stratum(row)
        ranks[stratum] += 1
        audit_rows.append(
            {
                "audit_id": f"TCT-V2-MAP-{index:03d}",
                "audit_stratum": stratum,
                "audit_rank_within_stratum": str(ranks[stratum]),
                **row,
                "audit_status": "pending",
                "auditor": "",
                "audit_date": "",
                "scope_boundary_verified": "",
                "coordinate_slice_verified": "",
                "cath_candidate_set_verified": "",
                "cath_selection_verified": "",
                "ecod_candidate_set_verified": "",
                "ecod_selection_verified": "",
                "hierarchy_fields_verified": "",
                "audit_disposition": "",
                "audit_notes": "",
            }
        )
    write_csv(args.output, audit_rows)

    summary = {
        "input_mapping": str(args.mapping),
        "sample_size_requested": args.sample_size,
        "sample_size_written": len(audit_rows),
        "seed": args.seed,
        "selection_method": "round-robin across source-boundary/CATH/ECOD status strata, then SHA-256 seeded order within stratum",
        "full_mapping_status_counts": {
            "source_boundary_status": dict(sorted(Counter(row.get("source_boundary_status", "") for row in mapping_rows).items())),
            "cath_mapping_status": dict(sorted(Counter(row.get("cath_mapping_status", "") for row in mapping_rows).items())),
            "ecod_mapping_status": dict(sorted(Counter(row.get("ecod_mapping_status", "") for row in mapping_rows).items())),
        },
        "sample_stratum_counts": dict(sorted(Counter(row["audit_stratum"] for row in audit_rows).items())),
        "completion_requirement": "Every row remains pending until reviewed against raw SCOPe, CATH, ECOD, and coordinate evidence.",
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(audit_rows)} pending audit rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

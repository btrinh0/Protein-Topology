from __future__ import annotations

"""Create a deterministic structure-held-out calibration split for V2."""

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare the locked V2 calibration split.")
    parser.add_argument(
        "--topology",
        type=Path,
        default=Path("data/processed/v2/domain_topology_alpha_beta_70_160_provisional.csv"),
    )
    parser.add_argument(
        "--mapping",
        type=Path,
        default=Path("data/processed/v2/domain_mapping_v2.csv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/processed/v2/calibration_split_v2.csv"),
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=Path("data/processed/v2/calibration_split_summary_v2.json"),
    )
    parser.add_argument("--calibration-fraction", type=float, default=0.30)
    parser.add_argument("--seed", default="tct-v2-calibration-split-2026-08-23")
    return parser.parse_args()


def split_hash(seed: str, structure_id: str) -> tuple[str, float]:
    digest = hashlib.sha256(f"{seed}:{structure_id}".encode("utf-8")).hexdigest()
    return digest, int(digest[:12], 16) / float(16**12)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    args = parse_args()
    if not 0.0 < args.calibration_fraction < 1.0:
        raise SystemExit("--calibration-fraction must be between 0 and 1.")
    topology_rows = read_csv(args.topology)
    mapping_rows = {row["domain_uid"]: row for row in read_csv(args.mapping)}
    if not topology_rows:
        raise SystemExit(f"No topology rows found in {args.topology}")

    rows: list[dict[str, str]] = []
    for topology in topology_rows:
        domain_uid = topology["domain_uid"]
        mapping = mapping_rows.get(domain_uid, {})
        digest, uniform_value = split_hash(args.seed, topology["structure_id"])
        eligible = topology["topology_status"] == "computed" and topology["analysis_candidate_status"] == "provisional_alpha_beta_candidate"
        split = "calibration" if uniform_value < args.calibration_fraction else "held_out"
        rows.append(
            {
                "domain_uid": domain_uid,
                "source_domain_id": topology["source_domain_id"],
                "structure_id": topology["structure_id"],
                "chain_id": topology["chain_id"],
                "split": split,
                "split_hash": digest,
                "split_uniform_value": f"{uniform_value:.12f}",
                "eligible_precalibration": str(eligible).lower(),
                "scope_fold_id": mapping.get("scope_fold_id", ""),
                "scope_superfamily_id": mapping.get("scope_superfamily_id", ""),
                "cath_topology_id": mapping.get("cath_topology_id", ""),
                "cath_homology_id": mapping.get("cath_homology_id", ""),
                "sse_count": topology.get("sse_count", ""),
                "sse_order": topology.get("sse_order", ""),
                "topology_signature_exact": topology.get("topology_signature_exact", ""),
                "sse_assignment_method": topology.get("sse_assignment_method", ""),
                "contact_definition_id": topology.get("contact_definition_id", ""),
                "quality_flags": topology.get("quality_flags", ""),
            }
        )

    rows.sort(key=lambda row: row["domain_uid"])
    fields = list(rows[0].keys())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    eligible_rows = [row for row in rows if row["eligible_precalibration"] == "true"]
    summary = {
        "topology_input": str(args.topology),
        "mapping_input": str(args.mapping),
        "seed": args.seed,
        "calibration_fraction": args.calibration_fraction,
        "assignment_unit": "structure_id",
        "rows": len(rows),
        "eligible_precalibration_rows": len(eligible_rows),
        "eligible_rows_by_split": dict(sorted(Counter(row["split"] for row in eligible_rows).items())),
        "all_rows_by_split": dict(sorted(Counter(row["split"] for row in rows).items())),
        "structures_by_split": {
            split: len({row["structure_id"] for row in rows if row["split"] == split})
            for split in ("calibration", "held_out")
        },
        "status": "locked_before_occupancy_or_target_selection",
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(rows)} split rows; {len(eligible_rows)} are eligible before calibration.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

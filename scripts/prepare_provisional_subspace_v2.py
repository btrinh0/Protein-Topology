from __future__ import annotations

"""Filter the V2 mapping table into a documented provisional subspace."""

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare a provisional alpha/beta V2 subspace.")
    parser.add_argument("--mapping", type=Path, default=Path("data/processed/v2/domain_mapping_v2.csv"))
    parser.add_argument("--astral-fasta", type=Path, default=Path("data/raw/classifications/scope_astral_40.fa"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/domain_mapping_alpha_beta_70_160_provisional.csv"))
    parser.add_argument("--min-length", type=int, default=70)
    parser.add_argument("--max-length", type=int, default=160)
    parser.add_argument("--scope-class", default="c")
    parser.add_argument("--max-structures", type=int, default=0)
    return parser.parse_args()


def read_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    current_id = ""
    current_length = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    lengths[current_id] = current_length
                current_id = line[1:].split()[0]
                current_length = 0
            else:
                current_length += len(line)
    if current_id:
        lengths[current_id] = current_length
    return lengths


def main() -> int:
    args = parse_args()
    if args.min_length < 1 or args.max_length < args.min_length:
        raise SystemExit("Invalid length window.")
    lengths = read_lengths(args.astral_fasta)
    with args.mapping.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    selected = [
        row for row in rows
        if row.get("scope_class_id") == args.scope_class
        and args.min_length <= lengths.get(row.get("source_domain_id", ""), 0) <= args.max_length
    ]
    if args.max_structures > 0:
        structure_ids: list[str] = []
        for row in selected:
            if row["structure_id"] not in structure_ids:
                structure_ids.append(row["structure_id"])
            if len(structure_ids) >= args.max_structures:
                break
        allowed = set(structure_ids)
        selected = [row for row in selected if row["structure_id"] in allowed]
    if not rows:
        raise SystemExit(f"No rows found in {args.mapping}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(selected)
    print(f"Wrote {len(selected)} domains across {len({row['structure_id'] for row in selected})} structures.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

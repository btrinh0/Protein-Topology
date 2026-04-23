from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable


DEFAULT_ATLAS_COLUMNS = [
    "domain_id",
    "structure_id",
    "chain_id",
    "residue_start",
    "residue_end",
    "domain_length",
    "kingdom",
    "scope_id",
    "cath_id",
    "ecod_id",
    "topology_representation",
    "ps_fraction",
    "pp_fraction",
    "x_fraction",
    "occupancy_tier",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a first normalized atlas table by joining domain metadata and chain topology summaries."
    )
    parser.add_argument(
        "--chain-topology",
        type=Path,
        default=Path("data/processed/chain_topology_manifest.csv"),
        help="Chain topology manifest produced by compute_circuit_topology.py.",
    )
    parser.add_argument(
        "--domain-mapping",
        type=Path,
        required=True,
        help="Normalized domain mapping CSV.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/processed/domain_atlas_table.csv"),
        help="CSV path for the normalized atlas table.",
    )
    parser.add_argument(
        "--failures-output",
        type=Path,
        default=Path("data/processed/atlas_join_failures.csv"),
        help="CSV path for domain rows that could not be joined cleanly.",
    )
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def normalized_topology_representation(row: dict[str, str]) -> str:
    s = row.get("s_fraction", "")
    p = row.get("p_fraction", "")
    x = row.get("x_fraction", "")
    return f"S={s};P={p};X={x}"


def main() -> int:
    args = parse_args()
    topology_rows = read_csv(args.chain_topology)
    domain_rows = read_csv(args.domain_mapping)

    topology_index = {
        (row["structure_id"], row["chain_id"]): row
        for row in topology_rows
    }

    atlas_rows: list[dict[str, str]] = []
    failure_rows: list[dict[str, str]] = []

    for domain in domain_rows:
        key = (domain.get("structure_id", ""), domain.get("chain_id", ""))
        topology = topology_index.get(key)
        if topology is None:
            failure_rows.append(
                {
                    "domain_id": domain.get("domain_id", ""),
                    "structure_id": domain.get("structure_id", ""),
                    "chain_id": domain.get("chain_id", ""),
                    "error": "Missing chain topology summary for domain join.",
                }
            )
            continue

        atlas_rows.append(
            {
                "domain_id": domain.get("domain_id", ""),
                "structure_id": domain.get("structure_id", ""),
                "chain_id": domain.get("chain_id", ""),
                "residue_start": domain.get("residue_start", ""),
                "residue_end": domain.get("residue_end", ""),
                "domain_length": domain.get("domain_length", ""),
                "kingdom": domain.get("kingdom", ""),
                "scope_id": domain.get("scope_id", ""),
                "cath_id": domain.get("cath_id", ""),
                "ecod_id": domain.get("ecod_id", ""),
                "topology_representation": normalized_topology_representation(topology),
                "ps_fraction": topology.get("s_fraction", ""),
                "pp_fraction": topology.get("p_fraction", ""),
                "x_fraction": topology.get("x_fraction", ""),
                "occupancy_tier": "",
            }
        )

    write_csv(args.output, DEFAULT_ATLAS_COLUMNS, atlas_rows)
    write_csv(
        args.failures_output,
        ["domain_id", "structure_id", "chain_id", "error"],
        failure_rows,
    )

    print(
        f"Wrote {len(atlas_rows)} atlas rows and {len(failure_rows)} join failures."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

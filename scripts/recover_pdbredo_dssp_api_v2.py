from __future__ import annotations

"""Recover missing PDB-REDO DSSP records through its documented API."""

import argparse
import csv
import gzip
import hashlib
import json
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


API_URL = "https://pdb-redo.eu/dssp/do"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Request DSSP annotations for entries missing from the PDB-REDO archive.")
    parser.add_argument("--errors", type=Path, default=Path("data/processed/v2/domain_topology_pdbredo_dssp_errors_v2.csv"))
    parser.add_argument("--split", type=Path, default=Path("data/processed/v2/calibration_split_v2.csv"))
    parser.add_argument("--topology", type=Path, default=Path("data/processed/v2/domain_topology_alpha_beta_70_160_calibration_pool.csv"))
    parser.add_argument("--cache-dir", type=Path, default=Path("data/raw/v2/pdbredo_dssp"))
    parser.add_argument("--max-structures", type=int, default=0)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def post_dssp(text: str) -> bytes:
    form = urllib.parse.urlencode({"data": text, "format": "dssp"}).encode("utf-8")
    request = urllib.request.Request(API_URL, data=form, method="POST", headers={"Content-Type": "application/x-www-form-urlencoded", "User-Agent": "TCT-v2/1.0"})
    with urllib.request.urlopen(request, timeout=180) as response:
        payload = response.read()
    if not payload or b"TOTAL NUMBER OF RESIDUES" not in payload:
        raise ValueError("PDB-REDO API response is not a legacy DSSP file")
    return payload


def main() -> int:
    args = parse_args()
    split_rows = {row["domain_uid"]: row for row in read_csv(args.split)}
    topology_rows = {row["domain_uid"]: row for row in read_csv(args.topology)}
    errors = read_csv(args.errors)
    structure_ids = sorted({row["structure_id"] for row in errors})
    if args.max_structures > 0:
        structure_ids = structure_ids[: args.max_structures]
    results: list[dict[str, object]] = []
    failures: list[dict[str, str]] = []
    for structure_id in structure_ids:
        source_rows = [topology_rows.get(row["domain_uid"], split_rows.get(row["domain_uid"], {})) for row in errors if row["structure_id"] == structure_id]
        source_path = Path(source_rows[0].get("source_path", "")) if source_rows else Path()
        try:
            if not source_path.exists():
                raise FileNotFoundError(f"missing local source structure: {source_path}")
            with gzip.open(source_path, "rt", encoding="utf-8") if source_path.name.endswith(".gz") else source_path.open("r", encoding="utf-8") as handle:
                source_text = handle.read()
            payload = post_dssp(source_text)
            destination = args.cache_dir / f"{structure_id.lower()}.dssp"
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(payload)
            metadata = {
                "structure_id": structure_id,
                "url": API_URL,
                "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
                "sha256": hashlib.sha256(payload).hexdigest(),
                "bytes": len(payload),
                "source_structure_path": str(source_path),
                "retrieval_mode": "api_post",
            }
            destination.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
            results.append(metadata)
        except (OSError, ValueError, urllib.error.URLError, urllib.error.HTTPError) as exc:
            failures.append({"structure_id": structure_id, "error": str(exc)})
    summary = {
        "api_url": API_URL,
        "structures_requested": len(structure_ids),
        "structures_recovered": len(results),
        "failures": failures,
        "retrieval_mode": "api_post",
        "notes": ["Recovered records are kept separate from archive downloads through sidecar metadata.", "The resulting DSSP states still require the same structural-similarity gate as archive records."],
    }
    output = args.cache_dir / "api_recovery_summary_v2.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Recovered {len(results)} API DSSP records with {len(failures)} failures.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

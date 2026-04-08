from __future__ import annotations

import argparse
import csv
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests
import urllib3

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

RCSB_MMCIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif.gz"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download mmCIF files from RCSB PDB for domains in the mapping."
    )
    parser.add_argument(
        "--domain-mapping",
        type=Path,
        default=Path("data/processed/domain_mapping.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/raw/structures"),
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=4,
    )
    return parser.parse_args()


def unique_pdb_ids(path: Path) -> list[str]:
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        return sorted({row["structure_id"].lower() for row in reader})


def download_one(pdb_id: str, output_dir: Path) -> tuple[str, bool, str]:
    dest = output_dir / f"{pdb_id}.cif.gz"
    if dest.exists() and dest.stat().st_size > 100:
        return pdb_id, True, "already exists"
    url = RCSB_MMCIF_URL.format(pdb_id=pdb_id.upper())
    try:
        resp = requests.get(url, timeout=60, verify=False)
        resp.raise_for_status()
        dest.write_bytes(resp.content)
        return pdb_id, True, ""
    except Exception as exc:
        return pdb_id, False, str(exc)


def main() -> int:
    args = parse_args()
    pdb_ids = unique_pdb_ids(args.domain_mapping)
    print(f"Need {len(pdb_ids)} unique PDB structures.")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    already = sum(1 for pid in pdb_ids
                  if (args.output_dir / f"{pid}.cif.gz").exists()
                  and (args.output_dir / f"{pid}.cif.gz").stat().st_size > 100)
    remaining = len(pdb_ids) - already
    print(f"Already downloaded: {already}, remaining: {remaining}")

    if remaining == 0:
        print("All structures already present.")
        return 0

    successes = already
    failures = []
    done_count = already

    with ThreadPoolExecutor(max_workers=args.max_workers) as pool:
        futures = {
            pool.submit(download_one, pid, args.output_dir): pid
            for pid in pdb_ids
            if not (args.output_dir / f"{pid}.cif.gz").exists()
            or (args.output_dir / f"{pid}.cif.gz").stat().st_size <= 100
        }
        for future in as_completed(futures):
            pdb_id, ok, err = future.result()
            done_count += 1
            if ok:
                successes += 1
            else:
                failures.append((pdb_id, err))
            if done_count % 200 == 0:
                print(f"  Progress: {done_count}/{len(pdb_ids)}")

    print(f"\nDownloaded {successes}/{len(pdb_ids)} structures.")
    if failures:
        fail_path = args.output_dir / "download_failures.csv"
        with fail_path.open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["pdb_id", "error"])
            w.writerows(failures)
        print(f"{len(failures)} failures recorded in {fail_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

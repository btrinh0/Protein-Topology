from __future__ import annotations

import argparse
from pathlib import Path

import requests
import urllib3

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

DOWNLOADS = [
    {
        "url": "https://scop.berkeley.edu/downloads/parse/dir.cla.scope.2.08-stable.txt",
        "dest": "scope_cla.txt",
        "label": "SCOPe 2.08 classification",
    },
    {
        "url": "https://scop.berkeley.edu/downloads/parse/dir.des.scope.2.08-stable.txt",
        "dest": "scope_des.txt",
        "label": "SCOPe 2.08 descriptions",
    },
    {
        "url": "https://scop.berkeley.edu/downloads/scopeseq-2.08/astral-scopedom-seqres-gd-sel-gs-bib-40-2.08.fa",
        "dest": "scope_astral_40.fa",
        "label": "SCOPe ASTRAL 40% domain representatives",
    },
    {
        "url": "https://www.cathdb.info/version/v4_3_0/api/rest/domain_summary",
        "dest": "cath_domain_list.txt",
        "label": "CATH v4.3 domain list",
    },
    {
        "url": "http://prodata.swmed.edu/ecod/distributions/ecod.latest.domains.txt",
        "dest": "ecod_domains.txt",
        "label": "ECOD domain list",
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download SCOPe, CATH, and ECOD classification files."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/raw/classifications"),
    )
    return parser.parse_args()


def download_file(url: str, dest: Path, label: str) -> bool:
    print(f"Downloading {label} ...")
    try:
        resp = requests.get(url, timeout=300, stream=True, verify=False)
        resp.raise_for_status()
    except Exception as exc:
        print(f"  FAILED: {exc}")
        return False
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("wb") as f:
        for chunk in resp.iter_content(chunk_size=1 << 16):
            f.write(chunk)
    size_mb = dest.stat().st_size / (1 << 20)
    print(f"  Saved {dest} ({size_mb:.2f} MB)")
    return True


def main() -> int:
    args = parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

    successes = 0
    for item in DOWNLOADS:
        ok = download_file(item["url"], out / item["dest"], item["label"])
        if ok:
            successes += 1

    print(f"\n{successes}/{len(DOWNLOADS)} downloads completed.")
    return 0 if successes >= 2 else 1


if __name__ == "__main__":
    raise SystemExit(main())

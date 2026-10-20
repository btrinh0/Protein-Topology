"""Route A, step A1: fetch and verify ProteinGym v1.3 substitutions, and write the download manifest.

Each file goes into its own directory under data/raw/proteingym/, which must be empty before the
first download. Downloads are checked against the size and SHA-256 recorded here, so a changed
upstream file stops the run instead of flowing into the analysis. Archive members are checked for
unsafe paths before extraction, and nothing in the archive is executed.

Sources, both from the ProteinGym repository (https://github.com/OATML-Markslab/ProteinGym):
  * DMS_ProteinGym_substitutions.zip, release v1.3, as listed in the repository README.
  * reference_files/DMS_substitutions.csv, the reference table on the main branch (blob 8d1ea9a1).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/proteingym"
ARCHIVE_DIR = RAW / "v1.3_archive"
EXTRACT_DIR = RAW / "v1.3_extracted"
REFERENCE_DIR = RAW / "reference_github_main"
MANIFEST = ROOT / "results/route_a_20261012/download_manifest.json"

ZIP_NAME = "DMS_ProteinGym_substitutions.zip"
ZIP_URL = f"https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/{ZIP_NAME}"
ZIP_BYTES = 43021128
ZIP_SHA256 = "3a83766254ac9ac9984ec25cb73c6e010ea4418f5e35f143933e6b6e6473b921"
REF_NAME = "DMS_substitutions.csv"
REF_URL = f"https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/main/reference_files/{REF_NAME}"
REF_BYTES = 208734
REF_SHA256 = "a8f498011532a74aa9fe556a50555a75e928c5837d19c06a87592ae04049b308"
REF_GIT_BLOB = "8d1ea9a19c0404b511cd24378b25c2a5f86f10e9"
EXPECTED_ASSAYS = 217


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch(url: str, target: Path) -> None:
    """Download into a directory that must be empty. Only called when the file is missing."""
    if target.parent.exists() and any(target.parent.iterdir()) and not target.exists():
        raise SystemExit(f"{target.parent} is not empty; the download needs its own empty directory")
    target.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "project6-route-a-fetch"})
    with urllib.request.urlopen(request, timeout=300) as response, target.open("wb") as handle:
        while True:
            block = response.read(1024 * 1024)
            if not block:
                break
            handle.write(block)


def verify(path: Path, expected_bytes: int, expected_sha: str) -> dict:
    size, digest = path.stat().st_size, sha256(path)
    if size != expected_bytes or digest != expected_sha:
        raise SystemExit(f"{path.name}: size {size} sha256 {digest}; expected {expected_bytes} {expected_sha}")
    return {"path": path.relative_to(ROOT).as_posix(), "bytes": size, "sha256": digest}


def safe_extract(archive: Path, destination: Path) -> int:
    """Extract only regular CSV members with plain relative names into an empty directory."""
    if destination.exists() and any(destination.iterdir()):
        return sum(1 for _ in destination.rglob("*.csv"))
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as bundle:
        for member in bundle.infolist():
            name = member.filename
            if member.is_dir():
                continue
            if name.startswith("/") or ".." in Path(name).parts or ":" in name:
                raise SystemExit(f"Unsafe archive member: {name}")
            if not name.endswith(".csv"):
                raise SystemExit(f"Unexpected archive member: {name}")
            bundle.extract(member, destination)
    return sum(1 for _ in destination.rglob("*.csv"))


def run() -> dict:
    zip_path = ARCHIVE_DIR / ZIP_NAME
    ref_path = REFERENCE_DIR / REF_NAME
    if not zip_path.exists():
        fetch(ZIP_URL, zip_path)
    if not ref_path.exists():
        fetch(REF_URL, ref_path)
    archive = verify(zip_path, ZIP_BYTES, ZIP_SHA256)
    reference = verify(ref_path, REF_BYTES, REF_SHA256)
    assays = safe_extract(zip_path, EXTRACT_DIR)
    if assays != EXPECTED_ASSAYS:
        raise SystemExit(f"Expected {EXPECTED_ASSAYS} assay files, found {assays}")
    manifest = {
        "step": "A1 acquire ProteinGym",
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "proteingym_release": "v1.3 (latest, per the repository README)",
        "archive": {"url": ZIP_URL, **archive, "extracted_to": EXTRACT_DIR.relative_to(ROOT).as_posix(),
                    "assay_files": assays},
        "reference_table": {"url": REF_URL, "github_blob_sha": REF_GIT_BLOB, **reference},
        "gitignore": "data/raw/ (already listed in .gitignore; checked with git check-ignore)",
        "outcome_note": "Downloaded and extracted only. No DMS_score column has been read in this step.",
    }
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    argparse.ArgumentParser(description=__doc__).parse_args()
    print(json.dumps(run(), indent=2))


if __name__ == "__main__":
    main()

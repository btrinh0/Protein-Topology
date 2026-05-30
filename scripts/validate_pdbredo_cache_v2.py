from __future__ import annotations

"""Verify cached PDB-REDO DSSP bytes and parseability."""

import argparse
import hashlib
import json
from pathlib import Path

from Bio.PDB.DSSP import make_dssp_dict


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate the local PDB-REDO DSSP cache.")
    parser.add_argument("--cache-dir", type=Path, default=Path("data/raw/v2/pdbredo_dssp"))
    parser.add_argument("--summary-output", type=Path, default=Path("data/processed/v2/pdbredo_cache_validation_v2.json"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    files = sorted(args.cache_dir.glob("*.dssp"))
    failures: list[dict[str, str]] = []
    parsed_residue_counts: list[int] = []
    for path in files:
        metadata_path = path.with_suffix(".json")
        try:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            payload = path.read_bytes()
            digest = hashlib.sha256(payload).hexdigest()
            if digest != metadata.get("sha256") or len(payload) != int(metadata.get("bytes", -1)):
                raise ValueError("sha256 or byte-count mismatch")
            dssp_values, _ = make_dssp_dict(str(path))
            if not dssp_values:
                raise ValueError("parsed DSSP record is empty")
            parsed_residue_counts.append(len(dssp_values))
        except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as exc:
            failures.append({"file": str(path), "error": str(exc)})
    summary = {
        "cache_dir": str(args.cache_dir),
        "dssp_files": len(files),
        "validated_files": len(parsed_residue_counts),
        "failures": failures,
        "residue_count_min": min(parsed_residue_counts) if parsed_residue_counts else None,
        "residue_count_max": max(parsed_residue_counts) if parsed_residue_counts else None,
        "status": "pass" if files and not failures and len(parsed_residue_counts) == len(files) else "fail",
    }
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Validated {len(parsed_residue_counts)}/{len(files)} cached DSSP files with {len(failures)} failures.")
    return 0 if summary["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())

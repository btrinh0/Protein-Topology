"""Fetch only the existing development assays' public ProteinGym predictions."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
from pathlib import Path
import time
import urllib.request
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
TABLE = ROOT / "results/beat_am_20261016/features/l1_table.npz"
OUT = ROOT / "data/external/predictor_comparison_20261026"
ARCHIVE_URL = "https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/zero_shot_substitutions_scores.zip"
REPO_API = "https://api.github.com/repos/OATML-Markslab/ProteinGym/commits/main"
REQUESTED = ("GEMME", "VESPA", "VESPAl", "VespaG")


def sha(data):
    return hashlib.sha256(data).hexdigest()


def dump(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def request(url, *, method="GET", headers=None):
    supplied = {"User-Agent": "protein-rule-audit-comparison/1.0", "Accept-Encoding": "identity"}
    supplied.update(headers or {})
    last = None
    for attempt in range(3):
        try:
            return urllib.request.urlopen(urllib.request.Request(url, method=method, headers=supplied), timeout=90)
        except Exception as exc:
            last = exc
            if attempt < 2:
                time.sleep(1 + attempt)
    raise last


class HTTPRangeFile(io.RawIOBase):
    def __init__(self, url):
        self.url, self.pos, self.bytes_read, self.requests = url, 0, 0, 0
        with request(url, method="HEAD") as response:
            self.length = int(response.headers["Content-Length"])
            self.etag = response.headers.get("ETag")
            self.last_modified = response.headers.get("Last-Modified")
        self.requests += 1

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=0):
        value = offset if whence == 0 else self.pos + offset if whence == 1 else self.length + offset
        if value < 0:
            raise ValueError("Negative seek")
        self.pos = value
        return self.pos

    def read(self, size=-1):
        if size < 0:
            size = self.length - self.pos
        size = min(size, self.length - self.pos)
        if size <= 0:
            return b""
        start, end = self.pos, self.pos + size - 1
        headers = {"Range": f"bytes={start}-{end}"}
        if self.etag:
            headers["If-Match"] = self.etag
        with request(self.url, headers=headers) as response:
            if response.status != 206:
                raise ValueError("Server did not honor byte range; refusing full archive download")
            content_range = response.headers.get("Content-Range", "")
            if content_range != f"bytes {start}-{end}/{self.length}":
                raise ValueError(f"Unexpected range response: {content_range}")
            data = response.read(size + 1)
        if len(data) != size:
            raise ValueError("Truncated range response")
        self.pos += size
        self.bytes_read += size
        self.requests += 1
        return data


def fetch_source_metadata(out):
    marker = out / "upstream_manifest.json"
    if marker.exists():
        existing = json.loads(marker.read_text(encoding="utf-8"))
        for entry in existing["files"]:
            if sha((out / entry["filename"]).read_bytes()) != entry["sha256"]:
                raise ValueError("Cached upstream metadata hash mismatch")
        commit, files = existing["github_commit"], existing["files"]
    else:
        with request(REPO_API) as response:
            commit = json.load(response)["sha"]
        files = []
    for name in ("config.json", "README.md", "proteingym/merge.py", "scripts/scoring_DMS_zero_shot/merge_all_scores.sh"):
        filename = "upstream_" + name.replace("/", "_")
        if any(item["filename"] == filename for item in files):
            continue
        url = f"https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/{commit}/{name}"
        with request(url) as response:
            data = response.read()
        (out / filename).write_bytes(data)
        files.append(dict(filename=filename, url=url, sha256=sha(data), bytes=len(data)))
    result = dict(retrieved_utc=datetime.now(timezone.utc).isoformat(), github_commit=commit, files=files)
    dump(marker, result)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--inspect-only", action="store_true")
    parser.add_argument("--metadata-only", action="store_true")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    with np.load(TABLE, allow_pickle=False) as table:
        assays = sorted(set(table["dms_id"].astype(str)))
    metadata = fetch_source_metadata(args.out)
    if args.metadata_only:
        print(json.dumps(metadata), flush=True)
        return
    manifest_path = args.out / "download_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists() else dict(
        study="external comparators on existing exposed development assays only",
        requested_models=list(REQUESTED), requested_assays=assays,
        source_table_sha256=sha(TABLE.read_bytes()), entries=[], errors=[])
    if manifest["requested_assays"] != assays or manifest["source_table_sha256"] != sha(TABLE.read_bytes()):
        raise ValueError("Source cohort changed")
    old = {entry["assay"]: entry for entry in manifest["entries"]}
    manifest["errors"] = []
    remote = HTTPRangeFile(ARCHIVE_URL)
    with zipfile.ZipFile(remote) as archive:
        candidates = {}
        for info in archive.infolist():
            basename = Path(info.filename).name
            if basename.endswith(".csv") and basename[:-4] in assays:
                if basename[:-4] in candidates:
                    raise ValueError("Duplicate assay path in remote archive")
                candidates[basename[:-4]] = info
        manifest["archive"] = dict(url=ARCHIVE_URL, bytes=remote.length, etag=remote.etag,
            last_modified=remote.last_modified, matching_entries=len(candidates))
        manifest["upstream_commit"] = metadata["github_commit"]
        dump(args.out / "archive_index.json", [dict(assay=a, member=x.filename, crc32=x.CRC,
            compressed_bytes=x.compress_size, uncompressed_bytes=x.file_size) for a,x in sorted(candidates.items())])
        print(f"Archive indexed: {len(candidates)}/{len(assays)} requested assays", flush=True)
        if args.inspect_only:
            selected = [a for a in assays if a.startswith(("ACE2_", "OPSD_", "PTEN_"))][:3]
        else:
            selected = assays
        for number, assay in enumerate(selected, 1):
            try:
                if assay in old:
                    item = old[assay]
                    path = args.out / item["retained_file"]
                    if path.exists() and sha(path.read_bytes()) == item["retained_sha256"]:
                        print(f"[{number}/{len(selected)}] cached {assay}", flush=True)
                        continue
                    raise ValueError("Cached retained-score hash mismatch")
                if assay not in candidates:
                    raise ValueError("Assay absent from official archive")
                info = candidates[assay]
                if info.file_size > 750_000_000:
                    raise ValueError("Unexpectedly large member")
                payload = archive.read(info)  # ZIP extraction validates each member's CRC.
                reader = csv.DictReader(io.StringIO(payload.decode("utf-8-sig")))
                header = reader.fieldnames or []
                print(f"[{number}/{len(selected)}] {assay}: score columns {[x for x in header if any(t.lower() in x.lower() for t in REQUESTED)]}", flush=True)
                if "mutant" not in header:
                    raise ValueError("Mutation identifier missing")
                keep = ["mutant"] + [name for name in REQUESTED if name in header]
                # Inspect only the three declared external score sources.
                rows = [{name: row[name] for name in keep} for row in reader]
                filename = f"{assay}.scores.csv.gz"
                path = args.out / filename
                with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=keep)
                    writer.writeheader()
                    writer.writerows(rows)
                entry = dict(assay=assay, archive_member=info.filename, archive_crc32=info.CRC,
                    source_uncompressed_sha256=sha(payload), source_uncompressed_bytes=len(payload),
                    source_columns=header, retained_columns=keep, rows=len(rows),
                    retained_file=filename, retained_sha256=sha(path.read_bytes()))
                old[assay] = entry
            except Exception as exc:
                manifest["errors"].append(dict(assay=assay, error=f"{type(exc).__name__}: {exc}"))
                print(f"FAILED {assay}: {exc}", flush=True)
            manifest["entries"] = [old[a] for a in sorted(old)]
            manifest["retrieved_bytes_this_run"] = remote.bytes_read
            manifest["http_requests_this_run"] = remote.requests
            manifest["updated_utc"] = datetime.now(timezone.utc).isoformat()
            dump(manifest_path, manifest)
    print(json.dumps(dict(assays_retained=len(old), errors=manifest["errors"],
        bytes_this_run=remote.bytes_read, http_requests=remote.requests)), flush=True)
    if manifest["errors"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

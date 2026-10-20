"""Resolve ProteinGym UniProt entry names to accessions (Route A, step A2).

ProteinGym's `UniProt_ID` column holds UniProt entry names such as `PAI1_HUMAN`, while the
topology census is keyed by accession such as `P05121`. This script asks the UniProt REST
search for every distinct entry name in the reference table and writes the answer, with the
query and the UniProt release, to a JSON file. Nothing here reads DMS scores.

Input:  data/raw/proteingym/reference_github_main/DMS_substitutions.csv
Output: data/raw/proteingym/uniprot_entry_map_v1.3.json (raw API answer, gitignored with data/raw/)
"""

from __future__ import annotations

import argparse
import csv
import json
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REFERENCE = ROOT / "data/raw/proteingym/reference_github_main/DMS_substitutions.csv"
DEFAULT_OUT = ROOT / "data/raw/proteingym/uniprot_entry_map_v1.3.json"
SEARCH = "https://rest.uniprot.org/uniprotkb/search"
FIELDS = "accession,id,gene_primary,length,reviewed,organism_name"
BATCH = 20


def distinct_entries(reference: Path) -> list[str]:
    with reference.open(encoding="utf-8", newline="") as handle:
        names = {row["UniProt_ID"].strip() for row in csv.DictReader(handle)}
    return sorted(name for name in names if name)


def query_batch(entries: list[str]) -> tuple[list[dict], str]:
    query = " OR ".join(f"(id:{entry})" for entry in entries)
    url = SEARCH + "?" + urllib.parse.urlencode(
        {"query": query, "fields": FIELDS, "format": "tsv", "size": 500})
    request = urllib.request.Request(url, headers={"User-Agent": "project6-route-a-resolver"})
    with urllib.request.urlopen(request, timeout=120) as response:
        release = response.headers.get("X-UniProt-Release", "")
        body = response.read().decode("utf-8")
    lines = body.splitlines()
    header = lines[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in lines[1:] if line]
    return rows, release


def run(reference: Path, out: Path) -> dict:
    entries = distinct_entries(reference)
    requested = set(entries)
    resolved: dict[str, dict] = {}
    returned_not_requested: dict[str, dict] = {}
    obsolete_skipped: dict[str, list[str]] = {}
    releases = set()
    for start in range(0, len(entries), BATCH):
        batch = entries[start:start + BATCH]
        rows, release = query_batch(batch)
        releases.add(release)
        for row in rows:
            name = row["Entry Name"]
            # The id search also returns retired accessions that keep the name, with blank fields.
            if not row["Length"]:
                obsolete_skipped.setdefault(name, []).append(row["Entry"])
                continue
            record = {
                "accession": row["Entry"],
                "gene": row["Gene Names (primary)"],
                "length": int(row["Length"]),
                "reviewed": row["Reviewed"] == "reviewed",
                "organism": row["Organism"],
            }
            # The id search also returns neighbouring entries; only requested names are kept.
            if name not in requested:
                returned_not_requested[name] = record
                continue
            if name in resolved:
                raise SystemExit(f"Entry name {name} matches more than one active accession")
            resolved[name] = record
    unresolved = [entry for entry in entries if entry not in resolved]
    payload = {
        "source": SEARCH,
        "fields": FIELDS,
        "reference_table": str(reference.relative_to(ROOT).as_posix()),
        "query_utc": datetime.now(timezone.utc).isoformat(),
        "uniprot_release": sorted(releases),
        "n_entry_names_requested": len(entries),
        "n_resolved": len(resolved),
        "unresolved": unresolved,
        "returned_not_requested_ignored": dict(sorted(returned_not_requested.items())),
        "retired_accessions_skipped": dict(sorted(obsolete_skipped.items())),
        "entries": dict(sorted(resolved.items())),
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    payload = run(args.reference, args.out)
    print(json.dumps({k: v for k, v in payload.items() if k != "entries"}, indent=2))


if __name__ == "__main__":
    main()

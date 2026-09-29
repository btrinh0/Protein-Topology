"""Complete the ortholog sweep by querying every MaveDB target name.

The score-set search endpoint returns at most 100 results and exposes no paging
parameter, so the first sweep screened only 100 deposits and found Kir2.1 only by
direct URN lookup. This queries each of the ~1,157 target gene names in turn,
caching every response, so the screen covers the whole catalogue.

Applies the ortholog rule frozen in PHASE_B_INCLUSION_RULES_FROZEN.md section 4a.
Counts only; no measured effect is read.
"""

from __future__ import annotations

import argparse
import json
import time
import urllib.error
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from Bio import Align

from extract_am_flank_scores import build_flank_map
from run_gate1 import read_topology
from search_ortholog_score_sets import (align_to_census, build_kmer_index,
                                        candidates_for, translate)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TOPOLOGY = ROOT / "data" / "raw" / "topology_viewer.html"
DEFAULT_CACHE = ROOT / "data" / "external" / "mavedb"
DEFAULT_OUT = ROOT / "results" / "ortholog_sweep_20260927"
API = "https://api.mavedb.org"


def search_name(name: str, cache: Path, pause: float) -> list[dict]:
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in name)[:80]
    path = cache / "name_search" / f"{safe}.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8")).get("scoreSets", [])
    path.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(
        API + "/api/v1/score-sets/search",
        data=json.dumps({"text": name}).encode(),
        headers={"Content-Type": "application/json"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                payload = json.loads(response.read().decode())
            path.write_text(json.dumps(payload), encoding="utf-8")
            time.sleep(pause)
            return payload.get("scoreSets", [])
        except (urllib.error.URLError, TimeoutError):
            time.sleep(2 * (attempt + 1))
    return []


def run(args: argparse.Namespace) -> dict:
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cache = Path(args.cache)
    proteins = read_topology(Path(args.topology))
    index = build_kmer_index(proteins)
    aligner = Align.PairwiseAligner(mode="global", open_gap_score=-11,
                                    extend_gap_score=-1, match_score=2, mismatch_score=-1)
    symbols = {alias for record in proteins.values() for alias in record["_aliases"]}
    names = json.loads((cache / "target_gene_names.json").read_text(encoding="utf-8"))

    seen_urns: set[str] = set()
    findings, qc = [], Counter()
    for position, name in enumerate(names, start=1):
        qc["names_queried"] += 1
        for score_set in search_name(name, cache, args.pause):
            urn = score_set.get("urn") or ""
            if urn in seen_urns:
                continue
            seen_urns.add(urn)
            qc["score_sets_screened"] += 1
            for target in score_set.get("targetGenes") or []:
                target_name = str(target.get("name") or "").strip()
                if target_name.upper() in symbols:
                    qc["matched_by_gene_symbol_already"] += 1
                    continue
                block = target.get("targetSequence") or {}
                raw = (block.get("sequence") or "").strip().upper()
                if len(raw) < 60:
                    continue
                kind = str(block.get("sequenceType") or "").lower()
                attempts = ([raw] if kind.startswith("protein")
                            else [translate(raw, frame) for frame in range(3)])
                best = None
                for sequence in attempts:
                    sequence = sequence.split("*")[0]
                    if len(sequence) < 60:
                        continue
                    for accession in candidates_for(sequence, index):
                        identity, mapping, gaps = align_to_census(
                            sequence, proteins[accession]["sequence"], aligner)
                        if best is None or identity > best[1]:
                            best = (accession, identity, mapping, gaps)
                if best is None:
                    continue
                accession, identity, mapping, gaps = best
                if identity < args.min_identity:
                    qc["below_identity_threshold"] += 1
                    continue
                record = proteins[accession]
                flank_positions = set(build_flank_map(record, 10))
                findings.append({
                    "urn": urn, "title": (score_set.get("title") or "")[:120],
                    "target_name": target_name, "census_uniprot": accession,
                    "census_gene": record.get("gene"), "identity": round(identity, 4),
                    "n_tmds": len(record["tmds"]),
                    "flank10_positions": len(flank_positions),
                    "flank10_identical_usable": len(flank_positions & set(mapping)),
                    "flank10_in_alignment_gap": len(flank_positions & gaps),
                    "topology_transfer_clean": len(flank_positions & gaps) == 0,
                    "cross_species_flag": True,
                })
                qc["ortholog_candidates"] += 1
        if position % 100 == 0:
            print(f"  ...{position}/{len(names)} names, {qc['score_sets_screened']} score sets, "
                  f"{qc['ortholog_candidates']} candidates", flush=True)

    findings.sort(key=lambda f: -f["flank10_identical_usable"])
    summary = {
        "analysis": "Complete ortholog sweep over all MaveDB target names",
        "reads_dms_effects": False,
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "rule_source": "PHASE_B_INCLUSION_RULES_FROZEN.md section 4a",
        "min_identity": args.min_identity,
        "quality_control": dict(qc),
        "candidates": findings,
    }
    (out_dir / "ortholog_sweep.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"quality_control": dict(qc)}, indent=2))
    for finding in findings[:40]:
        print(f"  {finding['target_name']:<22} -> {finding['census_gene']:<8} "
              f"identity={finding['identity']:.3f}  flank10 usable "
              f"{finding['flank10_identical_usable']}/{finding['flank10_positions']}  "
              f"clean={finding['topology_transfer_clean']}  {finding['urn']}")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", default=str(DEFAULT_TOPOLOGY))
    parser.add_argument("--cache", default=str(DEFAULT_CACHE))
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    parser.add_argument("--min-identity", type=float, default=0.95)
    parser.add_argument("--pause", type=float, default=0.35, help="seconds between API calls")
    run(parser.parse_args())


if __name__ == "__main__":
    main()

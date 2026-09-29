"""Ortholog-aware rerun of the Phase B search.

Phase B matched MaveDB deposits to the census by human gene symbol. That hid
`mKir2.1`, a mouse deposit of a protein whose human ortholog KCNJ2 is in the
census, and the same bug could hide other rodent scans.

This enumerates every MaveDB score set, aligns each target sequence to the census
by global alignment, and applies the ortholog rule frozen in
PHASE_B_INCLUSION_RULES_FROZEN.md section 4a:

* at least 95% identity over the mapped region
* topology conserved, checked here as: the alignment covers every transmembrane
  segment and its +/-10 flank with no gaps, so census topology transfers position
  by position
* only positions where the two reference residues are identical are usable
* cross-species sets flagged

Counts only. No measured effect is read.
"""

from __future__ import annotations

import argparse
import json
import re
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from Bio import Align

from extract_am_flank_scores import build_flank_map
from run_gate1 import read_topology

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TOPOLOGY = ROOT / "data" / "raw" / "topology_viewer.html"
DEFAULT_CACHE = ROOT / "data" / "external" / "mavedb"
DEFAULT_OUT = ROOT / "results" / "ortholog_search_20260927"
API = "https://api.mavedb.org"

MIN_IDENTITY = 0.95
KMER = 6
CODON_TABLE = {}
for _codon, _aa in zip(
    [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"],
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"):
    CODON_TABLE[_codon] = _aa


def translate(dna: str, frame: int) -> str:
    dna = re.sub(r"[^ACGT]", "", dna.upper())[frame:]
    return "".join(CODON_TABLE.get(dna[i:i + 3], "X") for i in range(0, len(dna) - 2, 3))


def fetch_all_score_sets(cache: Path) -> list[dict]:
    path = cache / "all_score_sets.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))["scoreSets"]
    request = urllib.request.Request(
        API + "/api/v1/score-sets/search", data=json.dumps({"text": ""}).encode(),
        headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(request, timeout=300) as response:
        payload = json.loads(response.read().decode())
    cache.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    return payload["scoreSets"]


def build_kmer_index(proteins: dict) -> dict[str, set[str]]:
    index: dict[str, set[str]] = defaultdict(set)
    for accession, record in proteins.items():
        sequence = record["sequence"]
        for i in range(len(sequence) - KMER + 1):
            index[sequence[i:i + KMER]].add(accession)
    return index


def candidates_for(sequence: str, index: dict[str, set[str]], top: int = 5) -> list[str]:
    hits: Counter = Counter()
    for i in range(0, max(1, len(sequence) - KMER + 1), 3):
        for accession in index.get(sequence[i:i + KMER], ()):
            hits[accession] += 1
    return [accession for accession, _ in hits.most_common(top)]


def align_to_census(target: str, census: str, aligner: Align.PairwiseAligner):
    alignment = aligner.align(target, census)[0]
    target_aligned, census_aligned = alignment[0], alignment[1]
    mapping: dict[int, int] = {}      # census position -> target position, identical only
    target_pos = census_pos = 0
    identical = compared = 0
    gap_census: set[int] = set()
    for a, b in zip(target_aligned, census_aligned):
        if a != "-":
            target_pos += 1
        if b != "-":
            census_pos += 1
        if a == "-" or b == "-":
            if b != "-":
                gap_census.add(census_pos)
            continue
        compared += 1
        if a == b:
            identical += 1
            mapping[census_pos] = target_pos
    identity = identical / compared if compared else 0.0
    return identity, mapping, gap_census


def run(args: argparse.Namespace) -> dict:
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    proteins = read_topology(Path(args.topology))
    index = build_kmer_index(proteins)
    aligner = Align.PairwiseAligner(mode="global", open_gap_score=-11,
                                    extend_gap_score=-1, substitution_matrix=None,
                                    match_score=2, mismatch_score=-1)

    score_sets = fetch_all_score_sets(Path(args.cache))
    symbols = {alias for record in proteins.values() for alias in record["_aliases"]}

    findings, qc = [], Counter()
    for score_set in score_sets:
        qc["score_sets_seen"] += 1
        for target in score_set.get("targetGenes") or []:
            name = str(target.get("name") or "").strip()
            block = target.get("targetSequence") or {}
            raw = (block.get("sequence") or "").strip().upper()
            if len(raw) < 60:
                continue
            kind = str(block.get("sequenceType") or "").lower()
            attempts = ([raw] if kind.startswith("protein")
                        else [translate(raw, frame) for frame in range(3)])
            if name.upper() in symbols:
                qc["matched_by_gene_symbol_already"] += 1
                continue
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
            flanks = build_flank_map(record, 10)
            flank_positions = set(flanks)
            usable = flank_positions & set(mapping)
            gapped_flank = flank_positions & gaps
            findings.append({
                "urn": score_set.get("urn"),
                "title": (score_set.get("title") or "")[:120],
                "target_name": name,
                "census_uniprot": accession,
                "census_gene": record.get("gene"),
                "identity": round(identity, 4),
                "census_length": record["length"],
                "n_tmds": len(record["tmds"]),
                "flank10_positions": len(flank_positions),
                "flank10_positions_identical_and_usable": len(usable),
                "flank10_positions_in_alignment_gap": len(gapped_flank),
                "topology_transfer_clean": len(gapped_flank) == 0,
                "cross_species_flag": True,
            })
            qc["ortholog_candidates"] += 1

    findings.sort(key=lambda f: -f["flank10_positions_identical_and_usable"])
    summary = {
        "analysis": "Ortholog-aware MaveDB search against the census",
        "reads_dms_effects": False,
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "rule_source": "PHASE_B_INCLUSION_RULES_FROZEN.md section 4a, frozen 2026-09-27",
        "min_identity": args.min_identity,
        "quality_control": dict(qc),
        "candidates": findings,
    }
    (out_dir / "ortholog_candidates.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"quality_control": dict(qc)}, indent=2))
    for finding in findings:
        print(f"  {finding['target_name']:<18} -> {finding['census_gene']:<8} "
              f"{finding['census_uniprot']}  identity={finding['identity']:.3f}  "
              f"flank10 usable {finding['flank10_positions_identical_and_usable']}"
              f"/{finding['flank10_positions']}  clean={finding['topology_transfer_clean']}")
        print(f"      {finding['urn']}  {finding['title']}")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", default=str(DEFAULT_TOPOLOGY))
    parser.add_argument("--cache", default=str(DEFAULT_CACHE))
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    parser.add_argument("--min-identity", type=float, default=MIN_IDENTITY)
    run(parser.parse_args())


if __name__ == "__main__":
    main()

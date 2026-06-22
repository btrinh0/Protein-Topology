"""Recover DMS score sets that Phase B could not place on a census sequence.

Phase B located a deposit by searching for its target sequence as an exact
substring of the census protein. That fails whenever the deposit uses a different
isoform, and eleven targets were lost that way.

Exact substring matching is not needed. Every scored variant already carries its
own reference residue and position, so the offset can be recovered from the
variants themselves: scan every offset that keeps the variants inside the census
sequence and keep the one that reconciles the most reference residues. This is the
rule already frozen in PHASE_B_INCLUSION_RULES_FROZEN.md section 2, applied
without the substring precondition.

Reads counts only. No measured effect is examined.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from run_phaseB_dms_inventory import (census_symbols, flank_map, parse_missense,
                                      readout_kind, score_set_rows, translate)
from run_gate1 import read_topology

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TOPOLOGY = ROOT / "data" / "raw" / "topology_viewer.html"
DEFAULT_CACHE = ROOT / "data" / "external" / "mavedb"
DEFAULT_OUT = ROOT / "results" / "phaseB_recovery_20260927"

MIN_AGREEMENT = 0.95
MIN_MISSENSE = 150
POSITIVE = set("KR")


def charge_class(ref: str, alt: str) -> str | None:
    if ref in POSITIVE and alt not in POSITIVE:
        return "positive_charge_loss"
    if ref not in POSITIVE and alt in POSITIVE:
        return "positive_charge_gain"
    return None


def scan_offsets(variants: list[tuple[str, int, str]], sequence: str) -> tuple[int, float]:
    """Offset maximising reference-residue agreement, over every feasible offset."""
    if not variants:
        return 0, 0.0
    positions = [position for _, position, _ in variants]
    low = 1 - max(positions)
    high = len(sequence) - min(positions)
    best, best_rate = 0, -1.0
    for offset in range(low, high + 1):
        agree = sum(
            1 for ref, position, _ in variants
            if 1 <= position + offset <= len(sequence)
            and sequence[position + offset - 1] == ref
        )
        rate = agree / len(variants)
        if rate > best_rate:
            best, best_rate = offset, rate
        if best_rate == 1.0:
            break
    return best, best_rate


def run(args: argparse.Namespace) -> dict:
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cache = Path(args.cache)
    proteins = read_topology(Path(args.topology))
    symbols = census_symbols(proteins)

    recovered, rejected = [], []
    for search_path in sorted(cache.glob("search_*.json")):
        gene = search_path.stem.replace("search_", "")
        candidates = symbols.get(gene.upper(), [])
        if not candidates:
            continue
        payload = json.loads(search_path.read_text(encoding="utf-8"))
        score_sets = payload if isinstance(payload, list) else payload.get("scoreSets", [])
        for score_set in score_sets:
            urn = score_set.get("urn", "")
            is_expression, _, readout = readout_kind(score_set)
            if not is_expression:
                continue
            rows = score_set_rows(urn, cache)
            if not rows:
                continue
            variants = []
            for row in rows:
                parsed = parse_missense(row.get("hgvs_pro", "") or "")
                if parsed:
                    variants.append(parsed)
            if len(variants) < MIN_MISSENSE:
                continue

            best = None
            for accession in candidates:
                sequence = proteins[accession]["sequence"]
                offset, agreement = scan_offsets(variants, sequence)
                if best is None or agreement > best[2]:
                    best = (accession, offset, agreement)
            accession, offset, agreement = best
            entry = {
                "gene": gene, "urn": urn, "uniprot": accession, "readout": readout,
                "offset": offset, "reference_agreement": round(agreement, 4),
                "n_missense_scored": len(variants),
                "title": score_set.get("title", ""),
            }
            if agreement < MIN_AGREEMENT:
                entry["status"] = "below_agreement_threshold"
                rejected.append(entry)
                continue

            record = proteins[accession]
            flanks = flank_map(record)
            counts = Counter()
            for ref, position, alt in variants:
                census_position = position + offset
                flank = flanks.get(census_position)
                if flank is None:
                    continue
                side, distance = flank
                if distance > 10:
                    continue
                kind = charge_class(ref, alt)
                if kind:
                    counts[f"flank_10_{side}_{kind}"] += 1
            entry.update({
                "num_tmds": len(record["tmds"]),
                "protein_length": record["length"],
                "flank_10_charge_change": sum(counts.values()),
                **{key: counts.get(key, 0) for key in (
                    "flank_10_cytosolic_positive_charge_loss",
                    "flank_10_cytosolic_positive_charge_gain",
                    "flank_10_outer_positive_charge_loss",
                    "flank_10_outer_positive_charge_gain")},
                "status": "recovered",
            })
            recovered.append(entry)

    recovered.sort(key=lambda e: -e["flank_10_charge_change"])
    summary = {
        "analysis": "Recovery of score sets Phase B could not locate by substring match",
        "reads_dms_effects": False,
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "rule_source": "PHASE_B_INCLUSION_RULES_FROZEN.md section 2, frozen 2026-09-27",
        "min_reference_agreement": MIN_AGREEMENT,
        "min_missense": MIN_MISSENSE,
        "n_recovered": len(recovered),
        "n_rejected_below_threshold": len(rejected),
        "recovered": recovered,
        "rejected": rejected,
    }
    (out_dir / "recovery_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    if recovered:
        with (out_dir / "recovered_score_sets.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(recovered[0]), extrasaction="ignore")
            writer.writeheader()
            writer.writerows(recovered)

    print(f"recovered: {len(recovered)}   rejected below {MIN_AGREEMENT}: {len(rejected)}")
    for entry in recovered:
        print(f"  {entry['gene']:<10} {entry['uniprot']:<8} agree={entry['reference_agreement']:.3f} "
              f"offset={entry['offset']:>5} missense={entry['n_missense_scored']:>6} "
              f"flank10_charge={entry['flank_10_charge_change']:>4}  "
              f"cyto {entry['flank_10_cytosolic_positive_charge_loss']}/"
              f"{entry['flank_10_cytosolic_positive_charge_gain']} "
              f"outer {entry['flank_10_outer_positive_charge_loss']}/"
              f"{entry['flank_10_outer_positive_charge_gain']}  [{entry['readout']}]")
    for entry in rejected[:10]:
        print(f"  REJECT {entry['gene']:<10} agree={entry['reference_agreement']:.3f} ({entry['readout']})")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", default=str(DEFAULT_TOPOLOGY))
    parser.add_argument("--cache", default=str(DEFAULT_CACHE))
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    run(parser.parse_args())


if __name__ == "__main__":
    main()

"""Export the Gate 1 flank cohort as a per-event table with no clinical label.

The ClinVar classification is used only as an inclusion filter, exactly as in
run_gate1.py: an event must carry a high-confidence P/LP or B/LB aggregate
classification at two stars or better. The label itself is never written out, so
this table can be used to plan the analysis (cluster sizes, power, covariate
structure) while the outcome stays blind.

Labels are attached later, after the preregistration is frozen, by joining on
variation_id.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from extract_am_flank_scores import build_flank_map
from run_gate1 import AA_CHANGE_RE, charge_change, norm_chrom, read_topology, load_clinvar_candidates, sha256

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TOPOLOGY = ROOT / "data" / "raw" / "topology_viewer.html"
DEFAULT_AM = ROOT / "data" / "raw" / "AlphaMissense_hg38.tsv.gz"
DEFAULT_CLINVAR = ROOT / "data" / "raw" / "clinvar_variant_summary.txt.gz"
DEFAULT_OUT = ROOT / "data" / "processed" / "gate0"

SIDE_NAME = {0: "cytosolic", 1: "outer"}


def run(args: argparse.Namespace) -> dict:
    topology_path, am_path = Path(args.topology), Path(args.alphamissense)
    clinvar_path, out_dir = Path(args.clinvar), Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    proteins = read_topology(topology_path)
    flank_maps = {accession: build_flank_map(record, args.max_distance)
                  for accession, record in proteins.items()}
    clinvar_by_key, qc_counts = load_clinvar_candidates(clinvar_path, proteins)
    qc = defaultdict(int, qc_counts)

    events: dict[tuple[str, str, str], dict] = {}
    with gzip.open(am_path, "rt", encoding="utf-8-sig", newline="") as handle:
        header = next((line for line in handle if line.startswith("#CHROM\t")), None)
        if header is None:
            raise ValueError("AlphaMissense file is missing the #CHROM header")
        reader = csv.DictReader(itertools.chain([header], handle), delimiter="\t")
        for row in reader:
            qc["alphamissense_rows_scanned"] += 1
            accession = (row.get("uniprot_id") or "").strip()
            if accession not in proteins:
                continue
            key = (norm_chrom(row.get("#CHROM", "")), int(row["POS"]),
                   (row.get("REF") or "").strip().upper(), (row.get("ALT") or "").strip().upper())
            clinvar_rows = clinvar_by_key.get(key)
            if not clinvar_rows:
                continue
            protein_variant = (row.get("protein_variant") or "").strip().upper()
            change = AA_CHANGE_RE.fullmatch(protein_variant)
            if not change:
                continue
            ref_aa, position_text, alt_aa = change.groups()
            position = int(position_text)
            record = proteins[accession]
            if position > len(record["sequence"]) or record["sequence"][position - 1] != ref_aa:
                continue
            flank = flank_maps[accession].get(position)
            if flank is None:
                continue
            for clinvar_row in clinvar_rows:
                if clinvar_row["symbols"] and not (clinvar_row["symbols"] & record["_aliases"]):
                    continue
                event_key = (clinvar_row["variation_id"], accession, protein_variant)
                event = events.setdefault(event_key, {
                    "variation_id": clinvar_row["variation_id"],
                    "uniprot": accession,
                    "gene": record.get("gene", ""),
                    "protein_variant": protein_variant,
                    "position": position,
                    "ref_aa": ref_aa,
                    "alt_aa": alt_aa,
                    "side": SIDE_NAME[flank[0]],
                    "distance_to_tmd": flank[1],
                    "tmd_index": flank[2],
                    "n_tmds": len(record["tmds"]),
                    "protein_length": record["length"],
                    "charge_change": charge_change(ref_aa, alt_aa) or "",
                    "am_pathogenicity": row.get("am_pathogenicity", ""),
                    "_labels": set(),
                })
                event["_labels"].add(clinvar_row["label"])

    qc["flank_events_before_label_conflict_filter"] = len(events)
    kept = [event for event in events.values() if len(event["_labels"]) == 1]
    qc["flank_events_excluded_for_conflicting_labels"] = len(events) - len(kept)
    qc["flank_events_exported"] = len(kept)

    fields = ["variation_id", "uniprot", "gene", "protein_variant", "position", "ref_aa",
              "alt_aa", "side", "distance_to_tmd", "tmd_index", "n_tmds", "protein_length",
              "charge_change", "am_pathogenicity"]
    out_path = out_dir / "blinded_flank_events.csv"
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for event in sorted(kept, key=lambda e: (e["uniprot"], e["position"], e["variation_id"])):
            writer.writerow(event)

    summary = {
        "analysis": "Blinded flank-event export (no clinical label column)",
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "max_distance_residues": args.max_distance,
        "label_use": "inclusion filter only; label not written",
        "quality_control": dict(qc),
        "output": {"path": str(out_path), "rows": len(kept)},
        "data_files": {
            "topology_viewer_html": {"sha256": sha256(topology_path)},
            "alphamissense_hg38": {"sha256": sha256(am_path)},
            "clinvar_variant_summary": {"sha256": sha256(clinvar_path)},
        },
    }
    (out_dir / "blinded_export_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k != "data_files"}, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", default=str(DEFAULT_TOPOLOGY))
    parser.add_argument("--alphamissense", default=str(DEFAULT_AM))
    parser.add_argument("--clinvar", default=str(DEFAULT_CLINVAR))
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    parser.add_argument("--max-distance", type=int, default=60)
    run(parser.parse_args())


if __name__ == "__main__":
    main()

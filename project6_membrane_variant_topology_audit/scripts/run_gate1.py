"""Compute outcome-blind feasibility counts for the membrane-topology study.

Clinical labels select eligible variants, but cell-level outputs combine benign
and pathogenic events so the analysis remains blinded before preregistration.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TOPOLOGY = ROOT / "data" / "raw" / "topology_viewer.html"
DEFAULT_AM = ROOT / "data" / "raw" / "AlphaMissense_hg38.tsv.gz"
DEFAULT_CLINVAR = ROOT / "data" / "raw" / "clinvar_variant_summary.txt.gz"
DEFAULT_OUT = ROOT / "results" / "gate1_20260923"

REVIEW_2_STAR_OR_BETTER = {
    "criteria provided, multiple submitters, no conflicts",
    "reviewed by expert panel",
    "practice guideline",
}
LABELS = {
    "pathogenic": "pathogenic",
    "likely pathogenic": "pathogenic",
    "pathogenic/likely pathogenic": "pathogenic",
    "likely pathogenic/pathogenic": "pathogenic",
    "benign": "benign",
    "likely benign": "benign",
    "benign/likely benign": "benign",
    "likely benign/benign": "benign",
}
AA_CHANGE_RE = re.compile(r"^([A-Z])([0-9]+)([A-Z])$")
WINDOWS_DEFAULT = (5, 10, 15)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def norm_chrom(value: str) -> str:
    value = (value or "").strip().upper()
    if value.startswith("CHR"):
        value = value[3:]
    if value in {"M", "MT", "MITO"}:
        return "M"
    return value


def one_base(value: str) -> bool:
    return len((value or "").strip()) == 1 and value.strip().upper() in "ACGT"


def read_topology(path: Path) -> dict[str, dict]:
    text = path.read_text(encoding="utf-8-sig")
    match = re.search(r"const PROTEINS_INLINE = (\[.*?\]);", text, flags=re.S)
    if not match:
        raise ValueError("Could not find the embedded Oxford census in the viewer HTML")
    records = json.loads(match.group(1))
    proteins = {}
    for record in records:
        if record["uniprot"] in proteins:
            raise ValueError(f"Duplicate UniProt accession in census: {record['uniprot']}")
        record["_aliases"] = {
            str(symbol).strip().upper()
            for symbol in [record.get("gene", ""), *record.get("synonyms", [])]
            if str(symbol).strip()
        }
        record["_flanks"] = build_flanks(record)
        proteins[record["uniprot"]] = record
    return proteins


def build_flanks(protein: dict) -> dict[int, dict[int, str]]:
    """Index loop residues by side and distance from the nearest adjacent TMD."""
    windows = {window: {} for window in WINDOWS_DEFAULT}
    tmds = protein["tmds"]
    loops = protein["loops"]
    for loop_id, loop in loops.items():
        index = int(loop_id)
        adjacent_distances = []
        if index > 1:
            previous_tmd = tmds[index - 2]
            adjacent_distances.append(("left", previous_tmd["end"]))
        if index <= len(tmds):
            next_tmd = tmds[index - 1]
            adjacent_distances.append(("right", next_tmd["start"]))
        side = loop.get("location", "").strip().lower()
        if side not in {"inside", "outside"}:
            continue
        side = "cytosolic" if side == "inside" else "outer"
        start, end = int(loop["start"]), int(loop["end"])
        for position in range(start, end + 1):
            distances = []
            for direction, boundary in adjacent_distances:
                distance = position - boundary if direction == "left" else boundary - position
                if distance > 0:
                    distances.append(distance)
            if not distances:
                continue
            distance = min(distances)
            for window in WINDOWS_DEFAULT:
                if distance <= window:
                    windows[window][position] = side
    return windows


def clinical_class(value: str) -> str | None:
    normalized = re.sub(r"\s+", " ", (value or "").strip().lower())
    return LABELS.get(normalized)


def parse_clinvar_key(row: dict[str, str]) -> tuple[str, int, str, str] | None:
    chrom = norm_chrom(row.get("Chromosome", ""))
    position_text = row.get("PositionVCF", "") or row.get("Start", "")
    try:
        position = int(position_text)
    except (TypeError, ValueError):
        return None
    ref = (row.get("ReferenceAlleleVCF", "") or row.get("ReferenceAllele", "")).strip().upper()
    alt = (row.get("AlternateAlleleVCF", "") or row.get("AlternateAllele", "")).strip().upper()
    if not chrom or position < 1 or not one_base(ref) or not one_base(alt):
        return None
    return chrom, position, ref, alt


def symbol_tokens(value: str) -> set[str]:
    return {
        part.strip().upper()
        for part in re.split(r"[;,|/]+", value or "")
        if part.strip()
    }


def load_clinvar_candidates(path: Path, proteins: dict[str, dict]) -> tuple[dict, dict]:
    aliases_to_proteins: dict[str, set[str]] = defaultdict(set)
    for accession, record in proteins.items():
        for alias in record["_aliases"]:
            aliases_to_proteins[alias].add(accession)

    by_key: dict[tuple[str, int, str, str], list[dict]] = defaultdict(list)
    qc = defaultdict(int)
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("ClinVar file is missing a header")
        required = {"Assembly", "Type", "ReviewStatus", "ClinicalSignificance", "GeneSymbol"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"ClinVar header missing required columns: {sorted(missing)}")
        for row in reader:
            qc["clinvar_rows_scanned"] += 1
            if (row.get("Assembly") or "").strip().upper() != "GRCH38":
                continue
            if (row.get("Type") or "").strip().lower() != "single nucleotide variant":
                continue
            review = re.sub(r"\s+", " ", (row.get("ReviewStatus") or "").strip().lower())
            if review not in REVIEW_2_STAR_OR_BETTER:
                continue
            label = clinical_class(row.get("ClinicalSignificance", ""))
            if label is None:
                continue
            symbols = symbol_tokens(row.get("GeneSymbol", ""))
            if symbols and not symbols.intersection(aliases_to_proteins):
                continue
            key = parse_clinvar_key(row)
            if key is None:
                qc["clinvar_confident_missense_without_snv_key"] += 1
                continue
            try:
                variation_id = str(int(row.get("VariationID", "")))
            except (TypeError, ValueError):
                qc["clinvar_missing_variation_id"] += 1
                continue
            # Retain labels for cohort selection and pooled totals, never for cell-level output.
            by_key[key].append({
                "variation_id": variation_id,
                "label": label,
                "symbols": symbols,
                "gene_symbol": (row.get("GeneSymbol") or "").strip(),
                "review_status": review,
            })
    qc["clinvar_candidate_rows"] = sum(len(rows) for rows in by_key.values())
    qc["clinvar_candidate_allele_keys"] = len(by_key)
    return by_key, qc


def charge_change(ref: str, alt: str) -> str | None:
    positive = {"K", "R"}
    if ref in positive and alt not in positive:
        return "positive_charge_loss"
    if ref not in positive and alt in positive:
        return "positive_charge_gain"
    return None


def run(args: argparse.Namespace) -> dict:
    topology_path = Path(args.topology)
    am_path = Path(args.alphamissense)
    clinvar_path = Path(args.clinvar)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    proteins = read_topology(topology_path)
    clinvar_by_key, clinvar_qc = load_clinvar_candidates(clinvar_path, proteins)
    qc = defaultdict(int, clinvar_qc)
    # Preserve conflicting labels until repeated mappings can be excluded as one event.
    mapped: dict[tuple[str, str, str], dict] = {}

    with gzip.open(am_path, "rt", encoding="utf-8-sig", newline="") as handle:
        header = next((line for line in handle if line.startswith("#CHROM\t")), None)
        if header is None:
            raise ValueError("AlphaMissense file is missing the #CHROM header")
        reader = csv.DictReader(itertools.chain([header], handle), delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("AlphaMissense file is missing a header")
        required = {"#CHROM", "POS", "REF", "ALT", "uniprot_id", "protein_variant"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"AlphaMissense header missing required columns: {sorted(missing)}")
        for row in reader:
            qc["alphamissense_rows_scanned"] += 1
            accession = (row.get("uniprot_id") or "").strip()
            if accession not in proteins:
                continue
            key = (
                norm_chrom(row.get("#CHROM", "")),
                int(row["POS"]),
                (row.get("REF") or "").strip().upper(),
                (row.get("ALT") or "").strip().upper(),
            )
            clinvar_rows = clinvar_by_key.get(key)
            if not clinvar_rows:
                continue
            protein_variant = (row.get("protein_variant") or "").strip().upper()
            change = AA_CHANGE_RE.fullmatch(protein_variant)
            if not change:
                qc["alphamissense_non_simple_protein_change"] += 1
                continue
            ref_aa, position_text, alt_aa = change.groups()
            position = int(position_text)
            record = proteins[accession]
            sequence = record["sequence"]
            if position > len(sequence) or sequence[position - 1] != ref_aa:
                qc["alphamissense_reference_sequence_mismatch"] += 1
                continue
            qc["alphamissense_clinvar_allele_mappings"] += 1
            for clinvar_row in clinvar_rows:
                if clinvar_row["symbols"] and not (clinvar_row["symbols"] & record["_aliases"]):
                    qc["gene_symbol_mapping_disagreement"] += 1
                    continue
                event_key = (clinvar_row["variation_id"], accession, protein_variant)
                event = mapped.setdefault(event_key, {
                    "labels": set(),
                    "position": position,
                    "ref": ref_aa,
                    "alt": alt_aa,
                    "uniprot": accession,
                })
                event["labels"].add(clinvar_row["label"])

    counts = {
        window: {(side, kind): 0 for side in ("cytosolic", "outer")
                 for kind in ("positive_charge_loss", "positive_charge_gain")}
        for window in WINDOWS_DEFAULT
    }
    flank_labeled = {window: 0 for window in WINDOWS_DEFAULT}
    flank_pathogenic = {window: 0 for window in WINDOWS_DEFAULT}
    flank_benign = {window: 0 for window in WINDOWS_DEFAULT}
    transition_counts = {
        window: {
            kind: defaultdict(lambda: {"cytosolic": 0, "outer": 0})
            for kind in ("positive_charge_loss", "positive_charge_gain")
        }
        for window in WINDOWS_DEFAULT
    }
    qc["mapped_protein_events_before_label_conflict_filter"] = len(mapped)
    qc["label_conflicted_protein_events_excluded"] = 0
    qc["mapped_protein_events_with_high_confidence_label"] = 0
    for event in mapped.values():
        if len(event["labels"]) != 1:
            qc["label_conflicted_protein_events_excluded"] += 1
            continue
        label = next(iter(event["labels"]))
        qc["mapped_protein_events_with_high_confidence_label"] += 1
        accession = event["uniprot"]
        position = event["position"]
        kind = charge_change(event["ref"], event["alt"])
        for window in WINDOWS_DEFAULT:
            side = proteins[accession]["_flanks"][window].get(position)
            if side is None:
                continue
            flank_labeled[window] += 1
            if label == "pathogenic":
                flank_pathogenic[window] += 1
            else:
                flank_benign[window] += 1
            if kind is not None:
                counts[window][(side, kind)] += 1
                transition_counts[window][kind][f"{event['ref']}>{event['alt']}"][side] += 1

    matched_substitution_support = {}
    for window in WINDOWS_DEFAULT:
        matched_substitution_support[str(window)] = {}
        for kind in ("positive_charge_loss", "positive_charge_gain"):
            by_transition = transition_counts[window][kind]
            shared = [values for values in by_transition.values()
                      if values["cytosolic"] > 0 and values["outer"] > 0]
            matched_substitution_support[str(window)][kind] = {
                "distinct_transition_types_on_cytosolic_side": sum(
                    values["cytosolic"] > 0 for values in by_transition.values()
                ),
                "distinct_transition_types_on_outer_side": sum(
                    values["outer"] > 0 for values in by_transition.values()
                ),
                "transition_types_present_on_both_sides": len(shared),
                "matched_support_n_sum_of_sidewise_minima": sum(
                    min(values["cytosolic"], values["outer"]) for values in shared
                ),
            }

    counts_path = out_dir / "counts_by_window.csv"
    with counts_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["window_residues", "side", "charge_change", "n_labeled_variants"])
        for window in WINDOWS_DEFAULT:
            for side in ("cytosolic", "outer"):
                for kind in ("positive_charge_loss", "positive_charge_gain"):
                    writer.writerow([window, side, kind, counts[window][(side, kind)]])

    primary_window = 10
    primary_cells = [counts[primary_window][(side, kind)]
                     for side in ("cytosolic", "outer")
                     for kind in ("positive_charge_loss", "positive_charge_gain")]
    summary = {
        "analysis": "Outcome-blind Gate 1 sample-size feasibility",
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "primary_window_residues": primary_window,
        "unit": "unique ClinVar VariationID × canonical UniProt accession × amino-acid substitution",
        "clinvar_labels_in_cells": "Pathogenic/Likely pathogenic and Benign/Likely benign combined; no per-cell label split",
        "clinvar_overall_pathogenic_count_scope": "All confidently labeled ClinVar missense protein events within the selected TMD-flank window, without side/change stratification",
        "primary_window_total_labeled_flank_events": flank_labeled[primary_window],
        "primary_window_total_pathogenic_flank_events": flank_pathogenic[primary_window],
        "primary_window_total_benign_flank_events": flank_benign[primary_window],
        "gate_thresholds": {
            "at_least_60_labeled_events_in_each_of_four_cells": 60,
            "at_least_150_pathogenic_flank_events_overall": 150,
        },
        "gate_pass": min(primary_cells) >= 60 and flank_pathogenic[primary_window] >= 150,
        "window_sensitivity": {
            str(window): {
                "total_labeled_flank_events": flank_labeled[window],
                "total_pathogenic_flank_events": flank_pathogenic[window],
                "total_benign_flank_events": flank_benign[window],
            }
            for window in WINDOWS_DEFAULT
        },
        "outcome_blind_matched_substitution_support": matched_substitution_support,
        "cell_counts_file": str(counts_path),
        "quality_control": dict(qc),
        "topology_proteins": len(proteins),
        "data_files": {
            "topology_viewer_html": {
                "path": str(topology_path), "bytes": topology_path.stat().st_size,
                "sha256": sha256(topology_path),
            },
            "alphamissense_hg38": {
                "path": str(am_path), "bytes": am_path.stat().st_size,
                "sha256": sha256(am_path),
            },
            "clinvar_variant_summary": {
                "path": str(clinvar_path), "bytes": clinvar_path.stat().st_size,
                "sha256": sha256(clinvar_path),
            },
        },
    }
    summary_path = out_dir / "gate_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "topology_proteins": len(proteins),
        "primary_window_cell_counts": {
            f"{side}_{kind}": counts[primary_window][(side, kind)]
            for side in ("cytosolic", "outer")
            for kind in ("positive_charge_loss", "positive_charge_gain")
        },
        "primary_window_total_pathogenic_flank_events": flank_pathogenic[primary_window],
        "primary_window_matched_substitution_support": matched_substitution_support[str(primary_window)],
        "gate_pass": summary["gate_pass"],
        "summary": str(summary_path),
    }, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", default=str(DEFAULT_TOPOLOGY))
    parser.add_argument("--alphamissense", default=str(DEFAULT_AM))
    parser.add_argument("--clinvar", default=str(DEFAULT_CLINVAR))
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()

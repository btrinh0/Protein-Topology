"""Phase A: does AlphaMissense already encode the positive-inside rule?

Label-free audit. Reads only the Oxford-MRC census topology and AlphaMissense
pathogenicity scores. No ClinVar classification is opened, so this cannot
unblind the Gate 1 side-by-charge-change cells.

The audit asks how AlphaMissense's score varies with (a) which side of the
membrane an aqueous residue sits on, (b) the signed change in formal charge the
substitution makes, and (c) how far the residue is from the nearest TMD edge.
If the predictor already encodes the rule, cytosolic positive-charge losses
should score higher than the same losses on the outer side, and the gap should
decay with distance from the membrane.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TOPOLOGY = ROOT / "data" / "raw" / "topology_viewer.html"
DEFAULT_AM = ROOT / "data" / "raw" / "AlphaMissense_hg38.tsv.gz"
DEFAULT_OUT = ROOT / "results" / "phaseA_20260924"

AA_CHANGE_RE = re.compile(r"^([A-Z])([0-9]+)([A-Z])$")
AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {aa: i for i, aa in enumerate(AA_ORDER)}

POSITIVE = {"K", "R"}
NEGATIVE = {"D", "E"}

# Upper edge of each distance band, in residues from the nearest TMD edge.
# The far bands act as an internal control: the rule should fade out there.
DISTANCE_BANDS = ((5, "1-5"), (10, "6-10"), (15, "11-15"),
                  (25, "16-25"), (50, "26-50"), (math.inf, "51+"))
PRIMARY_BAND_MAX = 10


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def norm_chrom(value: str) -> str:
    value = (value or "").strip().upper()
    return value[3:] if value.startswith("CHR") else value


def formal_charge(aa: str) -> int:
    if aa in POSITIVE:
        return 1
    if aa in NEGATIVE:
        return -1
    return 0


def band_for(distance: int) -> str:
    for upper, label in DISTANCE_BANDS:
        if distance <= upper:
            return label
    raise AssertionError("distance bands must cover every distance")


def gate1_class(ref: str, alt: str) -> str:
    """The coarse loss/gain categories used by the Gate 1 feasibility counts."""
    if ref in POSITIVE and alt not in POSITIVE:
        return "positive_charge_loss"
    if ref not in POSITIVE and alt in POSITIVE:
        return "positive_charge_gain"
    return "other"


def read_topology(path: Path) -> dict[str, dict]:
    text = path.read_text(encoding="utf-8-sig")
    match = re.search(r"const PROTEINS_INLINE = (\[.*?\]);", text, flags=re.S)
    if not match:
        raise ValueError("Could not find the embedded Oxford census in the viewer HTML")
    proteins = {}
    for record in json.loads(match.group(1)):
        record["_loop_residues"] = loop_residues(record)
        record["_seen"] = bytearray((len(record["sequence"]) * len(AA_ORDER) + 7) // 8)
        proteins[record["uniprot"]] = record
    return proteins


def loop_residues(protein: dict) -> dict[int, tuple[str, int]]:
    """Map every annotated aqueous-loop residue to (side, distance to nearest TMD edge).

    Unlike the Gate 1 flank builder this keeps residues at every distance, which
    is what lets the audit show how far from the membrane the rule still reaches.
    """
    residues: dict[int, tuple[str, int]] = {}
    tmds = protein["tmds"]
    for loop_id, loop in protein["loops"].items():
        index = int(loop_id)
        boundaries = []
        if index > 1:
            boundaries.append(("left", tmds[index - 2]["end"]))
        if index <= len(tmds):
            boundaries.append(("right", tmds[index - 1]["start"]))
        side = loop.get("location", "").strip().lower()
        if side not in {"inside", "outside"}:
            continue
        side = "cytosolic" if side == "inside" else "outer"
        for position in range(int(loop["start"]), int(loop["end"]) + 1):
            distances = [
                position - boundary if direction == "left" else boundary - position
                for direction, boundary in boundaries
            ]
            distances = [d for d in distances if d > 0]
            if distances:
                residues[position] = (side, min(distances))
    return residues


def mark_seen(protein: dict, position: int, alt: str) -> bool:
    """Record one (position, alt) substitution; return False if already counted.

    AlphaMissense repeats a protein substitution once per nucleotide change that
    produces it, so codon-degenerate substitutions would otherwise be upweighted.
    """
    bit = (position - 1) * len(AA_ORDER) + AA_INDEX[alt]
    byte, offset = divmod(bit, 8)
    mask = 1 << offset
    if protein["_seen"][byte] & mask:
        return False
    protein["_seen"][byte] |= mask
    return True


class Running:
    """Streaming count/mean/variance so no per-variant scores are held in memory."""

    __slots__ = ("n", "total", "total_sq")

    def __init__(self) -> None:
        self.n = 0
        self.total = 0.0
        self.total_sq = 0.0

    def add(self, value: float) -> None:
        self.n += 1
        self.total += value
        self.total_sq += value * value

    @property
    def mean(self) -> float | None:
        return self.total / self.n if self.n else None

    @property
    def sd(self) -> float | None:
        if self.n < 2:
            return None
        variance = (self.total_sq - self.total * self.total / self.n) / (self.n - 1)
        return math.sqrt(max(variance, 0.0))


def run(args: argparse.Namespace) -> dict:
    topology_path = Path(args.topology)
    am_path = Path(args.alphamissense)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    proteins = read_topology(topology_path)
    qc: dict[str, int] = defaultdict(int)

    by_delta: dict[tuple[str, str, int], Running] = defaultdict(Running)
    by_class: dict[tuple[str, str, str], Running] = defaultdict(Running)
    by_transition: dict[tuple[str, str], Running] = defaultdict(Running)

    with gzip.open(am_path, "rt", encoding="utf-8-sig", newline="") as handle:
        header = next((line for line in handle if line.startswith("#CHROM\t")), None)
        if header is None:
            raise ValueError("AlphaMissense file is missing the #CHROM header")
        reader = csv.DictReader(itertools.chain([header], handle), delimiter="\t")
        required = {"uniprot_id", "protein_variant", "am_pathogenicity"}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"AlphaMissense header missing required columns: {sorted(missing)}")
        for row in reader:
            qc["alphamissense_rows_scanned"] += 1
            protein = proteins.get((row.get("uniprot_id") or "").strip())
            if protein is None:
                continue
            qc["rows_in_census_proteins"] += 1
            change = AA_CHANGE_RE.fullmatch((row.get("protein_variant") or "").strip().upper())
            if not change:
                qc["non_simple_protein_change"] += 1
                continue
            ref_aa, position_text, alt_aa = change.groups()
            position = int(position_text)
            if ref_aa not in AA_INDEX or alt_aa not in AA_INDEX:
                qc["unknown_amino_acid_code"] += 1
                continue
            sequence = protein["sequence"]
            if position > len(sequence) or sequence[position - 1] != ref_aa:
                qc["reference_sequence_mismatch"] += 1
                continue
            placement = protein["_loop_residues"].get(position)
            if placement is None:
                qc["rows_outside_annotated_loops"] += 1
                continue
            if not mark_seen(protein, position, alt_aa):
                qc["duplicate_codon_paths_collapsed"] += 1
                continue
            try:
                score = float(row["am_pathogenicity"])
            except (TypeError, ValueError, KeyError):
                qc["unparsable_pathogenicity_score"] += 1
                continue

            side, distance = placement
            band = band_for(distance)
            delta = formal_charge(alt_aa) - formal_charge(ref_aa)
            qc["substitutions_counted"] += 1

            by_delta[(side, band, delta)].add(score)
            by_class[(side, band, gate1_class(ref_aa, alt_aa))].add(score)
            if distance <= PRIMARY_BAND_MAX:
                by_transition[(side, f"{ref_aa}>{alt_aa}")].add(score)

    band_labels = [label for _, label in DISTANCE_BANDS]
    band_rank = {label: i for i, label in enumerate(band_labels)}

    delta_path = out_dir / "am_score_by_side_charge_delta.csv"
    with delta_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["side", "distance_band", "charge_delta", "n_substitutions",
                         "mean_am_pathogenicity", "sd_am_pathogenicity"])
        for (side, band, delta), stats in sorted(
            by_delta.items(), key=lambda kv: (kv[0][0], band_rank[kv[0][1]], kv[0][2])
        ):
            writer.writerow([side, band, delta, stats.n,
                             _round(stats.mean), _round(stats.sd)])

    class_path = out_dir / "am_score_by_side_gate1_class.csv"
    with class_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["side", "distance_band", "charge_change", "n_substitutions",
                         "mean_am_pathogenicity", "sd_am_pathogenicity"])
        for (side, band, kind), stats in sorted(
            by_class.items(), key=lambda kv: (kv[0][0], band_rank[kv[0][1]], kv[0][2])
        ):
            writer.writerow([side, band, kind, stats.n,
                             _round(stats.mean), _round(stats.sd)])

    # Exact-transition contrast inside the primary window: same amino-acid change,
    # compared across the two sides, which removes the substitution-identity confound.
    transition_path = out_dir / "am_score_by_transition_primary_window.csv"
    transitions = sorted({transition for _, transition in by_transition})
    with transition_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["transition", "charge_delta", "n_cytosolic", "mean_cytosolic",
                         "n_outer", "mean_outer", "cytosolic_minus_outer"])
        for transition in transitions:
            ref_aa, alt_aa = transition.split(">")
            cyto = by_transition.get(("cytosolic", transition), Running())
            outer = by_transition.get(("outer", transition), Running())
            gap = (cyto.mean - outer.mean) if (cyto.mean is not None and outer.mean is not None) else None
            writer.writerow([transition, formal_charge(alt_aa) - formal_charge(ref_aa),
                             cyto.n, _round(cyto.mean), outer.n, _round(outer.mean), _round(gap)])

    summary = {
        "analysis": "Phase A label-free AlphaMissense topology-rule audit",
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "reads_clinvar": False,
        "unit": "unique canonical UniProt accession x protein substitution",
        "score_column": "am_pathogenicity",
        "primary_window_residues": PRIMARY_BAND_MAX,
        "distance_bands": band_labels,
        "charge_convention": "K/R = +1, D/E = -1, all other residues = 0; "
                             "charge_delta = charge(alt) - charge(ref)",
        "topology_proteins": len(proteins),
        "headline_contrasts": headline_contrasts(by_class, by_delta),
        "quality_control": dict(qc),
        "output_files": {
            "by_charge_delta": str(delta_path),
            "by_gate1_class": str(class_path),
            "by_transition_primary_window": str(transition_path),
        },
        "data_files": {
            "topology_viewer_html": {
                "path": str(topology_path), "bytes": topology_path.stat().st_size,
                "sha256": sha256(topology_path),
            },
            "alphamissense_hg38": {
                "path": str(am_path), "bytes": am_path.stat().st_size,
                "sha256": sha256(am_path),
            },
        },
    }
    summary_path = out_dir / "phaseA_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "substitutions_counted": qc["substitutions_counted"],
        "headline_contrasts": summary["headline_contrasts"],
        "summary": str(summary_path),
    }, indent=2))
    return summary


def _round(value: float | None, digits: int = 4) -> float | str:
    return "" if value is None else round(value, digits)


def headline_contrasts(by_class: dict, by_delta: dict) -> dict:
    """Cytosolic-minus-outer score gaps, per distance band, for the key categories."""
    contrasts: dict[str, dict] = {}
    for kind in ("positive_charge_loss", "positive_charge_gain"):
        contrasts[kind] = {}
        for _, band in DISTANCE_BANDS:
            cyto = by_class.get(("cytosolic", band, kind))
            outer = by_class.get(("outer", band, kind))
            if cyto and outer and cyto.mean is not None and outer.mean is not None:
                contrasts[kind][band] = {
                    "mean_cytosolic": round(cyto.mean, 4),
                    "mean_outer": round(outer.mean, 4),
                    "cytosolic_minus_outer": round(cyto.mean - outer.mean, 4),
                    "n_cytosolic": cyto.n,
                    "n_outer": outer.n,
                }
    return contrasts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", default=str(DEFAULT_TOPOLOGY))
    parser.add_argument("--alphamissense", default=str(DEFAULT_AM))
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    run(parser.parse_args())


if __name__ == "__main__":
    main()

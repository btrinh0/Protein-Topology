"""Shared loading, alignment-based mapping and flank geometry for Route A.

Route A asks whether the one-parameter correction helps on ProteinGym deep mutational scans.
Every ProteinGym assay has to be placed on the topology census before any score is looked up.
This module does that placement and nothing else: it never reads a DMS_score value.

Mapping is an alignment, not a global offset. The SERT cassette (scripts/load_sert.py) is the
reason: a best-fit offset of -24 mismaps 215 residues, while the two-segment alignment is exact.

Geometry (side and distance to the nearest TMD edge) comes from read_topology() in run_gate1.py
followed by build_flank_map() in extract_am_flank_scores.py, which is what produced
data/processed/gate0/am_full_flank_scores.npz. Using the same path keeps Route A on the geometry
that every earlier analysis used.
"""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from extract_am_flank_scores import build_flank_map  # noqa: E402
from run_gate1 import read_topology  # noqa: E402

REFERENCE = ROOT / "data/raw/proteingym/reference_github_main/DMS_substitutions.csv"
ASSAY_DIR = ROOT / "data/raw/proteingym/v1.3_extracted/DMS_ProteinGym_substitutions"
UNIPROT_MAP = ROOT / "data/raw/proteingym/uniprot_entry_map_v1.3.json"
CENSUS = ROOT / "data/raw/topology_viewer.html"
AM_ARCHIVE = ROOT / "data/processed/gate0/am_full_flank_scores.npz"
ESM_ARCHIVE = ROOT / "data/processed/gate0/esm1b_flank_scores.npz"

WINDOW = 10                 # flank sites within 10 residues of a TMD edge
MAX_DISTANCE = 60           # geometry depth used by the existing archives
NEUTRAL = "ACFGILMNPQSTVWY"  # native residues eligible for the K/R-gain contrast (validated panel)
KR = "KR"
MIN_AGREEMENT = 0.99        # reject an assay whose reference residues disagree with the census
CONTAMINATION_RADIUS = 10   # flank site this close to a construct insertion is flagged

MUTANT_TOKEN = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def load_reference(path: Path = REFERENCE) -> list[dict]:
    """ProteinGym reference table, one dict per assay. Metadata only."""
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def load_uniprot_map(path: Path = UNIPROT_MAP) -> dict:
    """Entry name -> {accession, gene, length, reviewed, organism}, from route_a_resolve_uniprot.py."""
    return json.loads(path.read_text(encoding="utf-8"))["entries"]


def load_census(path: Path = CENSUS) -> dict[str, dict]:
    """accession -> census record with `sequence`, `tmds`, `loops`, plus `_flanks` from run_gate1."""
    return read_topology(path)


def census_flanks(record: dict) -> dict[int, tuple[int, int, int]]:
    """census position -> (side index, distance to nearest TMD edge, TMD index), depth MAX_DISTANCE."""
    return build_flank_map(record, MAX_DISTANCE)


def parse_mutant(mutant: str) -> list[tuple[str, int, str]]:
    """'A1P:D2N' -> [('A', 1, 'P'), ('D', 2, 'N')]. Raises on anything that is not a substitution."""
    tokens = []
    for token in mutant.strip().split(":"):
        match = MUTANT_TOKEN.match(token)
        if not match:
            raise ValueError(f"Not a substitution token: {token!r} in {mutant!r}")
        wt, position, alt = match.groups()
        tokens.append((wt, int(position), alt))
    return tokens


def _aligner():
    from Bio import Align
    from Bio.Align import substitution_matrices

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    # Free end gaps: a construct that is a fragment of the census protein is not penalised at its ends.
    aligner.end_insertion_score = 0.0
    aligner.end_deletion_score = 0.0
    return aligner


_ALIGNER = None


def align_construct(census_seq: str, construct_seq: str) -> np.ndarray:
    """Map each construct position (1-based index into the returned array) to a census position.

    Returns an int array of length len(construct_seq)+1; entry 0 is unused, and 0 marks a construct
    position with no census counterpart (a tag, linker or other inserted residue).
    """
    global _ALIGNER
    if _ALIGNER is None:
        _ALIGNER = _aligner()
    alignment = _ALIGNER.align(census_seq, construct_seq)[0]
    census_blocks, construct_blocks = alignment.aligned[0], alignment.aligned[1]
    mapping = np.zeros(len(construct_seq) + 1, dtype=np.int64)
    for (c_start, c_end), (q_start, q_end) in zip(census_blocks, construct_blocks):
        length = c_end - c_start
        mapping[q_start + 1:q_end + 1] = np.arange(c_start + 1, c_start + 1 + length)
    return mapping


def segments_of(mapping: np.ndarray) -> list[dict]:
    """Runs of construct positions that map to consecutive census positions, with their offset."""
    segments = []
    positions = np.flatnonzero(mapping)
    if not len(positions):
        return segments
    start = previous = positions[0]
    for position in positions[1:]:
        consecutive = position == previous + 1 and mapping[position] == mapping[previous] + 1
        if not consecutive:
            segments.append((start, previous))
            start = position
        previous = position
    segments.append((start, previous))
    return [{"construct_start": int(a), "construct_end": int(b),
             "census_start": int(mapping[a]), "census_end": int(mapping[b]),
             "offset_census_minus_construct": int(mapping[a] - a)} for a, b in segments]


def insertions_of(mapping: np.ndarray, census_length: int) -> list[dict]:
    """Construct runs with no census counterpart, located by the census residues either side of them."""
    mapped = np.flatnonzero(mapping)
    if not len(mapped):
        return []
    out = []
    last = len(mapping) - 1
    # Leading and trailing runs are construct residues outside the census protein.
    if mapped[0] > 1:
        out.append({"construct_start": 1, "construct_end": int(mapped[0] - 1),
                    "census_before": 0, "census_after": int(mapping[mapped[0]])})
    for a, b in zip(mapped[:-1], mapped[1:]):
        if b - a > 1:
            out.append({"construct_start": int(a + 1), "construct_end": int(b - 1),
                        "census_before": int(mapping[a]), "census_after": int(mapping[b])})
    if mapped[-1] < last:
        out.append({"construct_start": int(mapped[-1] + 1), "construct_end": int(last),
                    "census_before": int(mapping[mapped[-1]]), "census_after": int(census_length + 1)})
    return out


def insertion_flanks(insertions: list[dict], flank_positions) -> list[int]:
    """Measured flank census positions within CONTAMINATION_RADIUS of any insertion point."""
    hits = set()
    for flank in flank_positions:
        for item in insertions:
            near = min(abs(flank - item["census_before"]), abs(flank - item["census_after"]))
            if near <= CONTAMINATION_RADIUS:
                hits.add(int(flank))
    return sorted(hits)


def map_assay(reference_row: dict, census_record: dict, mutants: list[list[tuple]]) -> dict:
    """Mapping manifest for one assay. `mutants` is the parsed mutant column (no scores involved)."""
    construct = reference_row["target_seq"]
    census_seq = census_record["sequence"]
    mapping = align_construct(census_seq, construct)
    aligned = mapping[1:] > 0
    aligned_pairs = int(aligned.sum())
    identity = float(np.mean([census_seq[mapping[i] - 1] == construct[i - 1]
                              for i in np.flatnonzero(mapping)])) if aligned_pairs else 0.0

    # Reference residues: the wild-type letter of every substitution, against the census residue it maps to.
    checked = agree = mapped_variants = dropped_variants = 0
    mismatch_to_target = 0
    for variant in mutants:
        variant_mapped = True
        for wt, position, _alt in variant:
            if position < 1 or position > len(construct):
                raise ValueError(f"{reference_row['DMS_id']}: position {position} outside target sequence")
            if construct[position - 1] != wt:
                mismatch_to_target += 1
            census_pos = mapping[position]
            if census_pos == 0:
                variant_mapped = False
                continue
            checked += 1
            agree += int(census_seq[census_pos - 1] == wt)
        if variant_mapped:
            mapped_variants += 1
        else:
            dropped_variants += 1

    agreement = agree / checked if checked else 0.0
    insertions = insertions_of(mapping, len(census_seq))
    dropped_positions = int(len(construct) - aligned_pairs)
    return {
        "dms_id": reference_row["DMS_id"],
        "uniprot_entry": reference_row["UniProt_ID"],
        "target_length": len(construct),
        "census_length": len(census_seq),
        "aligned_positions": aligned_pairs,
        "identity_over_aligned": identity,
        "reference_residues_checked": checked,
        "reference_residues_agree": agree,
        "reference_agreement": agreement,
        "wt_letter_vs_target_mismatches": mismatch_to_target,
        "positions_dropped_as_tag_or_construct": dropped_positions,
        "insertions": insertions,
        "segments": segments_of(mapping),
        "variants_mapped": mapped_variants,
        "variants_dropped": dropped_variants,
        "accepted": agreement >= MIN_AGREEMENT,
        "_mapping": mapping,
    }

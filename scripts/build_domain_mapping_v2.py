from __future__ import annotations

"""Build the domain-valid mapping table required by the Topological Constraint Test.

SCOPe provides the source-domain definitions; CATH and ECOD annotations are
evaluated as candidates and selected based on reciprocal residue overlap rather
than first-matching chain heuristics. Every candidate is written to a companion
table so the mapping can be audited and reproduced.
"""

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


SCOPE_SCCS_RE = re.compile(r"^([a-z])(?:\.(\d+))?(?:\.(\d+))?(?:\.(\d+))?$")
RESIDUE_RANGE_RE = re.compile(r"^(-?\d+)([A-Za-z]?)-(-?\d+)([A-Za-z]?)$")


@dataclass(frozen=True)
class Segment:
    chain_id: str
    start_resseq: int | None
    start_icode: str
    end_resseq: int | None
    end_icode: str
    is_full_chain: bool = False

    def as_dict(self) -> dict[str, object]:
        return {
            "chain_id": self.chain_id,
            "start_resseq": self.start_resseq,
            "start_icode": self.start_icode,
            "end_resseq": self.end_resseq,
            "end_icode": self.end_icode,
            "is_full_chain": self.is_full_chain,
        }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build V2 SCOPe-to-CATH/ECOD domain mappings using residue overlap. "
            "Existing outputs are not overwritten by default."
        )
    )
    parser.add_argument(
        "--classifications-dir",
        type=Path,
        default=Path("data/raw/classifications"),
        help="Directory containing the frozen classification source files.",
    )
    parser.add_argument(
        "--structures-dir",
        type=Path,
        default=Path("data/raw/structures"),
        help="Restrict the mapping to locally available experimental structures.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/processed/v2/domain_mapping_v2.csv"),
    )
    parser.add_argument(
        "--candidates-output",
        type=Path,
        default=Path("data/processed/v2/domain_mapping_candidates_v2.csv"),
    )
    parser.add_argument(
        "--diagnostics-output",
        type=Path,
        default=Path("data/processed/v2/domain_mapping_diagnostics_v2.json"),
    )
    parser.add_argument(
        "--cath-boundaries",
        type=Path,
        default=Path("data/raw/classifications/cath_domain_boundaries_v4_3_0.txt"),
        help="CATH Domall boundary file; required for CATH selection.",
    )
    parser.add_argument(
        "--min-reciprocal-overlap",
        type=float,
        default=0.80,
        help="Minimum coverage required for both source and candidate domains.",
    )
    parser.add_argument(
        "--ambiguity-margin",
        type=float,
        default=0.02,
        help="Do not select a candidate if its score is this close to the runner-up.",
    )
    parser.add_argument(
        "--include-non-astral40",
        action="store_true",
        help="Keep SCOPe domains absent from the ASTRAL 40% representative set.",
    )
    return parser.parse_args()


def normalise_chain(chain_id: str) -> str:
    chain_id = chain_id.strip()
    return "" if chain_id in {"", "-", "0", ".", "?"} else chain_id


def parse_segment_token(token: str) -> Segment:
    """Parse SCOPe/ECOD text such as ``A:1-40`` or ``A:1B-99A``."""

    token = token.strip()
    if ":" not in token:
        raise ValueError(f"Missing chain separator in segment {token!r}")
    chain_id, residue_range = token.split(":", 1)
    chain_id = normalise_chain(chain_id)
    residue_range = residue_range.strip()
    if not residue_range:
        return Segment(chain_id, None, "", None, "", is_full_chain=True)
    match = RESIDUE_RANGE_RE.fullmatch(residue_range)
    if not match:
        raise ValueError(f"Unsupported residue range {residue_range!r}")
    return Segment(
        chain_id=chain_id,
        start_resseq=int(match.group(1)),
        start_icode=match.group(2),
        end_resseq=int(match.group(3)),
        end_icode=match.group(4),
    )


def parse_segment_string(value: str) -> tuple[list[Segment], list[str]]:
    """Return parsed segments and any source tokens that need manual review."""

    if not value or value.strip() == "-":
        return [], [value or "<empty>"]
    segments: list[Segment] = []
    problems: list[str] = []
    for token in value.split(","):
        try:
            segments.append(parse_segment_token(token))
        except ValueError:
            problems.append(token.strip())
    return segments, problems


def parse_scope_hierarchy(sccs: str) -> dict[str, str]:
    match = SCOPE_SCCS_RE.fullmatch(sccs.strip())
    if not match:
        return {
            "scope_sccs": sccs,
            "scope_class_id": "",
            "scope_fold_id": "",
            "scope_superfamily_id": "",
            "scope_family_id": "",
        }
    values = [part for part in match.groups() if part is not None]
    return {
        "scope_sccs": sccs,
        "scope_class_id": values[0],
        "scope_fold_id": ".".join(values[:2]) if len(values) >= 2 else "",
        "scope_superfamily_id": ".".join(values[:3]) if len(values) >= 3 else "",
        "scope_family_id": ".".join(values[:4]) if len(values) >= 4 else "",
    }


def parse_scope_cla(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            source_domain_id, structure_id, source_range, sccs, sunid, hierarchy = parts[:6]
            segments, segment_problems = parse_segment_string(source_range)
            parsed_hierarchy = parse_scope_hierarchy(sccs)
            chain_ids = sorted({segment.chain_id for segment in segments})
            rows.append(
                {
                    "source_domain_id": source_domain_id,
                    "structure_id": structure_id.lower(),
                    "source_range": source_range,
                    "scope_sunid": sunid,
                    "scope_hierarchy_path": hierarchy,
                    "segments": segments,
                    "segment_problems": segment_problems,
                    "chain_ids": chain_ids,
                    **parsed_hierarchy,
                }
            )
    return rows


def parse_astral40_ids(path: Path) -> set[str]:
    domain_ids: set[str] = set()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                domain_ids.add(line[1:].split()[0].strip())
    return domain_ids


def parse_scope_species(path: Path) -> dict[str, str]:
    species: dict[str, str] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 5 and parts[1] == "sp":
                species[parts[0]] = parts[4]
    return species


def extract_species_description(scope_hierarchy_path: str, species_by_sunid: dict[str, str]) -> str:
    match = re.search(r"(?:^|,)sp=(\d+)(?:,|$)", scope_hierarchy_path)
    return species_by_sunid.get(match.group(1), "") if match else ""


def infer_kingdom(species_description: str) -> str:
    value = species_description.lower()
    if not value:
        return ""
    archaea = ("archaea", "archaeon", "methano", "halobacter", "sulfolobus", "pyrococcus", "thermococcus")
    bacteria = (
        "bacter", "escherichia", "bacillus", "staphyl", "strept", "salmonella",
        "mycobacter", "pseudomonas", "helicobacter", "clostridium", "vibrio",
        "neisseria", "campylobacter", "thermus", "thermotoga", "aquifex",
    )
    eukaryotes = (
        "human", "mouse", "rat", "homo sapiens", "mus musculus", "saccharomyces",
        "drosophila", "caenorhabditis", "arabidopsis", "xenopus", "danio", "gallus",
        "eukaryot", "fung", "plant", "mammal", "ciliate", "alga",
    )
    if any(term in value for term in archaea):
        return "archaea"
    if any(term in value for term in bacteria):
        return "bacteria"
    if any(term in value for term in eukaryotes):
        return "eukaryotes"
    return ""


def parse_cath_metadata(path: Path) -> dict[str, dict[str, str]]:
    """Parse the cached CATH REST JSON used by V0, without using it for mapping."""

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        first = handle.read(1)
        handle.seek(0)
        if first == "{":
            raw = json.load(handle)
            entries = raw.get("data", [])
            metadata: dict[str, dict[str, str]] = {}
            for entry in entries:
                domain_id = str(entry.get("domain_id", ""))
                if not domain_id:
                    continue
                superfamily_id = str(entry.get("superfamily_id", ""))
                fields = parse_cath_hierarchy(superfamily_id)
                metadata[domain_id] = {
                    "cath_domain_id": domain_id,
                    "cath_superfamily_id": superfamily_id,
                    **fields,
                }
            return metadata
    raise ValueError(
        f"{path} is not the expected cached CATH REST JSON. "
        "Supply a compatible metadata file before running V2 mapping."
    )


def parse_cath_hierarchy(superfamily_id: str) -> dict[str, str]:
    levels = superfamily_id.split(".") if superfamily_id else []
    return {
        "cath_class_id": levels[0] if len(levels) >= 1 else "",
        "cath_architecture_id": ".".join(levels[:2]) if len(levels) >= 2 else "",
        "cath_topology_id": ".".join(levels[:3]) if len(levels) >= 3 else "",
        "cath_homology_id": ".".join(levels[:4]) if len(levels) >= 4 else "",
    }


def parse_cath_boundaries(path: Path) -> dict[tuple[str, str], list[dict[str, Any]]]:
    """Parse CATH Domall (CDF) v2.0 into domain segments indexed by PDB/chain."""

    by_pdb_chain: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#") or not line.strip():
                continue
            tokens = line.split()
            if len(tokens) < 4 or not tokens[1].startswith("D"):
                continue
            chain_name = tokens[0]
            if len(chain_name) != 5:
                continue
            try:
                domain_count = int(tokens[1][1:])
                fragment_count = int(tokens[2][1:])
            except ValueError:
                continue
            position = 3
            domains: list[list[Segment]] = []
            try:
                for _ in range(domain_count):
                    segment_count = int(tokens[position])
                    position += 1
                    segments: list[Segment] = []
                    for _ in range(segment_count):
                        start_chain = normalise_chain(tokens[position])
                        start_resseq = int(tokens[position + 1])
                        start_icode = "" if tokens[position + 2] == "-" else tokens[position + 2]
                        end_chain = normalise_chain(tokens[position + 3])
                        end_resseq = int(tokens[position + 4])
                        end_icode = "" if tokens[position + 5] == "-" else tokens[position + 5]
                        if start_chain != end_chain:
                            raise ValueError("CATH segment spans two chain IDs")
                        segments.append(
                            Segment(start_chain, start_resseq, start_icode, end_resseq, end_icode)
                        )
                        position += 6
                    domains.append(segments)
                # Skip the trailing fragment descriptors, which are not part of a domain.
                for _ in range(fragment_count):
                    position += 7
            except (IndexError, ValueError) as exc:
                raise ValueError(f"CATH boundary parse failed at line {line_number}: {line!r}") from exc

            pdb_id = chain_name[:4].lower()
            chain_id = normalise_chain(chain_name[4])
            is_whole_chain = domain_count == 1 and fragment_count == 0
            for domain_index, segments in enumerate(domains, start=1):
                domain_id = f"{chain_name}{domain_index:02d}"
                candidate = {
                    "candidate_id": domain_id,
                    "segments": segments,
                    "is_whole_chain": is_whole_chain,
                    "source_line": line_number,
                }
                # A CATH domain can be discontinuous but its primary chain is the CDF chain name.
                by_pdb_chain[(pdb_id, chain_id)].append(candidate)
    return by_pdb_chain


def parse_ecod_domains(path: Path) -> dict[tuple[str, str], list[dict[str, Any]]]:
    by_pdb_chain: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        header: list[str] | None = None
        for line in handle:
            stripped = line.lstrip("#").strip()
            if header is None:
                if stripped.startswith("uid\t"):
                    header = stripped.split("\t")
                continue
            if line.startswith("#") or not line.strip():
                continue
            values = line.rstrip("\n").split("\t")
            if len(values) < len(header):
                continue
            row = dict(zip(header, values))
            pdb_id = row.get("pdb", "").lower()
            fallback_chain = normalise_chain(row.get("chain", ""))
            segments, parse_problems = parse_segment_string(row.get("pdb_range", ""))
            if not pdb_id:
                continue
            candidate = {
                "candidate_id": row.get("ecod_domain_id", ""),
                "segments": segments,
                "parse_problems": parse_problems,
                "ecod_f_id": row.get("f_id", ""),
                "ecod_architecture_name": row.get("architecture_name", ""),
                "ecod_x_name": row.get("x_name", ""),
                "ecod_h_name": row.get("h_name", ""),
                "ecod_t_name": row.get("t_name", ""),
                "ecod_f_name": row.get("f_name", ""),
            }
            candidate_chains = {segment.chain_id for segment in segments}
            if not candidate_chains and fallback_chain:
                candidate_chains = {fallback_chain}
            for chain_id in candidate_chains:
                by_pdb_chain[(pdb_id, chain_id)].append(candidate)
    return by_pdb_chain


def interval_length(segment: Segment) -> int:
    if segment.is_full_chain or segment.start_resseq is None or segment.end_resseq is None:
        return 0
    return abs(segment.end_resseq - segment.start_resseq) + 1


def segment_intersection_length(left: Segment, right: Segment) -> int:
    if left.chain_id != right.chain_id:
        return 0
    if left.is_full_chain or right.is_full_chain:
        return 0
    assert left.start_resseq is not None and left.end_resseq is not None
    assert right.start_resseq is not None and right.end_resseq is not None
    start = max(min(left.start_resseq, left.end_resseq), min(right.start_resseq, right.end_resseq))
    end = min(max(left.start_resseq, left.end_resseq), max(right.start_resseq, right.end_resseq))
    return max(0, end - start + 1)


def source_has_unbounded_segment(segments: Iterable[Segment]) -> bool:
    return any(segment.is_full_chain for segment in segments)


def overlap_metrics(source: list[Segment], candidate: list[Segment]) -> dict[str, float | str | None]:
    """Calculate reciprocal numeric overlap, refusing undefined full-chain comparisons."""

    if not source or not candidate:
        return {
            "score": None,
            "source_coverage": None,
            "candidate_coverage": None,
            "overlap_residue_span": 0.0,
            "basis": "missing_segments",
        }
    if source_has_unbounded_segment(source) or source_has_unbounded_segment(candidate):
        return {
            "score": None,
            "source_coverage": None,
            "candidate_coverage": None,
            "overlap_residue_span": 0.0,
            "basis": "unbounded_segment",
        }
    source_length = sum(interval_length(segment) for segment in source)
    candidate_length = sum(interval_length(segment) for segment in candidate)
    overlap = sum(segment_intersection_length(a, b) for a in source for b in candidate)
    if source_length <= 0 or candidate_length <= 0:
        return {
            "score": None,
            "source_coverage": None,
            "candidate_coverage": None,
            "overlap_residue_span": float(overlap),
            "basis": "zero_length_segment",
        }
    source_coverage = min(1.0, overlap / source_length)
    candidate_coverage = min(1.0, overlap / candidate_length)
    return {
        "score": min(source_coverage, candidate_coverage),
        "source_coverage": source_coverage,
        "candidate_coverage": candidate_coverage,
        "overlap_residue_span": float(overlap),
        "basis": "numeric_author_residue_overlap",
    }


def select_candidate(
    source_segments: list[Segment],
    candidates: list[dict[str, Any]],
    min_overlap: float,
    ambiguity_margin: float,
    source_is_whole_chain: bool,
) -> tuple[dict[str, Any] | None, str, list[dict[str, Any]]]:
    """Score every candidate and select one only when the match is auditable."""

    if not candidates:
        return None, "no_same_chain_candidate", []
    scored: list[dict[str, Any]] = []
    for candidate in candidates:
        metrics = overlap_metrics(source_segments, candidate["segments"])
        scored.append({**candidate, **metrics})

    # Whole-chain SCOPe records do not contain explicit end coordinates.  The one
    # safe exception is a one-domain, no-fragment CATH chain, which CATH itself
    # declares a whole-chain domain.  ECOD never gets this exception.
    if source_is_whole_chain:
        whole_chain_candidates = [row for row in scored if row.get("is_whole_chain")]
        if len(whole_chain_candidates) == 1:
            row = whole_chain_candidates[0]
            row.update(
                {
                    "score": 1.0,
                    "source_coverage": 1.0,
                    "candidate_coverage": 1.0,
                    "basis": "both_sources_whole_chain",
                }
            )
            return row, "matched_whole_chain", scored
        return None, "unresolved_full_chain_boundary", scored

    qualified = [
        row
        for row in scored
        if row["score"] is not None
        and row["source_coverage"] is not None
        and row["candidate_coverage"] is not None
        and row["source_coverage"] >= min_overlap
        and row["candidate_coverage"] >= min_overlap
    ]
    if not qualified:
        return None, "below_reciprocal_overlap_threshold", scored
    qualified.sort(key=lambda row: (float(row["score"]), float(row["overlap_residue_span"])), reverse=True)
    if len(qualified) > 1 and float(qualified[0]["score"]) - float(qualified[1]["score"]) < ambiguity_margin:
        return None, "ambiguous_overlap_match", scored
    return qualified[0], "matched", scored


def structure_ids_from_directory(path: Path) -> set[str]:
    suffixes = (".cif", ".mmcif", ".pdb", ".ent", ".cif.gz", ".mmcif.gz", ".pdb.gz", ".ent.gz")
    ids: set[str] = set()
    for candidate in path.rglob("*"):
        if not candidate.is_file() or not candidate.name.lower().endswith(suffixes):
            continue
        name = candidate.name.lower()
        if name.endswith(".gz"):
            name = name[:-3]
        ids.add(Path(name).stem)
    return ids


def format_float(value: object) -> str:
    return "" if value is None else f"{float(value):.6f}"


def stable_domain_uid(source_domain_id: str) -> str:
    digest = hashlib.sha256(f"scope-2.08:{source_domain_id}".encode("utf-8")).hexdigest()[:16]
    return f"scope208_{digest}"


def write_csv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


MAPPING_FIELDS = [
    "domain_uid", "source_domain_id", "structure_id", "chain_id", "domain_segments_json",
    "source_domain_range", "source_chain_count", "source_boundary_status", "source_boundary_parse_problems",
    "astral_cluster_id", "scope_sccs", "scope_class_id", "scope_fold_id", "scope_superfamily_id",
    "scope_family_id", "scope_sunid", "source_species", "kingdom", "mapping_status",
    "mapping_overlap_score", "cath_mapping_status", "cath_mapping_overlap_score",
    "cath_mapping_source_coverage", "cath_mapping_candidate_coverage", "cath_mapping_basis",
    "cath_domain_id", "cath_class_id", "cath_architecture_id", "cath_topology_id", "cath_homology_id",
    "ecod_mapping_status", "ecod_mapping_overlap_score", "ecod_mapping_source_coverage",
    "ecod_mapping_candidate_coverage", "ecod_mapping_basis", "ecod_domain_id",
    "ecod_architecture_id", "ecod_x_id", "ecod_h_id", "ecod_t_id", "ecod_f_id",
    "ecod_architecture_name", "ecod_x_name", "ecod_h_name", "ecod_t_name", "ecod_f_name",
    "ecod_hierarchy_status", "quality_flags",
]


CANDIDATE_FIELDS = [
    "domain_uid", "source_domain_id", "structure_id", "source_chain_id", "candidate_source",
    "candidate_id", "candidate_segments_json", "candidate_rank", "candidate_selected",
    "selection_status", "overlap_score", "source_coverage", "candidate_coverage",
    "overlap_residue_span", "overlap_basis", "candidate_parse_problems",
]


def main() -> int:
    args = parse_args()
    if not 0.0 < args.min_reciprocal_overlap <= 1.0:
        raise SystemExit("--min-reciprocal-overlap must be in (0, 1].")
    if not args.cath_boundaries.exists():
        raise SystemExit(
            f"Missing CATH boundary file: {args.cath_boundaries}. "
            "Download it before running V2 mapping; V0 positional mapping is forbidden."
        )

    source_dir = args.classifications_dir
    scope_domains = parse_scope_cla(source_dir / "scope_cla.txt")
    astral40_ids = parse_astral40_ids(source_dir / "scope_astral_40.fa")
    species_by_sunid = parse_scope_species(source_dir / "scope_des.txt")
    cath_metadata = parse_cath_metadata(source_dir / "cath_domain_list.txt")
    cath_by_chain = parse_cath_boundaries(args.cath_boundaries)
    ecod_by_chain = parse_ecod_domains(source_dir / "ecod_domains.txt")
    local_structure_ids = structure_ids_from_directory(args.structures_dir)

    mapping_rows: list[dict[str, str]] = []
    candidate_rows: list[dict[str, str]] = []
    skipped = Counter()

    for source in scope_domains:
        source_domain_id = source["source_domain_id"]
        structure_id = source["structure_id"]
        if not args.include_non_astral40 and source_domain_id not in astral40_ids:
            skipped["not_astral40"] += 1
            continue
        if structure_id not in local_structure_ids:
            skipped["structure_not_local"] += 1
            continue
        segments: list[Segment] = source["segments"]
        chain_ids: list[str] = source["chain_ids"]
        if source["segment_problems"]:
            boundary_status = "parse_failed"
        elif not segments:
            boundary_status = "missing"
        elif source_has_unbounded_segment(segments):
            boundary_status = "whole_chain"
        else:
            boundary_status = "explicit"
        source_is_whole_chain = boundary_status == "whole_chain" and len(chain_ids) == 1
        source_chain_id = chain_ids[0] if len(chain_ids) == 1 else ""
        domain_uid = stable_domain_uid(source_domain_id)
        species = extract_species_description(source["scope_hierarchy_path"], species_by_sunid)
        flags: list[str] = []
        if len(chain_ids) != 1:
            flags.append("not_single_chain_source_domain")
        if boundary_status in {"parse_failed", "missing"}:
            flags.append("source_boundary_unresolved")

        cath_selected: dict[str, Any] | None = None
        ecod_selected: dict[str, Any] | None = None
        cath_status = "not_attempted"
        ecod_status = "not_attempted"
        cath_scored: list[dict[str, Any]] = []
        ecod_scored: list[dict[str, Any]] = []
        if len(chain_ids) == 1 and boundary_status in {"explicit", "whole_chain"}:
            cath_selected, cath_status, cath_scored = select_candidate(
                segments,
                cath_by_chain.get((structure_id, source_chain_id), []),
                args.min_reciprocal_overlap,
                args.ambiguity_margin,
                source_is_whole_chain,
            )
            # ECOD cannot use the whole-chain exception: the current ECOD file
            # does not expose an equivalent whole-chain declaration.
            ecod_selected, ecod_status, ecod_scored = select_candidate(
                segments,
                ecod_by_chain.get((structure_id, source_chain_id), []),
                args.min_reciprocal_overlap,
                args.ambiguity_margin,
                False,
            )
        else:
            cath_status = "source_domain_not_single_chain_or_unresolved"
            ecod_status = "source_domain_not_single_chain_or_unresolved"

        for candidate_source, scored, selected, status in (
            ("CATH", cath_scored, cath_selected, cath_status),
            ("ECOD", ecod_scored, ecod_selected, ecod_status),
        ):
            ordered = sorted(
                scored,
                key=lambda row: (
                    -1.0 if row.get("score") is None else -float(row["score"]),
                    -float(row.get("overlap_residue_span", 0.0)),
                    str(row.get("candidate_id", "")),
                ),
            )
            for rank, candidate in enumerate(ordered, start=1):
                candidate_rows.append(
                    {
                        "domain_uid": domain_uid,
                        "source_domain_id": source_domain_id,
                        "structure_id": structure_id,
                        "source_chain_id": source_chain_id,
                        "candidate_source": candidate_source,
                        "candidate_id": str(candidate.get("candidate_id", "")),
                        "candidate_segments_json": json.dumps(
                            [segment.as_dict() for segment in candidate.get("segments", [])], separators=(",", ":")
                        ),
                        "candidate_rank": str(rank),
                        "candidate_selected": str(candidate is selected).lower(),
                        "selection_status": status,
                        "overlap_score": format_float(candidate.get("score")),
                        "source_coverage": format_float(candidate.get("source_coverage")),
                        "candidate_coverage": format_float(candidate.get("candidate_coverage")),
                        "overlap_residue_span": format_float(candidate.get("overlap_residue_span")),
                        "overlap_basis": str(candidate.get("basis", "")),
                        "candidate_parse_problems": ";".join(candidate.get("parse_problems", [])),
                    }
                )

        if cath_status.startswith("matched"):
            mapping_status = "mapped"
        elif boundary_status in {"parse_failed", "missing"}:
            mapping_status = "source_boundary_unresolved"
        elif len(chain_ids) != 1:
            mapping_status = "not_single_chain_source_domain"
        else:
            mapping_status = "needs_mapping_review"
        if not cath_status.startswith("matched"):
            flags.append(f"cath_{cath_status}")
        if not ecod_status.startswith("matched"):
            flags.append(f"ecod_{ecod_status}")
        if ecod_selected is not None:
            flags.append("ecod_t_group_id_not_available_in_source_file")

        cath_fields = cath_metadata.get(cath_selected["candidate_id"], {}) if cath_selected else {}
        ecod_fields = ecod_selected or {}
        scores = [row.get("score") for row in (cath_selected, ecod_selected) if row and row.get("score") is not None]
        mapping_rows.append(
            {
                "domain_uid": domain_uid,
                "source_domain_id": source_domain_id,
                "structure_id": structure_id,
                "chain_id": source_chain_id,
                "domain_segments_json": json.dumps([segment.as_dict() for segment in segments], separators=(",", ":")),
                "source_domain_range": source["source_range"],
                "source_chain_count": str(len(chain_ids)),
                "source_boundary_status": boundary_status,
                "source_boundary_parse_problems": ";".join(source["segment_problems"]),
                "astral_cluster_id": "ASTRAL40" if source_domain_id in astral40_ids else "",
                "scope_sccs": source["scope_sccs"],
                "scope_class_id": source["scope_class_id"],
                "scope_fold_id": source["scope_fold_id"],
                "scope_superfamily_id": source["scope_superfamily_id"],
                "scope_family_id": source["scope_family_id"],
                "scope_sunid": source["scope_sunid"],
                "source_species": species,
                "kingdom": infer_kingdom(species),
                "mapping_status": mapping_status,
                "mapping_overlap_score": format_float(min(scores) if scores else None),
                "cath_mapping_status": cath_status,
                "cath_mapping_overlap_score": format_float(cath_selected.get("score") if cath_selected else None),
                "cath_mapping_source_coverage": format_float(cath_selected.get("source_coverage") if cath_selected else None),
                "cath_mapping_candidate_coverage": format_float(cath_selected.get("candidate_coverage") if cath_selected else None),
                "cath_mapping_basis": str(cath_selected.get("basis", "") if cath_selected else ""),
                "cath_domain_id": cath_fields.get("cath_domain_id", ""),
                "cath_class_id": cath_fields.get("cath_class_id", ""),
                "cath_architecture_id": cath_fields.get("cath_architecture_id", ""),
                "cath_topology_id": cath_fields.get("cath_topology_id", ""),
                "cath_homology_id": cath_fields.get("cath_homology_id", ""),
                "ecod_mapping_status": ecod_status,
                "ecod_mapping_overlap_score": format_float(ecod_selected.get("score") if ecod_selected else None),
                "ecod_mapping_source_coverage": format_float(ecod_selected.get("source_coverage") if ecod_selected else None),
                "ecod_mapping_candidate_coverage": format_float(ecod_selected.get("candidate_coverage") if ecod_selected else None),
                "ecod_mapping_basis": str(ecod_selected.get("basis", "") if ecod_selected else ""),
                "ecod_domain_id": ecod_fields.get("candidate_id", ""),
                # The ECOD domain list provides F-group IDs and group names but
                # does not expose distinct machine-readable A/X/H/T IDs.  Leave
                # those blank rather than inventing a topology ID from F-group.
                "ecod_architecture_id": "",
                "ecod_x_id": "",
                "ecod_h_id": "",
                "ecod_t_id": "",
                "ecod_f_id": ecod_fields.get("ecod_f_id", ""),
                "ecod_architecture_name": ecod_fields.get("ecod_architecture_name", ""),
                "ecod_x_name": ecod_fields.get("ecod_x_name", ""),
                "ecod_h_name": ecod_fields.get("ecod_h_name", ""),
                "ecod_t_name": ecod_fields.get("ecod_t_name", ""),
                "ecod_f_name": ecod_fields.get("ecod_f_name", ""),
                "ecod_hierarchy_status": "needs_machine_readable_t_group_source" if ecod_selected else "not_mapped",
                "quality_flags": ";".join(flags),
            }
        )

    mapping_rows.sort(key=lambda row: row["source_domain_id"])
    candidate_rows.sort(key=lambda row: (row["source_domain_id"], row["candidate_source"], int(row["candidate_rank"])))
    write_csv(args.output, MAPPING_FIELDS, mapping_rows)
    write_csv(args.candidates_output, CANDIDATE_FIELDS, candidate_rows)

    diagnostics = {
        "configuration": {
            "scope_release": "2.08",
            "cath_boundary_file": str(args.cath_boundaries),
            "minimum_reciprocal_overlap": args.min_reciprocal_overlap,
            "ambiguity_margin": args.ambiguity_margin,
            "astral40_only": not args.include_non_astral40,
        },
        "counts": {
            "mapping_rows": len(mapping_rows),
            "candidate_rows": len(candidate_rows),
            "skipped": dict(sorted(skipped.items())),
            "mapping_status": dict(sorted(Counter(row["mapping_status"] for row in mapping_rows).items())),
            "source_boundary_status": dict(sorted(Counter(row["source_boundary_status"] for row in mapping_rows).items())),
            "cath_mapping_status": dict(sorted(Counter(row["cath_mapping_status"] for row in mapping_rows).items())),
            "ecod_mapping_status": dict(sorted(Counter(row["ecod_mapping_status"] for row in mapping_rows).items())),
        },
        "integrity_notes": [
            "No V0 first-on-chain mapping was used.",
            "All CATH and ECOD candidates are retained in the companion candidate table.",
            "ECOD T-group IDs are intentionally blank until a source exposing them is added.",
        ],
    }
    args.diagnostics_output.parent.mkdir(parents=True, exist_ok=True)
    args.diagnostics_output.write_text(json.dumps(diagnostics, indent=2) + "\n", encoding="utf-8")

    print(f"Wrote {len(mapping_rows):,} V2 domain mappings to {args.output}")
    print(f"Wrote {len(candidate_rows):,} mapping candidates to {args.candidates_output}")
    print(json.dumps(diagnostics["counts"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

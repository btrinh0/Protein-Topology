from __future__ import annotations

"""Compute domain-sliced contacts and provisional SSE contact signatures."""

import argparse
import csv
import gzip
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO

from Bio.PDB import MMCIFParser, PDBParser, Polypeptide
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


REPRESENTATIVE_ATOMS = ("CB", "CA")
CONTACT_DEFINITION_ID = "residue_cb_ca_8A_sep4_v2_domain_slice"


@dataclass(frozen=True)
class ResidueKey:
    number: int
    insertion: str

    def order(self) -> tuple[int, int]:
        return (self.number, 0 if not self.insertion else ord(self.insertion[0]))


@dataclass(frozen=True)
class Annotation:
    state: str
    chain_id: str
    start: ResidueKey
    end: ResidueKey
    sheet_id: str = ""
    range_id: str = ""


@dataclass(frozen=True)
class SSESegment:
    state: str
    start_position: int
    end_position: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute domain-specific V2 contact features and SSE signatures."
    )
    parser.add_argument(
        "--mapping",
        type=Path,
        default=Path("data/processed/v2/domain_mapping_v2.csv"),
    )
    parser.add_argument(
        "--chain-manifest",
        type=Path,
        default=Path("data/processed/chain_manifest.csv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/processed/v2/domain_topology_v2.csv"),
    )
    parser.add_argument(
        "--errors-output",
        type=Path,
        default=Path("data/processed/v2/domain_topology_errors_v2.csv"),
    )
    parser.add_argument(
        "--config-output",
        type=Path,
        default=Path("data/processed/v2/domain_topology_config_v2.json"),
    )
    parser.add_argument("--distance-threshold", type=float, default=8.0)
    parser.add_argument("--min-sequence-separation", type=int, default=4)
    parser.add_argument("--min-sse-length", type=int, default=3)
    parser.add_argument("--min-sse-contact-support", type=int, default=2)
    return parser.parse_args()


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def detect_format(path: Path) -> str:
    name = path.name.lower()
    if name.endswith((".cif", ".mmcif", ".cif.gz", ".mmcif.gz")):
        return "mmcif"
    if name.endswith((".pdb", ".ent", ".pdb.gz", ".ent.gz")):
        return "pdb"
    raise ValueError(f"Unsupported structure format: {path}")


def load_structure(path: Path, file_format: str, structure_id: str):
    parser = MMCIFParser(QUIET=True) if file_format == "mmcif" else PDBParser(QUIET=True)
    with open_text(path) as handle:
        return parser.get_structure(structure_id, handle)


def normalise_insertion(value: object) -> str:
    text = str(value or "").strip()
    return "" if text in {"", "?", ".", "-"} else text


def parse_residue_key(number: object, insertion: object = "") -> ResidueKey:
    return ResidueKey(int(str(number).strip()), normalise_insertion(insertion))


def residue_key(residue) -> ResidueKey:
    return ResidueKey(int(residue.id[1]), normalise_insertion(residue.id[2]))


def key_in_range(key: ResidueKey, start: ResidueKey, end: ResidueKey) -> bool:
    left, right = sorted((start.order(), end.order()))
    return left <= key.order() <= right


def value_at(values: object, index: int, default: str = "") -> str:
    if isinstance(values, list) and index < len(values):
        value = values[index]
    elif values is None:
        value = default
    else:
        value = values
    text = str(value)
    return default if text in {"?", ".", "None"} else text


def author_annotations(path: Path) -> list[Annotation]:
    if detect_format(path) != "mmcif":
        return []
    with open_text(path) as handle:
        data = MMCIF2Dict(handle)
    annotations: list[Annotation] = []
    conf_types = data.get("_struct_conf.conf_type_id", [])
    conf_count = len(conf_types) if isinstance(conf_types, list) else 0
    for index in range(conf_count):
        kind = value_at(conf_types, index).upper()
        state = "H" if kind.startswith("HELX") else ""
        if not state:
            continue
        try:
            annotations.append(
                Annotation(
                    state=state,
                    chain_id=value_at(data.get("_struct_conf.beg_auth_asym_id", []), index),
                    start=parse_residue_key(
                        value_at(data.get("_struct_conf.beg_auth_seq_id", []), index),
                        value_at(data.get("_struct_conf.pdbx_beg_PDB_ins_code", []), index),
                    ),
                    end=parse_residue_key(
                        value_at(data.get("_struct_conf.end_auth_seq_id", []), index),
                        value_at(data.get("_struct_conf.pdbx_end_PDB_ins_code", []), index),
                    ),
                )
            )
        except (TypeError, ValueError):
            continue
    sheet_ids = data.get("_struct_sheet_range.id", [])
    sheet_count = len(sheet_ids) if isinstance(sheet_ids, list) else 0
    for index in range(sheet_count):
        try:
            annotations.append(
                Annotation(
                    state="E",
                    chain_id=value_at(data.get("_struct_sheet_range.beg_auth_asym_id", []), index),
                    start=parse_residue_key(
                        value_at(data.get("_struct_sheet_range.beg_auth_seq_id", []), index),
                        value_at(data.get("_struct_sheet_range.pdbx_beg_PDB_ins_code", []), index),
                    ),
                    end=parse_residue_key(
                        value_at(data.get("_struct_sheet_range.end_auth_seq_id", []), index),
                        value_at(data.get("_struct_sheet_range.pdbx_end_PDB_ins_code", []), index),
                    ),
                    sheet_id=value_at(data.get("_struct_sheet_range.sheet_id", []), index),
                    range_id=value_at(data.get("_struct_sheet_range.id", []), index),
                )
            )
        except (TypeError, ValueError):
            continue
    return annotations


def representative_atom(residue):
    for atom_name in REPRESENTATIVE_ATOMS:
        if atom_name == "CA" and residue.get_resname() != "GLY" and "CB" in residue:
            continue
        if atom_name in residue:
            return residue[atom_name]
    return None


def amino_acid_residues(chain) -> list[tuple[object, ResidueKey]]:
    result: list[tuple[object, ResidueKey]] = []
    for residue in chain:
        if residue.id[0] == " " and Polypeptide.is_aa(residue, standard=False):
            result.append((residue, residue_key(residue)))
    return result


def segment_contains(key: ResidueKey, segment: dict[str, object]) -> bool:
    if segment.get("is_full_chain"):
        return True
    start = parse_residue_key(segment.get("start_resseq"), segment.get("start_icode", ""))
    end = parse_residue_key(segment.get("end_resseq"), segment.get("end_icode", ""))
    return key_in_range(key, start, end)


def select_domain_residues(
    chain_residues: list[tuple[object, ResidueKey]],
    segments: list[dict[str, object]],
    chain_id: str,
) -> list[tuple[object, ResidueKey]]:
    selected: list[tuple[object, ResidueKey]] = []
    chain_segments = [segment for segment in segments if str(segment.get("chain_id", "")) == chain_id]
    for residue, key in chain_residues:
        if any(segment_contains(key, segment) for segment in chain_segments):
            selected.append((residue, key))
    return selected


def contact_pairs(
    residues: list[tuple[object, ResidueKey]],
    distance_threshold: float,
    min_sequence_separation: int,
) -> tuple[list[tuple[int, int]], list[tuple[object, ResidueKey]]]:
    eligible = [(residue, key) for residue, key in residues if representative_atom(residue) is not None]
    contacts: list[tuple[int, int]] = []
    threshold_sq = distance_threshold * distance_threshold
    for left_index, (left_residue, _) in enumerate(eligible):
        left_atom = representative_atom(left_residue)
        for right_index in range(left_index + 1, len(eligible)):
            right_residue, _ = eligible[right_index]
            if right_index - left_index < min_sequence_separation:
                continue
            right_atom = representative_atom(right_residue)
            delta = left_atom.coord - right_atom.coord
            distance_sq = float(delta[0] ** 2 + delta[1] ** 2 + delta[2] ** 2)
            if distance_sq <= threshold_sq:
                contacts.append((left_index, right_index))
    return contacts, eligible


def circuit_relations(contacts: list[tuple[int, int]]) -> tuple[int, int, int, int]:
    series = parallel = cross = pair_count = 0
    for index, (left_start, left_end) in enumerate(contacts):
        for right_start, right_end in contacts[index + 1 :]:
            if len({left_start, left_end, right_start, right_end}) < 4:
                continue
            if left_end < right_start or right_end < left_start:
                relation = "series"
            elif (left_start < right_start < right_end < left_end) or (right_start < left_start < left_end < right_end):
                relation = "parallel"
            elif (left_start < right_start < left_end < right_end) or (right_start < left_start < right_end < left_end):
                relation = "cross"
            else:
                continue
            pair_count += 1
            if relation == "series":
                series += 1
            elif relation == "parallel":
                parallel += 1
            else:
                cross += 1
    return series, parallel, cross, pair_count


def fraction(numerator: int, denominator: int) -> str:
    return "" if denominator == 0 else f"{numerator / denominator:.6f}"


def states_for_residues(
    residues: list[tuple[object, ResidueKey]],
    annotations: list[Annotation],
    chain_id: str,
) -> list[str]:
    states = ["L"] * len(residues)
    for index, (_, key) in enumerate(residues):
        matching = [
            annotation.state
            for annotation in annotations
            if annotation.chain_id == chain_id and key_in_range(key, annotation.start, annotation.end)
        ]
        if "H" in matching:
            states[index] = "H"
        elif "E" in matching:
            states[index] = "E"
    return states


def build_sse_segments(states: list[str], min_sse_length: int) -> tuple[list[SSESegment], dict[int, int]]:
    segments: list[SSESegment] = []
    index = 0
    while index < len(states):
        state = states[index]
        if state not in {"H", "E"}:
            index += 1
            continue
        end = index
        while end + 1 < len(states) and states[end + 1] == state:
            end += 1
        if end - index + 1 >= min_sse_length:
            segments.append(SSESegment(state, index, end))
        index = end + 1
    position_to_sse: dict[int, int] = {}
    for sse_index, segment in enumerate(segments):
        for position in range(segment.start_position, segment.end_position + 1):
            position_to_sse[position] = sse_index
    return segments, position_to_sse


def topology_signature(
    segments: list[SSESegment],
    eligible_indices: list[int],
    contacts: list[tuple[int, int]],
    min_support: int,
) -> tuple[str, int, str]:
    if not segments:
        return "", 0, "no_sse_segments"
    position_to_sse = {segment_position: sse_index for sse_index, segment in enumerate(segments) for segment_position in range(segment.start_position, segment.end_position + 1)}
    edge_counts: dict[tuple[int, int], int] = defaultdict(int)
    for left, right in contacts:
        left_sse = position_to_sse.get(eligible_indices[left])
        right_sse = position_to_sse.get(eligible_indices[right])
        if left_sse is None or right_sse is None or left_sse == right_sse:
            continue
        edge_counts[tuple(sorted((left_sse, right_sse)))] += 1
    edges: list[str] = []
    for (left_sse, right_sse), count in sorted(edge_counts.items()):
        if count < min_support:
            continue
        left_state = segments[left_sse].state
        right_state = segments[right_sse].state
        edges.append(f"{left_sse + 1}-{right_sse + 1}:{left_state}{right_state}")
    nodes = "-".join(segment.state for segment in segments)
    return f"{nodes} | {';'.join(edges)}", len(edges), "computed"


def safe_float(value: float | None) -> str:
    return "" if value is None or not math.isfinite(value) else f"{value:.6f}"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


OUTPUT_FIELDS = [
    "domain_uid", "source_domain_id", "structure_id", "chain_id", "source_path", "domain_segments_json",
    "topology_status", "analysis_candidate_status", "sse_assignment_method", "sse_assignment_status",
    "sse_count", "sse_order", "sse_segments_json", "topology_signature_exact", "topology_macroclass",
    "macroclass_status", "sse_contact_edge_count", "sse_orientation_status", "resolved_residue_count",
    "eligible_residue_count", "source_span_residue_count", "coordinate_coverage_fraction",
    "missing_coordinate_fraction", "contact_definition_id", "contact_distance_threshold_angstrom",
    "min_sequence_separation", "min_sse_contact_support", "contact_count", "contact_density", "contact_order",
    "contact_pair_count", "series_pair_count", "parallel_pair_count", "cross_pair_count", "series_fraction",
    "parallel_fraction", "cross_fraction", "quality_flags",
]


def main() -> int:
    args = parse_args()
    mapping_rows = read_csv(args.mapping)
    chain_rows = read_csv(args.chain_manifest)
    path_by_chain = {(row["structure_id"], row["chain_id"]): row for row in chain_rows}
    domains_by_structure: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in mapping_rows:
        domains_by_structure[row["structure_id"]].append(row)

    output_rows: list[dict[str, str]] = []
    error_rows: list[dict[str, str]] = []
    for structure_id in sorted(domains_by_structure):
        structure_domains = domains_by_structure[structure_id]
        structure_path = ""
        structure_format = ""
        for domain in structure_domains:
            chain_manifest_row = path_by_chain.get((structure_id, domain["chain_id"]))
            if chain_manifest_row:
                structure_path = chain_manifest_row["source_path"]
                structure_format = chain_manifest_row["file_format"]
                break
        if not structure_path:
            for domain in structure_domains:
                error_rows.append({"domain_uid": domain["domain_uid"], "source_domain_id": domain["source_domain_id"], "structure_id": structure_id, "error": "structure_not_in_chain_manifest"})
            continue
        try:
            structure = load_structure(Path(structure_path), structure_format or detect_format(Path(structure_path)), structure_id)
            annotations = author_annotations(Path(structure_path))
            first_model = next(structure.get_models())
        except Exception as exc:
            for domain in structure_domains:
                error_rows.append({"domain_uid": domain["domain_uid"], "source_domain_id": domain["source_domain_id"], "structure_id": structure_id, "error": str(exc)})
            continue
        chain_cache = {chain.id: amino_acid_residues(chain) for chain in first_model}
        for domain in structure_domains:
            flags: list[str] = []
            try:
                segments = json.loads(domain["domain_segments_json"])
                chain_id = domain["chain_id"]
                chain_residue_rows = chain_cache.get(chain_id, [])
                selected_residues = select_domain_residues(chain_residue_rows, segments, chain_id)
                domain_residue_count = len(selected_residues)
                eligible_contacts, eligible_residues = contact_pairs(
                    selected_residues,
                    args.distance_threshold,
                    args.min_sequence_separation,
                )
                eligible_positions = [selected_residues.index(item) for item in eligible_residues]
                states = states_for_residues(selected_residues, annotations, chain_id)
                sse_segments, position_to_sse = build_sse_segments(states, args.min_sse_length)
                signature, edge_count, signature_status = topology_signature(
                    sse_segments,
                    eligible_positions,
                    eligible_contacts,
                    args.min_sse_contact_support,
                )
                series, parallel, cross, relation_count = circuit_relations(eligible_contacts)
                span_count = sum(
                    abs(int(segment.get("end_resseq")) - int(segment.get("start_resseq"))) + 1
                    for segment in segments
                    if not segment.get("is_full_chain")
                )
                if any(segment.get("is_full_chain") for segment in segments):
                    span_count = len(chain_residue_rows)
                coverage = domain_residue_count / span_count if span_count else None
                contact_order = (
                    sum(right - left for left, right in eligible_contacts)
                    / (len(eligible_contacts) * max(len(eligible_residues) - 1, 1))
                    if eligible_contacts else None
                )
                contact_density = len(eligible_contacts) / max(len(eligible_residues), 1)
                if len(domain["chain_id"]) == 0:
                    flags.append("chain_id_unresolved")
                if not annotations:
                    flags.append("author_sse_annotations_unavailable")
                if signature_status != "computed":
                    flags.append(signature_status)
                if domain.get("source_boundary_status") == "whole_chain":
                    flags.append("whole_chain_source_boundary")
                flags.append("sse_assignment_requires_dssp_or_stride_calibration")
                if domain_residue_count < 70 or domain_residue_count > 160:
                    candidate_status = "outside_provisional_70_160_window"
                elif not any(segment.state == "H" for segment in sse_segments) or not any(segment.state == "E" for segment in sse_segments):
                    candidate_status = "not_alpha_beta_by_author_annotation"
                else:
                    candidate_status = "provisional_alpha_beta_candidate"
                output_rows.append(
                    {
                        **{field: domain.get(field, "") for field in ("domain_uid", "source_domain_id", "structure_id", "chain_id")},
                        "source_path": structure_path,
                        "domain_segments_json": domain["domain_segments_json"],
                        "topology_status": "computed",
                        "analysis_candidate_status": candidate_status,
                        "sse_assignment_method": "author_annotation",
                        "sse_assignment_status": "provisional",
                        "sse_count": str(len(sse_segments)),
                        "sse_order": "-".join(segment.state for segment in sse_segments),
                        "sse_segments_json": json.dumps([segment.__dict__ for segment in sse_segments], separators=(",", ":")),
                        "topology_signature_exact": signature,
                        "topology_macroclass": "",
                        "macroclass_status": "awaiting_calibration",
                        "sse_contact_edge_count": str(edge_count),
                        "sse_orientation_status": "not_assigned",
                        "resolved_residue_count": str(domain_residue_count),
                        "eligible_residue_count": str(len(eligible_residues)),
                        "source_span_residue_count": str(span_count),
                        "coordinate_coverage_fraction": safe_float(coverage),
                        "missing_coordinate_fraction": safe_float(1.0 - coverage if coverage is not None else None),
                        "contact_definition_id": CONTACT_DEFINITION_ID,
                        "contact_distance_threshold_angstrom": f"{args.distance_threshold:.3f}",
                        "min_sequence_separation": str(args.min_sequence_separation),
                        "min_sse_contact_support": str(args.min_sse_contact_support),
                        "contact_count": str(len(eligible_contacts)),
                        "contact_density": f"{contact_density:.6f}",
                        "contact_order": safe_float(contact_order),
                        "contact_pair_count": str(relation_count),
                        "series_pair_count": str(series),
                        "parallel_pair_count": str(parallel),
                        "cross_pair_count": str(cross),
                        "series_fraction": fraction(series, relation_count),
                        "parallel_fraction": fraction(parallel, relation_count),
                        "cross_fraction": fraction(cross, relation_count),
                        "quality_flags": ";".join(sorted(set(flags))),
                    }
                )
            except Exception as exc:
                error_rows.append({"domain_uid": domain["domain_uid"], "source_domain_id": domain["source_domain_id"], "structure_id": structure_id, "error": str(exc)})

    output_rows.sort(key=lambda row: row["domain_uid"])
    write_csv(args.output, OUTPUT_FIELDS, output_rows)
    write_csv(args.errors_output, ["domain_uid", "source_domain_id", "structure_id", "error"], error_rows)
    config = {
        "contact_definition_id": CONTACT_DEFINITION_ID,
        "distance_threshold_angstrom": args.distance_threshold,
        "min_sequence_separation": args.min_sequence_separation,
        "min_sse_length": args.min_sse_length,
        "min_sse_contact_support": args.min_sse_contact_support,
        "sse_assignment_method": "author_annotation",
        "sse_assignment_status": "provisional; replace or compare with DSSP/STRIDE during calibration",
        "macroclass_status": "not assigned before calibration",
    }
    args.config_output.parent.mkdir(parents=True, exist_ok=True)
    args.config_output.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote topology rows for {len(output_rows):,} domains and {len(error_rows):,} errors.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

# Copyright (c) 2026, btrinh0.
# SPDX-License-Identifier: LicenseRef-Biopython OR BSD-3-Clause
"""Audit positive-inside score patterns in annotated Biopython records.

Biopython is an optional dependency. Install membrane-rule-audit[biopython]
to use this module. SeqIO's UniProt XML parser supplies the feature names used
here; other sources must provide the same explicit annotation contract.

"""

from __future__ import annotations

import operator

import numpy as np

try:
    from Bio.SeqFeature import ExactPosition
    from Bio.SeqRecord import SeqRecord
except ModuleNotFoundError as exc:
    if exc.name == "Bio":
        raise ImportError(
            "Install the optional dependency with: "
            'python -m pip install "membrane-rule-audit[biopython]"'
        ) from exc
    raise

from . import panel
from .audit import flank_map

_SIDES = {
    "cytoplasmic": "inside",
    "cytosolic": "inside",
    "inside": "inside",
    "extracellular": "outside",
    "lumenal": "outside",
    "luminal": "outside",
    "periplasmic": "outside",
    "outside": "outside",
}


def _integer(value, name, minimum):
    """Validate an integer option without accepting Boolean values."""
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be an integer >= {minimum}")
    try:
        value = operator.index(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be an integer >= {minimum}") from exc
    if value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def topology_from_seqrecord(record):
    """Convert exact UniProt-style membrane features to the audit schema.

    Parameters
    ----------
    record : Bio.SeqRecord.SeqRecord
        A protein record with an explicit ID and sequence. If present, the
        molecule_type annotation must be "protein". Relevant features
        have type "transmembrane region" or "topological domain".
        Their locations must be exact, local, contiguous, non-overlapping
        intervals with strand None or +1.

    Returns
    -------
    dict
        Sequence, ID, transmembrane intervals, and topological domains in
        one-based inclusive coordinates. Domain descriptions Cytoplasmic,
        Cytosolic, and Inside map to the cytosolic side; Extracellular,
        Lumenal, Luminal, Periplasmic, and Outside map to the other side.
        Matching ignores case and surrounding whitespace. Unrecognized or
        absent descriptions remain unknown and cannot enter the comparison.
        A description may be a string or a singleton list.

    Raises
    ------
    ValueError
        If an ID, sequence, transmembrane annotation, or relevant feature
        location is missing or inconsistent.

    Notes
    -----
    Biopython uses zero-based, end-exclusive locations. A feature [5:15]
    becomes residues 6 through 15. The input record is left unchanged.
    Unknown sides are never inferred from alternating membrane crossings.
    """
    if not isinstance(record, SeqRecord):
        raise ValueError("Provide Biopython SeqRecord objects")
    if not isinstance(record.id, str) or not record.id or record.id == "<unknown id>":
        raise ValueError("Every record needs a unique, explicit sequence ID")
    molecule_type = record.annotations.get("molecule_type")
    if molecule_type is not None and (
        not isinstance(molecule_type, str) or molecule_type.lower() != "protein"
    ):
        raise ValueError("The molecule_type annotation must identify a protein")
    if record.seq is None or not len(record.seq):
        raise ValueError("The protein sequence is empty")
    sequence = str(record.seq).upper()
    segments, domains, intervals = [], [], []
    for feature in record.features:
        if feature.type not in {"transmembrane region", "topological domain"}:
            continue
        location = feature.location
        if (
            location is None
            or len(location.parts) != 1
            or type(location.start) is not ExactPosition
            or type(location.end) is not ExactPosition
            or location.ref is not None
            or location.ref_db is not None
            or location.strand not in {None, 1}
        ):
            raise ValueError(
                "Relevant features require exact local contiguous protein locations"
            )
        start, end = int(location.start) + 1, int(location.end)
        if not 1 <= start <= end <= len(sequence):
            raise ValueError("A feature extends outside the protein sequence")
        intervals.append((start, end))
        if feature.type == "transmembrane region":
            segments.append({"start": start, "end": end})
            continue
        description = feature.qualifiers.get("description", "")
        if isinstance(description, list):
            if len(description) != 1:
                raise ValueError("A domain description list must contain one string")
            description = description[0]
        if not isinstance(description, str):
            raise ValueError("Side descriptions must be strings")
        side = _SIDES.get(description.strip().lower(), "unknown")
        domains.append({"start": start, "end": end, "location": side})
    if not segments:
        raise ValueError("A record needs annotated transmembrane regions")
    intervals.sort()
    if any(left[1] >= right[0] for left, right in zip(intervals, intervals[1:])):
        raise ValueError("Relevant feature intervals overlap")
    return {
        "uniprot": record.id,
        "sequence": sequence,
        "tmds": sorted(segments, key=lambda item: item["start"]),
        "loops": {
            str(i): domain
            for i, domain in enumerate(
                sorted(domains, key=lambda item: item["start"]), 1
            )
        },
    }


def audit_seqrecords(
    records,
    variant_scores,
    *,
    score_direction,
    window=10,
    draws=1999,
    seed=20261003,
):
    """Measure a matched positive-inside contrast from existing variant scores.

    Parameters
    ----------
    records : iterable of Bio.SeqRecord.SeqRecord
        Annotated protein records satisfying topology_from_seqrecord.
        IDs must be unique and match the score rows exactly.
    variant_scores : iterable of tuple
        Five-tuples (record_id, position, reference, alternate, score).
        Positions are one-based integers in the supplied sequence. Residues
        are uppercase standard amino acids and must differ. Scores may be
        finite numbers or NaN. Missing substitution rows are also allowed.
        Conflicting finite duplicates in retained flanks raise ValueError;
        a finite duplicate replaces a missing value.
    score_direction : str
        Required: "higher-is-damaging" or "higher-is-tolerated".
        The latter is negated before the calculation.
    window : int, default 10
        Maximum distance in residues from a transmembrane edge. Only domains
        directly adjacent to an annotated transmembrane interval contribute.
    draws : int, default 1999
        Number of whole-protein bootstrap draws; 0 disables the interval.
        An interval also requires at least two contributing proteins.
    seed : int, default 20261003
        Nonnegative NumPy random-generator seed.

    Returns
    -------
    dict
        Input settings, annotation counts, and a result dictionary with
        estimate_raw, ci95_raw, and contributing site/stratum/protein
        counts. A negative estimate indicates a lower relative penalty for
        K/R introduction on the cytosolic side. No matched strata yields
        estimate_raw=None and ci95_raw=None.

    Notes
    -----
    The existing panel estimator first subtracts each site's mean uncharged
    alternative score from its mean K/R score. It then compares membrane
    sides within protein/reference-residue strata, weighting each difference
    by n_inside*n_outside/(n_inside+n_outside). A site needs an uncharged
    reference, at least one finite K/R score, and one finite uncharged score.
    The uncharged comparator is ACFGILMNPQSTVWY.

    This reports an output-score pattern on the supplied annotations and
    comparator. It does not fit a predictor or return a clinical verdict.
    """
    if score_direction not in {"higher-is-damaging", "higher-is-tolerated"}:
        raise ValueError("Declare higher-is-damaging or higher-is-tolerated")
    window = _integer(window, "window", 1)
    draws = _integer(draws, "draws", 0)
    seed = _integer(seed, "seed", 0)
    topology = {}
    for record in records:
        converted = topology_from_seqrecord(record)
        accession = converted["uniprot"]
        if accession in topology:
            raise ValueError("Duplicate sequence IDs")
        topology[accession] = converted
    if not topology:
        raise ValueError("Provide at least one annotated protein record")
    accessions = sorted(topology)
    protein_index = {accession: i for i, accession in enumerate(accessions)}
    flanks = {
        accession: flank_map(value, window) for accession, value in topology.items()
    }
    positions = {}
    n_input_rows = 0
    for accession, position, reference, alternate, score in variant_scores:
        n_input_rows += 1
        if accession not in topology:
            raise ValueError("A score refers to an unknown sequence ID")
        position = _integer(position, "position", 1)
        sequence = topology[accession]["sequence"]
        if position > len(sequence):
            raise ValueError("A score position is outside its protein sequence")
        if (
            not isinstance(reference, str)
            or not isinstance(alternate, str)
            or len(reference) != 1
            or len(alternate) != 1
            or reference not in panel.ALPHABET
            or alternate not in panel.ALPHABET
            or reference == alternate
        ):
            raise ValueError("Scores must describe standard missense substitutions")
        if sequence[position - 1] != reference:
            raise ValueError("A score's reference residue disagrees with its sequence")
        try:
            score = float(score)
        except (TypeError, ValueError) as exc:
            raise ValueError("Scores must be numbers or NaN") from exc
        if np.isinf(score):
            raise ValueError("Scores may be finite or missing, but not infinite")
        if position not in flanks[accession]:
            continue
        site = positions.setdefault((accession, position), np.full(20, np.nan))
        alternate_index = panel.ALPHABET.index(alternate)
        previous = site[alternate_index]
        if np.isfinite(previous):
            if np.isfinite(score) and abs(previous - score) > 1e-12:
                raise ValueError("Conflicting duplicate substitution scores")
        else:
            site[alternate_index] = score
    keys = sorted(positions)
    score_matrix = np.array([positions[key] for key in keys], dtype=float).reshape(
        -1, 20
    )
    if score_direction == "higher-is-tolerated":
        score_matrix = -score_matrix
    mapped = {
        "protein": np.array([protein_index[a] for a, _ in keys], dtype=int),
        "ref_aa": np.array(
            [panel.ALPHABET.index(topology[a]["sequence"][p - 1]) for a, p in keys],
            dtype=int,
        ),
        "side": np.array([flanks[a][p][0] for a, p in keys], dtype=int),
        "scores": score_matrix,
    }
    contrast = panel.gain_contrast(mapped)
    parts = panel.stratum_parts(mapped, contrast["penalty"], contrast["eligible"])
    return {
        "rule": "positive_charge_gain",
        "comparator": panel.NEUTRAL,
        "input_type": "Biopython SeqRecord with exact UniProt-style features",
        "records": len(accessions),
        "input_score_rows": n_input_rows,
        "annotated_scored_flank_sites": len(keys),
        "ignored_unknown_side_domains": sum(
            domain["location"] == "unknown"
            for value in topology.values()
            for domain in value["loops"].values()
        ),
        "input_score_direction": score_direction,
        "output_score_direction": "higher-is-damaging",
        "window": window,
        "draws": draws,
        "seed": seed,
        "sign_convention": "cytosolic minus noncytosolic relative K/R penalty",
        "result": panel.summarize_parts(parts, draws, seed),
    }

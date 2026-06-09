"""Phase B: inventory deep mutational scans that map onto the census topology.

The gate asks whether at least four census membrane proteins have an abundance,
surface-expression or trafficking scan with enough charge changes in the aqueous
flanks to fit the Phase C model. Function readouts (currents, signalling, binding,
catalysis) are reported but flagged, because mixing them with expression readouts
would muddy the model.

Numbering is resolved against the target sequence that MaveDB stores with each
score set, not against an accession, so an isoform offset cannot silently shift
every residue. A score set is only counted when its target sequence can be
located exactly inside the census sequence.

Label-free with respect to ClinVar: no clinical classification is read here.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import re
import time
import urllib.error
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TOPOLOGY = ROOT / "data" / "raw" / "topology_viewer.html"
DEFAULT_CACHE = ROOT / "data" / "external" / "mavedb"
DEFAULT_OUT = ROOT / "results" / "phaseB_20260924"

API = "https://api.mavedb.org"

# Readout keywords. EXPRESSION_RE marks the assays the Phase C model can use;
# FUNCTION_RE marks readouts that must not be pooled with them.
EXPRESSION_RE = re.compile(
    r"abundance|surface expression|surface-expression|cell.surface|surface level|"
    r"surface presentation|trafficking|traffick|localization|localisation|"
    r"VAMP.?seq|protein level|expression level|steady.state|stability|"
    r"degradation|maturation|folding", re.I)
FUNCTION_RE = re.compile(
    r"current|conductance|electrophysiolog|patch.clamp|signalling|signaling|"
    r"cAMP|binding|transport activity|catalytic|activity|uptake|"
    r"growth|proliferat|fitness|luciferase", re.I)

MIN_VARIANTS = 150
# A score set is only used when this fraction of its reference residues match the census.
MIN_REFERENCE_AGREEMENT = 0.95

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
HGVS_PRO_RE = re.compile(r"^p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})$")
POSITIVE = {"K", "R"}
WINDOWS = (5, 10, 15)

# NCBI genetic code table 1, in the conventional TCAG x TCAG x TCAG codon order.
_BASES = "TCAG"
_AMINO_ACIDS = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
CODON_TABLE = {
    first + second + third: _AMINO_ACIDS[index]
    for index, (first, second, third) in enumerate(
        (f, s, t) for f in _BASES for s in _BASES for t in _BASES)
}


def translate(dna: str, frame: int) -> str:
    """Translate one reading frame, stopping at the first stop codon."""
    protein = []
    for start in range(frame, len(dna) - 2, 3):
        residue = CODON_TABLE.get(dna[start:start + 3])
        if residue is None or residue == "*":
            break
        protein.append(residue)
    return "".join(protein)


def fetch(path: str, cache: Path, body: dict | None = None, retries: int = 3) -> object:
    """GET or POST the MaveDB API, caching the raw response on disk."""
    cache.parent.mkdir(parents=True, exist_ok=True)
    if cache.exists():
        text = cache.read_text(encoding="utf-8")
    else:
        data = json.dumps(body).encode() if body is not None else None
        headers = {"Content-Type": "application/json"} if body is not None else {}
        request = urllib.request.Request(API + path, data=data, headers=headers)
        last: Exception | None = None
        for attempt in range(retries):
            try:
                with urllib.request.urlopen(request, timeout=180) as response:
                    text = response.read().decode("utf-8")
                break
            except (urllib.error.URLError, TimeoutError) as exc:
                last = exc
                time.sleep(2 * (attempt + 1))
        else:
            raise RuntimeError(f"MaveDB request failed for {path}: {last}")
        cache.write_text(text, encoding="utf-8")
    return text if cache.suffix == ".csv" else json.loads(text)


def read_topology(path: Path) -> dict[str, dict]:
    text = path.read_text(encoding="utf-8-sig")
    match = re.search(r"const PROTEINS_INLINE = (\[.*?\]);", text, flags=re.S)
    if not match:
        raise ValueError("Could not find the embedded Oxford census in the viewer HTML")
    proteins = {}
    for record in json.loads(match.group(1)):
        record["_flanks"] = flank_map(record)
        proteins[record["uniprot"]] = record
    return proteins


def flank_map(protein: dict) -> dict[int, tuple[str, int]]:
    """Every annotated aqueous-loop residue mapped to (side, distance to TMD edge)."""
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


def census_symbols(proteins: dict[str, dict]) -> dict[str, list[str]]:
    """Gene symbol (and synonym) -> census accessions. Single characters are dropped
    because one-letter synonyms collide with unrelated MaveDB target names."""
    index: dict[str, list[str]] = defaultdict(list)
    for accession, record in proteins.items():
        names = [record.get("gene", ""), *(record.get("synonyms") or [])]
        for name in names:
            name = str(name).strip().upper()
            if len(name) > 1:
                index[name].append(accession)
    return index


def readout_kind(score_set: dict) -> tuple[bool, bool, str]:
    experiment = score_set.get("experiment") or {}
    parts = [score_set.get("title"), score_set.get("shortDescription"),
             experiment.get("title"), experiment.get("shortDescription"),
             experiment.get("abstractText"), experiment.get("methodText")]
    # The score set's own title decides the readout; the experiment abstract often
    # mentions every assay in the paper, so it is used only as a fallback signal.
    own = " | ".join(str(p) for p in parts[:2] if p)
    full = " | ".join(str(p) for p in parts if p)
    expression = bool(EXPRESSION_RE.search(own)) or bool(EXPRESSION_RE.search(full))
    function = bool(FUNCTION_RE.search(own))
    label = ("expression" if bool(EXPRESSION_RE.search(own)) and not function
             else "expression_or_function" if expression else "function")
    return expression, function, label


def locate_in_census(target_sequence: str, proteins: dict[str, dict],
                     accessions: list[str]) -> tuple[str | None, int | None]:
    """Find the census protein containing this target sequence; return (accession, offset).

    offset is added to a MaveDB residue number to get the census residue number.
    """
    for accession in accessions:
        census = proteins[accession]["sequence"]
        if target_sequence == census:
            return accession, 0
        index = census.find(target_sequence)
        if index >= 0:
            return accession, index
    return None, None


def best_offset(variants: list[tuple[str, int, str]], sequence: str,
                located_offset: int | None) -> tuple[int, float]:
    """Pick the residue offset that best reconciles DMS reference residues.

    Returns (offset, fraction of variants whose reference residue then agrees).
    """
    candidates = {0}
    if located_offset is not None:
        candidates.update(located_offset + delta for delta in range(-3, 4))
    best, best_rate = 0, -1.0
    for offset in sorted(candidates):
        agree = sum(
            1 for ref, position, _ in variants
            if 1 <= position + offset <= len(sequence)
            and sequence[position + offset - 1] == ref
        )
        rate = agree / len(variants)
        if rate > best_rate:
            best, best_rate = offset, rate
    return best, best_rate


def parse_missense(hgvs: str) -> tuple[str, int, str] | None:
    match = HGVS_PRO_RE.fullmatch((hgvs or "").strip())
    if not match:
        return None
    ref3, position, alt3 = match.groups()
    ref, alt = AA3.get(ref3.upper()), AA3.get(alt3.upper())
    if ref is None or alt is None or ref == alt:
        return None
    return ref, int(position), alt


def charge_change(ref: str, alt: str) -> str | None:
    if ref in POSITIVE and alt not in POSITIVE:
        return "positive_charge_loss"
    if ref not in POSITIVE and alt in POSITIVE:
        return "positive_charge_gain"
    return None


def score_set_rows(urn: str, cache_dir: Path) -> list[dict]:
    text = fetch(f"/api/v1/score-sets/{urn}/scores",
                 cache_dir / f"{urn.replace(':', '_')}.scores.csv")
    return list(csv.DictReader(io.StringIO(text)))


def audit_score_set(score_set: dict, proteins: dict[str, dict],
                    symbols: dict[str, list[str]], cache_dir: Path) -> dict | None:
    targets = score_set.get("targetGenes") or []
    for target in targets:
        sequence_block = target.get("targetSequence") or {}
        sequence_type = (sequence_block.get("sequenceType") or "").lower()
        raw_sequence = (sequence_block.get("sequence") or "").strip().upper()
        if not raw_sequence:
            continue
        # Most MaveDB targets store the coding DNA, so translate it. The deposited
        # region is not always in frame 0, so every frame is tried and the one that
        # lands inside a census sequence wins.
        if sequence_type == "protein":
            attempts = [(raw_sequence, "protein")]
        elif sequence_type == "dna":
            attempts = [(translate(raw_sequence, frame), f"dna_frame{frame}")
                        for frame in (0, 1, 2)]
        else:
            continue

        name = str(target.get("name") or "").strip().upper()
        candidates = symbols.get(name, [])
        if not candidates:
            candidates = list(proteins)  # fall back to a sequence search over the census
        accession = offset = None
        for target_sequence, how in attempts:
            if len(target_sequence) < 20:
                continue
            accession, offset = locate_in_census(target_sequence, proteins, candidates)
            if accession is not None:
                resolved_via = how
                break
        if accession is None:
            return {"urn": score_set["urn"], "target_name": name,
                    "status": f"target_sequence_not_found_in_census ({sequence_type})"}

        record = proteins[accession]
        sequence = record["sequence"]

        variants = []
        for row in score_set_rows(score_set["urn"], cache_dir):
            parsed = parse_missense(row.get("hgvs_pro", ""))
            if parsed is None or row.get("score") in (None, "", "NA"):
                continue
            variants.append(parsed)
        if not variants:
            return {"urn": score_set["urn"], "status": "no_scored_simple_missense"}

        # Locating the target sequence gives a candidate offset, but some deposits
        # omit the start codon while still numbering hgvs_pro in full-protein
        # coordinates. Choose the offset that actually reconciles the reference
        # residues, and keep the agreement rate as the evidence that it is right.
        offset, agreement = best_offset(variants, sequence, offset)
        if agreement < MIN_REFERENCE_AGREEMENT:
            return {"urn": score_set["urn"], "target_name": name,
                    "status": f"reference_residues_disagree (best {agreement:.3f} "
                              f"at offset {offset})"}

        flanks = record["_flanks"]
        counts: dict[str, int] = defaultdict(int)
        scored = ref_mismatch = 0
        for ref, position, alt in variants:
            census_position = position + offset
            if census_position < 1 or census_position > len(sequence) or \
                    sequence[census_position - 1] != ref:
                ref_mismatch += 1
                continue
            scored += 1
            placement = flanks.get(census_position)
            if placement is None:
                continue
            side, distance = placement
            kind = charge_change(ref, alt)
            for window in WINDOWS:
                if distance <= window:
                    counts[f"flank_{window}"] += 1
                    if kind:
                        counts[f"flank_{window}_{side}_{kind}"] += 1
                        counts[f"flank_{window}_charge_change"] += 1

        expression, function, label = readout_kind(score_set)
        return {
            "urn": score_set["urn"], "status": "ok", "target_name": name,
            "uniprot": accession, "gene": record["gene"], "offset": offset,
            "resolved_via": resolved_via,
            "reference_agreement": round(agreement, 4),
            "num_tmds": record["numTMDs"], "protein_length": record["length"],
            "readout": label, "expression_readout": expression,
            "title": score_set.get("title") or "",
            "n_variants_reported": score_set.get("numVariants"),
            "n_missense_scored": scored, "n_reference_mismatch": ref_mismatch,
            **{key: counts[key] for key in sorted(counts)},
        }
    return {"urn": score_set["urn"], "status": "no_protein_target_sequence"}


def run(args: argparse.Namespace) -> dict:
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = Path(args.cache)
    proteins = read_topology(Path(args.topology))
    symbols = census_symbols(proteins)

    names = fetch("/api/v1/target-genes/names", cache_dir / "target_gene_names.json")
    matching = sorted({str(n) for n in names
                       if isinstance(n, str) and str(n).strip().upper() in symbols})
    print(f"MaveDB target names matching a census gene: {len(matching)}")

    score_sets: dict[str, dict] = {}
    for name in matching:
        result = fetch("/api/v1/score-sets/search",
                       cache_dir / f"search_{name.upper()}.json",
                       body={"targets": [name], "published": True})
        found = result.get("scoreSets", []) if isinstance(result, dict) else result
        for record in found:
            score_sets.setdefault(record["urn"], record)
    print(f"published score sets on those targets: {len(score_sets)}")

    audited, skipped = [], defaultdict(int)
    for urn, score_set in sorted(score_sets.items()):
        expression, _, _ = readout_kind(score_set)
        if not expression:
            skipped["no_expression_readout"] += 1
            continue
        if (score_set.get("numVariants") or 0) < MIN_VARIANTS:
            skipped["too_few_variants"] += 1
            continue
        row = audit_score_set(score_set, proteins, symbols, cache_dir)
        if row is None or row.get("status") != "ok":
            skipped[(row or {}).get("status", "unknown")] += 1
            print(f"  skip {urn}: {(row or {}).get('status')}")
            continue
        audited.append(row)
        print(f"  {row['gene']:<8} {row['urn']:<26} flank10={row.get('flank_10', 0):>5} "
              f"charge10={row.get('flank_10_charge_change', 0):>4}  {row['readout']}")

    fields = ["gene", "uniprot", "urn", "readout", "expression_readout", "resolved_via", "reference_agreement", "num_tmds",
              "protein_length", "offset", "n_variants_reported", "n_missense_scored",
              "n_reference_mismatch", "flank_5", "flank_10", "flank_15",
              "flank_5_charge_change", "flank_10_charge_change", "flank_15_charge_change",
              "flank_10_cytosolic_positive_charge_loss", "flank_10_cytosolic_positive_charge_gain",
              "flank_10_outer_positive_charge_loss", "flank_10_outer_positive_charge_gain",
              "title"]
    inventory_path = out_dir / "dms_inventory.csv"
    with inventory_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in sorted(audited, key=lambda r: -r.get("flank_10_charge_change", 0)):
            writer.writerow({key: row.get(key, 0) for key in fields})

    # Gate: proteins whose best expression score set has usable flank charge changes.
    best: dict[str, dict] = {}
    for row in audited:
        if row["readout"] != "expression":
            continue
        key = row["uniprot"]
        if row.get("flank_10_charge_change", 0) > best.get(key, {}).get("flank_10_charge_change", -1):
            best[key] = row
    qualifying = {k: v for k, v in best.items() if v.get("flank_10_charge_change", 0) >= 10}

    summary = {
        "analysis": "Phase B DMS inventory against the ER membrane census",
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "reads_clinvar": False,
        "expression_readouts_only": True,
        "min_variants_per_score_set": MIN_VARIANTS,
        "census_matching_target_names": len(matching),
        "published_score_sets_considered": len(score_sets),
        "score_sets_audited": len(audited),
        "skipped": dict(skipped),
        "proteins_with_expression_readout": len(best),
        "proteins_with_at_least_10_flank_charge_changes": len(qualifying),
        "gate_threshold_proteins": 4,
        "gate_pass": len(qualifying) >= 4,
        "qualifying_proteins": {
            row["gene"]: {
                "uniprot": row["uniprot"], "urn": row["urn"], "num_tmds": row["num_tmds"],
                "n_missense_scored": row["n_missense_scored"],
                "flank_10": row.get("flank_10", 0),
                "flank_10_charge_change": row.get("flank_10_charge_change", 0),
                "title": row["title"],
            }
            for row in sorted(qualifying.values(), key=lambda r: -r.get("flank_10_charge_change", 0))
        },
        "inventory_file": str(inventory_path),
    }
    summary_path = out_dir / "phaseB_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: summary[k] for k in (
        "score_sets_audited", "proteins_with_expression_readout",
        "proteins_with_at_least_10_flank_charge_changes", "gate_pass")}, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", default=str(DEFAULT_TOPOLOGY))
    parser.add_argument("--cache", default=str(DEFAULT_CACHE))
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    run(parser.parse_args())


if __name__ == "__main__":
    main()

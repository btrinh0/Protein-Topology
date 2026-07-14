"""Audit mutation-score contrasts near membrane boundaries.

For each charge and neutral-placebo rule, the audit reports:

1. A within-site side difference that removes additive position effects.
2. A descriptive, cross-fitted composition baseline.
3. A distance profile relative to the nearest membrane boundary.

Each rule has a declared direction; placebo rules expect zero. Differences against
a placebo use shared protein resamples to preserve correlation. The composition
calibration is descriptive and does not identify an internal model mechanism.

Usage
-----
    python -m membrane_rule_audit --topology census.html --scores scores.tsv \\
        --score-direction higher-is-damaging --out results/

See README.md for the score-file formats and for how to add a new predictor.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
POSITIVE, NEGATIVE, AMBIGUOUS = set("KR"), set("DE"), set("H")
CHARGED = POSITIVE | NEGATIVE | AMBIGUOUS
DEFAULT_BINS = ((1, 5), (6, 10), (11, 15), (16, 30), (31, 60))
PLACEBO_VI, PLACEBO_LM = set("VI"), set("LM")
PLACEBO_NQ, PLACEBO_Q = set("NQ"), set("Q")
RULES = {
    "positive_charge_loss": ("loss", POSITIVE),
    "positive_charge_gain": ("gain", POSITIVE),
    "negative_charge_loss": ("loss", NEGATIVE),
    "negative_charge_gain": ("gain", NEGATIVE),
    # Uncharged, conservative, and with no side asymmetry of their own: whatever
    # these pick up is a side baseline the charge rules also inherit. Two sets, so
    # a placebo reading cannot rest on one choice of residues.
    "neutral_placebo_gain_vi": ("gain", PLACEBO_VI),
    "neutral_placebo_gain_lm": ("gain", PLACEBO_LM),
    # Polar probes rather than placebos. If a side effect tracks polarity instead of
    # charge, these carry it and the hydrophobic pair above does not. N is the
    # N-glycosylation residue, so gaining it on the outer side can create a sequon;
    # Q is its non-glycosylatable isostere and isolates polarity on its own.
    "neutral_placebo_gain_nq": ("gain", PLACEBO_NQ),
    "neutral_placebo_gain_q": ("gain", PLACEBO_Q),
}
# What the positive-inside rule predicts for each cytosolic-minus-outer contrast.
# Losing a cytosolic K/R should hurt more inside; gaining a K/R should hurt more
# outside. None means this tool declares no direction, so no verdict is issued.
EXPECTED_SIGN = {
    "positive_charge_loss": +1,
    "positive_charge_gain": -1,
    "negative_charge_loss": None,
    "negative_charge_gain": None,
    "neutral_placebo_gain_vi": 0,
    "neutral_placebo_gain_lm": 0,
    "neutral_placebo_gain_nq": 0,
    "neutral_placebo_gain_q": 0,
}
CONTRASTS = {
    "positive_loss_minus_gain": ("positive_charge_loss", "positive_charge_gain"),
    "positive_gain_minus_placebo_vi": ("positive_charge_gain", "neutral_placebo_gain_vi"),
    "positive_gain_minus_placebo_lm": ("positive_charge_gain", "neutral_placebo_gain_lm"),
    "negative_loss_minus_gain": ("negative_charge_loss", "negative_charge_gain"),
    "negative_gain_minus_placebo_vi": ("negative_charge_gain", "neutral_placebo_gain_vi"),
    "positive_gain_minus_placebo_nq": ("positive_charge_gain", "neutral_placebo_gain_nq"),
    "positive_gain_minus_placebo_q": ("positive_charge_gain", "neutral_placebo_gain_q"),
    # Does a polar swap behave differently by side than a hydrophobic one? A large
    # value here means the side effect tracks polarity, not charge.
    "placebo_nq_minus_placebo_vi": ("neutral_placebo_gain_nq", "neutral_placebo_gain_vi"),
    "placebo_q_minus_placebo_vi": ("neutral_placebo_gain_q", "neutral_placebo_gain_vi"),
    # Two hydrophobic pairs should agree; this says whether they do.
    "placebo_lm_minus_placebo_vi": ("neutral_placebo_gain_lm", "neutral_placebo_gain_vi"),
}
# Every gain rule's group held out at once, so one comparator serves them all and a
# contrast between two gain rules differs only in its numerator. Narrower than the
# per-rule comparator, so it is reported beside the 14-residue definition, not instead.
COMMON_COMPARATOR = frozenset(
    set(AMINO_ACIDS) - CHARGED - PLACEBO_VI - PLACEBO_LM - PLACEBO_NQ - PLACEBO_Q)
VARIANT = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def load_topology(path: Path) -> dict:
    """Read a topology census: {accession: {sequence, tmds, loops, gene}}.

    Accepts the Oxford-MRC viewer HTML with an embedded PROTEINS_INLINE array, or
    a plain JSON file with the same record shape.
    """
    text = Path(path).read_text(encoding="utf-8-sig")
    if path.suffix.lower() in {".json"}:
        records = json.loads(text)
    else:
        match = re.search(r"const PROTEINS_INLINE = (\[.*?\]);", text, flags=re.S)
        if not match:
            raise ValueError("No embedded PROTEINS_INLINE array found in the topology file")
        records = json.loads(match.group(1))
    return {r["uniprot"]: r for r in records}


def flank_map(record: dict, max_distance: int) -> dict[int, tuple[int, int, int]]:
    """position -> (side, distance to nearest TMD edge, 1-based TMD index).

    Boundaries come from transmembrane coordinates rather than loop numbering,
    because censuses commonly skip loop indices where two segments abut.
    """
    ends = {t["end"]: i for i, t in enumerate(record["tmds"], start=1)}
    starts = {t["start"]: i for i, t in enumerate(record["tmds"], start=1)}
    out: dict[int, tuple[int, int, int]] = {}
    for loop in record["loops"].values():
        side = str(loop.get("location", "")).strip().lower()
        if side not in {"inside", "outside"}:
            continue
        side_index = 0 if side == "inside" else 1
        start, end = int(loop["start"]), int(loop["end"])
        borders = []
        if start - 1 in ends:
            borders.append((start - 1, ends[start - 1], -1))
        if end + 1 in starts:
            borders.append((end + 1, starts[end + 1], 1))
        for position in range(start, end + 1):
            best = None
            for boundary, tmd_index, direction in borders:
                distance = (position - boundary) if direction == -1 else (boundary - position)
                if 0 < distance <= max_distance and (best is None or distance < best[0]):
                    best = (distance, tmd_index)
            if best:
                out[position] = (side_index, best[0], best[1])
    return out


def load_scores(path: Path, columns: dict[str, str] | None = None) -> list[tuple]:
    """Read predictor scores as (accession, position, reference, alternate, score).

    Supported: .npz with arrays protein/position/ref_aa/alt_aa/score plus an
    `accessions` index; or a delimited text file, optionally gzipped. A text file
    needs an accession column, a score column, and either a `protein_variant`
    column like `R123W` or separate reference/position/alternate columns.
    """
    path = Path(path)
    if path.suffix == ".npz":
        blob = np.load(path, allow_pickle=False)
        accessions = list(blob["accessions"])
        alphabet = np.array(list("".join(blob["aa_alphabet"].tolist()))
                            if "aa_alphabet" in blob else list(AMINO_ACIDS))
        score_key = "score" if "score" in blob else (
            "am_pathogenicity" if "am_pathogenicity" in blob else None)
        if score_key is None:
            raise ValueError("No score array found in the .npz")
        return [(accessions[p], int(pos), str(alphabet[r]), str(alphabet[a]), float(s))
                for p, pos, r, a, s in zip(blob["protein"], blob["position"],
                                           blob["ref_aa"], blob["alt_aa"], blob[score_key])]

    columns = columns or {}
    opener = gzip.open if path.suffix == ".gz" else open
    rows = []
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        sample = handle.readline()
        handle.seek(0)
        delimiter = "\t" if sample.count("\t") >= sample.count(",") else ","
        reader = csv.DictReader(handle, delimiter=delimiter)
        accession_col = columns.get("accession", "uniprot_id")
        score_col = columns.get("score", "score")
        variant_col = columns.get("variant", "protein_variant")
        for row in reader:
            accession = (row.get(accession_col) or "").strip()
            try:
                score = float(row.get(score_col))
            except (TypeError, ValueError):
                continue
            if variant_col in row and row.get(variant_col):
                match = VARIANT.fullmatch((row[variant_col] or "").strip().upper())
                if not match:
                    continue
                reference, position, alternate = match.group(1), int(match.group(2)), match.group(3)
            else:
                try:
                    position = int(row[columns.get("position", "position")])
                except (TypeError, ValueError, KeyError):
                    continue
                reference = (row.get(columns.get("reference", "ref")) or "").strip().upper()
                alternate = (row.get(columns.get("alternate", "alt")) or "").strip().upper()
            if reference in AMINO_ACIDS and alternate in AMINO_ACIDS and reference != alternate:
                rows.append((accession, position, reference, alternate, score))
    return rows


def build_sites(scores, proteins: dict, max_distance: int, flip: bool) -> dict:
    """Group scores by residue, keeping only flank positions that match the census."""
    maps = {a: flank_map(r, max_distance) for a, r in proteins.items()}
    sites: dict = defaultdict(dict)
    context: dict = {}
    qc = defaultdict(int)
    for accession, position, reference, alternate, score in scores:
        record = proteins.get(accession)
        if record is None:
            qc["accession_not_in_census"] += 1
            continue
        sequence = record["sequence"]
        if position > len(sequence) or sequence[position - 1] != reference:
            qc["reference_mismatch"] += 1
            continue
        flank = maps[accession].get(position)
        if flank is None:
            continue
        key = (accession, position, reference)
        sites[key][alternate] = -score if flip else score
        context[key] = {"accession": accession, "gene": record.get("gene", accession),
                        "side": flank[0], "distance": flank[1], "tmd_index": flank[2]}
        qc["kept"] += 1
    return {"sites": sites, "context": context, "qc": dict(qc)}


def penalties(table: dict, kind: str, group: set[str],
              comparator: frozenset[str] | None = None) -> list[dict]:
    """Within-site penalty for a residue class, with charge reversals excluded from every arm.

    gain: reference uncharged, numerator introduces the class, comparator stays uncharged.
    loss: reference carries the class, numerator removes it, comparator swaps within the group.

    The group is held out of the opposing arm so the two arms stay disjoint. For a
    charged group that is already true, since the group sits outside the uncharged
    comparator; for an uncharged placebo group it is what keeps the arms from
    sharing terms and diluting the contrast toward zero.

    Holding the group out also shrinks the comparator, so two gain rules do not share
    one. Pass `comparator` to fix an explicit pool for every rule instead, which makes
    a contrast between two gain rules differ only in its numerator. Gain only: the loss
    estimand's comparator is a property-preserving swap inside the group.
    """
    if comparator is not None and kind != "gain":
        raise ValueError("an explicit comparator applies to gain rules only")
    out = []
    for key, alternates in table["sites"].items():
        reference = key[2]
        if kind == "gain":
            if reference in CHARGED:
                continue
            numerator = [v for a, v in alternates.items() if a in group]
            if comparator is None:
                denominator = [v for a, v in alternates.items()
                               if a not in CHARGED and a not in group]
            else:
                denominator = [v for a, v in alternates.items()
                               if a in comparator and a not in group]
        else:
            if reference not in group:
                continue
            numerator = [v for a, v in alternates.items()
                         if a not in CHARGED and a not in group]
            denominator = [v for a, v in alternates.items()
                           if a in group and a != reference]
        if not numerator or not denominator:
            continue
        info = table["context"][key]
        out.append({**info, "position": key[1], "reference": reference,
                    "penalty": float(np.mean(numerator)) - float(np.mean(denominator))})
    return out


def strata(rows: list[dict]) -> dict:
    """Per-(gene, reference residue) side contrasts, dropping single-sided cells."""
    buckets: dict = defaultdict(lambda: {"c": [], "o": []})
    for row in rows:
        buckets[(row["gene"], row["reference"])]["c" if row["side"] == 0 else "o"].append(
            row["penalty"])
    genes, weight, diff, used, dropped = [], [], [], 0, 0
    for (gene, _), arms in sorted(buckets.items()):
        n1, n0 = len(arms["c"]), len(arms["o"])
        if n1 == 0 or n0 == 0:
            dropped += n1 + n0
            continue
        genes.append(gene)
        weight.append(n1 * n0 / (n1 + n0))
        diff.append(float(np.mean(arms["c"]) - np.mean(arms["o"])))
        used += n1 + n0
    names = sorted(set(genes))
    index = {g: i for i, g in enumerate(names)}
    return {"cluster": np.array([index[g] for g in genes], dtype=int),
            "weight": np.array(weight), "diff": np.array(diff), "genes": names,
            "n_contributing_sites": used, "n_dropped_single_sided": dropped}


def point_and_interval(block: dict, draws: int, seed: int) -> dict:
    if block["diff"].size < 2:
        return {"estimate": None, "ci95": [None, None], "n_strata": int(block["diff"].size)}
    weight, diff = block["weight"], block["diff"]
    estimate = float((weight * diff).sum() / weight.sum())
    rng = np.random.default_rng(seed)
    n = len(block["genes"])
    samples = np.empty(draws)
    for i in range(draws):
        counts = np.bincount(rng.integers(0, n, n), minlength=n).astype(float)
        m = counts[block["cluster"]]
        total = (weight * m).sum()
        samples[i] = float((weight * m * diff).sum() / total) if total > 0 else np.nan
    finite = samples[np.isfinite(samples)]
    return {"estimate": estimate,
            "ci95": [float(np.percentile(finite, 2.5)), float(np.percentile(finite, 97.5))]
            if finite.size > 1 else [None, None],
            "n_strata": int(diff.size), "n_genes": n,
            "n_contributing_sites": block["n_contributing_sites"],
            "n_dropped_single_sided": block["n_dropped_single_sided"]}


def shared_sites(rows_by_rule: dict) -> dict:
    """Restrict every rule to the sites eligible for all of them.

    A singleton group has no numerator where its own residue is native, so that rule
    loses those sites while its neighbours keep them. Intersecting first makes every
    contrast a comparison on one population instead of across two.
    """
    identify = lambda r: (r["accession"], r["position"], r["reference"])
    common: set | None = None
    for rows in rows_by_rule.values():
        keys = {identify(r) for r in rows}
        common = keys if common is None else common & keys
    common = common or set()
    return {rule: [r for r in rows if identify(r) in common]
            for rule, rows in rows_by_rule.items()}


def verdict(estimate: float, ci95: list, expected, tol: float = 1e-12) -> str:
    """Read a rule's estimate against its declared prediction.

    `tol` stops accumulated floating-point error from deciding a verdict. A generator
    that plants no effect lands near 1e-17 with an interval that excludes zero
    numerically, which is not a finding. The tolerance sits many orders of magnitude
    below any substantive effect on any plausible score scale.
    """
    if expected is None:
        return "no declared prediction"
    decided = ci95[0] is not None and (ci95[0] > tol or ci95[1] < -tol)
    if expected == 0:
        return "placebo excludes zero" if decided else "placebo consistent with zero"
    if not decided:
        return "indeterminate: interval includes zero"
    return "consistent" if (estimate > 0) == (expected > 0) else "inconsistent"


def paired_contrasts(blocks: dict, contrasts: dict, draws: int, seed: int) -> dict:
    """Per-rule estimates and differences in differences from one shared protein draw.

    Every rule is resampled under the same draw over the union of its contributing
    proteins, so the two rules in a contrast keep their correlation. Resampling each
    rule on its own draw would treat them as independent and widen the contrast
    interval. A rule whose protein set is a strict subset of the union contributes a
    varying number of clusters per draw, which is the usual behaviour of a shared
    cluster bootstrap over a common population.
    """
    usable = {name: b for name, b in blocks.items() if b["diff"].size >= 2}
    if not usable:
        return {"rules": {}, "contrasts": {}, "note": "no rule had two or more strata"}
    universe = sorted(set().union(*(set(b["genes"]) for b in usable.values())))
    position = {gene: i for i, gene in enumerate(universe)}
    lift = {name: np.array([position[g] for g in b["genes"]], dtype=int)
            for name, b in usable.items()}
    rng = np.random.default_rng(seed)
    samples = {name: np.empty(draws) for name in usable}
    for i in range(draws):
        counts = np.bincount(rng.integers(0, len(universe), len(universe)),
                             minlength=len(universe)).astype(float)
        for name, block in usable.items():
            m = counts[lift[name]][block["cluster"]]
            total = (block["weight"] * m).sum()
            samples[name][i] = (float((block["weight"] * m * block["diff"]).sum() / total)
                                if total > 0 else np.nan)

    def interval(draw):
        finite = draw[np.isfinite(draw)]
        return ([float(np.percentile(finite, 2.5)), float(np.percentile(finite, 97.5))]
                if finite.size > 1 else [None, None])

    def point(block):
        return float((block["weight"] * block["diff"]).sum() / block["weight"].sum())

    rules = {}
    for name, block in usable.items():
        estimate, ci95 = point(block), interval(samples[name])
        rules[name] = {"estimate": estimate, "ci95": ci95,
                       "expected_sign": EXPECTED_SIGN.get(name),
                       "verdict": verdict(estimate, ci95, EXPECTED_SIGN.get(name)),
                       "n_strata": int(block["diff"].size),
                       "n_genes": len(block["genes"]),
                       "n_contributing_sites": block["n_contributing_sites"]}
    pairs = {}
    for name, (left, right) in contrasts.items():
        if left not in usable or right not in usable:
            continue
        ci95 = interval(samples[left] - samples[right])
        pairs[name] = {"estimate": point(usable[left]) - point(usable[right]),
                       "ci95": ci95, "minuend": left, "subtrahend": right,
                       "excludes_zero": bool(ci95[0] is not None
                                             and (ci95[0] > 0 or ci95[1] < 0))}
    return {"draws": draws, "shared_universe_genes": len(universe),
            "rules": rules, "contrasts": pairs}


def composition_baseline(table: dict, bins, folds: int, seed: int) -> dict:
    """Frequency baseline B1: log f(reference | side, bin) - log f(alternate | side, bin).

    Fitted on training families only and scored out of fold. Families are the gene
    symbols present in the census; pass a richer clustering if you have one.
    """
    keys = sorted(table["sites"])
    families = sorted({table["context"][k]["gene"] for k in keys})
    rng = np.random.default_rng(seed)
    fold_of = {f: int(rng.integers(0, folds)) for f in families}

    def bin_of(distance: int) -> int:
        for i, (low, high) in enumerate(bins):
            if low <= distance <= high:
                return i
        return -1

    index = {a: i for i, a in enumerate(AMINO_ACIDS)}
    scored: dict = defaultdict(dict)
    for held in range(folds):
        counts = np.full((2, len(bins), 20), 0.5)
        for key in keys:
            info = table["context"][key]
            if fold_of[info["gene"]] == held:
                continue
            bucket = bin_of(info["distance"])
            if bucket >= 0:
                counts[info["side"], bucket, index[key[2]]] += 1.0
        logf = np.log(counts / counts.sum(axis=2, keepdims=True))
        for key in keys:
            info = table["context"][key]
            if fold_of[info["gene"]] != held:
                continue
            bucket = bin_of(info["distance"])
            if bucket < 0:
                continue
            for alternate in table["sites"][key]:
                scored[key][alternate] = float(
                    logf[info["side"], bucket, index[key[2]]]
                    - logf[info["side"], bucket, index[alternate]])
    return {"sites": scored, "context": table["context"], "qc": {}}


def calibrated_residual(predictor_rows, baseline_rows, draws: int, seed: int) -> dict:
    """How much of the side difference the predictor's own use of composition implies.

    Regress predictor penalties on baseline penalties over all shared sites, then
    compare the observed side difference with slope x the baseline's. The residual
    is side-specific over- or under-use of the footprint. Predictor noise does not
    bias this, and the baseline being near-deterministic is an asset in a regressor.
    """
    key = lambda r: (r["accession"], r["position"], r["reference"])
    left = {key(r): r for r in predictor_rows}
    right = {key(r): r for r in baseline_rows}
    shared = sorted(set(left) & set(right))
    if len(shared) < 50:
        return {"note": "too few shared sites", "n_shared": len(shared)}
    y = np.array([left[k]["penalty"] for k in shared])
    x = np.array([right[k]["penalty"] for k in shared])
    design = np.column_stack([np.ones(len(x)), x])
    coefficients, *_ = np.linalg.lstsq(design, y, rcond=None)
    slope = float(coefficients[1])

    predictor_block = strata([left[k] for k in shared])
    baseline_block = strata([right[k] for k in shared])
    if predictor_block["diff"].size != baseline_block["diff"].size:
        return {"note": "stratum sets differ", "slope": slope}
    residual = predictor_block["diff"] - slope * baseline_block["diff"]
    weight = predictor_block["weight"]
    rng = np.random.default_rng(seed)
    n = len(predictor_block["genes"])
    samples = np.empty(draws)
    for i in range(draws):
        counts = np.bincount(rng.integers(0, n, n), minlength=n).astype(float)
        m = counts[predictor_block["cluster"]]
        total = (weight * m).sum()
        samples[i] = float((weight * m * residual).sum() / total) if total > 0 else np.nan
    finite = samples[np.isfinite(samples)]
    return {
        "slope_on_composition": slope,
        "observed_side_difference": float((weight * predictor_block["diff"]).sum() / weight.sum()),
        "implied_by_calibration": float(slope * (weight * baseline_block["diff"]).sum() / weight.sum()),
        "residual": float((weight * residual).sum() / weight.sum()),
        "residual_ci95": [float(np.percentile(finite, 2.5)), float(np.percentile(finite, 97.5))]
        if finite.size > 1 else [None, None],
        "n_shared_sites": len(shared), "n_strata": int(residual.size), "n_genes": n,
    }


def run_audit(topology: Path, scores_path: Path, flip: bool, out_dir: Path,
              max_distance: int = 60, bins=DEFAULT_BINS, folds: int = 5,
              draws: int = 1000, seed: int = 0, columns=None) -> dict:
    proteins = load_topology(Path(topology))
    raw = load_scores(Path(scores_path), columns)
    table = build_sites(raw, proteins, max_distance, flip)
    baseline = composition_baseline(table, bins, folds, seed)

    results: dict = {}
    rows_by_rule: dict = {}
    for rule, (kind, group) in RULES.items():
        rows = penalties(table, kind, group)
        rows_by_rule[rule] = rows
        base_rows = penalties(baseline, kind, group)
        near = [r for r in rows if r["distance"] <= bins[1][1]]
        results[rule] = {
            "side_difference_all": point_and_interval(strata(rows), draws, seed),
            "side_difference_near": point_and_interval(strata(near), draws, seed),
            "by_distance_bin": {
                f"{low}-{high}": point_and_interval(
                    strata([r for r in rows if low <= r["distance"] <= high]), draws, seed)
                for low, high in bins
            },
            "cytosolic_mean_penalty_by_bin": {
                f"{low}-{high}": (float(np.mean([r["penalty"] for r in rows
                                                 if r["side"] == 0 and low <= r["distance"] <= high]))
                                  if any(r["side"] == 0 and low <= r["distance"] <= high for r in rows)
                                  else None)
                for low, high in bins
            },
            "calibrated_vs_composition": calibrated_residual(rows, base_rows, draws, seed),
        }

    near_limit = bins[1][1]
    gain_rules = {r: g for r, (k, g) in RULES.items() if k == "gain"}
    gain_contrasts = {n: pair for n, pair in CONTRASTS.items()
                      if pair[0] in gain_rules and pair[1] in gain_rules}
    matched_rows = shared_sites({r: penalties(table, "gain", g, COMMON_COMPARATOR)
                                 for r, g in gain_rules.items()})
    matched_gain_audit = {
        "comparator": "".join(sorted(COMMON_COMPARATOR)),
        "note": "Every gain rule on one comparator and one site population, so a "
                "contrast differs only in its numerator. The comparator is narrower "
                "than the per-rule one, so these estimates are not the headline "
                "definition; read them as the matched contrast.",
        "all": paired_contrasts({r: strata(rw) for r, rw in matched_rows.items()},
                                gain_contrasts, draws, seed),
        f"within_{near_limit}_residues": paired_contrasts(
            {r: strata([x for x in rw if x["distance"] <= near_limit])
             for r, rw in matched_rows.items()}, gain_contrasts, draws, seed),
    }
    rule_audit = {
        "all": paired_contrasts({r: strata(rw) for r, rw in rows_by_rule.items()},
                                CONTRASTS, draws, seed),
        f"within_{near_limit}_residues": paired_contrasts(
            {r: strata([x for x in rw if x["distance"] <= near_limit])
             for r, rw in rows_by_rule.items()}, CONTRASTS, draws, seed),
    }

    summary = {
        "tool": "membrane_rule_audit",
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "topology": str(topology), "scores": str(scores_path),
        "score_direction": "higher is more damaging"
        + (" (input was flipped)" if flip else " (input taken as given)"),
        "max_distance": max_distance, "bins": [list(b) for b in bins],
        "folds": folds, "bootstrap_draws": draws, "seed": seed,
        "input_quality_control": table["qc"],
        "sign_convention": "Side difference is cytosolic minus outer. Positive means the "
                           "predictor calls the change more damaging on the cytosolic side.",
        "declared_predictions": {k: v for k, v in EXPECTED_SIGN.items()},
        "results": results,
        "rule_audit": rule_audit,
        "matched_gain_audit": matched_gain_audit,
    }
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "membrane_rule_audit.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--topology", required=True, help="census HTML or JSON")
    parser.add_argument("--scores", required=True, help=".npz, .tsv, .csv, optionally .gz")
    parser.add_argument("--score-direction", required=True,
                        choices=("higher-is-damaging", "higher-is-tolerated"),
                        help="declare it; the tool will not guess")
    parser.add_argument("--out", required=True)
    parser.add_argument("--max-distance", type=int, default=60)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--draws", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--accession-column", default="uniprot_id")
    parser.add_argument("--score-column", default="score")
    parser.add_argument("--variant-column", default="protein_variant")
    args = parser.parse_args()

    summary = run_audit(
        Path(args.topology), Path(args.scores),
        flip=(args.score_direction == "higher-is-tolerated"),
        out_dir=Path(args.out), max_distance=args.max_distance,
        folds=args.folds, draws=args.draws, seed=args.seed,
        columns={"accession": args.accession_column, "score": args.score_column,
                 "variant": args.variant_column})

    print(f"sites kept: {summary['input_quality_control'].get('kept', 0):,}")
    print(f"{'rule':<24}{'side diff':>12}{'95% CI':>22}{'slope':>9}{'residual':>11}  genes")
    for rule, block in summary["results"].items():
        side = block["side_difference_all"]
        cal = block["calibrated_vs_composition"]
        if side["estimate"] is None:
            print(f"{rule:<24}{'n/a':>12}")
            continue
        ci = f"[{side['ci95'][0]:+.3f},{side['ci95'][1]:+.3f}]"
        slope = cal.get("slope_on_composition")
        residual = cal.get("residual")
        slope_text = f"{slope:+9.3f}" if slope is not None else f"{'n/a':>9}"
        residual_text = f"{residual:+11.4f}" if residual is not None else f"{'n/a':>11}"
        print(f"{rule:<24}{side['estimate']:+12.4f}{ci:>22}{slope_text}{residual_text}"
              f"  {side.get('n_genes', 0)}")

    blocks = [("audit", scope, block)
              for scope, block in summary["rule_audit"].items()]
    matched = summary["matched_gain_audit"]
    blocks += [(f"matched-comparator audit [{matched['comparator']}]", scope, block)
               for scope, block in matched.items() if isinstance(block, dict)]
    for label, scope, block in blocks:
        if not block.get("rules"):
            continue
        print()
        print(f"{label}, {scope.replace('_', ' ')} "
              f"({block['draws']} shared draws over {block['shared_universe_genes']} proteins)")
        print(f"  {'rule':<24}{'estimate':>11}{'95% CI':>22}{'want':>6}  verdict")
        for rule, row in block["rules"].items():
            ci = f"[{row['ci95'][0]:+.4f},{row['ci95'][1]:+.4f}]"
            want = {1: "+", -1: "-", 0: "0", None: "."}[row["expected_sign"]]
            print(f"  {rule:<24}{row['estimate']:+11.4f}{ci:>22}{want:>6}  {row['verdict']}")
        if block["contrasts"]:
            print(f"  {'contrast':<24}{'estimate':>11}{'95% CI':>22}{'':>6}  excludes zero")
            for name, row in block["contrasts"].items():
                ci = f"[{row['ci95'][0]:+.4f},{row['ci95'][1]:+.4f}]"
                print(f"  {name:<24}{row['estimate']:+11.4f}{ci:>22}{'':>6}  "
                      f"{'yes' if row['excludes_zero'] else 'no'}")


if __name__ == "__main__":
    main()

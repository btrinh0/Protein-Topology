"""Step 2: does a predictor follow the compositional footprint, fall short, or go beyond it?

Frozen design: docs/STEP2_COMPOSITION_CONTROL_PREREG.md
SHA-256 06cdc2e8449ab03ced939b6cf1599d7d754a5fcec36c4a2ca992c5e9d1de5385

Composition by side and distance is the evolutionary footprint of the topology
rules, not a rival explanation to them. So the baselines here are not controls to
be subtracted; they are alternative "predictors" scored through the same
within-site estimator, and the question is whether a real predictor reproduces
them, falls short, or exceeds them.

Label-free: no clinical label and no deep-mutational-scanning effect is read.

Within-site penalties use the audit's residue definitions (run_gate0_within_site),
so the predictor numbers here reproduce the ones disclosed in section 2 of the
frozen design rather than a differently-defined variant of them.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from run_gate0 import load_frame
from run_gate0_within_site import site_gain_penalties, site_penalties
from run_gate1 import read_topology

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "results" / "step2_20260929"
ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
BINS = [(1, 5), (6, 10), (11, 15), (16, 30), (31, 60)]
KMER, SIGNATURE, BANDS = 5, 64, 16
PSEUDOCOUNT = 0.5
MINHASH_HASH_METHOD = "blake2b-64 ASCII k-mers, masked to 63 bits (v1)"


def minhash_clusters(proteins: dict, seed: int) -> dict[str, int]:
    """Sequence-similarity clusters via MinHash signatures and banded LSH.

    Keeps related proteins together so a family cannot straddle the train/test
    split. All-against-all alignment over 4,863 proteins is not affordable; this
    is the standard cheap substitute.
    """
    rng = np.random.default_rng(seed)
    coefficients = rng.integers(1, 2 ** 61 - 1, size=(SIGNATURE, 2))
    prime = (1 << 61) - 1
    accessions = sorted(proteins)
    signatures = np.zeros((len(accessions), SIGNATURE), dtype=np.int64)

    for row, accession in enumerate(accessions):
        sequence = proteins[accession]["sequence"]
        kmers = {sequence[i:i + KMER] for i in range(len(sequence) - KMER + 1)}
        if not kmers:
            signatures[row] = np.arange(SIGNATURE)
            continue
        # Python's hash(str) changes between processes, independently of the
        # analysis seed. A stable digest makes cross-validation families repeat.
        hashed = np.array([
            int.from_bytes(hashlib.blake2b(k.encode("ascii"), digest_size=8,
                                         person=b"p6-kmers-v1").digest(), "little")
            & 0x7FFFFFFFFFFFFFFF for k in kmers
        ], dtype=np.int64)
        for slot in range(SIGNATURE):
            a, b = coefficients[slot]
            signatures[row, slot] = int(((a * hashed + b) % prime).min())

    parent = list(range(len(accessions)))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x: int, y: int) -> None:
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[max(rx, ry)] = min(rx, ry)

    per_band = SIGNATURE // BANDS
    for band in range(BANDS):
        buckets: dict = defaultdict(list)
        block = signatures[:, band * per_band:(band + 1) * per_band]
        for row in range(len(accessions)):
            buckets[tuple(block[row].tolist())].append(row)
        for members in buckets.values():
            for other in members[1:]:
                union(members[0], other)
    return {accessions[i]: find(i) for i in range(len(accessions))}


def symbol_families(proteins: dict) -> dict[str, str]:
    """Alternative split: gene-symbol family root, e.g. SLC6A1 and SLC6A2 together."""
    import re
    out = {}
    for accession, record in proteins.items():
        gene = str(record.get("gene") or accession).upper()
        match = re.match(r"^([A-Z]+\d*[A-Z]*?)(?:\d+)?$", gene)
        root = re.sub(r"\d+$", "", gene) or gene
        out[accession] = root
    return out


def fit_b1(frame: pd.DataFrame, train_mask: np.ndarray) -> np.ndarray:
    """log f(ref | side, bin) - log f(alt | side, bin), frequencies from training rows only.

    Frequencies are over *reference residues*, which is the composition of the
    proteome at that side and distance - the footprint itself.
    """
    train = frame[train_mask]
    counts = np.full((2, len(BINS), 20), PSEUDOCOUNT)
    bin_of_train = train["bin"].to_numpy()
    side_of_train = train["side"].to_numpy()
    ref_of_train = train["ref_aa"].to_numpy()
    # One count per residue position, not per substitution, so the alphabet is
    # not weighted by how many alternates each codon reaches.
    seen = set()
    for side, bucket, ref, protein, position in zip(
            side_of_train, bin_of_train, ref_of_train,
            train["protein"].to_numpy(), train["position"].to_numpy()):
        key = (int(protein), int(position))
        if key in seen:
            continue
        seen.add(key)
        counts[int(side), int(bucket), int(ref)] += 1.0
    frequencies = counts / counts.sum(axis=2, keepdims=True)
    log_frequencies = np.log(frequencies)
    side = frame["side"].to_numpy().astype(int)
    bucket = frame["bin"].to_numpy().astype(int)
    return (log_frequencies[side, bucket, frame["ref_aa"].to_numpy()]
            - log_frequencies[side, bucket, frame["alt_aa"].to_numpy()])


def fit_b0(frame: pd.DataFrame, train_mask: np.ndarray) -> np.ndarray:
    """Substitution identity only: the same score wherever that pair occurs."""
    train = frame[train_mask]
    table = np.zeros((20, 20))
    counts = np.zeros((20, 20))
    b1_train = fit_b1(frame, train_mask)[train_mask]
    np.add.at(table, (train["ref_aa"].to_numpy(), train["alt_aa"].to_numpy()), b1_train)
    np.add.at(counts, (train["ref_aa"].to_numpy(), train["alt_aa"].to_numpy()), 1.0)
    mean = np.divide(table, counts, out=np.zeros_like(table), where=counts > 0)
    return mean[frame["ref_aa"].to_numpy(), frame["alt_aa"].to_numpy()]


def b2_features(frame: pd.DataFrame, proteins: dict, accessions: list[str]) -> np.ndarray:
    """Protein-specific topological context: charge margin, first-segment flag, loop length."""
    margin, first, length = {}, {}, {}
    for index, accession in enumerate(accessions):
        record = proteins.get(accession)
        if record is None:
            continue
        sequence = record["sequence"]
        for loop in record["loops"].values():
            side = str(loop.get("location", "")).lower()
            if side not in {"inside", "outside"}:
                continue
            start, end = int(loop["start"]), int(loop["end"])
            window = sequence[max(0, start - 1):end]
            charge = sum(1 for a in window if a in "KR") - sum(1 for a in window if a in "DE")
            signed = charge if side == "inside" else -charge
            for position in range(start, end + 1):
                margin[(index, position)] = signed
                length[(index, position)] = np.log1p(end - start + 1)
    keys = list(zip(frame["protein"].to_numpy(), frame["position"].to_numpy()))
    return np.column_stack([
        np.array([margin.get(k, 0.0) for k in keys], dtype=float),
        (frame["tmd_index"].to_numpy() == 1).astype(float),
        np.array([length.get(k, 0.0) for k in keys], dtype=float),
    ])


def fit_b2(frame: pd.DataFrame, train_mask: np.ndarray, features: np.ndarray) -> np.ndarray:
    """B1 conditioned on protein-specific context instead of pooled across proteins.

    The frozen design says B2 is "B1 plus protein-specific context ... entered as a
    linear adjustment fitted on training families only", but never says what that
    adjustment is fitted *to*. There is no target that is not circular: fitting it
    to a predictor would make B2 a model of that predictor rather than of
    composition. Deviation logged. B2 is instead built the same way as B1 with the
    context entering as extra conditioning variables, which is what "goes beyond
    averages" has to mean for a frequency model.

    Conditioning adds the sign of the segment's inside-minus-outside charge margin
    and whether the site sits on the first transmembrane segment.
    """
    margin_sign = np.sign(features[:, 0]).astype(int) + 1
    first_segment = features[:, 1].astype(int)
    context = margin_sign * 2 + first_segment
    counts = np.full((2, len(BINS), 6, 20), PSEUDOCOUNT)
    train = frame[train_mask]
    seen = set()
    for side, bucket, ctx, ref, protein, position in zip(
            train["side"].to_numpy(), train["bin"].to_numpy(), context[train_mask],
            train["ref_aa"].to_numpy(), train["protein"].to_numpy(),
            train["position"].to_numpy()):
        key = (int(protein), int(position))
        if key in seen:
            continue
        seen.add(key)
        counts[int(side), int(bucket), int(ctx), int(ref)] += 1.0
    log_frequencies = np.log(counts / counts.sum(axis=3, keepdims=True))
    side = frame["side"].to_numpy().astype(int)
    bucket = frame["bin"].to_numpy().astype(int)
    return (log_frequencies[side, bucket, context, frame["ref_aa"].to_numpy()]
            - log_frequencies[side, bucket, context, frame["alt_aa"].to_numpy()])


def side_difference(sites: pd.DataFrame) -> tuple[float, int]:
    if sites.empty:
        return float("nan"), 0
    buckets: dict = defaultdict(lambda: {"c": [], "o": []})
    for protein, substitution, cyto, penalty in zip(
            sites["protein"].to_numpy(), sites["substitution"].to_numpy(),
            sites["is_cytosolic"].to_numpy(), sites["penalty"].to_numpy()):
        buckets[(int(protein), str(substitution))]["c" if cyto else "o"].append(float(penalty))
    numerator = denominator = 0.0
    used = 0
    for arms in buckets.values():
        n1, n0 = len(arms["c"]), len(arms["o"])
        if n1 == 0 or n0 == 0:
            continue
        weight = n1 * n0 / (n1 + n0)
        numerator += weight * (np.mean(arms["c"]) - np.mean(arms["o"]))
        denominator += weight
        used += n1 + n0
    return (numerator / denominator if denominator > 0 else float("nan")), used


def rank_side(sites: pd.DataFrame) -> float | None:
    """P(cytosolic penalty > outer penalty) within stratum, ties at 0.5.

    Scale-free, so a near-deterministic baseline and a noisy predictor are on the
    same footing. The standardised difference is not comparable between them:
    dividing by each source's own spread rewards a baseline that has almost no
    within-site variance, which is exactly what B1 is.
    """
    if sites.empty:
        return None
    wins = pairs = 0.0
    for _, group in sites.groupby(["protein", "substitution"], sort=False):
        cyto = group.loc[group["is_cytosolic"], "penalty"].to_numpy()
        outer = group.loc[~group["is_cytosolic"], "penalty"].to_numpy()
        if cyto.size == 0 or outer.size == 0:
            continue
        comparison = cyto[:, None] - outer[None, :]
        wins += float((comparison > 0).sum() + 0.5 * (comparison == 0).sum())
        pairs += float(comparison.size)
    return wins / pairs if pairs else None


def profile(frame: pd.DataFrame, column: str, kind: str, group: set[str]) -> dict:
    """Standardised side difference by distance bin for one score column."""
    work = frame.rename(columns={column: "score"})
    sites = (site_gain_penalties(work, group) if kind == "gain"
             else site_penalties(work, group))
    if sites.empty:
        return {}
    spread = float(sites["penalty"].std(ddof=1)) or 1.0
    out = {}
    for low, high in BINS:
        subset = sites[(sites["distance"] >= low) & (sites["distance"] <= high)]
        value, used = side_difference(subset)
        out[f"{low}-{high}"] = {"standardised": value / spread if np.isfinite(value) else None,
                                "rank": rank_side(subset), "n_sites": used}
    subset = sites[sites["distance"] <= 10]
    value, used = side_difference(subset)
    out["window_10"] = {"standardised": value / spread if np.isfinite(value) else None,
                        "rank": rank_side(subset), "n_sites": used}
    return out


def paired_gap(frame: pd.DataFrame, left: str, right: str, kind: str, group: set[str],
               draws: int, seed: int) -> dict:
    """P2: left minus right on identical sites, one shared protein-clustered draw."""
    def stats(column):
        work = frame.rename(columns={column: "score"})
        sites = (site_gain_penalties(work, group) if kind == "gain"
                 else site_penalties(work, group))
        sites = sites[sites["distance"] <= 10]
        spread = float(sites["penalty"].std(ddof=1)) or 1.0
        buckets: dict = defaultdict(lambda: {"c": [], "o": []})
        for protein, substitution, cyto, penalty in zip(
                sites["protein"].to_numpy(), sites["substitution"].to_numpy(),
                sites["is_cytosolic"].to_numpy(), sites["penalty"].to_numpy()):
            buckets[(int(protein), str(substitution))]["c" if cyto else "o"].append(
                float(penalty) / spread)
        keys, weights, diffs, proteins_ = [], [], [], []
        for key, arms in sorted(buckets.items()):
            n1, n0 = len(arms["c"]), len(arms["o"])
            if n1 == 0 or n0 == 0:
                continue
            keys.append(key); weights.append(n1 * n0 / (n1 + n0))
            diffs.append(float(np.mean(arms["c"]) - np.mean(arms["o"])))
            proteins_.append(key[0])
        return keys, np.array(weights), np.array(diffs), np.array(proteins_)

    lk, lw, ld, lp = stats(left)
    rk, rw, rd, _ = stats(right)
    shared = sorted(set(lk) & set(rk))
    if len(shared) < 3:
        return {"n_shared": len(shared)}
    li = {k: i for i, k in enumerate(lk)}; ri = {k: i for i, k in enumerate(rk)}
    idx_l = np.array([li[k] for k in shared]); idx_r = np.array([ri[k] for k in shared])
    weight = lw[idx_l]
    gap = ld[idx_l] - rd[idx_r]
    names = sorted({k[0] for k in shared})
    cluster = np.array([names.index(k[0]) for k in shared])
    point = float((weight * gap).sum() / weight.sum())
    rng = np.random.default_rng(seed)
    samples = np.empty(draws)
    for i in range(draws):
        counts = np.bincount(rng.integers(0, len(names), len(names)),
                             minlength=len(names)).astype(float)
        m = counts[cluster]
        total = (weight * m).sum()
        samples[i] = float((weight * m * gap).sum() / total) if total > 0 else np.nan
    finite = samples[np.isfinite(samples)]
    return {
        "left": left, "right": right, "n_shared_strata": len(shared),
        "n_proteins": len(names),
        "left_standardised": float((weight * ld[idx_l]).sum() / weight.sum()),
        "right_standardised": float((weight * rd[idx_r]).sum() / weight.sum()),
        "gap": point,
        "ci95": [float(np.percentile(finite, 2.5)), float(np.percentile(finite, 97.5))]
        if finite.size > 1 else [None, None],
    }


def run(args) -> dict:
    proteins = read_topology(Path(args.topology))
    frame, _ = load_frame(Path(args.scores))
    blob = np.load(args.scores, allow_pickle=False)
    accessions = list(blob["accessions"])

    bin_index = np.full(len(frame), -1, dtype=int)
    distance = frame["distance"].to_numpy()
    for i, (low, high) in enumerate(BINS):
        bin_index[(distance >= low) & (distance <= high)] = i
    frame = frame[bin_index >= 0].copy()
    frame["bin"] = bin_index[bin_index >= 0]

    clusters = (minhash_clusters(proteins, args.seed) if args.split == "sequence"
                else symbol_families(proteins))
    cluster_of = np.array([clusters.get(accessions[p], -1) for p in frame["protein"].to_numpy()])
    unique = sorted(set(cluster_of.tolist()))
    rng = np.random.default_rng(args.seed)
    fold_of_cluster = {c: int(rng.integers(0, args.folds)) for c in unique}
    fold = np.array([fold_of_cluster[c] for c in cluster_of])

    features = b2_features(frame, proteins, accessions)
    frame["b0"] = np.nan; frame["b1"] = np.nan; frame["b2"] = np.nan
    for f in range(args.folds):
        train = fold != f
        held = ~train
        if train.sum() == 0 or held.sum() == 0:
            continue
        b1 = fit_b1(frame, train)
        b0 = fit_b0(frame, train)
        b2 = fit_b2(frame, train, features)
        frame.loc[held, "b1"] = b1[held]
        frame.loc[held, "b0"] = b0[held]
        frame.loc[held, "b2"] = b2[held]
    frame = frame.dropna(subset=["b0", "b1", "b2"])

    esm = None
    if Path(args.esm).exists():
        esm_frame, _ = load_frame(Path(args.esm))
        key = ["protein", "position", "alt_aa"]
        esm = esm_frame[key + ["score"]].rename(columns={"score": "esm1b"})
        frame = frame.merge(esm, on=key, how="left")
    frame = frame.rename(columns={"score": "alphamissense"})

    sources = ["b0", "b1", "b2", "alphamissense"] + (["esm1b"] if esm is not None else [])
    rules = {"positive_charge_gain": ("gain", set("KR")),
             "positive_charge_loss": ("loss", set("KR")),
             "negative_charge_gain": ("gain", set("DE")),
             "negative_charge_loss": ("loss", set("DE"))}

    profiles: dict = {}
    for rule, (kind, group) in rules.items():
        profiles[rule] = {}
        for source in sources:
            work = frame.dropna(subset=[source])
            profiles[rule][source] = profile(work, source, kind, group)

    gaps: dict = {}
    for rule, (kind, group) in rules.items():
        gaps[rule] = {}
        for predictor in ("alphamissense", "esm1b"):
            if predictor not in sources:
                continue
            work = frame.dropna(subset=[predictor, "b1"])
            gaps[rule][f"{predictor}_minus_b1"] = paired_gap(
                work, predictor, "b1", kind, group, args.boot_draws, args.seed)

    summary = {
        "analysis": "Step 2 composition control: footprint reproduction on held-out families",
        "prereg_sha256": "06cdc2e8449ab03ced939b6cf1599d7d754a5fcec36c4a2ca992c5e9d1de5385",
        "clinical_labels_used": False, "dms_effects_used": False,
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "split": args.split, "folds": args.folds,
        "n_clusters": len(unique),
        "n_rows_scored_held_out": int(len(frame)),
        "profiles_by_rule": profiles,
        "primary_gaps_window_10": gaps,
    }
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    (out / f"step2_{args.split}.json").write_text(json.dumps(summary, indent=2) + "\n",
                                                  encoding="utf-8")

    print(f"split={args.split}  clusters={len(unique)}  held-out rows={len(frame):,}")
    for rule in rules:
        print(f"\n=== {rule} : standardised side difference, window <=10 ===")
        for source in sources:
            value = profiles[rule][source].get("window_10", {}).get("standardised")
            n = profiles[rule][source].get("window_10", {}).get("n_sites", 0)
            print(f"   {source:<14} {value:+.4f}  n={n}" if value is not None
                  else f"   {source:<14}    n/a")
        for name, g in gaps[rule].items():
            if "gap" in g:
                print(f"   * {name:<22} {g['gap']:+.4f}  CI [{g['ci95'][0]:+.4f},"
                      f"{g['ci95'][1]:+.4f}]  strata={g['n_shared_strata']}")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", default=str(ROOT / "data/raw/topology_viewer.html"))
    parser.add_argument("--scores", default=str(ROOT / "data/processed/gate0/am_full_flank_scores.npz"))
    parser.add_argument("--esm", default=str(ROOT / "data/processed/gate0/esm1b_flank_scores.npz"))
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    parser.add_argument("--split", choices=("sequence", "symbol"), default="sequence")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--boot-draws", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260929)
    run(parser.parse_args())


if __name__ == "__main__":
    main()

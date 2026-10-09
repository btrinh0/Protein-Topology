"""Compare each predictor with the compositional footprint after calibration.

The preregistered self-normalized metric is unsuitable here because B1 is a
deterministic function of side, distance bin, reference residue, and alternate
residue. Within a stratum it has no measurement noise, so comparisons with noisy
predictors primarily reflect noise rather than footprint use. The deviation is
recorded in the decision log.

Primary analysis: fit the predictor's within-site penalties on B0 and B1 using
training families. On held-out families, compare the observed side difference
with beta times the B1 side difference. The residual measures side-specific over-
or under-use of the footprint in the predictor's own units.

Secondary analysis: average penalties within reference, alternate, and distance
strata, then correlate the predictor and B1 side differences.

Label-free throughout.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from run_gate0 import load_frame
from run_gate0_within_site import site_gain_penalties, site_penalties
from run_gate1 import read_topology
from run_step2_composition_control import (BINS, PSEUDOCOUNT, b2_features,
                                           MINHASH_HASH_METHOD, minhash_clusters, symbol_families)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "results" / "step2_20260929"
# Treat each B2 context cell as if it had 50 pooled B1 observations.
B2_SHRINKAGE = 50.0


def unique_positions(frame: pd.DataFrame, mask: np.ndarray) -> pd.DataFrame:
    """One row per residue position, so the alphabet is not weighted by codon reach."""
    return frame[mask].drop_duplicates(subset=["protein", "position"])


def frequency_tables(frame: pd.DataFrame, mask: np.ndarray, context: np.ndarray):
    """B1 (side x bin) and B2 (side x bin x context), B2 shrunk toward B1."""
    rows = unique_positions(frame, mask)
    side = rows["side"].to_numpy().astype(int)
    bucket = rows["bin"].to_numpy().astype(int)
    ref = rows["ref_aa"].to_numpy().astype(int)

    b1_counts = np.full((2, len(BINS), 20), PSEUDOCOUNT)
    np.add.at(b1_counts, (side, bucket, ref), 1.0)
    b1_freq = b1_counts / b1_counts.sum(axis=2, keepdims=True)

    # Preserve original row indices when applying context after deduplication.
    ctx = context[rows.index.to_numpy()]
    b2_counts = np.zeros((2, len(BINS), 6, 20))
    np.add.at(b2_counts, (side, bucket, ctx.astype(int), ref), 1.0)
    # Hierarchical shrinkage toward B1: sparse context cells fall back to the
    # pooled footprint rather than to noise.
    b2_freq = ((b2_counts + B2_SHRINKAGE * b1_freq[:, :, None, :])
               / (b2_counts.sum(axis=3, keepdims=True) + B2_SHRINKAGE))
    return np.log(b1_freq), np.log(b2_freq)


def build_scored_frame(args) -> tuple[pd.DataFrame, dict]:
    proteins = read_topology(Path(args.topology))
    frame, _ = load_frame(Path(args.scores))
    blob = np.load(args.scores, allow_pickle=False)
    accessions = list(blob["accessions"])

    bin_index = np.full(len(frame), -1, dtype=int)
    distance = frame["distance"].to_numpy()
    for i, (low, high) in enumerate(BINS):
        bin_index[(distance >= low) & (distance <= high)] = i
    frame = frame[bin_index >= 0].reset_index(drop=True)
    frame["bin"] = bin_index[bin_index >= 0]

    features = b2_features(frame, proteins, accessions)
    context = (np.sign(features[:, 0]).astype(int) + 1) * 2 + features[:, 1].astype(int)

    clusters = (minhash_clusters(proteins, args.seed) if args.split == "sequence"
                else symbol_families(proteins))
    cluster_of = np.array([clusters.get(accessions[p], -1)
                           for p in frame["protein"].to_numpy()])
    unique = sorted(set(cluster_of.tolist()), key=str)
    rng = np.random.default_rng(args.seed)
    fold_of = {c: int(rng.integers(0, args.folds)) for c in unique}
    frame["fold"] = np.array([fold_of[c] for c in cluster_of])

    frame["b0"] = np.nan; frame["b1"] = np.nan; frame["b2"] = np.nan
    ref = frame["ref_aa"].to_numpy().astype(int)
    alt = frame["alt_aa"].to_numpy().astype(int)
    side = frame["side"].to_numpy().astype(int)
    bucket = frame["bin"].to_numpy().astype(int)

    for f in range(args.folds):
        train = (frame["fold"].to_numpy() != f)
        held = ~train
        if train.sum() == 0 or held.sum() == 0:
            continue
        log_b1, log_b2 = frequency_tables(frame, train, context)
        b1 = log_b1[side, bucket, ref] - log_b1[side, bucket, alt]
        b2 = log_b2[side, bucket, context.astype(int), ref] - log_b2[side, bucket, context.astype(int), alt]
        table = np.zeros((20, 20)); counts = np.zeros((20, 20))
        np.add.at(table, (ref[train], alt[train]), b1[train])
        np.add.at(counts, (ref[train], alt[train]), 1.0)
        mean = np.divide(table, counts, out=np.zeros_like(table), where=counts > 0)
        frame.loc[held, "b0"] = mean[ref[held], alt[held]]
        frame.loc[held, "b1"] = b1[held]
        frame.loc[held, "b2"] = b2[held]

    frame = frame.dropna(subset=["b0", "b1", "b2"]).reset_index(drop=True)
    frame = frame.rename(columns={"score": "alphamissense"})
    if Path(args.esm).exists():
        esm_frame, _ = load_frame(Path(args.esm))
        key = ["protein", "position", "alt_aa"]
        frame = frame.merge(esm_frame[key + ["score"]].rename(columns={"score": "esm1b"}),
                            on=key, how="left")
    return frame, {"n_clusters": len(unique), "n_rows": int(len(frame)),
                   "family_split": args.split, "folds": args.folds,
                   "seed": args.seed,
                   "minhash_hash_method": MINHASH_HASH_METHOD if args.split == "sequence" else None}


def penalties(frame: pd.DataFrame, column: str, kind: str, group: set[str]) -> pd.DataFrame:
    if kind == "loss":
        # The older Gate 0 helper defines loss as every alternate outside the
        # reference charge group, including reversals and histidine. Step 2 must
        # use the clean, uncharged comparator used by the project's later rules.
        alphabet = "ACDEFGHIKLMNPQRSTVWY"
        allowed = {alphabet.index(aa) for aa in (set("ACFGILMNPQSTVWY") | group)}
        frame = frame[frame["alt_aa"].isin(allowed)]
    work = frame.rename(columns={column: "score"})
    sites = site_gain_penalties(work, group) if kind == "gain" else site_penalties(work, group)
    return sites[sites["distance"] <= 10]


def side_difference(sites: pd.DataFrame, values: np.ndarray | None = None):
    """Efficiency-weighted side difference, returned with its per-protein parts."""
    penalty = sites["penalty"].to_numpy() if values is None else values
    buckets: dict = defaultdict(lambda: {"c": [], "o": []})
    for protein, substitution, cyto, value in zip(
            sites["protein"].to_numpy(), sites["substitution"].to_numpy(),
            sites["is_cytosolic"].to_numpy(), penalty):
        buckets[(int(protein), str(substitution))]["c" if cyto else "o"].append(float(value))
    proteins, weights, diffs = [], [], []
    for (protein, _), arms in buckets.items():
        n1, n0 = len(arms["c"]), len(arms["o"])
        if n1 == 0 or n0 == 0:
            continue
        proteins.append(protein)
        weights.append(n1 * n0 / (n1 + n0))
        diffs.append(float(np.mean(arms["c"]) - np.mean(arms["o"])))
    return np.array(proteins), np.array(weights), np.array(diffs)


def calibrated(frame: pd.DataFrame, predictor: str, kind: str, group: set[str],
               draws: int, seed: int, folds: int = 5) -> dict:
    """Primary: observed side difference minus beta x B1's, cross-fitted over all folds.

    For each fold the slope is fitted on the other folds and applied to that fold's
    held-out strata, so every quantity is evaluated out of fold. Folds are
    cluster-disjoint by protein, so a protein sits in exactly one fold and the
    pooled strata do not overlap.
    """
    work = frame.dropna(subset=[predictor, "b0", "b1"])
    key = ["protein", "position", "substitution"]
    betas, proteins_all, weights_all, residual_all = [], [], [], []
    observed_num = implied_num = weight_sum = 0.0
    per_fold = []

    for f in range(folds):
        train = work[work["fold"] != f]
        held = work[work["fold"] == f]
        if len(train) < 1000 or len(held) < 500:
            continue
        merged = (penalties(train, predictor, kind, group)[key + ["penalty"]]
                  .rename(columns={"penalty": "y"})
                  .merge(penalties(train, "b1", kind, group)[key + ["penalty"]]
                         .rename(columns={"penalty": "x1"}), on=key)
                  .merge(penalties(train, "b0", kind, group)[key + ["penalty"]]
                         .rename(columns={"penalty": "x0"}), on=key))
        if len(merged) < 100:
            continue
        design = np.column_stack([np.ones(len(merged)), merged["x0"].to_numpy(),
                                  merged["x1"].to_numpy()])
        coefficients, *_ = np.linalg.lstsq(design, merged["y"].to_numpy(), rcond=None)
        beta = float(coefficients[2])

        pair = (penalties(held, predictor, kind, group)[key + ["penalty", "is_cytosolic", "distance"]]
                .rename(columns={"penalty": "y"})
                .merge(penalties(held, "b1", kind, group)[key + ["penalty"]]
                       .rename(columns={"penalty": "x1"}), on=key))
        if len(pair) < 30:
            continue
        proteins, weights, dy = side_difference(pair, pair["y"].to_numpy())
        _, _, dx = side_difference(pair, pair["x1"].to_numpy())
        if dy.size == 0:
            continue
        betas.append(beta)
        proteins_all.append(proteins); weights_all.append(weights)
        residual_all.append(dy - beta * dx)
        observed_num += float((weights * dy).sum())
        implied_num += float((weights * beta * dx).sum())
        weight_sum += float(weights.sum())
        per_fold.append({"fold": f, "beta": beta, "n_strata": int(dy.size),
                         "n_held_out_sites": int(len(pair))})

    if not residual_all or weight_sum <= 0:
        return {"note": "no fold produced usable held-out strata"}

    proteins = np.concatenate(proteins_all)
    weights = np.concatenate(weights_all)
    residual = np.concatenate(residual_all)
    names = sorted(set(proteins.tolist()))
    index = {p: i for i, p in enumerate(names)}
    cluster = np.array([index[p] for p in proteins])

    rng = np.random.default_rng(seed)
    samples = np.empty(draws)
    for i in range(draws):
        counts = np.bincount(rng.integers(0, len(names), len(names)),
                             minlength=len(names)).astype(float)
        m = counts[cluster]
        total = (weights * m).sum()
        samples[i] = float((weights * m * residual).sum() / total) if total > 0 else np.nan
    finite = samples[np.isfinite(samples)]
    return {
        "folds_used": len(per_fold),
        "beta_on_b1_mean": float(np.mean(betas)),
        "beta_on_b1_by_fold": betas,
        "observed_side_difference": observed_num / weight_sum,
        "implied_by_calibration": implied_num / weight_sum,
        "residual": float((weights * residual).sum() / weights.sum()),
        "residual_ci95": [float(np.percentile(finite, 2.5)),
                          float(np.percentile(finite, 97.5))] if finite.size > 1 else [None, None],
        "n_strata": int(residual.size), "n_proteins": len(names),
        "per_fold": per_fold,
    }


def shape_agreement(frame: pd.DataFrame, predictor: str, kind: str, group: set[str]) -> dict:
    """Secondary: correlation of stratum-level side differences, noise averaged down."""
    work = frame.dropna(subset=[predictor, "b1"])
    sites_p = penalties(work, predictor, kind, group)
    sites_1 = penalties(work, "b1", kind, group)
    key = ["protein", "position", "substitution"]
    pair = (sites_p[key + ["penalty", "is_cytosolic", "distance"]].rename(columns={"penalty": "y"})
            .merge(sites_1[key + ["penalty"]].rename(columns={"penalty": "x1"}), on=key))
    if len(pair) < 100:
        return {"note": "too few shared sites"}
    pair["bin"] = pd.cut(pair["distance"], bins=[0, 5, 10, 15, 30, 60], labels=False)
    rows = []
    for (substitution, bucket), group_rows in pair.groupby(["substitution", "bin"], observed=True):
        cyto = group_rows[group_rows["is_cytosolic"]]
        outer = group_rows[~group_rows["is_cytosolic"]]
        if len(cyto) < 3 or len(outer) < 3:
            continue
        rows.append((float(cyto["y"].mean() - outer["y"].mean()),
                     float(cyto["x1"].mean() - outer["x1"].mean()), len(group_rows)))
    if len(rows) < 5:
        return {"note": "too few strata", "n_strata": len(rows)}
    y = np.array([r[0] for r in rows]); x = np.array([r[1] for r in rows])
    return {"pearson_r": float(np.corrcoef(y, x)[0, 1]),
            "spearman_r": float(np.corrcoef(np.argsort(np.argsort(y)),
                                            np.argsort(np.argsort(x)))[0, 1]),
            "n_strata": len(rows)}


def run(args) -> dict:
    frame, meta = build_scored_frame(args)
    rules = {"positive_charge_gain": ("gain", set("KR")),
             "positive_charge_loss": ("loss", set("KR")),
             "negative_charge_gain": ("gain", set("DE")),
             "negative_charge_loss": ("loss", set("DE"))}
    rules = {rule: spec for rule, spec in rules.items()
             if not args.rules or rule in args.rules}
    predictors = [p for p in ("alphamissense", "esm1b", "b2")
                  if p in frame.columns and (not args.predictors or p in args.predictors)]

    results: dict = {}
    for rule, (kind, group) in rules.items():
        results[rule] = {}
        for predictor in predictors:
            results[rule][predictor] = {
                "calibrated": calibrated(frame, predictor, kind, group,
                                         args.boot_draws, args.seed, args.folds),
                "shape": shape_agreement(frame, predictor, kind, group),
            }

    summary = {
        "analysis": "Step 2 calibrated comparison (replaces the frozen self-normalised metric)",
        "prereg_sha256": "06cdc2e8449ab03ced939b6cf1599d7d754a5fcec36c4a2ca992c5e9d1de5385",
        "deviation": "Frozen metric divides each source by its own within-site spread, which "
                     "mechanically favours a noiseless baseline. Replaced with a calibrated "
                     "residual and a shape correlation. See the decision log.",
        "clinical_labels_used": False, "dms_effects_used": False,
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "split": args.split, **meta, "b2_shrinkage": B2_SHRINKAGE,
        "charge_loss_comparator": "uncharged alternates only; charge reversals and histidine excluded (2026-10-03 correction)",
        "selected_rules": list(rules), "selected_predictors": predictors,
        "results": results,
    }
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    (out / f"step2_calibrated_{args.split}.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    for rule in rules:
        print(f"\n=== {rule} ===")
        for predictor in predictors:
            c = results[rule][predictor]["calibrated"]
            s = results[rule][predictor]["shape"]
            if "residual" not in c:
                print(f"  {predictor:<14} {c.get('note')}")
                continue
            lo, hi = c["residual_ci95"]
            print(f"  {predictor:<14} beta={c['beta_on_b1_mean']:+.4f} "
                  f"({c['folds_used']} folds)  observed "
                  f"{c['observed_side_difference']:+.4f}  implied {c['implied_by_calibration']:+.4f}")
            print(f"   {'':<13} residual {c['residual']:+.4f} CI [{lo:+.4f},{hi:+.4f}]"
                  f"   shape r={s.get('pearson_r', float('nan')):+.3f} "
                  f"(strata {s.get('n_strata', 0)})")
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
    parser.add_argument("--rules", nargs="+", choices=("positive_charge_gain", "positive_charge_loss",
                                                       "negative_charge_gain", "negative_charge_loss"))
    parser.add_argument("--predictors", nargs="+", choices=("alphamissense", "esm1b", "b2"))
    run(parser.parse_args())


if __name__ == "__main__":
    main()

"""Gate 0 robustness check: within-site charge contrasts in AlphaMissense.

Label-free. Gate 0 stratified by protein and exact substitution, which controls
protein identity and substitution identity but *not* how constrained the
individual residue is. If cytosolic K/R residues sit in more constrained sites
than outer ones, Gate 0's side difference is still confounded.

This removes that. At each K or R site, AlphaMissense scores both the
charge-preserving substitution (K<->R) and charge-losing ones. Their difference
is a charge-loss penalty computed entirely inside one site, so any constraint
shared by the whole residue cancels:

    penalty(site) = mean(score | charge lost) - mean(score | charge preserved)

The side comparison is then a second difference, taken across sites within the
same protein and the same reference residue. The same construction runs at D/E
sites using D<->E as the preserving substitution.

If the positive-inside rule is what drives Gate 0, the penalty should be larger
on the cytosolic side and should shrink with distance from the membrane.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from run_gate0 import (DISTANCE_BINS, PRIMARY_WINDOW, WINDOWS, combine, interval,
                       load_frame, sign_p_value, stratum_statistics)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCORES = ROOT / "data" / "processed" / "gate0" / "am_flank_scores.npz"
DEFAULT_OUT = ROOT / "results" / "gate0_20260924"

GROUPS = {"positive": set("KR"), "negative": set("DE")}


def site_penalties(frame: pd.DataFrame, group: set[str]) -> pd.DataFrame:
    """One charge-loss penalty per site, from substitutions at that site alone."""
    alphabet = sorted({*"ACDEFGHIKLMNPQRSTVWY"})
    codes = {aa: i for i, aa in enumerate(alphabet)}
    group_codes = {codes[aa] for aa in group}

    work = frame[frame["ref_aa"].isin(group_codes)].copy()
    if work.empty:
        return pd.DataFrame()
    work["preserves"] = work["alt_aa"].isin(group_codes)

    keys = ["protein", "position"]
    grouped = work.groupby(keys, sort=False)
    stats = grouped.agg(
        n_preserving=("preserves", "sum"),
        n_total=("preserves", "size"),
        side=("side", "first"),
        distance=("distance", "first"),
        tmd_index=("tmd_index", "first"),
        ref_aa=("ref_aa", "first"),
    )
    preserved = work[work["preserves"]].groupby(keys, sort=False)["score"].mean()
    lost = work[~work["preserves"]].groupby(keys, sort=False)["score"].mean()
    stats = stats.join(preserved.rename("score_preserved")).join(lost.rename("score_lost"))
    stats = stats[(stats["n_preserving"] > 0)
                  & (stats["n_total"] - stats["n_preserving"] > 0)].reset_index()
    if stats.empty:
        return stats
    stats["penalty"] = stats["score_lost"] - stats["score_preserved"]
    stats["is_cytosolic"] = stats["side"].to_numpy() == 0
    # Stratify on the reference residue: a K site is only compared with K sites.
    stats["substitution"] = np.array(alphabet)[stats["ref_aa"].to_numpy()]
    return stats


def site_gain_penalties(frame: pd.DataFrame, group: set[str]) -> pd.DataFrame:
    """One charge-gain penalty per uncharged site.

    At a residue that carries no charge, AlphaMissense scores both substitutions
    that introduce a charge and substitutions that do not. Their difference is the
    cost of introducing that charge, measured inside one site:

        penalty(site) = mean(score | charge introduced) - mean(score | still uncharged)

    The positive-inside rule predicts introducing K/R is worse on the *outer* side,
    so the side difference should be NEGATIVE. Gate 0's between-site estimator found
    it positive, the wrong sign; this construction checks whether that survives once
    site constraint is removed.
    """
    alphabet = sorted({*"ACDEFGHIKLMNPQRSTVWY"})
    codes = {aa: i for i, aa in enumerate(alphabet)}
    charged = {codes[aa] for aa in set("KRDEH")}
    group_codes = {codes[aa] for aa in group}

    work = frame[~frame["ref_aa"].isin(charged)].copy()
    work = work[work["alt_aa"].isin(group_codes) | ~work["alt_aa"].isin(charged)]
    if work.empty:
        return pd.DataFrame()
    work["gains"] = work["alt_aa"].isin(group_codes)

    keys = ["protein", "position"]
    stats = work.groupby(keys, sort=False).agg(
        n_gaining=("gains", "sum"),
        n_total=("gains", "size"),
        side=("side", "first"),
        distance=("distance", "first"),
        tmd_index=("tmd_index", "first"),
        ref_aa=("ref_aa", "first"),
    )
    gaining = work[work["gains"]].groupby(keys, sort=False)["score"].mean()
    neutral = work[~work["gains"]].groupby(keys, sort=False)["score"].mean()
    stats = stats.join(gaining.rename("score_gaining")).join(neutral.rename("score_neutral"))
    stats = stats[(stats["n_gaining"] > 0)
                  & (stats["n_total"] - stats["n_gaining"] > 0)].reset_index()
    if stats.empty:
        return stats
    stats["penalty"] = stats["score_gaining"] - stats["score_neutral"]
    stats["is_cytosolic"] = stats["side"].to_numpy() == 0
    stats["substitution"] = np.array(alphabet)[stats["ref_aa"].to_numpy()]
    return stats


def estimate(sites: pd.DataFrame, column: str, n_proteins: int,
             draws: int, seed: int) -> dict:
    if sites.empty:
        return {"n_sites": 0}
    frame = sites.rename(columns={column: "score"})[
        ["protein", "substitution", "is_cytosolic", "score"]]
    stats = stratum_statistics(frame)
    difference, auc = combine(stats)
    rng = np.random.default_rng(seed)
    samples = np.full(draws, np.nan)
    if stats["protein"].size:
        for i in range(draws):
            counts = np.bincount(rng.integers(0, n_proteins, n_proteins),
                                 minlength=n_proteins).astype(float)
            samples[i], _ = combine(stats, counts[stats["protein"]])
    return {
        "n_sites": int(len(sites)),
        "n_sites_cytosolic": int(sites["is_cytosolic"].sum()),
        "n_sites_outer": int((~sites["is_cytosolic"]).sum()),
        "n_strata": int(stats["protein"].size),
        "cytosolic_mean": float(sites.loc[sites["is_cytosolic"], column].mean()),
        "outer_mean": float(sites.loc[~sites["is_cytosolic"], column].mean()),
        "side_difference": difference,
        "ci95": interval(samples),
        "bootstrap_two_sided_p": sign_p_value(samples),
        "p_cytosolic_higher": auc,
    }


def run(args: argparse.Namespace) -> dict:
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    frame, qc = load_frame(Path(args.scores))
    n_proteins = int(frame["protein"].max()) + 1

    results: dict = {}
    for name, group in GROUPS.items():
        sites = site_penalties(frame, group)
        block: dict = {
            "n_sites_with_both_comparators": int(len(sites)),
            "by_window": {}, "distance_dose_response": {}, "first_tmd_only": {},
            "diagnostic_preserving_substitution_only": {},
        }
        for window in WINDOWS:
            subset = sites[sites["distance"] <= window]
            block["by_window"][str(window)] = estimate(
                subset, "penalty", n_proteins, args.draws, args.seed)
        for low, high in DISTANCE_BINS:
            subset = sites[(sites["distance"] >= low) & (sites["distance"] <= high)]
            block["distance_dose_response"][f"{low}-{high}"] = estimate(
                subset, "penalty", n_proteins, args.draws, args.seed)
        first = sites[(sites["distance"] <= PRIMARY_WINDOW) & (sites["tmd_index"] == 1)]
        block["first_tmd_only"] = estimate(first, "penalty", n_proteins, args.draws, args.seed)

        # Diagnostic: the charge-preserving substitution on its own. A side
        # difference here is residual site-level constraint that the penalty
        # construction is designed to cancel.
        primary = sites[sites["distance"] <= PRIMARY_WINDOW]
        block["diagnostic_preserving_substitution_only"] = estimate(
            primary, "score_preserved", n_proteins, args.draws, args.seed)
        block["diagnostic_losing_substitutions_only"] = estimate(
            primary, "score_lost", n_proteins, args.draws, args.seed)

        gain_sites = site_gain_penalties(frame, group)
        gain_block: dict = {"n_sites_with_both_comparators": int(len(gain_sites)),
                            "by_window": {}, "distance_dose_response": {}}
        for window in WINDOWS:
            subset = gain_sites[gain_sites["distance"] <= window]
            gain_block["by_window"][str(window)] = estimate(
                subset, "penalty", n_proteins, args.draws, args.seed)
        for low, high in DISTANCE_BINS:
            subset = gain_sites[(gain_sites["distance"] >= low) & (gain_sites["distance"] <= high)]
            gain_block["distance_dose_response"][f"{low}-{high}"] = estimate(
                subset, "penalty", n_proteins, args.draws, args.seed)
        first = gain_sites[(gain_sites["distance"] <= PRIMARY_WINDOW)
                           & (gain_sites["tmd_index"] == 1)]
        gain_block["first_tmd_only"] = estimate(first, "penalty", n_proteins,
                                                args.draws, args.seed)
        block["gain_at_uncharged_sites"] = gain_block
        results[name] = block

    summary = {
        "analysis": "Gate 0 robustness: within-site charge-loss penalty by membrane side",
        "clinical_labels_used": False,
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "construction": "penalty(site) = mean(score | charge lost) - mean(score | charge preserved), "
                        "computed within a single residue; side comparison stratified by "
                        "(protein, reference residue) with a protein-clustered bootstrap",
        "sign_convention": "Positive side difference = the charge-loss penalty is larger on "
                           "the cytosolic side",
        "quality_control": qc,
        "results": results,
    }
    (out_dir / "gate0_within_site.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    for name, block in results.items():
        print(f"\n=== {name} charge, within-site penalty ===")
        print(f"  sites with both comparators: {block['n_sites_with_both_comparators']}")
        primary = block["by_window"][str(PRIMARY_WINDOW)]
        print(f"  window <=10: cyto {primary['cytosolic_mean']:+.4f} vs outer "
              f"{primary['outer_mean']:+.4f}  diff {primary['side_difference']:+.4f} "
              f"{primary['ci95']}  n_sites {primary['n_sites']}  strata {primary['n_strata']}")
        for key, value in block["distance_dose_response"].items():
            if value.get("n_strata"):
                print(f"    {key:>6}: diff {value['side_difference']:+.4f} "
                      f"[{value['ci95'][0]:+.4f},{value['ci95'][1]:+.4f}]  "
                      f"sites {value['n_sites']:>6}")
        diagnostic = block["diagnostic_preserving_substitution_only"]
        print(f"  diagnostic, preserving substitution alone: diff "
              f"{diagnostic['side_difference']:+.4f} {diagnostic['ci95']}")
        gain = block["gain_at_uncharged_sites"]
        gain_primary = gain["by_window"][str(PRIMARY_WINDOW)]
        print(f"  GAIN at uncharged sites (rule predicts NEGATIVE): "
              f"window<=10 diff {gain_primary['side_difference']:+.4f} "
              f"[{gain_primary['ci95'][0]:+.4f},{gain_primary['ci95'][1]:+.4f}]  "
              f"sites {gain_primary['n_sites']}")
        for key, value in gain["distance_dose_response"].items():
            if value.get("n_strata"):
                print(f"    {key:>6}: diff {value['side_difference']:+.4f} "
                      f"[{value['ci95'][0]:+.4f},{value['ci95'][1]:+.4f}]  "
                      f"sites {value['n_sites']:>6}")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", default=str(DEFAULT_SCORES))
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    parser.add_argument("--draws", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260924)
    run(parser.parse_args())


if __name__ == "__main__":
    main()

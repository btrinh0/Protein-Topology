"""Estimate slope uncertainty and decompose the gain penalty by arm.

The five folds share training data, so their slopes are not independent estimates.
This script computes a protein-clustered bootstrap interval for the slope.

The gain penalty is score(X -> K/R) - score(X -> uncharged). The analysis regresses
the K/R and comparator arms separately against their B1 counterparts, then repeats
the comparator analysis for polar-only and hydrophobic-only residue sets. This
distinguishes a charge-specific signal from constraints in the comparator arm.

Label-free throughout.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from run_step2_calibrated import build_scored_frame

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "results" / "step2_20260929"
ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
CODE = {a: i for i, a in enumerate(ALPHABET)}

POSITIVE = {CODE["K"], CODE["R"]}
CHARGED = {CODE[a] for a in "KRDEH"}
POLAR = {CODE[a] for a in "QNST"}
HYDROPHOBIC = {CODE[a] for a in "AVLIMFW"}


def arm_table(frame: pd.DataFrame, columns: list[str], comparator: set[int]) -> pd.DataFrame:
    """Per-site mean score in the K/R arm and in the comparator arm, for each column.

    Sites are uncharged references within 10 residues, matching the gain penalty.
    """
    work = frame[(frame["distance"] <= 10) & (~frame["ref_aa"].isin(CHARGED))]
    alt = work["alt_aa"].to_numpy()
    in_numerator = np.isin(alt, list(POSITIVE))
    in_comparator = np.isin(alt, list(comparator))
    work = work.assign(_num=in_numerator, _cmp=in_comparator)
    work = work[in_numerator | in_comparator]
    if work.empty:
        return pd.DataFrame()

    grouped = work.groupby(["protein", "position", "side"], sort=False)
    out = grouped.agg(n_num=("_num", "sum"), n_cmp=("_cmp", "sum")).reset_index()
    for column in columns:
        numerator = (work[work["_num"]].groupby(["protein", "position", "side"], sort=False)[column]
                     .mean().rename(f"{column}_num"))
        comparator_mean = (work[work["_cmp"]].groupby(["protein", "position", "side"], sort=False)[column]
                           .mean().rename(f"{column}_cmp"))
        out = out.merge(numerator.reset_index(), on=["protein", "position", "side"], how="left")
        out = out.merge(comparator_mean.reset_index(), on=["protein", "position", "side"], how="left")
    return out[(out["n_num"] > 0) & (out["n_cmp"] > 0)].dropna()


def clustered_slope(y: np.ndarray, x0: np.ndarray, x1: np.ndarray,
                    protein: np.ndarray, draws: int, seed: int) -> dict:
    """Slope on x1 controlling x0, with a protein-clustered bootstrap interval."""
    def fit(mask_rows: np.ndarray) -> float:
        design = np.column_stack([np.ones(mask_rows.sum()), x0[mask_rows], x1[mask_rows]])
        if design.shape[0] < 10:
            return float("nan")
        coefficients, *_ = np.linalg.lstsq(design, y[mask_rows], rcond=None)
        return float(coefficients[2])

    point = fit(np.ones(len(y), dtype=bool))
    names = np.unique(protein)
    index = {p: i for i, p in enumerate(names)}
    cluster = np.array([index[p] for p in protein])
    by_cluster = [np.flatnonzero(cluster == i) for i in range(len(names))]

    rng = np.random.default_rng(seed)
    samples = np.full(draws, np.nan)
    for i in range(draws):
        picked = rng.integers(0, len(names), len(names))
        rows = np.concatenate([by_cluster[p] for p in picked])
        design = np.column_stack([np.ones(rows.size), x0[rows], x1[rows]])
        try:
            coefficients, *_ = np.linalg.lstsq(design, y[rows], rcond=None)
            samples[i] = float(coefficients[2])
        except np.linalg.LinAlgError:
            continue
    finite = samples[np.isfinite(samples)]
    low, high = ((float(np.percentile(finite, 2.5)), float(np.percentile(finite, 97.5)))
                 if finite.size > 1 else (float("nan"), float("nan")))
    return {"slope": point, "ci95": [low, high],
            "excludes_zero": bool(np.isfinite(low) and (low > 0 or high < 0)),
            "n_sites": int(len(y)), "n_proteins": int(len(names))}


def run(args) -> dict:
    frame, meta = build_scored_frame(args)
    predictors = [p for p in ("alphamissense", "esm1b") if p in frame.columns]
    columns = predictors + ["b0", "b1"]

    results: dict = {"whole_penalty": {}, "by_arm": {}, "comparator_sensitivity": {}}

    comparators = {"all_uncharged": {c for c in range(20)} - CHARGED,
                   "polar_only": POLAR, "hydrophobic_only": HYDROPHOBIC}

    for name, comparator in comparators.items():
        table = arm_table(frame, columns, comparator)
        if table.empty:
            results["comparator_sensitivity"][name] = {"note": "no sites"}
            continue
        protein = table["protein"].to_numpy()
        b1_pen = table["b1_num"].to_numpy() - table["b1_cmp"].to_numpy()
        b0_pen = table["b0_num"].to_numpy() - table["b0_cmp"].to_numpy()

        block = {}
        for predictor in predictors:
            y_pen = table[f"{predictor}_num"].to_numpy() - table[f"{predictor}_cmp"].to_numpy()
            block[predictor] = clustered_slope(y_pen, b0_pen, b1_pen, protein,
                                               args.boot_draws, args.seed)
        if name == "all_uncharged":
            results["whole_penalty"] = block
            for predictor in predictors:
                results["by_arm"][predictor] = {
                    "kr_arm": clustered_slope(
                        table[f"{predictor}_num"].to_numpy(), table["b0_num"].to_numpy(),
                        table["b1_num"].to_numpy(), protein, args.boot_draws, args.seed),
                    "comparator_arm": clustered_slope(
                        table[f"{predictor}_cmp"].to_numpy(), table["b0_cmp"].to_numpy(),
                        table["b1_cmp"].to_numpy(), protein, args.boot_draws, args.seed),
                }
        else:
            results["comparator_sensitivity"][name] = block

    summary = {
        "analysis": "Step 2 follow-ups: clustered interval on the slope, and arm decomposition",
        "rule": "positive_charge_gain, window <=10",
        "clinical_labels_used": False, "dms_effects_used": False,
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        **meta, "results": results,
    }
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    (out / "step2_arms.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    def show(label, block):
        lo, hi = block["ci95"]
        flag = "EXCLUDES 0" if block["excludes_zero"] else "includes 0"
        print(f"   {label:<22} {block['slope']:+.5f}  CI [{lo:+.5f}, {hi:+.5f}]  {flag}"
              f"   sites={block['n_sites']}")

    print("=== whole gain penalty: slope on B1, protein-clustered ===")
    for predictor in predictors:
        show(predictor, results["whole_penalty"][predictor])
    print("\n=== decomposed by arm ===")
    for predictor in predictors:
        print(f"  {predictor}")
        show("K/R arm", results["by_arm"][predictor]["kr_arm"])
        show("comparator arm", results["by_arm"][predictor]["comparator_arm"])
    print("\n=== restricted comparator sets ===")
    for name, block in results["comparator_sensitivity"].items():
        print(f"  {name}")
        for predictor in predictors:
            if predictor in block:
                show(predictor, block[predictor])
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", default=str(ROOT / "data/raw/topology_viewer.html"))
    parser.add_argument("--scores", default=str(ROOT / "data/processed/gate0/am_full_flank_scores.npz"))
    parser.add_argument("--esm", default=str(ROOT / "data/processed/gate0/esm1b_flank_scores.npz"))
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    parser.add_argument("--split", choices=("sequence", "symbol"), default="sequence")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--boot-draws", type=int, default=400)
    parser.add_argument("--seed", type=int, default=20260929)
    run(parser.parse_args())


if __name__ == "__main__":
    main()

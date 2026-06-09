"""Gate 0: does AlphaMissense itself encode the positive-inside rule?

Label-free. This script never reads ClinVar. It asks whether AlphaMissense scores
the *same* amino-acid substitution in the *same* protein differently depending on
which side of the membrane the flank faces.

Side difference (per substitution class, per scope):
  Stratify by (protein, exact substitution). Keep strata with at least one scored
  variant on each side. Within a stratum take the cytosolic-minus-outer difference
  in AlphaMissense pathogenicity, then combine strata with efficient weights
  n_cyto*n_outer/(n_cyto+n_outer). Also report the stratified rank statistic
  P(cytosolic score > outer score), which assumes no scale.

  Positive = AlphaMissense calls the substitution more damaging on the cytosolic side.

A raw side difference is not evidence for the positive-inside rule on its own: if
cytosolic flanks are simply more constrained, every substitution class shifts the
same way. The rule makes a *directional* prediction, so the primary quantities are
differences in differences:

  positive_loss_minus_gain     losing K/R should hurt more inside while gaining K/R
                               should hurt more outside, so the contrast is signed
                               and cancels any side baseline shared by all classes.
  positive_loss_minus_placebo  excess over substitutions that change no charge.

Uncertainty is a protein-clustered bootstrap. Every class in a scope is resampled
under the same protein draw, so contrasts keep the correlation between classes.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCORES = ROOT / "data" / "processed" / "gate0" / "am_flank_scores.npz"
DEFAULT_OUT = ROOT / "results" / "gate0_20260924"

POSITIVE = set("KR")
NEGATIVE = set("DE")
AMBIGUOUS = set("H")  # histidine is only partly protonated at physiological pH

CLASSES = ["positive_charge_loss", "positive_charge_gain", "negative_charge_loss",
           "negative_charge_gain", "charge_neutral_placebo"]
CONTRASTS = {
    "positive_loss_minus_gain": ("positive_charge_loss", "positive_charge_gain"),
    "positive_loss_minus_placebo": ("positive_charge_loss", "charge_neutral_placebo"),
    "positive_gain_minus_placebo": ("positive_charge_gain", "charge_neutral_placebo"),
    "negative_loss_minus_gain": ("negative_charge_loss", "negative_charge_gain"),
    "negative_loss_minus_placebo": ("negative_charge_loss", "charge_neutral_placebo"),
    "negative_gain_minus_placebo": ("negative_charge_gain", "charge_neutral_placebo"),
}
DISTANCE_BINS = [(1, 5), (6, 10), (11, 15), (16, 30), (31, 60)]
WINDOWS = (5, 10, 15)
PRIMARY_WINDOW = 10


def substitution_class(ref: str, alt: str) -> str | None:
    if ref in POSITIVE and alt not in POSITIVE:
        return "positive_charge_loss"
    if ref not in POSITIVE and alt in POSITIVE:
        return "positive_charge_gain"
    if ref in NEGATIVE and alt not in NEGATIVE:
        return "negative_charge_loss"
    if ref not in NEGATIVE and alt in NEGATIVE:
        return "negative_charge_gain"
    if not ({ref, alt} & (POSITIVE | NEGATIVE | AMBIGUOUS)):
        return "charge_neutral_placebo"
    return None


def load_frame(path: Path) -> tuple[pd.DataFrame, dict]:
    blob = np.load(path, allow_pickle=False)
    alphabet = np.array(blob["aa_alphabet"].tolist())
    frame = pd.DataFrame({
        "protein": blob["protein"],
        "position": blob["position"],
        "side": blob["side"],
        "distance": blob["distance"].astype(np.int16),
        "tmd_index": blob["tmd_index"].astype(np.int16),
        "ref_aa": blob["ref_aa"],
        "alt_aa": blob["alt_aa"],
        "score": blob["am_pathogenicity"].astype(np.float64),
    })
    qc = {"rows_loaded": int(len(frame))}

    # One genomic SNV per row; two SNVs in a codon can encode the same substitution.
    qc["rows_in_duplicate_protein_variants"] = int(
        frame.duplicated(subset=["protein", "position", "alt_aa"], keep=False).sum()
    )
    frame = frame.groupby(
        ["protein", "position", "alt_aa"], as_index=False, sort=False
    ).agg({"side": "first", "distance": "first", "tmd_index": "first",
           "ref_aa": "first", "score": "mean"})
    qc["unique_protein_variants"] = int(len(frame))

    frame["substitution"] = (pd.Series(alphabet[frame["ref_aa"].to_numpy()], index=frame.index)
                             + ">" + pd.Series(alphabet[frame["alt_aa"].to_numpy()], index=frame.index))
    classes = {f"{ref}>{alt}": substitution_class(ref, alt)
               for ref in alphabet for alt in alphabet if ref != alt}
    frame["class"] = frame["substitution"].map(classes)
    frame["is_cytosolic"] = frame["side"].to_numpy() == 0
    qc["unique_protein_variants_in_a_named_class"] = int(frame["class"].notna().sum())
    return frame, qc


def stratum_statistics(frame: pd.DataFrame) -> dict[str, np.ndarray]:
    """Per-(protein, substitution) sufficient statistics for both estimators."""
    empty = {key: np.array([], dtype=float) for key in ("protein", "weight", "diff", "u", "pairs", "n")}
    if frame.empty:
        return empty
    work = frame[["protein", "substitution", "is_cytosolic", "score"]].copy()
    keys = ["protein", "substitution"]
    groups = work.groupby(keys, sort=False)

    # Rank-based U: within-stratum ranks of cytosolic scores, ties averaged.
    work["rank"] = groups["score"].rank(method="average")
    work["cyto_rank"] = np.where(work["is_cytosolic"], work["rank"], 0.0)
    work["cyto_score"] = np.where(work["is_cytosolic"], work["score"], 0.0)
    work["outer_score"] = np.where(work["is_cytosolic"], 0.0, work["score"])

    stats = work.groupby(keys, sort=False).agg(
        n1=("is_cytosolic", "sum"),
        n=("is_cytosolic", "size"),
        rank_sum=("cyto_rank", "sum"),
        s1=("cyto_score", "sum"),
        s0=("outer_score", "sum"),
    ).reset_index()

    n1 = stats["n1"].to_numpy(dtype=float)
    n0 = stats["n"].to_numpy(dtype=float) - n1
    keep = (n1 > 0) & (n0 > 0)
    if not keep.any():
        return empty
    n1, n0 = n1[keep], n0[keep]
    return {
        "protein": stats["protein"].to_numpy()[keep],
        "weight": n1 * n0 / (n1 + n0),
        "diff": stats["s1"].to_numpy()[keep] / n1 - stats["s0"].to_numpy()[keep] / n0,
        "u": stats["rank_sum"].to_numpy()[keep] - n1 * (n1 + 1.0) / 2.0,
        "pairs": n1 * n0,
        "n": n1 + n0,
    }


def combine(stats: dict[str, np.ndarray], multiplier: np.ndarray | None = None) -> tuple[float, float]:
    if stats["protein"].size == 0:
        return float("nan"), float("nan")
    m = 1.0 if multiplier is None else multiplier
    weight = stats["weight"] * m
    pairs = stats["pairs"] * m
    total_weight, total_pairs = weight.sum(), pairs.sum()
    if total_weight <= 0 or total_pairs <= 0:
        return float("nan"), float("nan")
    return float((weight * stats["diff"]).sum() / total_weight), float((stats["u"] * m).sum() / total_pairs)


def interval(samples: np.ndarray) -> list[float | None]:
    finite = samples[np.isfinite(samples)]
    if finite.size < 2:
        return [None, None]
    return [float(np.percentile(finite, 2.5)), float(np.percentile(finite, 97.5))]


def sign_p_value(samples: np.ndarray) -> float | None:
    finite = samples[np.isfinite(samples)]
    if finite.size < 2:
        return None
    below = float((finite <= 0).mean())
    above = float((finite >= 0).mean())
    return float(min(1.0, 2.0 * min(below, above)))


def estimate_block(frames: dict[str, pd.DataFrame], n_proteins: int,
                   draws: int, seed: int) -> dict:
    stats = {name: stratum_statistics(frame) for name, frame in frames.items()}
    rng = np.random.default_rng(seed)
    diff_boot = {name: np.full(draws, np.nan) for name in stats}
    auc_boot = {name: np.full(draws, np.nan) for name in stats}
    for i in range(draws):
        counts = np.bincount(rng.integers(0, n_proteins, n_proteins),
                             minlength=n_proteins).astype(float)
        for name, stat in stats.items():
            if stat["protein"].size:
                diff_boot[name][i], auc_boot[name][i] = combine(stat, counts[stat["protein"]])

    block: dict = {"classes": {}, "contrasts": {}}
    for name, stat in stats.items():
        difference, auc = combine(stat)
        block["classes"][name] = {
            "n_strata": int(stat["protein"].size),
            "n_proteins_with_both_sides": int(np.unique(stat["protein"]).size),
            "n_scored_variants_in_strata": int(stat["n"].sum()) if stat["n"].size else 0,
            "mean_difference_cytosolic_minus_outer": difference,
            "mean_difference_ci95": interval(diff_boot[name]),
            "p_cytosolic_higher": auc,
            "p_cytosolic_higher_ci95": interval(auc_boot[name]),
        }
    for label, (left, right) in CONTRASTS.items():
        if left not in stats or right not in stats:
            continue
        samples = diff_boot[left] - diff_boot[right]
        point = (block["classes"][left]["mean_difference_cytosolic_minus_outer"]
                 - block["classes"][right]["mean_difference_cytosolic_minus_outer"])
        block["contrasts"][label] = {
            "difference_in_differences": point,
            "ci95": interval(samples),
            "bootstrap_two_sided_p": sign_p_value(samples),
        }
    return block


def split_classes(frame: pd.DataFrame) -> dict[str, pd.DataFrame]:
    return {name: frame[frame["class"] == name] for name in CLASSES}


def run(args: argparse.Namespace) -> dict:
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    frame, qc = load_frame(Path(args.scores))
    n_proteins = int(frame["protein"].max()) + 1
    named = frame[frame["class"].notna()]

    results: dict = {"by_window": {}, "distance_dose_response": {}}

    for window in WINDOWS:
        subset = named[named["distance"] <= window]
        results["by_window"][str(window)] = estimate_block(
            split_classes(subset), n_proteins, args.draws, args.seed)

    for low, high in DISTANCE_BINS:
        subset = named[(named["distance"] >= low) & (named["distance"] <= high)]
        results["distance_dose_response"][f"{low}-{high}"] = estimate_block(
            split_classes(subset), n_proteins, args.draws, args.seed)

    # Distance-matched: strata must also share a 5-residue distance bin.
    matched = named[named["distance"] <= PRIMARY_WINDOW].copy()
    matched["substitution"] = (matched["substitution"] + "@"
                               + np.where(matched["distance"].to_numpy() <= 5, "1-5", "6-10"))
    results["distance_matched_strata"] = estimate_block(
        split_classes(matched), n_proteins, args.draws, args.seed)

    first = named[(named["distance"] <= PRIMARY_WINDOW) & (named["tmd_index"] == 1)]
    results["first_tmd_only"] = estimate_block(
        split_classes(first), n_proteins, args.draws, args.seed)

    primary = named[(named["distance"] <= PRIMARY_WINDOW)
                    & (named["class"] == "positive_charge_loss")]
    per_transition = {}
    for substitution, group in primary.groupby("substitution", sort=True):
        block = estimate_block({substitution: group}, n_proteins,
                               max(args.draws // 4, 200), args.seed)
        entry = block["classes"][substitution]
        if entry["n_strata"] >= args.min_strata:
            per_transition[substitution] = entry
    results["per_transition_primary"] = per_transition

    window_frame = named[named["distance"] <= PRIMARY_WINDOW]
    results["unstratified_means"] = {
        name: {
            "cytosolic_mean": float(subset.loc[subset["is_cytosolic"], "score"].mean()),
            "outer_mean": float(subset.loc[~subset["is_cytosolic"], "score"].mean()),
            "cytosolic_n": int(subset["is_cytosolic"].sum()),
            "outer_n": int((~subset["is_cytosolic"]).sum()),
        }
        for name, subset in split_classes(window_frame).items()
    }

    summary = {
        "analysis": "Gate 0: is the positive-inside rule present in AlphaMissense scores?",
        "clinical_labels_used": False,
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "primary_window_residues": PRIMARY_WINDOW,
        "primary_contrast": "positive_loss_minus_gain",
        "estimator": "Stratified by (protein, exact substitution); weights n1*n0/(n1+n0); "
                     f"protein-clustered bootstrap, {args.draws} draws, shared across "
                     "classes within a scope",
        "sign_convention": "Positive = more damaging on the cytosolic side. The "
                           "positive-inside rule predicts a positive difference for "
                           "K/R loss and a negative one for K/R gain.",
        "quality_control": qc,
        "results": results,
    }
    (out_dir / "gate0_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    rows = []
    scopes = [("window", results["by_window"]), ("distance_bin", results["distance_dose_response"]),
              ("distance_matched", {"10": results["distance_matched_strata"]}),
              ("first_tmd", {"10": results["first_tmd_only"]})]
    for scope, block in scopes:
        for key, content in block.items():
            for name, values in content["classes"].items():
                rows.append({
                    "scope": scope, "key": key, "quantity": "side_difference", "term": name,
                    "n_strata": values["n_strata"],
                    "n_variants": values["n_scored_variants_in_strata"],
                    "estimate": values["mean_difference_cytosolic_minus_outer"],
                    "ci_low": values["mean_difference_ci95"][0],
                    "ci_high": values["mean_difference_ci95"][1],
                    "p_cytosolic_higher": values["p_cytosolic_higher"],
                })
            for name, values in content["contrasts"].items():
                rows.append({
                    "scope": scope, "key": key, "quantity": "difference_in_differences",
                    "term": name, "n_strata": None, "n_variants": None,
                    "estimate": values["difference_in_differences"],
                    "ci_low": values["ci95"][0], "ci_high": values["ci95"][1],
                    "p_cytosolic_higher": None,
                })
    pd.DataFrame(rows).to_csv(out_dir / "gate0_estimates.csv", index=False)

    primary_block = results["by_window"][str(PRIMARY_WINDOW)]
    print(json.dumps({
        "primary_window_classes": {
            name: {"d": value["mean_difference_cytosolic_minus_outer"],
                   "ci": value["mean_difference_ci95"]}
            for name, value in primary_block["classes"].items()
        },
        "primary_window_contrasts": primary_block["contrasts"],
        "first_tmd_contrasts": results["first_tmd_only"]["contrasts"],
    }, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", default=str(DEFAULT_SCORES))
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    parser.add_argument("--draws", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260924)
    parser.add_argument("--min-strata", type=int, default=25)
    run(parser.parse_args())


if __name__ == "__main__":
    main()

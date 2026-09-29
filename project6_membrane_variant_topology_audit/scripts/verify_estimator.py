"""Verification of the primary estimator: independent implementation plus simulation.

Two checks that use no outcome data of any kind.

1. **Independent implementation.** The primary estimator is a stratified weighted
   mean difference with weights n1*n0/(n1+n0). This recomputes it by a different
   algorithm — a within-stratum (fixed-effects) least-squares regression of the
   penalty on the side indicator — and requires the two to agree. The two routes
   share no code beyond loading the data. They are algebraically equivalent, so
   any disagreement is a bug in one of them.

2. **Simulation on the real design.** The real sites, sides, strata, proteins and
   distances are kept exactly as they are; only the outcome is replaced by
   simulated values with a known planted effect. This measures whether the
   estimator is unbiased and whether the protein-clustered bootstrap interval
   actually covers at 95% under this design's clustering — which a dozen large
   proteins could easily break.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from run_gate0 import combine, load_frame, stratum_statistics
from run_gate0_within_site import site_penalties

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCORES = ROOT / "data" / "processed" / "gate0" / "am_flank_scores.npz"
DEFAULT_OUT = ROOT / "results" / "verification_20260927"

PRIMARY_WINDOW = 10


def fixed_effects_estimate(sites: pd.DataFrame, outcome: str) -> float:
    """Independent route: within-stratum demeaning, then least squares on the side.

    For a binary regressor this is the textbook fixed-effects estimator,
    beta = sum(x~ * y~) / sum(x~^2) over stratum-demeaned values. It equals the
    weighted mean difference with weights n1*n0/(n1+n0), but gets there without
    ever forming a per-stratum difference.
    """
    frame = sites[["protein", "substitution", "is_cytosolic", outcome]].copy()
    frame["x"] = frame["is_cytosolic"].astype(float)
    keys = ["protein", "substitution"]
    grouped = frame.groupby(keys, sort=False)
    # A stratum with one side only contributes nothing after demeaning.
    frame["x_within"] = frame["x"] - grouped["x"].transform("mean")
    frame["y_within"] = frame[outcome] - grouped[outcome].transform("mean")
    denominator = float((frame["x_within"] ** 2).sum())
    if denominator == 0:
        return float("nan")
    return float((frame["x_within"] * frame["y_within"]).sum() / denominator)


def stratified_estimate(sites: pd.DataFrame, outcome: str) -> float:
    frame = sites.rename(columns={outcome: "score"})[
        ["protein", "substitution", "is_cytosolic", "score"]]
    return combine(stratum_statistics(frame))[0]


def simulate(sites: pd.DataFrame, rng: np.random.Generator, effect: float,
             decay_lambda: float | None, protein_sd: float, site_sd: float) -> np.ndarray:
    """Planted effect on the real design.

    The effect applies to cytosolic sites only, optionally damped by distance, on
    top of a protein-level random intercept and site noise. The protein intercept
    is what makes naive intervals too narrow, so it is the thing the coverage
    check has to survive.
    """
    n_proteins = int(sites["protein"].max()) + 1
    protein_effect = rng.normal(0.0, protein_sd, n_proteins)[sites["protein"].to_numpy()]
    damping = (1.0 if decay_lambda is None
               else np.exp(-sites["distance"].to_numpy() / decay_lambda))
    signal = effect * sites["is_cytosolic"].to_numpy().astype(float) * damping
    return signal + protein_effect + rng.normal(0.0, site_sd, len(sites))


def bootstrap_interval(sites: pd.DataFrame, outcome: str, n_proteins: int,
                       draws: int, rng: np.random.Generator) -> tuple[float, float]:
    frame = sites.rename(columns={outcome: "score"})[
        ["protein", "substitution", "is_cytosolic", "score"]]
    stats = stratum_statistics(frame)
    if stats["protein"].size == 0:
        return float("nan"), float("nan")
    samples = np.full(draws, np.nan)
    for i in range(draws):
        weights = np.bincount(rng.integers(0, n_proteins, n_proteins),
                              minlength=n_proteins).astype(float)
        samples[i], _ = combine(stats, weights[stats["protein"]])
    finite = samples[np.isfinite(samples)]
    if finite.size < 2:
        return float("nan"), float("nan")
    return float(np.percentile(finite, 2.5)), float(np.percentile(finite, 97.5))


def run(args: argparse.Namespace) -> dict:
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    frame, _ = load_frame(Path(args.scores))
    n_proteins = int(frame["protein"].max()) + 1
    sites = site_penalties(frame, set("KR"))
    sites = sites[sites["distance"] <= PRIMARY_WINDOW].reset_index(drop=True)

    # --- Check 1: two implementations on the real AlphaMissense penalties ---
    stratified = stratified_estimate(sites, "penalty")
    fixed_effects = fixed_effects_estimate(sites, "penalty")
    agreement = {
        "stratified_weighted_mean_difference": stratified,
        "within_stratum_fixed_effects": fixed_effects,
        "absolute_difference": abs(stratified - fixed_effects),
        "agree_to_1e-10": bool(abs(stratified - fixed_effects) < 1e-10),
        "n_sites": int(len(sites)),
    }

    # --- Check 2: simulation on the real design ---
    scenarios = {
        "null_no_effect": {"effect": 0.0, "decay_lambda": None},
        "flat_effect_0.05": {"effect": 0.05, "decay_lambda": None},
        "decaying_effect_0.10_lambda8": {"effect": 0.10, "decay_lambda": 8.0},
    }
    simulation: dict = {}
    for name, spec in scenarios.items():
        rng = np.random.default_rng(args.seed)
        estimates, covered, rejected = [], 0, 0
        for _ in range(args.replicates):
            sites["simulated"] = simulate(sites, rng, spec["effect"], spec["decay_lambda"],
                                          args.protein_sd, args.site_sd)
            estimate = stratified_estimate(sites, "simulated")
            estimates.append(estimate)
            low, high = bootstrap_interval(sites, "simulated", n_proteins,
                                           args.boot_draws, rng)
            # The estimand under damping is the design-weighted mean of the
            # damped signal over cytosolic sites, not the raw effect size.
            if spec["decay_lambda"] is None:
                truth = spec["effect"]
            else:
                cytosolic = sites["is_cytosolic"].to_numpy()
                truth = float(spec["effect"] * np.exp(
                    -sites["distance"].to_numpy()[cytosolic] / spec["decay_lambda"]).mean())
            if np.isfinite(low) and low <= truth <= high:
                covered += 1
            if np.isfinite(low) and not (low <= 0.0 <= high):
                rejected += 1
        estimates = np.array(estimates)
        simulation[name] = {
            "planted_effect": spec["effect"],
            "decay_lambda": spec["decay_lambda"],
            "target_estimand": truth,
            "mean_estimate": float(estimates.mean()),
            "bias": float(estimates.mean() - truth),
            "sd_of_estimates": float(estimates.std(ddof=1)),
            "coverage_of_95pct_interval": covered / args.replicates,
            "rejection_rate_against_zero": rejected / args.replicates,
            "replicates": args.replicates,
        }

    summary = {
        "analysis": "Estimator verification: independent implementation and simulation",
        "clinical_labels_used": False,
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "design_source": "real K/R sites within 10 residues of a boundary; only the "
                         "outcome is simulated",
        "simulation_settings": {"protein_sd": args.protein_sd, "site_sd": args.site_sd,
                                "bootstrap_draws": args.boot_draws},
        "implementation_agreement": agreement,
        "simulation": simulation,
    }
    (out_dir / "verification_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print("=== implementation agreement ===")
    print(f"  stratified weighted mean : {stratified:.12f}")
    print(f"  fixed-effects regression : {fixed_effects:.12f}")
    print(f"  agree to 1e-10           : {agreement['agree_to_1e-10']}")
    print("\n=== simulation on the real design ===")
    for name, value in simulation.items():
        print(f"  {name}")
        print(f"    target {value['target_estimand']:+.5f}  mean estimate "
              f"{value['mean_estimate']:+.5f}  bias {value['bias']:+.5f}")
        print(f"    95% interval coverage {value['coverage_of_95pct_interval']:.3f}   "
              f"rejects zero {value['rejection_rate_against_zero']:.3f}")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", default=str(DEFAULT_SCORES))
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--boot-draws", type=int, default=300)
    parser.add_argument("--protein-sd", type=float, default=0.08)
    parser.add_argument("--site-sd", type=float, default=0.25)
    parser.add_argument("--seed", type=int, default=20260927)
    run(parser.parse_args())


if __name__ == "__main__":
    main()

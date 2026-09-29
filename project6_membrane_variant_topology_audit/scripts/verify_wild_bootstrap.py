"""Compare inference methods on the real design: percentile bootstrap vs wild cluster bootstrap-t.

The audit currently uses a protein-clustered percentile bootstrap, which the
simulation showed rejects at 6.0% against a nominal 5%. That is typical with a
few dozen clusters, but the restricted wild cluster bootstrap-t usually controls
size better, because it imposes the null when generating replicates and
studentises the statistic.

Implementation, for the fixed-effects form of the estimator:

    y~ and x~ are stratum-demeaned outcome and side indicator
    beta = sum(x~ y~) / sum(x~^2)
    cluster-robust se from per-protein score sums

Restricted (null-imposed) replicates take the residual under H0: beta = 0, which
is simply y~, and flip it by a protein-level Rademacher weight. The bootstrap t
distribution is then compared with the observed t.

Uses no outcome data: the outcome is simulated on the real design.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from run_gate0 import load_frame
from run_gate0_within_site import site_penalties
from verify_estimator import bootstrap_interval, simulate, stratified_estimate

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCORES = ROOT / "data" / "processed" / "gate0" / "am_full_flank_scores.npz"
DEFAULT_OUT = ROOT / "results" / "verification_20260927"
PRIMARY_WINDOW = 10


def demean(sites: pd.DataFrame, outcome: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    frame = sites[["protein", "substitution", outcome]].copy()
    frame["x"] = sites["is_cytosolic"].to_numpy().astype(float)
    grouped = frame.groupby(["protein", "substitution"], sort=False)
    x_within = frame["x"].to_numpy() - grouped["x"].transform("mean").to_numpy()
    y_within = frame[outcome].to_numpy() - grouped[outcome].transform("mean").to_numpy()
    return x_within, y_within, frame["protein"].to_numpy()


def fixed_effects_t(x: np.ndarray, y: np.ndarray, cluster: np.ndarray,
                    n_clusters: int) -> tuple[float, float]:
    denominator = float((x * x).sum())
    if denominator == 0:
        return float("nan"), float("nan")
    beta = float((x * y).sum() / denominator)
    residual = y - beta * x
    scores = np.bincount(cluster, weights=x * residual, minlength=n_clusters)
    variance = float((scores ** 2).sum()) / (denominator ** 2)
    if variance <= 0:
        return beta, float("nan")
    return beta, beta / np.sqrt(variance)


def wild_cluster_p_value(x: np.ndarray, y: np.ndarray, cluster: np.ndarray,
                         n_clusters: int, draws: int,
                         rng: np.random.Generator) -> tuple[float, float]:
    """Restricted wild cluster bootstrap-t p-value against beta = 0."""
    beta, t_observed = fixed_effects_t(x, y, cluster, n_clusters)
    if not np.isfinite(t_observed):
        return beta, float("nan")
    # Under H0 the fitted value is zero, so the restricted residual is y itself.
    extreme = 0
    for _ in range(draws):
        flips = rng.choice((-1.0, 1.0), size=n_clusters)
        _, t_star = fixed_effects_t(x, y * flips[cluster], cluster, n_clusters)
        if np.isfinite(t_star) and abs(t_star) >= abs(t_observed):
            extreme += 1
    return beta, (extreme + 1) / (draws + 1)


def run(args: argparse.Namespace) -> dict:
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    frame, _ = load_frame(Path(args.scores))
    n_proteins = int(frame["protein"].max()) + 1
    sites = site_penalties(frame, set("KR"))
    sites = sites[sites["distance"] <= PRIMARY_WINDOW].reset_index(drop=True)

    # Observed result under both inference methods.
    x, y, cluster = demean(sites, "penalty")
    beta, t_observed = fixed_effects_t(x, y, cluster, n_proteins)
    rng = np.random.default_rng(args.seed)
    _, wild_p = wild_cluster_p_value(x, y, cluster, n_proteins, args.wild_draws, rng)
    low, high = bootstrap_interval(sites, "penalty", n_proteins, args.boot_draws, rng)
    observed = {
        "estimate": beta,
        "t_statistic": t_observed,
        "percentile_bootstrap_ci95": [low, high],
        "wild_cluster_bootstrap_t_p_value": wild_p,
        "n_sites": int(len(sites)),
        "n_clusters_with_data": int(np.unique(cluster).size),
    }

    # Size and power on the real design.
    scenarios = {"null_no_effect": 0.0, "flat_effect_0.05": 0.05}
    simulation: dict = {}
    for name, effect in scenarios.items():
        rng = np.random.default_rng(args.seed)
        reject_percentile = reject_wild = 0
        for _ in range(args.replicates):
            sites["simulated"] = simulate(sites, rng, effect, None,
                                          args.protein_sd, args.site_sd)
            low, high = bootstrap_interval(sites, "simulated", n_proteins,
                                           args.boot_draws, rng)
            if np.isfinite(low) and not (low <= 0.0 <= high):
                reject_percentile += 1
            xs, ys, cs = demean(sites, "simulated")
            _, p_value = wild_cluster_p_value(xs, ys, cs, n_proteins,
                                              args.wild_draws, rng)
            if np.isfinite(p_value) and p_value < 0.05:
                reject_wild += 1
        simulation[name] = {
            "planted_effect": effect,
            "percentile_bootstrap_rejection_rate": reject_percentile / args.replicates,
            "wild_cluster_bootstrap_t_rejection_rate": reject_wild / args.replicates,
            "replicates": args.replicates,
        }

    summary = {
        "analysis": "Inference comparison: percentile bootstrap vs restricted wild cluster bootstrap-t",
        "clinical_labels_used": False,
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "source_scores": str(args.scores),
        "nominal_level": 0.05,
        "observed": observed,
        "simulation": simulation,
    }
    (out_dir / "wild_bootstrap_comparison.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", default=str(DEFAULT_SCORES))
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--boot-draws", type=int, default=300)
    parser.add_argument("--wild-draws", type=int, default=399)
    parser.add_argument("--protein-sd", type=float, default=0.08)
    parser.add_argument("--site-sd", type=float, default=0.25)
    parser.add_argument("--seed", type=int, default=20260927)
    run(parser.parse_args())


if __name__ == "__main__":
    main()

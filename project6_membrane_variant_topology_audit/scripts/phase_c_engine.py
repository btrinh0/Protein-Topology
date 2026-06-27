"""Shared estimator for Phase C: one implementation, used by the analysis and the simulation.

Fix F3 of the frozen correction plan. The original analysis formed per-(gene,
reference residue) stratum differences and took an efficiency-weighted mean of
them, while the design simulation demeaned site-level values within gene and
regressed on the side indicator. Those are different estimators, so the reported
size and power did not describe the procedure that produced the result. Both now
call the functions here.

Also carries fix F2 — the weighting scheme reaches the point estimate, the
percentile bootstrap and the wild bootstrap alike — and fix F4, an interval
obtained by inverting the preregistered test rather than from the percentile
bootstrap the frozen protocol demotes.
"""

from __future__ import annotations

import math
from collections import defaultdict

import numpy as np

WEBB = np.array([-math.sqrt(1.5), -1.0, -math.sqrt(0.5),
                 math.sqrt(0.5), 1.0, math.sqrt(1.5)])


def build_strata(rows: list[dict], equal_gene_weight: bool = False) -> dict:
    """Per-(gene, reference residue) stratum statistics, with fix F5 counts.

    A stratum with observations on one side only carries no information about a
    side difference. Those rows are counted separately rather than folded into the
    reported sample size.
    """
    buckets: dict = defaultdict(lambda: {"cyto": [], "outer": []})
    for row in rows:
        key = (row["gene"], row["ref"])
        buckets[key]["cyto" if row["cytosolic"] else "outer"].append(row["penalty"])

    genes, weights, diffs, n_cyto, n_outer = [], [], [], [], []
    dropped_sites = 0
    for (gene, _), arms in sorted(buckets.items()):
        n1, n0 = len(arms["cyto"]), len(arms["outer"])
        if n1 == 0 or n0 == 0:
            dropped_sites += n1 + n0
            continue
        genes.append(gene)
        weights.append(n1 * n0 / (n1 + n0))
        diffs.append(float(np.mean(arms["cyto"]) - np.mean(arms["outer"])))
        n_cyto.append(n1)
        n_outer.append(n0)

    names = sorted(set(genes))
    index = {g: i for i, g in enumerate(names)}
    cluster = np.array([index[g] for g in genes], dtype=int)
    weight = np.array(weights, dtype=float)

    if equal_gene_weight and len(names):
        # Each gene contributes the same total weight, spread over its strata.
        per_gene = np.bincount(cluster, weights=weight, minlength=len(names))
        scale = np.where(per_gene > 0, 1.0 / np.maximum(per_gene, 1e-12), 0.0)
        weight = weight * scale[cluster]

    return {
        "cluster": cluster, "weight": weight, "diff": np.array(diffs, dtype=float),
        "genes": names, "n_strata": len(diffs),
        "n_cytosolic_sites": int(sum(n_cyto)), "n_outer_sites": int(sum(n_outer)),
        "n_contributing_sites": int(sum(n_cyto) + sum(n_outer)),
        "n_dropped_single_sided_sites": dropped_sites,
    }


def statistic(strata: dict, values: np.ndarray | None = None) -> tuple[float, float]:
    """Weighted mean side difference and its cluster-robust t."""
    weight = strata["weight"]
    diff = strata["diff"] if values is None else values
    total = weight.sum()
    if total <= 0 or diff.size == 0:
        return float("nan"), float("nan")
    beta = float((weight * diff).sum() / total)
    contribution = np.bincount(strata["cluster"], weights=weight * (diff - beta),
                               minlength=len(strata["genes"])) / total
    variance = float((contribution ** 2).sum())
    return beta, (beta / math.sqrt(variance) if variance > 0 else float("nan"))


def wild_p_value(strata: dict, null_value: float, draws: int, weights_kind: str,
                 seed: int) -> float:
    """Restricted wild cluster bootstrap-t against H0: beta = null_value."""
    n_clusters = len(strata["genes"])
    if n_clusters < 2 or strata["n_strata"] < 2:
        return float("nan")
    shifted = strata["diff"] - null_value
    beta, _ = statistic(strata, shifted)
    _, t_observed = statistic(strata, shifted)
    if not np.isfinite(t_observed):
        return float("nan")
    # Center once so every bootstrap draw is generated under the null.
    centred = shifted - beta
    rng = np.random.default_rng(seed)
    extreme = 0
    for _ in range(draws):
        flips = (rng.choice((-1.0, 1.0), size=n_clusters) if weights_kind == "rademacher"
                 else rng.choice(WEBB, size=n_clusters))
        _, t_star = statistic(strata, centred * flips[strata["cluster"]])
        if np.isfinite(t_star) and abs(t_star) >= abs(t_observed):
            extreme += 1
    return (extreme + 1) / (draws + 1)


def wild_interval(strata: dict, draws: int, weights_kind: str, seed: int,
                  level: float = 0.05, span: float = 4.0,
                  tolerance: float = 1e-3) -> list[float]:
    """Interval by inverting the wild bootstrap test: the null values it does not reject.

    Fix F4. The frozen protocol makes this test primary but the decision table asks
    for an interval, and the original analysis filled that gap with the percentile
    bootstrap the same protocol demotes.
    """
    beta, _ = statistic(strata)
    if not np.isfinite(beta):
        return [float("nan"), float("nan")]
    scale = float(np.std(strata["diff"], ddof=1)) if strata["n_strata"] > 1 else 1.0
    reach = max(span * max(scale, abs(beta), 1e-3), 1e-3)

    def rejected(value: float) -> bool:
        p = wild_p_value(strata, value, draws, weights_kind, seed)
        return np.isfinite(p) and p < level

    bounds = []
    for direction in (-1.0, 1.0):
        inside, outside = beta, beta + direction * reach
        if not rejected(outside):
            # Return the search boundary when inversion finds no rejection.
            bounds.append(float(outside))
            continue
        while abs(outside - inside) > tolerance:
            middle = (inside + outside) / 2.0
            if rejected(middle):
                outside = middle
            else:
                inside = middle
        bounds.append(float(inside))
    return [min(bounds), max(bounds)]


def percentile_interval(strata: dict, draws: int, seed: int) -> list[float]:
    """Gene-clustered percentile interval. Superseded by wild_interval; kept for comparison."""
    n_clusters = len(strata["genes"])
    if n_clusters < 2:
        return [float("nan"), float("nan")]
    rng = np.random.default_rng(seed)
    samples = np.full(draws, np.nan)
    for i in range(draws):
        counts = np.bincount(rng.integers(0, n_clusters, n_clusters),
                             minlength=n_clusters).astype(float)
        multiplier = counts[strata["cluster"]]
        total = (strata["weight"] * multiplier).sum()
        if total > 0:
            samples[i] = float((strata["weight"] * multiplier * strata["diff"]).sum() / total)
    finite = samples[np.isfinite(samples)]
    return ([float(np.percentile(finite, 2.5)), float(np.percentile(finite, 97.5))]
            if finite.size > 1 else [float("nan"), float("nan")])


def fit(rows: list[dict], label: str, wild_draws: int, boot_draws: int, seed: int,
        equal_gene_weight: bool = False, interval_draws: int | None = None) -> dict:
    """Full corrected fit: estimate, inverted-test interval, both weightings, honest counts."""
    strata = build_strata(rows, equal_gene_weight)
    if strata["n_strata"] < 2:
        return {"scope": label, "n_contributing_sites": strata["n_contributing_sites"],
                "n_strata": strata["n_strata"], "estimate": float("nan")}
    beta, t_stat = statistic(strata)
    result = {
        "scope": label,
        "estimate": beta,
        "t_statistic": t_stat,
        "n_contributing_sites": strata["n_contributing_sites"],
        "n_cytosolic_sites": strata["n_cytosolic_sites"],
        "n_outer_sites": strata["n_outer_sites"],
        "n_dropped_single_sided_sites": strata["n_dropped_single_sided_sites"],
        "n_strata": strata["n_strata"],
        "n_genes": len(strata["genes"]),
        "equal_gene_weight": equal_gene_weight,
        "percentile_ci95_superseded": percentile_interval(strata, boot_draws, seed),
    }
    for kind in ("webb", "rademacher"):
        result[f"wild_{kind}"] = {
            "p_value": wild_p_value(strata, 0.0, wild_draws, kind, seed),
            "ci95_inverted": wild_interval(strata, interval_draws or max(wild_draws // 10, 199),
                                           kind, seed),
        }
    return result

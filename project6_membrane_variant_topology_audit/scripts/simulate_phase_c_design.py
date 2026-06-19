"""Simulation on Phase C's own design, not the audit's.

The verification already on file was run on the AlphaMissense audit, which has
thousands of protein clusters. Phase C has sixteen genes, and twelve once GPCRs
are dropped. That is the regime where cluster-robust inference misbehaves, so the
size check has to be run there.

Compares three inference methods at the real gene sizes:

* percentile bootstrap over genes
* restricted wild cluster bootstrap-t, Rademacher weights (2-point)
* restricted wild cluster bootstrap-t, Webb weights (6-point)

Rademacher weights admit only 2^G distinct draws. At G = 12 that is 4,096, fewer
than the 9,999 replicates planned, so the bootstrap distribution is coarse and the
attainable p-values are granular. Webb's 6-point distribution is the standard fix
below about 12 clusters.

No measured data of any kind is used: gene sizes come from blinded counts, and
every outcome is simulated.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "results" / "verification_phase_c_20260927"

# Blinded K/R gain cell counts, window <=10: gene -> (cytosolic sites, outer sites)
GAIN_CELL = {
    "RHO": (92, 90), "CXCR4": (90, 98), "GPR68": (86, 92), "CCR5": (80, 90),
    "KCNH2": (66, 94), "VKORC1": (32, 36), "KCNJ2": (26, 28), "LDLR": (16, 14),
    "CYP2C9": (15, 0), "INSR": (0, 15), "MPL": (14, 18), "SGCB": (14, 18),
    "KCNE1": (13, 16), "CYP2C19": (12, 2), "SGCA": (6, 4), "CD86": (2, 4),
}
GPCR = {"CCR5", "CXCR4", "GPR68", "RHO"}
WEBB = np.array([-math.sqrt(1.5), -1.0, -math.sqrt(0.5),
                 math.sqrt(0.5), 1.0, math.sqrt(1.5)])


def build_design(cells: dict[str, tuple[int, int]]):
    """Site-level design: one row per site, carrying its gene and side."""
    gene_index, side, gene_names = [], [], []
    for index, (gene, (n_cytosolic, n_outer)) in enumerate(sorted(cells.items())):
        gene_names.append(gene)
        gene_index.extend([index] * (n_cytosolic + n_outer))
        side.extend([1.0] * n_cytosolic + [0.0] * n_outer)
    return np.array(gene_index), np.array(side), gene_names


def demean(values: np.ndarray, side: np.ndarray, gene: np.ndarray, n_genes: int):
    """Stratum-demean by gene, which is the Phase C stratifier."""
    counts = np.bincount(gene, minlength=n_genes)
    counts[counts == 0] = 1
    value_mean = np.bincount(gene, weights=values, minlength=n_genes) / counts
    side_mean = np.bincount(gene, weights=side, minlength=n_genes) / counts
    return values - value_mean[gene], side - side_mean[gene]


def estimate_t(values: np.ndarray, side: np.ndarray, gene: np.ndarray,
               n_genes: int) -> tuple[float, float]:
    y, x = demean(values, side, gene, n_genes)
    denominator = float((x * x).sum())
    if denominator == 0:
        return float("nan"), float("nan")
    beta = float((x * y).sum() / denominator)
    residual = y - beta * x
    scores = np.bincount(gene, weights=x * residual, minlength=n_genes)
    variance = float((scores ** 2).sum()) / denominator ** 2
    return beta, (beta / math.sqrt(variance) if variance > 0 else float("nan"))


def wild_p_value(values, side, gene, n_genes, draws, rng, weights) -> float:
    _, t_observed = estimate_t(values, side, gene, n_genes)
    if not np.isfinite(t_observed):
        return float("nan")
    y, x = demean(values, side, gene, n_genes)   # restricted: null residual is y
    extreme = 0
    for _ in range(draws):
        flips = (rng.choice((-1.0, 1.0), size=n_genes) if weights == "rademacher"
                 else rng.choice(WEBB, size=n_genes))
        _, t_star = estimate_t(y * flips[gene], side, gene, n_genes)
        if np.isfinite(t_star) and abs(t_star) >= abs(t_observed):
            extreme += 1
    return (extreme + 1) / (draws + 1)


def percentile_rejects(values, side, gene, n_genes, draws, rng) -> bool:
    beta, _ = estimate_t(values, side, gene, n_genes)
    samples = np.empty(draws)
    for i in range(draws):
        picked = rng.integers(0, n_genes, n_genes)
        mask = np.concatenate([np.flatnonzero(gene == g) for g in picked])
        relabel = np.concatenate([np.full((gene == g).sum(), i)
                                  for i, g in enumerate(picked)])
        samples[i], _ = estimate_t(values[mask], side[mask], relabel, n_genes)
    finite = samples[np.isfinite(samples)]
    if finite.size < 2:
        return False
    low, high = np.percentile(finite, [2.5, 97.5])
    return not (low <= 0.0 <= high)


def run_scenario(cells, effect, args, rng) -> dict:
    gene, side, names = build_design(cells)
    n_genes = len(names)
    counts = {"percentile": 0, "rademacher": 0, "webb": 0}
    for _ in range(args.replicates):
        gene_effect = rng.normal(0.0, args.gene_sd, n_genes)[gene]
        values = effect * side + gene_effect + rng.normal(0.0, args.site_sd, len(side))
        if percentile_rejects(values, side, gene, n_genes, args.boot_draws, rng):
            counts["percentile"] += 1
        for label in ("rademacher", "webb"):
            p = wild_p_value(values, side, gene, n_genes, args.wild_draws, rng, label)
            if np.isfinite(p) and p < 0.05:
                counts[label] += 1
    return {
        "planted_effect": effect,
        "n_genes": n_genes,
        "n_sites": int(len(side)),
        "distinct_rademacher_draws": 2 ** n_genes,
        "rejection_rates": {k: v / args.replicates for k, v in counts.items()},
    }


def run(args: argparse.Namespace) -> dict:
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    no_gpcr = {k: v for k, v in GAIN_CELL.items() if k not in GPCR}

    results = {}
    for scope, cells in (("all_16_genes", GAIN_CELL), ("no_gpcr_12_genes", no_gpcr)):
        results[scope] = {}
        for label, effect in (("null", 0.0), ("effect_0.30_sd", 0.30 * args.site_sd)):
            rng = np.random.default_rng(args.seed)
            results[scope][label] = run_scenario(cells, effect, args, rng)

    summary = {
        "analysis": "Size and power of Phase C inference on Phase C's own design",
        "uses_measured_data": False,
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "nominal_level": 0.05,
        "settings": {"gene_sd": args.gene_sd, "site_sd": args.site_sd,
                     "replicates": args.replicates, "wild_draws": args.wild_draws,
                     "boot_draws": args.boot_draws},
        "results": results,
    }
    (out_dir / "phase_c_inference_simulation.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    for scope, block in results.items():
        print(f"\n=== {scope} ===")
        for label, value in block.items():
            rates = value["rejection_rates"]
            print(f"  {label:<16} genes={value['n_genes']} sites={value['n_sites']} "
                  f"2^G={value['distinct_rademacher_draws']}")
            print(f"      percentile {rates['percentile']:.3f}   "
                  f"rademacher {rates['rademacher']:.3f}   webb {rates['webb']:.3f}")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    parser.add_argument("--replicates", type=int, default=400)
    parser.add_argument("--wild-draws", type=int, default=999)
    parser.add_argument("--boot-draws", type=int, default=499)
    parser.add_argument("--gene-sd", type=float, default=0.30)
    parser.add_argument("--site-sd", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=20260927)
    run(parser.parse_args())


if __name__ == "__main__":
    main()

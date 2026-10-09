"""Simulate Phase C using the realized design and shared analysis engine.

Uses the observed per-(gene, reference residue) cell sizes from the corrected run.
The analysis and simulation share the estimator implementation so the reported
size and power apply to the fitted procedure.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from phase_c_engine import build_strata, statistic, wild_p_value
from run_phase_c_corrected import assemble_variants, penalties_from
from run_phase_c import POSITIVE, PRIMARY_WINDOW


def realized_shape(rows: list[dict]) -> list[tuple[str, str, int, int]]:
    cells: dict = defaultdict(lambda: [0, 0])
    for row in rows:
        cells[(row["gene"], row["ref"])][0 if row["cytosolic"] else 1] += 1
    return [(g, r, n1, n0) for (g, r), (n1, n0) in cells.items() if n1 and n0]


def simulate_once(shape, effect, gene_sd, site_sd, rng) -> list[dict]:
    genes = sorted({g for g, _, _, _ in shape})
    offset = {g: rng.normal(0.0, gene_sd) for g in genes}
    rows = []
    for gene, ref, n1, n0 in shape:
        for _ in range(n1):
            rows.append({"gene": gene, "ref": ref, "cytosolic": True,
                         "penalty": effect + offset[gene] + rng.normal(0.0, site_sd)})
        for _ in range(n0):
            rows.append({"gene": gene, "ref": ref, "cytosolic": False,
                         "penalty": offset[gene] + rng.normal(0.0, site_sd)})
    return rows


def run(args) -> dict:
    lab, context, meta = assemble_variants(args)
    gain = [r for r in penalties_from(lab, context, "gain", POSITIVE, "positive_charge_gain")
            if r["distance"] <= PRIMARY_WINDOW]
    shape = realized_shape(gain)
    observed = build_strata(gain)
    site_sd = float(np.std([r["penalty"] for r in gain], ddof=1))

    results = {}
    for label, effect in (("null", 0.0), ("effect_0.30_sd", 0.30 * site_sd),
                          ("effect_0.50_sd", 0.50 * site_sd)):
        counts = {"webb": 0, "rademacher": 0}
        estimates = []
        rng = np.random.default_rng(args.seed)
        for _ in range(args.replicates):
            rows = simulate_once(shape, effect, args.gene_sd * site_sd, site_sd, rng)
            strata = build_strata(rows)
            estimates.append(statistic(strata)[0])
            for kind in ("webb", "rademacher"):
                p = wild_p_value(strata, 0.0, args.wild_draws, kind, args.seed)
                if np.isfinite(p) and p < 0.05:
                    counts[kind] += 1
        results[label] = {
            "planted_effect_sd_units": effect / site_sd,
            "mean_estimate": float(np.mean(estimates)),
            "rejection_rates": {k: v / args.replicates for k, v in counts.items()},
        }

    summary = {
        "analysis": "Size and power on the realized Phase C design, via the analysis engine",
        "fix": "F3 of the frozen correction plan",
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "realized_design": {
            "contributing_sites": observed["n_contributing_sites"],
            "cytosolic": observed["n_cytosolic_sites"], "outer": observed["n_outer_sites"],
            "strata": observed["n_strata"], "genes": len(observed["genes"]),
            "dropped_single_sided_sites": observed["n_dropped_single_sided_sites"],
        },
        "settings": {"gene_sd_in_site_sd_units": args.gene_sd,
                     "site_sd_observed": site_sd, "replicates": args.replicates,
                     "wild_draws": args.wild_draws},
        "results": results,
    }
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    (out / "realized_design_simulation.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", default=str(Path(__file__).resolve().parents[1] / "data/raw/topology_viewer.html"))
    parser.add_argument("--intramem", default=str(Path(__file__).resolve().parents[1] / "data/external/uniprot/human_intramem.tsv"))
    parser.add_argument("--outdir", default=str(Path(__file__).resolve().parents[1] / "results" / "phase_c_corrected_20260929"))
    parser.add_argument("--replicates", type=int, default=400)
    parser.add_argument("--wild-draws", type=int, default=399)
    parser.add_argument("--gene-sd", type=float, default=0.30)
    parser.add_argument("--seed", type=int, default=20260929)
    run(parser.parse_args())


if __name__ == "__main__":
    main()

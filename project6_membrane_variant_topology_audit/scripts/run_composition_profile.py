"""Positive control: are the charge-topology rules visible in the census sequences?

Label-free and model-free. This reads only the amino-acid sequences and the
topology annotations, and measures how often K/R and D/E actually occur at each
distance from a transmembrane boundary, on each side.

If the positive-inside rule is real in this proteome, cytosolic flanks should be
K/R-enriched relative to outer flanks, and the enrichment should fade with
distance. That gives an evolutionary decay length to compare against the decay
length in AlphaMissense's scores: a predictor that learned the rule from sequence
should show a similar profile.

Nothing here depends on a predictor or a clinical label, so it is the reference
both of those are measured against.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from extract_am_flank_scores import build_flank_map
from run_gate1 import read_topology

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TOPOLOGY = ROOT / "data" / "raw" / "topology_viewer.html"
DEFAULT_OUT = ROOT / "results" / "composition_20260927"

MAX_DISTANCE = 60
DISTANCE_BINS = [(1, 5), (6, 10), (11, 15), (16, 30), (31, 60)]
GROUPS = {"positive_KR": set("KR"), "negative_DE": set("DE")}


def collect(proteins: dict) -> dict[str, np.ndarray]:
    """One row per flank residue: protein index, side, distance, and residue identity."""
    protein_idx, sides, distances, residues = [], [], [], []
    for index, (accession, record) in enumerate(sorted(proteins.items())):
        sequence = record["sequence"]
        for position, (side, distance, _tmd) in build_flank_map(record, MAX_DISTANCE).items():
            if position > len(sequence):
                continue
            protein_idx.append(index)
            sides.append(side)
            distances.append(distance)
            residues.append(sequence[position - 1])
    return {
        "protein": np.array(protein_idx, dtype=np.int32),
        "side": np.array(sides, dtype=np.int8),
        "distance": np.array(distances, dtype=np.int16),
        "residue": np.array(residues, dtype="<U1"),
    }


def enrichment(data: dict, member: np.ndarray, mask: np.ndarray,
               n_proteins: int, draws: int, seed: int) -> dict:
    """Cytosolic minus outer frequency of a residue group, protein-clustered."""
    protein = data["protein"][mask]
    cytosolic = data["side"][mask] == 0
    hit = member[mask]
    if cytosolic.sum() == 0 or (~cytosolic).sum() == 0:
        return {}

    # Per-protein sufficient statistics keep the bootstrap cheap.
    counts = np.zeros((n_proteins, 4))
    np.add.at(counts, protein, np.column_stack([
        cytosolic & hit, cytosolic, (~cytosolic) & hit, ~cytosolic]).astype(float))

    def point(weights: np.ndarray | None = None) -> float:
        block = counts if weights is None else counts * weights[:, None]
        inside, inside_n, outside, outside_n = block.sum(axis=0)
        if inside_n == 0 or outside_n == 0:
            return float("nan")
        return inside / inside_n - outside / outside_n

    rng = np.random.default_rng(seed)
    samples = np.full(draws, np.nan)
    for i in range(draws):
        weights = np.bincount(rng.integers(0, n_proteins, n_proteins),
                              minlength=n_proteins).astype(float)
        samples[i] = point(weights)
    finite = samples[np.isfinite(samples)]

    inside, inside_n, outside, outside_n = counts.sum(axis=0)
    return {
        "cytosolic_frequency": float(inside / inside_n),
        "outer_frequency": float(outside / outside_n),
        "difference": float(inside / inside_n - outside / outside_n),
        "ci95": ([float(np.percentile(finite, 2.5)), float(np.percentile(finite, 97.5))]
                 if finite.size > 1 else [None, None]),
        "n_cytosolic_residues": int(inside_n),
        "n_outer_residues": int(outside_n),
    }


def fit_decay(distances: list[float], values: list[float]) -> dict:
    """Least-squares exponential A*exp(-d/lambda) on positive values, via log-linear fit."""
    pairs = [(d, v) for d, v in zip(distances, values) if v > 0]
    if len(pairs) < 3:
        return {"lambda_residues": None, "note": "too few positive points to fit"}
    x = np.array([d for d, _ in pairs], dtype=float)
    y = np.log(np.array([v for _, v in pairs], dtype=float))
    slope, intercept = np.polyfit(x, y, 1)
    if slope >= 0:
        return {"lambda_residues": None, "note": "no decay: fitted slope is non-negative"}
    predicted = slope * x + intercept
    ss_res = float(((y - predicted) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    return {
        "lambda_residues": round(float(-1.0 / slope), 2),
        "amplitude_at_zero": round(float(math.exp(intercept)), 5),
        "log_scale_r_squared": round(1 - ss_res / ss_tot, 3) if ss_tot > 0 else None,
    }


def run(args: argparse.Namespace) -> dict:
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    proteins = read_topology(Path(args.topology))
    data = collect(proteins)
    n_proteins = len(proteins)

    results: dict = {"by_distance_bin": {}, "by_residue_position": {}, "decay_fit": {}}
    for name, group in GROUPS.items():
        member = np.isin(data["residue"], list(group))
        results["by_distance_bin"][name] = {}
        for low, high in DISTANCE_BINS:
            mask = (data["distance"] >= low) & (data["distance"] <= high)
            results["by_distance_bin"][name][f"{low}-{high}"] = enrichment(
                data, member, mask, n_proteins, args.draws, args.seed)

        per_position = {}
        for distance in range(1, MAX_DISTANCE + 1):
            mask = data["distance"] == distance
            value = enrichment(data, member, mask, n_proteins, 0, args.seed)
            if value:
                per_position[str(distance)] = {
                    "difference": value["difference"],
                    "cytosolic_frequency": value["cytosolic_frequency"],
                    "outer_frequency": value["outer_frequency"],
                }
        results["by_residue_position"][name] = per_position
        results["decay_fit"][name] = fit_decay(
            [float(d) for d in per_position],
            [v["difference"] for v in per_position.values()])

    summary = {
        "analysis": "Amino-acid composition of TMD flanks by side and distance",
        "clinical_labels_used": False,
        "predictor_scores_used": False,
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "sign_convention": "Positive difference = the residue group is more frequent on the "
                           "cytosolic side",
        "topology_proteins": n_proteins,
        "flank_residues_examined": int(data["residue"].size),
        "results": results,
    }
    (out_dir / "composition_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    for name in GROUPS:
        print(f"\n=== {name}: cytosolic minus outer frequency ===")
        for key, value in results["by_distance_bin"][name].items():
            print(f"  {key:>6}: {value['difference']:+.4f} "
                  f"[{value['ci95'][0]:+.4f},{value['ci95'][1]:+.4f}]   "
                  f"cyto {value['cytosolic_frequency']:.4f} vs outer {value['outer_frequency']:.4f}")
        print(f"  decay fit: {results['decay_fit'][name]}")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", default=str(DEFAULT_TOPOLOGY))
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    parser.add_argument("--draws", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260927)
    run(parser.parse_args())


if __name__ == "__main__":
    main()

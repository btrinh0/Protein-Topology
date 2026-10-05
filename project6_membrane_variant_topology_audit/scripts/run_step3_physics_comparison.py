"""Step 3: compare the cytosolic K/R gain profile of each source against ER insertion physics.

Frozen design: docs/STEP3_PHYSICAL_ANCHOR_PREREG.md
SHA-256 40a411e3968bb1c0f0ad1365bc4cf47878411a82894371c1814e75aa701021af

The physics (Lerch-Bader et al. 2008) measures the free-energy cost or benefit of
adding a positive charge at a given distance on the CYTOSOLIC side. So the matching
model quantity is the cytosolic-side within-site gain penalty by distance, not the
side difference used elsewhere in this project.

Both raw quantities describe cost: negative insertion ddG and negative model
penalty each indicate favourability. Negate BOTH for comparisons on the
favourability axis. These descriptive comparisons do not establish model failure.

Label-free.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from run_gate0_within_site import site_gain_penalties
from run_step2_calibrated import build_scored_frame

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "results" / "step3_20260930"
PHYSICS = ROOT / "data" / "external" / "physics" / "lerch_bader_2008_table1.json"
BINS = [(1, 5), (6, 10), (11, 15), (16, 30)]
NEAR, FAR = [(1, 5), (6, 10)], [(11, 15), (16, 30)]
PREREG = ROOT / "docs" / "STEP3_PHYSICAL_ANCHOR_PREREG.md"
PREREG_SHA256 = "40a411e3968bb1c0f0ad1365bc4cf47878411a82894371c1814e75aa701021af"


def clustered_difference(near_values: np.ndarray, near_proteins: np.ndarray,
                         far_values: np.ndarray, far_proteins: np.ndarray,
                         draws: int, seed: int) -> dict:
    """Resample whole proteins jointly, preserving covariance between arms."""
    if draws < 1:
        raise ValueError("Bootstrap draws must be positive")
    if not near_values.size or not far_values.size:
        return {"difference": None, "ci95": [None, None], "n_proteins": 0}
    names = sorted(set(near_proteins.tolist()) | set(far_proteins.tolist()))
    index = {p: i for i, p in enumerate(names)}
    near_cluster = np.array([index[p] for p in near_proteins])
    far_cluster = np.array([index[p] for p in far_proteins])
    near_sum = np.bincount(near_cluster, weights=near_values, minlength=len(names))
    far_sum = np.bincount(far_cluster, weights=far_values, minlength=len(names))
    near_n = np.bincount(near_cluster, minlength=len(names))
    far_n = np.bincount(far_cluster, minlength=len(names))
    rng = np.random.default_rng(seed)
    samples = []
    for _ in range(draws):
        counts = np.bincount(rng.integers(0, len(names), len(names)), minlength=len(names))
        n_near, n_far = counts @ near_n, counts @ far_n
        if n_near > 0 and n_far > 0:
            samples.append(float(counts @ near_sum / n_near - counts @ far_sum / n_far))
    return {"difference": float(near_values.mean() - far_values.mean()),
            "ci95": [float(np.percentile(samples, 2.5)), float(np.percentile(samples, 97.5))]
            if len(samples) > 1 else [None, None],
            "n_proteins": len(names), "valid_bootstrap_draws": len(samples)}


def clustered_mean(values: np.ndarray, proteins: np.ndarray, draws: int,
                   seed: int) -> dict:
    if draws < 1:
        raise ValueError("Bootstrap draws must be positive")
    if values.size == 0:
        return {"mean": None, "ci95": [None, None], "n": 0}
    names = sorted(set(proteins.tolist()))
    index = {p: i for i, p in enumerate(names)}
    cluster = np.array([index[p] for p in proteins])
    rng = np.random.default_rng(seed)
    samples = np.empty(draws)
    for i in range(draws):
        counts = np.bincount(rng.integers(0, len(names), len(names)),
                             minlength=len(names)).astype(float)
        weight = counts[cluster]
        total = weight.sum()
        samples[i] = float((weight * values).sum() / total) if total > 0 else np.nan
    finite = samples[np.isfinite(samples)]
    return {"mean": float(values.mean()),
            "ci95": [float(np.percentile(finite, 2.5)),
                     float(np.percentile(finite, 97.5))] if finite.size > 1 else [None, None],
            "n": int(values.size), "n_proteins": len(names)}


def cytosolic_profile(frame, column: str, draws: int, seed: int) -> dict:
    work = frame.dropna(subset=[column]).rename(columns={column: "score"})
    sites = site_gain_penalties(work, set("KR"))
    if sites.empty:
        return {}
    spread = float(sites["penalty"].std(ddof=1))
    if not np.isfinite(spread) or spread <= 0:
        spread = 1.0
    out: dict = {"penalty_sd": spread, "by_bin": {}}
    for side_name, want in (("cytosolic", True), ("outer", False)):
        block = {}
        subset_side = sites[sites["is_cytosolic"] == want]
        for low, high in BINS:
            chunk = subset_side[(subset_side["distance"] >= low)
                                & (subset_side["distance"] <= high)]
            stats = clustered_mean(chunk["penalty"].to_numpy() / spread,
                                   chunk["protein"].to_numpy(), draws, seed)
            block[f"{low}-{high}"] = stats
        near = subset_side[(subset_side["distance"] >= 1) & (subset_side["distance"] <= 10)]
        far = subset_side[(subset_side["distance"] >= 11) & (subset_side["distance"] <= 30)]
        near_stats = clustered_mean(near["penalty"].to_numpy() / spread,
                                    near["protein"].to_numpy(), draws, seed)
        far_stats = clustered_mean(far["penalty"].to_numpy() / spread,
                                   far["protein"].to_numpy(), draws, seed)
        block["near_minus_far"] = {
            "near": near_stats["mean"], "far": far_stats["mean"],
            **clustered_difference(near["penalty"].to_numpy() / spread,
                                   near["protein"].to_numpy(),
                                   far["penalty"].to_numpy() / spread,
                                   far["protein"].to_numpy(), draws, seed),
        }
        out["by_bin"][side_name] = block
    return out


def physics_profile() -> dict:
    data = json.loads(PHYSICS.read_text(encoding="utf-8"))
    points = [p for p in data["distance_series"]["points"]
              if p["n_lys"] == 1 and not p["censored"]]
    mapping = {5: "1-5", 8: "6-10", 13: "11-15", 19: "16-30", 24: "16-30"}
    buckets: dict = defaultdict(list)
    for p in points:
        label = mapping.get(p["distance"])
        if label:
            buckets[label].append(p["ddGapp"])
    # Physical favourability. Compare with NEGATED model penalties below.
    return {label: -float(np.mean(values)) for label, values in buckets.items()}


def compare_profiles(physics: dict, profiles: dict) -> dict:
    """Compare both sources as favourability, with correlations descriptive only."""
    labels = [f"{a}-{b}" for a, b in BINS]
    physics_vector = np.array([physics.get(label, np.nan) for label in labels], dtype=float)
    comparison = {}
    for source, profile in profiles.items():
        block = profile.get("by_bin", {}).get("cytosolic", {})
        model_vector = np.array([block.get(label, {}).get("mean")
                                 if block.get(label, {}).get("mean") is not None else np.nan
                                 for label in labels], dtype=float)
        favourability = -model_vector
        usable = np.isfinite(physics_vector) & np.isfinite(favourability)
        correlation = None
        if usable.sum() >= 3 and np.std(physics_vector[usable]) > 0 and np.std(favourability[usable]) > 0:
            correlation = float(np.corrcoef(physics_vector[usable], favourability[usable])[0, 1])
        sign_agreement = [
            {"bin": label, "physics_negated": float(physics_vector[i]),
             "model": float(model_vector[i]) if np.isfinite(model_vector[i]) else None,
             "model_favourability": float(favourability[i]) if np.isfinite(favourability[i]) else None,
             "same_sign": bool(np.sign(physics_vector[i]) == np.sign(favourability[i]))
             if usable[i] else None}
            for i, label in enumerate(labels) if np.isfinite(physics_vector[i])
        ]
        comparison[source] = {
            "comparison_axis": "favourability: negative ddG versus negative within-site penalty",
            "sign_agreement_by_bin": sign_agreement,
            "n_bins_same_sign": sum(1 for item in sign_agreement if item["same_sign"]),
            "descriptive_correlation": correlation,
            "near_minus_far": block.get("near_minus_far"),
        }
    return comparison


def run(args) -> dict:
    if hashlib.sha256(PREREG.read_bytes()).hexdigest() != PREREG_SHA256:
        raise ValueError("Frozen Step 3 preregistration fingerprint changed")
    frame, meta = build_scored_frame(args)
    sources = [c for c in ("alphamissense", "esm1b", "b1", "b2") if c in frame.columns]
    profiles = {s: cytosolic_profile(frame, s, args.boot_draws, args.seed) for s in sources}
    physics = physics_profile()

    comparison = compare_profiles(physics, profiles)
    labels = [f"{a}-{b}" for a, b in BINS]

    summary = {
        "analysis": "Step 3: cytosolic K/R gain profile against ER insertion physics",
        "prereg_sha256": PREREG_SHA256,
        "clinical_labels_used": False, "dms_effects_used": False,
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "physics_source": "Lerch-Bader et al. 2008, Table 1, 1-Lys series, censored points excluded",
        "physics_negated_by_bin": physics,
        "sign_note": "Comparisons use physical favourability (-ddG) and model favourability "
                     "(-within-site penalty). The model field preserves the raw penalty; "
                     "model_favourability is its negation. A within-site difference can be "
                     "negative even when absolute pathogenicity scores are non-negative. "
                     "Different scaffolds and endpoints prevent interpreting sign alone as model failure.",
        "provenance": {"script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                       "physics_input_sha256": hashlib.sha256(PHYSICS.read_bytes()).hexdigest(),
                       "bootstrap_draws": args.boot_draws},
        **meta,
        "model_profiles": profiles,
        "comparison": comparison,
    }
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    (out / "step3_physics_comparison.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    print("physics (negated ddG, 1-Lys):", {k: round(v, 3) for k, v in physics.items()})
    print("\n=== cytosolic K/R gain penalty, standardised (negative = adding charge looks favourable) ===")
    header = "  {:<14}" + "".join(f"{l:>16}" for l in labels) + "{:>18}"
    print(header.format("source", "near-far"))
    for source in sources:
        block = profiles[source].get("by_bin", {}).get("cytosolic", {})
        row = f"  {source:<14}"
        for l in labels:
            stats = block.get(l, {})
            m = stats.get("mean")
            row += f"{m:+.4f}".rjust(16) if m is not None else "n/a".rjust(16)
        nf = block.get("near_minus_far", {}).get("difference")
        row += (f"{nf:+.4f}".rjust(18) if nf is not None else "n/a".rjust(18))
        print(row)
    print("\n=== outer side, for reference ===")
    for source in sources:
        block = profiles[source].get("by_bin", {}).get("outer", {})
        row = f"  {source:<14}"
        for l in labels:
            m = block.get(l, {}).get("mean")
            row += f"{m:+.4f}".rjust(16) if m is not None else "n/a".rjust(16)
        print(row)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", default=str(ROOT / "data/raw/topology_viewer.html"))
    parser.add_argument("--scores", default=str(ROOT / "data/processed/gate0/am_full_flank_scores.npz"))
    parser.add_argument("--esm", default=str(ROOT / "data/processed/gate0/esm1b_flank_scores.npz"))
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    parser.add_argument("--split", choices=("sequence", "symbol"), default="sequence")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--boot-draws", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260930)
    run(parser.parse_args())


if __name__ == "__main__":
    main()

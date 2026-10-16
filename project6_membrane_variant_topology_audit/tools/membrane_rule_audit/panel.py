"""Reusable raw-score panel analysis using the earlier uncharged comparator.

This module packages the existing within-site estimator. It does not run
alternative comparators, transformed scales, subgroups, or calibration studies.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
NEUTRAL = "ACFGILMNPQSTVWY"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_panel(path: Path, window: int = 10) -> dict:
    """Load the project's NPZ schema and reject conflicting variant metadata."""
    with np.load(path, allow_pickle=False) as blob:
        alphabet = "".join(blob["aa_alphabet"].tolist())
        if alphabet != ALPHABET:
            raise ValueError(f"Unexpected archive alphabet: {alphabet}")
        keep = (blob["distance"] > 0) & (blob["distance"] <= window)
        arrays = {name: blob[name][keep] for name in
                  ("protein", "position", "side", "distance", "tmd_index", "ref_aa", "alt_aa")}
        score_key = "score" if "score" in blob else "am_pathogenicity"
        values = blob[score_key][keep].astype(float)
        accessions = blob["accessions"].copy()
    if not len(values):
        raise ValueError("Empty score archive within the selected window")
    if (np.any(arrays["position"] < 1) or np.any(arrays["position"] >= 2**32)
            or np.any(arrays["protein"] < 0) or np.any(arrays["protein"] >= len(accessions))):
        raise ValueError("Invalid protein/position")
    if np.any(~np.isin(arrays["side"], [0, 1])):
        raise ValueError("Unknown membrane side")
    if any(np.any((arrays[name] < 0) | (arrays[name] >= 20)) for name in ("ref_aa", "alt_aa")):
        raise ValueError("Invalid residue index")
    code = (arrays["protein"].astype(np.int64) << 32) + arrays["position"].astype(np.int64)
    keys, first, inverse = np.unique(code, return_index=True, return_inverse=True)
    panel = {name: arrays[name][first] for name in
             ("protein", "position", "side", "distance", "tmd_index", "ref_aa")}
    for name in ("side", "distance", "tmd_index", "ref_aa"):
        if np.any(arrays[name] != panel[name][inverse]):
            raise ValueError(f"Conflicting metadata at a site: {name}")
    if np.any(arrays["alt_aa"] == arrays["ref_aa"]):
        raise ValueError("Synonymous scores must not be treated as missense alternatives")
    cells = inverse * 20 + arrays["alt_aa"].astype(int)
    finite = np.isfinite(values)
    sums = np.bincount(cells[finite], weights=values[finite], minlength=len(keys)*20)
    counts = np.bincount(cells[finite], minlength=len(keys)*20)
    scores = np.divide(sums, counts, out=np.full(len(keys)*20, np.nan), where=counts > 0)
    repeated = finite & (counts[cells] > 1)
    if np.any(np.abs(values[repeated] - scores[cells[repeated]]) > 1e-12):
        raise ValueError("Conflicting scores for a duplicate variant")
    panel.update(scores=scores.reshape(-1, 20), keys=keys, accessions=accessions,
                 n_input_rows=int(len(values)),
                 n_duplicate_rows=int(finite.sum()-np.count_nonzero(counts)))
    return panel


def arm_mean(scores: np.ndarray, residues: str) -> tuple:
    selected = scores[:, [ALPHABET.index(a) for a in residues]]
    count = np.isfinite(selected).sum(axis=1)
    mean = np.divide(np.nansum(selected, axis=1), count,
                     out=np.full(len(scores), np.nan), where=count > 0)
    return mean, count


def gain_contrast(panel: dict) -> dict:
    """K/R introductions minus uncharged alternatives at the same native site."""
    introduced, n_introduced = arm_mean(panel["scores"], "KR")
    neutral, n_neutral = arm_mean(panel["scores"], NEUTRAL)
    eligible = np.isin(panel["ref_aa"], [ALPHABET.index(a) for a in NEUTRAL])
    eligible &= (n_introduced > 0) & (n_neutral > 0)
    return {"penalty": introduced-neutral, "eligible": eligible}


def stratum_parts(panel: dict, values: np.ndarray, eligible: np.ndarray) -> dict:
    """Harmonic-size-weighted contrasts in protein/reference-residue strata."""
    keep = eligible & np.isfinite(values)
    protein = panel["protein"][keep].astype(int)
    ref = panel["ref_aa"][keep].astype(int)
    side = panel["side"][keep].astype(int)
    if not len(protein):
        return {"numerator": np.array([]), "denominator": np.array([]),
                "n_sites": 0, "n_strata": 0, "n_contributing_sites": 0,
                "n_cytosolic": 0, "n_outer": 0, "penalty_sd": None}
    group = protein*20 + ref
    cyto = side == 0
    n1 = np.bincount(group[cyto], minlength=int(group.max())+1)
    n0 = np.bincount(group[~cyto], minlength=len(n1))
    y = values[keep]
    sum1 = np.bincount(group[cyto], weights=y[cyto], minlength=len(n1))
    sum0 = np.bincount(group[~cyto], weights=y[~cyto], minlength=len(n1))
    shared = (n1 > 0) & (n0 > 0)
    ids = np.flatnonzero(shared)
    weights = n1[shared]*n0[shared]/(n1[shared]+n0[shared])
    differences = sum1[shared]/n1[shared] - sum0[shared]/n0[shared]
    _, index = np.unique(ids//20, return_inverse=True)
    return {"numerator": np.bincount(index, weights=weights*differences),
            "denominator": np.bincount(index, weights=weights),
            "n_sites": int(keep.sum()), "n_strata": int(shared.sum()),
            "n_contributing_sites": int((n1[shared]+n0[shared]).sum()),
            "n_cytosolic": int(cyto.sum()), "n_outer": int((~cyto).sum()),
            "penalty_sd": float(np.std(y, ddof=1)) if len(y)>1 else None}


def summarize_parts(parts: dict, draws: int, seed: int) -> dict:
    numerator, denominator = parts["numerator"], parts["denominator"]
    n = len(denominator)
    summary = {k: v for k, v in parts.items() if k not in {"numerator", "denominator"}}
    summary.update(n_proteins=n, estimate_raw=None, ci95_raw=None)
    if not n or denominator.sum() <= 0:
        return summary
    point = float(numerator.sum()/denominator.sum())
    summary["estimate_raw"] = point
    if n < 2 or draws < 2:
        summary["interval_status"] = "No across-protein interval computed"
        return summary
    rng = np.random.default_rng(seed)
    samples = []
    for start in range(0, draws, 64):
        counts = rng.multinomial(n, np.full(n, 1/n), size=min(64, draws-start))
        samples.extend((counts@numerator/(counts@denominator)).tolist())
    summary["ci95_raw"] = np.percentile(samples, [2.5, 97.5]).tolist()
    summary["interval_status"] = "Percentile protein-resampling stability interval"
    return summary


def run_panel_audit(scores_path: Path, out_dir: Path, model: str,
                    score_direction: str, window: int = 10,
                    draws: int = 1999, seed: int = 20261003) -> dict:
    """Apply the fixed earlier raw-score comparison to a named score source."""
    if score_direction not in {"higher-is-damaging", "higher-is-tolerated"}:
        raise ValueError("Declare a supported score direction")
    if window < 1 or draws < 0:
        raise ValueError("Window must be positive and draws nonnegative")
    panel = load_panel(scores_path, window)
    if score_direction == "higher-is-tolerated":
        panel["scores"] = -panel["scores"]
    contrast = gain_contrast(panel)
    result = summarize_parts(stratum_parts(panel, contrast["penalty"], contrast["eligible"]), draws, seed)
    summary = {
        "tool": "membrane_rule_audit", "version": "0.1.0", "model": model,
        "scope": "Software demonstration of the existing raw-score comparison",
        "scores": str(Path(scores_path).resolve()), "input_sha256": sha256(scores_path),
        "implementation_sha256": sha256(Path(__file__)),
        "input_score_direction": score_direction,
        "output_score_direction": "higher-is-damaging",
        "rule": "positive_charge_gain", "comparator": NEUTRAL,
        "scale": "raw", "window_residues": window, "draws": draws, "seed": seed,
        "robustness_analyses_run": False, "clinical_labels_used": False,
        "sign_convention": "Cytosolic minus outer within-site charge penalty",
        "input_quality_control": {"n_rows_within_window": panel["n_input_rows"],
                                  "n_duplicate_rows": panel["n_duplicate_rows"],
                                  "n_archive_sites_within_window": len(panel["keys"])},
        "result": result,
    }
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir/"membrane_rule_panel.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, required=True, help="Precomputed flank-score NPZ")
    parser.add_argument("--model", required=True, help="Label for this score source")
    parser.add_argument("--score-direction", required=True,
                        choices=("higher-is-damaging", "higher-is-tolerated"))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--window", type=int, default=10)
    parser.add_argument("--draws", type=int, default=1999)
    parser.add_argument("--seed", type=int, default=20261003)
    args = parser.parse_args()
    summary = run_panel_audit(args.scores, args.out, args.model, args.score_direction,
                              args.window, args.draws, args.seed)
    print(json.dumps({"model": args.model, **summary["result"]}, indent=2))


if __name__ == "__main__":
    main()

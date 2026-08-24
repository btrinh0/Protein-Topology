"""Assess outcome-blind pool coverage, prior exposure and confirmation feasibility."""

from __future__ import annotations

import csv
import hashlib
import json
import statistics as st
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "beat_am_20261016" / "inventory"

A0 = ROOT / "results" / "route_a_20261012" / "a0_feasibility.json"
CMP = ROOT / "results" / "outperformance_research_20261012" / "ALL_COMPARISONS.csv"
PHASE_C = ROOT / "scripts" / "run_phase_c.py"

# Track prior exposure by accession because gene symbols differ between sources.
OPENED_FOUR = {  # Exposure was recorded in results/outperformance_research_20261012.
    "Q12809": "KCNH2 full-length trafficking (Oct 12 benchmark)",
    "P08100": "RHO/OPSD 2 endpoints (Oct 12 benchmark)",
    "P51787": "KCNQ1 4 endpoints (Oct 12 benchmark)",
    "P30518": "AVPR2 1 endpoint (Oct 12 benchmark)",
}
SERT = {"P31645": "SERT/SC6A4 trafficking (scripts/load_sert.py, docs/SERT_TRAFFICKING_*)"}


DELTA_GRID = [0.01, 0.02, 0.03, 0.04, 0.05]

# The SD scenarios come from paired correlations on the four exposed development proteins.
SD_FULL = {"optimistic": 0.017, "central": 0.025, "pessimistic": 0.036}
SD_FLANK = {"optimistic": 0.028, "central": 0.050, "pessimistic": 0.134}

N_REPS = 2000
N_BOOT = 1000
SEED = 20261016


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def phase_c_accessions() -> dict[str, str]:
    """Parse SCORE_SETS out of run_phase_c.py rather than trusting a prose list."""
    text = PHASE_C.read_text(encoding="utf-8")
    start = text.index("SCORE_SETS = {")
    end = text.index("}", start)
    block = text[start:end]
    out: dict[str, str] = {}
    for line in block.splitlines()[1:]:
        line = line.strip()
        if not line.startswith('"'):
            continue
        gene = line.split('"')[1]
        parts = [p.strip().strip('"') for p in line.split("(")[1].split(")")[0].split(",")]
        out[parts[1]] = gene
    return out




def build_pool(a0: dict, phase_c: dict[str, str]) -> list[dict]:
    """One row per accepted assay, with exclusion reasons attached."""
    rows = []
    for a in a0["per_assay"]:
        if not a.get("in_census") or not a.get("accepted"):
            continue
        acc = a["accession"]
        reasons = []
        if acc in phase_c:
            reasons.append(f"Phase C score set ({phase_c[acc]}); docs/PHASE_C_RESULTS_20260928.md")
        if acc in OPENED_FOUR:
            reasons.append(OPENED_FOUR[acc])
        if acc in SERT:
            reasons.append(SERT[acc])
        rows.append(
            {
                "dms_id": a["dms_id"],
                "protein": a["uniprot_entry"].split("_")[0],
                "accession": acc,
                "endpoint_type": endpoint_type(a["dms_id"]),
                "reference_agreement": a["reference_agreement"],
                "census_length": a.get("census_length"),
                "n_single": a["n_single"],
                "flank_sites_measured": a.get("flank_sites_measured_in_window", 0),
                "kr_gain_flank_single": a.get("kr_gain_flank_single", 0),
                "kr_gain_sites_eligible": a.get("kr_gain_sites_eligible", 0),
                "confirmation_eligible": not reasons,
                "exclusion_reasons": " | ".join(reasons),
            }
        )
    return sorted(rows, key=lambda r: (r["protein"], r["dms_id"]))


def endpoint_type(dms_id: str) -> str:
    low = dms_id.lower()
    for key, label in (
        ("abundance", "abundance"),
        ("expression", "abundance"),
        ("activity", "function"),
        ("function", "function"),
        ("trafficking", "trafficking"),
    ):
        if key in low:
            return label
    return "unspecified"


def collapse_to_proteins(pool: list[dict]) -> list[dict]:
    """The unit of evidence is the protein, so collapse assays within protein."""
    by: dict[str, list[dict]] = {}
    for r in pool:
        by.setdefault(r["protein"], []).append(r)
    out = []
    for prot, rs in sorted(by.items()):
        flank = max(r["flank_sites_measured"] for r in rs)
        out.append(
            {
                "protein": prot,
                "accession": rs[0]["accession"],
                "n_assays": len(rs),
                "endpoint_types": ",".join(sorted({r["endpoint_type"] for r in rs})),
                "n_single_max": max(r["n_single"] for r in rs),
                "flank_sites_measured_max": flank,
                "kr_gain_flank_max": max(r["kr_gain_flank_single"] for r in rs),
                "flank_tier": flank_tier(flank),
                "confirmation_eligible": all(r["confirmation_eligible"] for r in rs),
                "exclusion_reasons": " | ".join(
                    sorted({r["exclusion_reasons"] for r in rs if r["exclusion_reasons"]})
                ),
            }
        )
    return out


def flank_tier(n_flank_sites: int) -> str:
    """How usable a within-protein flank Spearman is at this coverage."""
    if n_flank_sites == 0:
        return "none"
    if n_flank_sites < 25:
        return "token"      # Correlations from this small sample would be dominated by sampling noise.
    return "substantial"




def exact_signflip_p(diffs: np.ndarray, signs: np.ndarray) -> np.ndarray:
    """Two-sided exact sign-flip p for each column of `diffs`.

    diffs: (n, R). signs: (M, n) of +-1 with M = 2**n (all flips).
    Smallest attainable value is 2 / 2**n, matching the prompt's arithmetic.
    """
    n, _ = diffs.shape
    perm_means = (signs @ diffs) / n
    observed = diffs.mean(axis=0)
    hits = (np.abs(perm_means) >= np.abs(observed) - 1e-15).sum(axis=0)
    return hits / signs.shape[0]


def sign_matrix(n: int) -> np.ndarray:
    grid = np.arange(1 << n, dtype=np.uint32)[:, None]
    bits = (grid >> np.arange(n, dtype=np.uint32)[None, :]) & 1
    return np.where(bits == 1, 1.0, -1.0)


def power_at_d(n: int, d: float, rng: np.random.Generator,
               signs: np.ndarray) -> dict:
    """Power at standardized effect d = delta / sd_between_proteins.

    All three win conditions are invariant to a positive rescaling of the
    paired differences -- the sign-flip p compares |perm mean| to |observed
    mean|, the percentile bootstrap bound only has its sign tested, and
    mean > 0 is scale-free. So power is a function of (n, d) alone, and the
    minimum detectable improvement for any SD is just d* * SD. That is what
    makes this grid cheap enough to run exactly rather than by approximation.
    """
    diffs = rng.normal(d, 1.0, size=(n, N_REPS))

    mean_pos = diffs.mean(axis=0) > 0
    p_ok = exact_signflip_p(diffs, signs) < 0.05

    # Accumulate resamples iteratively to bound memory at B by R rather than B by n by R.
    boot = np.empty((N_BOOT, N_REPS))
    for b in range(N_BOOT):
        idx = rng.integers(0, n, size=(n, N_REPS))
        boot[b] = np.take_along_axis(diffs, idx, axis=0).mean(axis=0)
    ci_ok = np.quantile(boot, 0.025, axis=0) > 0

    win = mean_pos & p_ok & ci_ok
    return {
        "n_proteins": n,
        "standardized_effect_d": round(d, 4),
        "min_attainable_exact_p": 2 / (1 << n),
        "power_all_three_conditions": float(win.mean()),
        "power_signflip_only": float(p_ok.mean()),
        "power_bootstrap_only": float(ci_ok.mean()),
    }


def d_star(n: int, rng: np.random.Generator, signs: np.ndarray,
           target: float = 0.80) -> float | None:
    """Standardized effect reaching `target` power, by bisection.

    Power is monotone in d, so bisect on [0, 6] to 0.01 resolution.
    Returns None when a decisive exact p is unreachable at this n.
    """
    if 2 / (1 << n) >= 0.05:
        return None
    lo, hi = 0.0, 6.0
    if power_at_d(n, hi, rng, signs)["power_all_three_conditions"] < target:
        return None
    while hi - lo > 0.01:
        mid = (lo + hi) / 2
        if power_at_d(n, mid, rng, signs)["power_all_three_conditions"] >= target:
            hi = mid
        else:
            lo = mid
    return round(hi, 3)




def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    a0 = json.loads(A0.read_text(encoding="utf-8"))
    phase_c = phase_c_accessions()

    pool = build_pool(a0, phase_c)
    proteins = collapse_to_proteins(pool)

    write_csv(OUT / "pool_table_assays.csv", pool)
    write_csv(OUT / "pool_table_proteins.csv", proteins)

    eligible = [p for p in proteins if p["confirmation_eligible"]]
    n_elig = len(eligible)
    flank_sub = [p for p in eligible if p["flank_tier"] == "substantial"]
    flank_any = [p for p in eligible if p["flank_tier"] != "none"]

    rng = np.random.default_rng(SEED)
    signs_cache: dict[int, np.ndarray] = {}

    def signs_for(n: int) -> np.ndarray:
        if n not in signs_cache:
            signs_cache[n] = sign_matrix(n)
        return signs_cache[n]


    ns = list(range(3, n_elig + 1))
    curve, dstars = [], {}
    for n in ns:
        signs = signs_for(n)
        dstars[n] = d_star(n, rng, signs)
        for d in (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0):
            curve.append(power_at_d(n, d, rng, signs))
    write_csv(OUT / "standardized_power_curve.csv", curve)

    # Use delta divided by SD as the standardized effect for each power scenario.
    grid, mdi = [], []
    for scope, sds, pool_n in (
        ("full_scan", SD_FULL, n_elig),
        ("flank_subset", SD_FLANK, len(flank_any)),
    ):
        for n in ns:
            if scope == "flank_subset" and n > max(len(flank_any), 3):
                continue
            signs = signs_for(n)
            for label, sd in sds.items():
                for delta in DELTA_GRID:
                    row = power_at_d(n, delta / sd, rng, signs)
                    row.update(
                        scope=scope, sd_label=label, sd_between_proteins=sd,
                        delta=delta,
                    )
                    grid.append(row)
                mdi.append(
                    {
                        "scope": scope,
                        "n_confirmation_proteins": n,
                        "n_left_for_all_fitted_stages": pool_n - n,
                        "sd_label": label,
                        "sd_between_proteins": sd,
                        "min_attainable_exact_p": 2 / (1 << n),
                        "decisive_p_reachable": 2 / (1 << n) < 0.05,
                        "d_star_80pct": dstars[n],
                        "min_detectable_improvement_80pct": (
                            round(dstars[n] * sd, 4) if dstars[n] else None
                        ),
                    }
                )
    write_csv(OUT / "confirmation_power_grid.csv", grid)
    write_csv(OUT / "minimum_detectable_improvement.csv", mdi)

    summary = {
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "outcome_blind": "no confirmation outcome label read; A0 metadata + "
                         "development paired differences on opened proteins only",
        "inputs": {
            "a0_feasibility": {"path": str(A0.relative_to(ROOT)), "sha256": sha256(A0)},
            "all_comparisons": {"path": str(CMP.relative_to(ROOT)), "sha256": sha256(CMP)},
            "run_phase_c": {"path": str(PHASE_C.relative_to(ROOT)), "sha256": sha256(PHASE_C)},
        },
        "pool": {
            "assays_accepted": len(pool),
            "proteins_total": len(proteins),
            "proteins_confirmation_eligible": n_elig,
            "eligible": [p["protein"] for p in eligible],
            "excluded": {
                p["protein"]: p["exclusion_reasons"]
                for p in proteins if not p["confirmation_eligible"]
            },
            "phase_c_score_set_parsed": phase_c,
        },
        "flank_coverage_among_eligible": {
            "substantial_ge25_sites": [p["protein"] for p in flank_sub],
            "token_1_to_24_sites": [
                p["protein"] for p in flank_any if p["flank_tier"] == "token"
            ],
            "none": [p["protein"] for p in proteins
                     if p["confirmation_eligible"] and p["flank_tier"] == "none"],
        },
        "assay_averaging_rule": "arithmetic mean of per-assay Spearman within a "
                                "protein, so a protein with 2 assays carries the "
                                "same weight as one with 1",
        "signflip_arithmetic": {
            "smallest_exact_two_sided_p": "2 / 2**n",
            "n_required_for_p_below_0.05": 6,
            "n_eligible_available": n_elig,
        },
        "seed": SEED,
        "n_reps": N_REPS,
        "n_boot": N_BOOT,
    }
    (OUT / "confirmation_feasibility.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    print(f"accepted assays: {len(pool)}")
    print(f"proteins: {len(proteins)}  confirmation-eligible: {n_elig}")
    print(f"eligible: {[p['protein'] for p in eligible]}")
    print(f"flank substantial among eligible: {[p['protein'] for p in flank_sub]}")
    print(f"flank token among eligible: "
          f"{[p['protein'] for p in flank_any if p['flank_tier'] == 'token']}")
    print(f"wrote {len(list(OUT.glob('*')))} files to {OUT.relative_to(ROOT)}")


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()

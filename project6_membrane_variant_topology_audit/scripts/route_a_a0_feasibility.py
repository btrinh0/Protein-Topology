"""Route A, step A0: the feasibility count, which is the decision gate.

Three numbers, then a power simulation, then the decision floor. The DMS_score column is never
read; every count uses variant identity and position only.

Counts
  1. ProteinGym assays whose protein maps onto the topology census (accepted mapping).
  2. Of those, membrane proteins with at least one measured flank site inside the 10-residue window.
  3. K/R-gain substitutions (single mutants, alternative K or R, native residue in ACFGILMNPQSTVWY)
     at in-window flank sites.

Floor, fixed in the build plan before looking: at least 8 assays with at least 30 K/R-gain flank
substitutions each. Below the floor, Route A's capability test is too thin to answer.

The power simulation uses the real K/R-gain flank variants and their real uncorrected AlphaMissense
scores, with a hypothetical DMS outcome whose side effect is planted in SD units. It answers "what
paired difference in per-assay Spearman can this design detect?" Every outcome in it is simulated.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))

import route_a_common as common  # noqa: E402
import form_ab_lib as lib  # noqa: E402

FLOOR_ASSAYS = 8
FLOOR_SITES = 30
BETA = 0.02971251410524106  # Form B value fitted on the sequence-split training buckets
DECAY_LENGTH = 8.0
OUT = common.ROOT / "results/route_a_20261012/a0_feasibility.json"
ALPHA_INDEX = {aa: i for i, aa in enumerate(lib.ALPHABET)}
NEUTRAL_SET = set(common.NEUTRAL)


def single_variants(mutants: list[list[tuple]]) -> list[tuple]:
    return [tokens[0] for tokens in mutants if len(tokens) == 1]


def assay_counts(row: dict, census: dict, uniprot: dict, mutants: list[list[tuple]]) -> dict:
    """Counts for one assay. Reads no score."""
    entry = uniprot.get(row["UniProt_ID"])
    out = {"dms_id": row["DMS_id"], "uniprot_entry": row["UniProt_ID"], "taxon": row["taxon"],
           "accession": entry["accession"] if entry else None,
           "in_census": bool(entry) and entry["accession"] in census,
           "n_variants": len(mutants), "n_single": len(single_variants(mutants)),
           "n_multiple": sum(1 for t in mutants if len(t) > 1)}
    if not out["in_census"]:
        return out
    record = census[entry["accession"]]
    mapping = common.map_assay(row, record, mutants)
    flanks = common.census_flanks(record)
    mapping_array = mapping.pop("_mapping")
    out.update({k: mapping[k] for k in ("reference_agreement", "positions_dropped_as_tag_or_construct",
                                        "accepted", "identity_over_aligned", "census_length",
                                        "target_length")})
    if not mapping["accepted"]:
        out["reject_reason"] = f"reference agreement {mapping['reference_agreement']:.4f} below {common.MIN_AGREEMENT}"
        return out

    sequence = record["sequence"]
    flank_sites_measured = set()
    kr_gain = []   # (census_position, wt, alt, side, distance)
    neutral_alt_sites = defaultdict(set)   # census position -> neutral alternatives measured
    kr_sites = set()
    for wt, position, alt in single_variants(mutants):
        census_pos = int(mapping_array[position])
        if census_pos == 0 or sequence[census_pos - 1] != wt:
            continue
        flank = flanks.get(census_pos)
        if flank is None or flank[1] > common.WINDOW:
            continue
        side, distance, _ = flank
        flank_sites_measured.add(census_pos)
        if wt in NEUTRAL_SET and alt in common.KR:
            kr_gain.append((census_pos, wt, alt, side, distance))
            kr_sites.add(census_pos)
        if wt in NEUTRAL_SET and alt in NEUTRAL_SET and alt != wt:
            neutral_alt_sites[census_pos].add(alt)

    # Eligible for the contrast: a site with both a K/R alternative and a neutral alternative measured.
    eligible_sites = {s for s in kr_sites if neutral_alt_sites.get(s)}
    eligible_gain = [v for v in kr_gain if v[0] in eligible_sites]
    # Strata of the estimator: protein x native residue, with both sides present among eligible sites.
    strata = defaultdict(lambda: {0: 0, 1: 0})
    for census_pos, wt, _alt, side, _d in eligible_gain:
        strata[wt][side] += 1
    two_sided = {wt for wt, sides in strata.items() if sides[0] and sides[1]}
    out.update({
        "flank_sites_measured_in_window": len(flank_sites_measured),
        "kr_gain_flank_single": len(kr_gain),
        "kr_gain_flank_cytosolic": sum(1 for v in kr_gain if v[3] == 0),
        "kr_gain_flank_outer": sum(1 for v in kr_gain if v[3] == 1),
        "kr_gain_sites_eligible": len(eligible_sites),
        "kr_gain_variants_in_two_sided_strata": sum(1 for v in eligible_gain if v[1] in two_sided),
        "two_sided_strata": len(two_sided),
    })
    out["_kr_variants"] = kr_gain
    out["_eligible_variants"] = eligible_gain
    return out


def load_am_lookup() -> dict:
    """Sorted keys for the uncorrected AlphaMissense flank scores: (protein index, position, ref, alt)."""
    with np.load(common.AM_ARCHIVE, allow_pickle=False) as blob:
        accessions = blob["accessions"].astype(str)
        protein = blob["protein"].astype(np.int64)
        position = blob["position"].astype(np.int64)
        ref = blob["ref_aa"].astype(np.int64)
        alt = blob["alt_aa"].astype(np.int64)
        score = blob["am_pathogenicity"].astype(float)
    key = (protein << 40) + (position << 16) + (ref << 8) + alt
    order = np.argsort(key)
    return {"index": {a: i for i, a in enumerate(accessions)}, "key": key[order], "score": score[order]}


def am_scores(lookup: dict, accession: str, positions, refs, alts) -> np.ndarray:
    """Uncorrected AlphaMissense for each (position, ref, alt); NaN where the archive has no row."""
    positions = np.asarray(positions, dtype=np.int64)
    if accession not in lookup["index"]:
        return np.full(len(positions), np.nan)
    protein = lookup["index"][accession]
    keys = (np.int64(protein) << 40) + (positions << 16) + (np.asarray(refs, dtype=np.int64) << 8) \
        + np.asarray(alts, dtype=np.int64)
    found = np.minimum(np.searchsorted(lookup["key"], keys), len(lookup["key"]) - 1)
    ok = lookup["key"][found] == keys
    return np.where(ok, lookup["score"][found], np.nan)


def corrected(scores: np.ndarray, side: np.ndarray, distance: np.ndarray, ref: np.ndarray, alt: np.ndarray) -> np.ndarray:
    """Apply the selected Form B correction to one row per variant."""
    matrix = np.full((len(scores), 20), np.nan)
    matrix[np.arange(len(scores)), alt] = scores
    panel = {"protein": np.zeros(len(scores), dtype=np.int64), "ref_aa": ref, "side": side,
             "distance": distance, "scores": matrix}
    shifted = lib.apply_correction(panel, BETA, DECAY_LENGTH)
    return shifted[np.arange(len(scores)), alt]


def simulate(assays: list[dict], lookup: dict, gammas, reps: int, boot: int, seed: int) -> dict:
    """Power for the paired mean Spearman difference (corrected minus uncorrected), protein-clustered.

    Assumption-labelled outcome model, per assay a and variant v, all simulated:
        y_v = -( r * z_v + sqrt(1 - r^2) * e_v + gamma * outer_v * decay_v ),
    z_v the within-assay standardised uncorrected AlphaMissense score (higher = more damaging),
    e_v ~ N(0, 1), r = 0.4 the baseline correlation, gamma the planted side effect in SD units.
    A study is one draw of e for every assay; the bootstrap resamples proteins (clusters).
    """
    rng = np.random.default_rng(seed)
    prepared = []
    for a in assays:
        s = a["am"]
        keep = np.isfinite(s)
        if keep.sum() < 10 or np.std(s[keep]) == 0:
            continue
        z = (s[keep] - s[keep].mean()) / s[keep].std()
        outer = a["side"][keep] == 1
        decay = np.exp(-(a["distance"][keep] - 1) / DECAY_LENGTH)
        corr = a["corrected"][keep]
        prepared.append({"accession": a["accession"], "z": z, "s": s[keep], "c": corr,
                         "outer": outer.astype(float), "decay": decay})
    clusters = sorted({p["accession"] for p in prepared})
    cluster_of = np.array([clusters.index(p["accession"]) for p in prepared])
    result = {"n_assays": len(prepared), "n_clusters": len(clusters),
              "baseline_r": 0.4, "reps": reps, "bootstrap_draws": boot, "gammas": {}}
    for gamma in gammas:
        d_mean, positives = [], 0
        for _ in range(reps):
            diffs = []
            for p in prepared:
                e = rng.standard_normal(len(p["z"]))
                y = -(0.4 * p["z"] + np.sqrt(1 - 0.16) * e + gamma * p["outer"] * p["decay"])
                rho_u = spearmanr(p["s"], -y).statistic
                rho_c = spearmanr(p["c"], -y).statistic
                diffs.append(rho_c - rho_u)
            diffs = np.asarray(diffs)
            d_mean.append(diffs.mean())
            # Cluster bootstrap over proteins: resample clusters with replacement, average within each.
            idx = rng.integers(0, len(clusters), size=(boot, len(clusters)))
            cluster_sums = np.bincount(cluster_of, weights=diffs, minlength=len(clusters))
            cluster_counts = np.bincount(cluster_of, minlength=len(clusters))
            boot_means = (cluster_sums[idx].sum(axis=1)) / (cluster_counts[idx].sum(axis=1))
            low = np.percentile(boot_means, 2.5)
            if low > 0:
                positives += 1
        result["gammas"][str(gamma)] = {"mean_paired_spearman_difference": float(np.mean(d_mean)),
                                        "sd_of_study_mean": float(np.std(d_mean, ddof=1)),
                                        "power_one_sided_lower_ci_above_zero": positives / reps}
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reps", type=int, default=300)
    parser.add_argument("--boot", type=int, default=999)
    parser.add_argument("--seed", type=int, default=20261012)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()

    reference = common.load_reference()
    uniprot = common.load_uniprot_map()
    census = common.load_census()
    rows = []
    for row in reference:
        path = common.ASSAY_DIR / row["DMS_filename"]
        mutants = []
        with path.open(encoding="utf-8") as handle:
            header = handle.readline().rstrip("\n").split(",")
            if header[0] != "mutant":
                raise SystemExit(f"{path.name}: expected mutant in the first column, found {header[0]!r}")
            for line in handle:
                # Only the first field is split off; mutated_sequence and DMS_score are never parsed.
                token = line.split(",", 1)[0].strip()
                if token:
                    mutants.append(common.parse_mutant(token))
        rows.append(assay_counts(row, census, uniprot, mutants))

    counts = {
        "assays_total": len(rows),
        "assays_human": sum(r["taxon"] == "Human" for r in rows),
        "assays_mapped_to_census": sum(bool(r.get("accepted")) for r in rows),
        "assays_membrane_with_window_flank": sum(bool(r.get("accepted")) and r.get("flank_sites_measured_in_window", 0) > 0 for r in rows),
        "kr_gain_flank_single_total": sum(r.get("kr_gain_flank_single", 0) for r in rows),
    }
    passing = [r for r in rows if r.get("kr_gain_flank_single", 0) >= FLOOR_SITES]
    gate = {"floor_assays": FLOOR_ASSAYS, "floor_kr_gain_per_assay": FLOOR_SITES,
            "assays_meeting_per_assay_floor": len(passing),
            "passes": len(passing) >= FLOOR_ASSAYS,
            "assays_meeting_per_assay_floor_ids": sorted(r["dms_id"] for r in passing)}

    lookup = load_am_lookup()
    candidates = []
    for r in rows:
        if not r.get("accepted") or not r.get("_kr_variants"):
            continue
        variants = r["_kr_variants"]
        positions = np.array([v[0] for v in variants]); wt = np.array([ALPHA_INDEX[v[1]] for v in variants])
        alt = np.array([ALPHA_INDEX[v[2]] for v in variants]); side = np.array([v[3] for v in variants])
        dist = np.array([v[4] for v in variants])
        am = am_scores(lookup, r["accession"], positions, wt, alt)
        cor = corrected(am, side, dist, wt, alt)
        candidates.append({"dms_id": r["dms_id"], "accession": r["accession"], "am": am, "corrected": cor,
                           "side": side, "distance": dist, "n_kr": len(variants), "n_am_missing": int(np.isnan(am).sum())})
    # Predictor-only perturbation: how much the correction reorders the capability variants (no outcome used).
    perturbation = []
    for c in candidates:
        keep = np.isfinite(c["am"])
        if keep.sum() >= 3:
            perturbation.append({"dms_id": c["dms_id"], "n": int(keep.sum()),
                                 "spearman_uncorrected_vs_corrected": float(spearmanr(c["am"][keep], c["corrected"][keep]).statistic),
                                 "fraction_rows_changed": float(np.mean(c["am"][keep] != c["corrected"][keep]))})
    simulation = simulate(candidates, lookup, gammas=(0.0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0),
                          reps=args.reps, boot=args.boot, seed=args.seed)

    out = {
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "outcome_blind": "DMS_score not read; mutant column only",
        "reference_table_sha256": lib.P.sha256(common.REFERENCE),
        "uniprot_map_sha256": lib.P.sha256(common.UNIPROT_MAP),
        "correction": {"beta": BETA, "decay_length": DECAY_LENGTH, "source": "results/form_b_20261008/alphamissense_sequence.json"},
        "counts": counts,
        "gate": gate,
        "per_assay": [{k: v for k, v in r.items() if not k.startswith("_")} for r in rows],
        "predictor_only_perturbation": perturbation,
        "simulation": simulation,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o)) + "\n", encoding="utf-8")
    print(json.dumps({"counts": counts, "gate": gate, "simulation": simulation}, indent=2, default=str))


if __name__ == "__main__":
    main()

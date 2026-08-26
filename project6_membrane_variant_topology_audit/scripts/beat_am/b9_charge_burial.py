"""Audit charge-introduction scores by membrane depth and side."""

from __future__ import annotations

import gzip
import json
import sys
import time
import warnings
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from run_gate1 import read_topology  # noqa: E402

AM_GZ = ROOT / "data" / "raw" / "AlphaMissense_aa_substitutions.tsv.gz"
ESM_ZIP = ROOT / "data" / "raw" / "ESM1b_ALL_hum_isoforms_LLR.zip"
ESM_PREFIX = "content/ALL_hum_isoforms_ESM1b_LLR/"
OUT = ROOT / "results" / "beat_am_20261016" / "charge_burial"

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
AIDX = {a: i for i, a in enumerate(ALPHABET)}
UNCHARGED = "ACFGILMNPQSTVWY"          # Exclude charged residues and histidine from the uncharged comparator.
CHARGED = "DEKR"
POSITIVE = "KR"
NEGATIVE = "DE"

MIN_TM_LEN = 10
N_BINS = 5
LOOP_GROUP = N_BINS
N_GROUPS = N_BINS + 1
MAX_LOOP_PER_STRATUM = 10
MIN_REF_ARM = 2                         # Require at least two reference residues as in the frozen B8 stratum definition.
MIN_UNCHARGED_ALTS = 10
MIN_ALTS_FOR_SPREAD = 15
N_BOOT = 1999
SEED = 20261016

_LOG_LINES: list[str] = []


def log(m: str) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {m}"
    _LOG_LINES.append(line)
    print(line, flush=True)



def depth_bin(p: int, s: int, seglen: int) -> int:
    """Fractional depth d in [0,1] (0 = interface, 1 = centre) -> bin 0..4."""
    half = (seglen - 1) / 2.0
    d = 1.0 - abs((p - s) - half) / half
    return min(int(d * N_BINS), N_BINS - 1)


def enumerate_sites(census: dict) -> dict:
    """Eligible in-TM positions by depth bin, plus the capped loop anchor pool.

    Returns parallel arrays. `group` is 0..4 for depth bins and 5 for loops.
    `side` is 0 unknown, 1 cytosolic, 2 non-cytosolic (only set for loops).
    """
    rng = np.random.default_rng(SEED)
    accs: list[str] = []
    prot_names: list[str] = []
    pos_l: list[int] = []
    nat_l: list[int] = []
    grp_l: list[int] = []
    prot_l: list[int] = []
    side_l: list[int] = []
    unch_idx = {a: i for i, a in enumerate(UNCHARGED)}

    n_tm_elig = 0
    for acc in sorted(census):
        rec = census[acc]
        seq = rec.get("sequence") or ""
        if not seq:
            continue
        L = len(seq)
        pidx = len(prot_names)
        rows_before = len(pos_l)

        # Include every TM position when excluding membrane segments from the loop-control pool.
        tm_cover: set[int] = set()
        for t in rec.get("tmds") or []:
            s, e = int(t["start"]), int(t["end"])
            tm_cover.update(range(s, e + 1))
        for t in rec.get("tmds") or []:
            s, e = int(t["start"]), int(t["end"])
            seglen = e - s + 1
            if seglen < MIN_TM_LEN:
                continue
            n_tm_elig += 1
            for p in range(s, e + 1):
                if p > L:
                    continue
                nat = seq[p - 1]
                if nat not in unch_idx:
                    continue
                pos_l.append(p)
                nat_l.append(unch_idx[nat])
                grp_l.append(depth_bin(p, s, seglen))
                prot_l.append(pidx)
                side_l.append(0)


        pool: dict[tuple[int, int], list[int]] = defaultdict(list)
        for _lid, lp in (rec.get("loops") or {}).items():
            loc = (lp.get("location") or "").strip().lower()
            if loc not in ("inside", "outside"):
                continue
            sd = 1 if loc == "inside" else 2
            for p in range(int(lp["start"]), int(lp["end"]) + 1):
                if p > L or p in tm_cover:
                    continue
                nat = seq[p - 1]
                if nat not in unch_idx:
                    continue
                pool[(unch_idx[nat], sd)].append(p)
        for (ni, sd), plist in pool.items():
            if len(plist) > MAX_LOOP_PER_STRATUM:
                take = rng.choice(len(plist), size=MAX_LOOP_PER_STRATUM,
                                  replace=False)
                plist = [plist[j] for j in take]
            for p in plist:
                pos_l.append(p)
                nat_l.append(ni)
                grp_l.append(LOOP_GROUP)
                prot_l.append(pidx)
                side_l.append(sd)

        if len(pos_l) > rows_before:
            prot_names.append(acc)
            accs.extend([acc] * (len(pos_l) - rows_before))

    sites = {
        "acc": np.array(accs),
        "pos": np.asarray(pos_l, dtype=np.int32),
        "nat": np.asarray(nat_l, dtype=np.int8),
        "grp": np.asarray(grp_l, dtype=np.int8),
        "prot": np.asarray(prot_l, dtype=np.int32),
        "side": np.asarray(side_l, dtype=np.int8),
        "prot_names": prot_names,
        "n_tm_eligible_segments": n_tm_elig,
    }
    g = sites["grp"]
    log(f"  eligible TM segments (len>={MIN_TM_LEN}): {n_tm_elig:,}")
    for b in range(N_BINS):
        log(f"    depth bin {b}: {int((g == b).sum()):,} uncharged-native positions")
    log(f"    loop anchor (capped at {MAX_LOOP_PER_STRATUM} per protein x residue): "
        f"{int((g == LOOP_GROUP).sum()):,}")
    log(f"  positions to score: {len(g):,} across {len(prot_names):,} proteins")
    return sites



def build_index(sites: dict) -> dict:
    # Use plain string lookup keys for the streaming AlphaMissense pass.
    return {(a, p): i for i, (a, p) in
            enumerate(zip(sites["acc"].tolist(), sites["pos"].tolist()))}


def extract_am(index: dict, n_rows: int) -> np.ndarray:
    """Single streaming pass over the 216M-line AlphaMissense table."""
    S = np.full((n_rows, 20), np.nan, dtype=np.float32)
    seen = filled = 0
    t0 = time.time()
    with gzip.open(AM_GZ, "rt") as fh:
        for line in fh:
            if line[0] == "#" or line.startswith("uniprot_id"):
                continue
            seen += 1
            if seen % 40_000_000 == 0:
                log(f"    {seen:,} lines, {time.time() - t0:.0f}s, {filled:,} cells")
            t1 = line.find("\t")
            t2 = line.find("\t", t1 + 1)
            var = line[t1 + 1:t2]
            pos_s = var[1:-1]
            if not pos_s.isdigit():
                continue
            row = index.get((line[:t1], int(pos_s)))
            if row is None:
                continue
            col = AIDX.get(var[-1])
            if col is None:
                continue
            t3 = line.find("\t", t2 + 1)
            S[row, col] = float(line[t2 + 1:t3])
            filled += 1
    log(f"  AM pass: {seen:,} lines in {time.time() - t0:.0f}s, {filled:,} scores")
    return S


def extract_esm(sites: dict, index: dict, n_rows: int) -> np.ndarray:
    """Per-accession CSVs out of the zip. Writes damage = -LLR (negated)."""
    S = np.full((n_rows, 20), np.nan, dtype=np.float32)
    wanted_by_acc: dict[str, set[int]] = defaultdict(set)
    for (a, p) in index:
        wanted_by_acc[a].add(p)
    found = missing = filled = 0
    t0 = time.time()
    with zipfile.ZipFile(ESM_ZIP) as z:
        available = set(z.namelist())
        for n_done, acc in enumerate(sorted(wanted_by_acc), start=1):
            if n_done % 1000 == 0:
                log(f"    {n_done:,} accessions, {time.time() - t0:.0f}s, "
                    f"{filled:,} cells")
            name = f"{ESM_PREFIX}{acc}_LLR.csv"
            if name not in available:
                missing += 1
                continue
            with z.open(name) as fh:
                rows = fh.read().decode("utf-8").splitlines()
            if len(rows) < 2:
                missing += 1
                continue
            found += 1
            want = wanted_by_acc[acc]
            header = rows[0].split(",")
            # The first header cell is empty and remaining labels encode native residue and position.
            cols: list[tuple[int, int]] = []
            for ci, cell in enumerate(header[1:], start=1):
                parts = cell.strip().split()
                if len(parts) != 2:
                    continue
                p = int(parts[1])
                if p not in want:
                    continue
                r = index.get((acc, p))
                if r is not None:
                    cols.append((ci, r))
            if not cols:
                continue
            for line in rows[1:]:
                cells = line.split(",")
                alt = cells[0].strip()
                col = AIDX.get(alt)
                if col is None:
                    continue
                for ci, r in cols:
                    v = cells[ci]
                    if not v or v == "nan":
                        continue
                    try:
                        # Negate fitness-oriented LLR values to obtain damage-oriented scores.
                        S[r, col] = -float(v)
                    except ValueError:
                        continue
                    filled += 1
    log(f"  ESM1b: {found:,} accessions read, {missing:,} absent from the zip, "
        f"{filled:,} scores in {time.time() - t0:.0f}s")
    return S


def load_or_extract(tag: str, sites: dict, index: dict, fn) -> np.ndarray:
    cache = OUT / f"position_scores_{tag}.npz"
    n = len(sites["pos"])
    if cache.exists():
        z = np.load(cache, allow_pickle=False)
        if (z["pos"].shape == sites["pos"].shape
                and bool((z["pos"] == sites["pos"]).all())
                and bool((z["acc"] == sites["acc"]).all())):
            log(f"  reusing cached {tag} scores from {cache.name}")
            return z["scores"]
        log(f"  cached {tag} scores do not match the current site table; re-extracting")
    S = fn(index, n) if tag == "am" else fn(sites, index, n)
    np.savez_compressed(cache, scores=S, acc=sites["acc"], pos=sites["pos"])
    log(f"  cached {tag} scores -> {cache.name}")
    return S



def contrasts(S: np.ndarray, nat: np.ndarray) -> dict:
    """z, z_KR, z_DE and the within-position spread, vectorised over sites.

    Reference arm = the uncharged residues other than the native one. Histidine
    is in neither arm. A site is usable only when all four charged alternatives
    are present and at least MIN_UNCHARGED_ALTS of the reference arm are.
    """
    ch = [AIDX[a] for a in CHARGED]
    po = [AIDX[a] for a in POSITIVE]
    ne = [AIDX[a] for a in NEGATIVE]
    unch = np.array([AIDX[a] for a in UNCHARGED])


    ref = S[:, unch]
    is_native = (np.arange(len(UNCHARGED))[None, :] == nat[:, None])
    ref = np.where(is_native, np.nan, ref)
    n_ref = np.sum(~np.isnan(ref), axis=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        ref_mean = np.nanmean(ref, axis=1)

    ok_ch = ~np.isnan(S[:, ch]).any(axis=1)
    usable = ok_ch & (n_ref >= MIN_UNCHARGED_ALTS)

    z = S[:, ch].mean(axis=1) - ref_mean
    z_kr = S[:, po].mean(axis=1) - ref_mean
    z_de = S[:, ne].mean(axis=1) - ref_mean

    # The S4 scale uses each position's spread across all 19 non-native alternatives.
    allalt = S.copy()
    native_col = np.array([AIDX[a] for a in UNCHARGED])[nat]
    allalt[np.arange(len(nat)), native_col] = np.nan
    n_alt = np.sum(~np.isnan(allalt), axis=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        spread = np.nanstd(allalt, axis=1, ddof=1)
    ok_spread = usable & (n_alt >= MIN_ALTS_FOR_SPREAD) & (spread > 0)
    z_scaled = np.where(ok_spread, z / np.where(spread > 0, spread, np.nan), np.nan)

    return {"z": z, "z_KR": z_kr, "z_DE": z_de, "z_scaled": z_scaled,
            "usable": usable, "usable_scaled": ok_spread, "spread": spread}



def stratum_table(prot, nat, grp, z, n_prot, n_groups):
    """Dense per-stratum (protein x native residue) counts and means by group."""
    key = (prot.astype(np.int64) * 16 + nat) * n_groups + grp
    uniq, inv = np.unique(key, return_inverse=True)
    cnt = np.bincount(inv).astype(np.float64)
    tot = np.bincount(inv, weights=z.astype(np.float64))
    strat = uniq // n_groups
    grp_of = (uniq % n_groups).astype(np.int64)
    sids, sinv = np.unique(strat, return_inverse=True)
    N = np.zeros((len(sids), n_groups))
    M = np.zeros((len(sids), n_groups))
    N[sinv, grp_of] = cnt
    M[sinv, grp_of] = np.divide(tot, cnt, out=np.zeros_like(tot), where=cnt > 0)
    prot_of = (sids // 16).astype(np.int64)
    return N, M, prot_of, n_prot


def raw_arm(N, M, prot_of, n_prot, g):
    """Per-protein numerator/denominator for the plain site mean within group g."""
    num = np.bincount(prot_of, weights=N[:, g] * M[:, g], minlength=n_prot)
    den = np.bincount(prot_of, weights=N[:, g], minlength=n_prot)
    return num, den


def contrast_arm(N, M, prot_of, n_prot, g1, g0):
    """Per-protein numerator/denominator for the Algorithm 1 weighted difference.

    Strata need >= 1 site in the test arm and >= MIN_REF_ARM in the reference
    arm; weight w = n1*n0/(n1+n0), exactly as B8.
    """
    m = (N[:, g1] >= 1) & (N[:, g0] >= MIN_REF_ARM)
    if not m.any():
        return (np.zeros(n_prot), np.zeros(n_prot), 0, 0, 0)
    n1, n0 = N[m, g1], N[m, g0]
    w = n1 * n0 / (n1 + n0)
    d = M[m, g1] - M[m, g0]
    num = np.bincount(prot_of[m], weights=w * d, minlength=n_prot)
    den = np.bincount(prot_of[m], weights=w, minlength=n_prot)
    return num, den, int(m.sum()), int(n1.sum()), int(n0.sum())


def spearman_vs_index(V: np.ndarray) -> np.ndarray:
    """Spearman rho between the column index 0..k-1 and each row of V."""
    n = V.shape[1]
    order = np.argsort(V, axis=1, kind="stable")
    ranks = np.empty(V.shape, dtype=np.float64)
    np.put_along_axis(ranks, order, np.arange(n, dtype=np.float64)[None, :], axis=1)
    d = ranks - np.arange(n, dtype=np.float64)[None, :]
    return 1.0 - 6.0 * (d * d).sum(axis=1) / (n * (n * n - 1))


def ci(draws: np.ndarray) -> list[float]:
    lo, hi = np.nanquantile(draws, [0.025, 0.975])
    return [float(lo), float(hi)]


def analyse(prot, nat, grp, z, n_prot, label):
    """The whole depth analysis for one z-variant of one scorer.

    One shared whole-protein resample per draw, so the five bins, the
    monotonicity statistic and the endpoint gap are consistent within a draw.
    """
    N, M, prot_of, n_prot = stratum_table(prot, nat, grp, z, n_prot, N_GROUPS)

    quantities: list[tuple[str, np.ndarray, np.ndarray]] = []
    meta: dict = {}
    for b in range(N_BINS):
        num, den = raw_arm(N, M, prot_of, n_prot, b)
        quantities.append((f"raw_bin{b}", num, den))
    num, den = raw_arm(N, M, prot_of, n_prot, LOOP_GROUP)
    quantities.append(("raw_loop", num, den))
    for b in range(1, N_BINS):
        num, den, ns, n1, n0 = contrast_arm(N, M, prot_of, n_prot, b, 0)
        quantities.append((f"adj_bin{b}_vs_bin0", num, den))
        meta[f"adj_bin{b}_vs_bin0"] = {"n_strata": ns, "n_test_sites": n1,
                                       "n_reference_sites": n0}
    for name, g1, g0 in (("centre_vs_loop", N_BINS - 1, LOOP_GROUP),
                         ("interface_vs_loop", 0, LOOP_GROUP)):
        num, den, ns, n1, n0 = contrast_arm(N, M, prot_of, n_prot, g1, g0)
        quantities.append((name, num, den))
        meta[name] = {"n_strata": ns, "n_test_sites": n1, "n_reference_sites": n0}

    NUM = np.vstack([q[1] for q in quantities])
    DEN = np.vstack([q[2] for q in quantities])
    names = [q[0] for q in quantities]
    ones = np.ones(n_prot)
    with np.errstate(invalid="ignore", divide="ignore"):
        point = (NUM @ ones) / (DEN @ ones)

    rng = np.random.default_rng(SEED)
    draws = np.empty((N_BOOT, len(names)))
    chunk = 250
    for start in range(0, N_BOOT, chunk):
        k = min(chunk, N_BOOT - start)
        pick = rng.integers(0, n_prot, size=(k, n_prot))
        C = np.zeros((k, n_prot))
        for i in range(k):
            C[i] = np.bincount(pick[i], minlength=n_prot)
        with np.errstate(invalid="ignore", divide="ignore"):
            draws[start:start + k] = (C @ NUM.T) / (C @ DEN.T)

    col = {n: i for i, n in enumerate(names)}
    raw_cols = [col[f"raw_bin{b}"] for b in range(N_BINS)]
    adj_curve_draws = np.column_stack(
        [np.zeros(N_BOOT)] + [draws[:, col[f"adj_bin{b}_vs_bin0"]]
                              for b in range(1, N_BINS)])
    adj_curve_point = np.concatenate(
        [[0.0], [point[col[f"adj_bin{b}_vs_bin0"]] for b in range(1, N_BINS)]])

    rho_raw = spearman_vs_index(draws[:, raw_cols])
    rho_adj = spearman_vs_index(adj_curve_draws)
    mono_raw = float(np.mean(np.all(np.diff(draws[:, raw_cols], axis=1) > 0, axis=1)))
    mono_adj = float(np.mean(np.all(np.diff(adj_curve_draws, axis=1) > 0, axis=1)))

    sd_z = float(np.nanstd(z, ddof=1))
    n_by_grp = {int(g): int((grp == g).sum()) for g in range(N_GROUPS)}

    out = {
        "variant": label,
        "n_sites_used": int(len(z)),
        "n_proteins": int(n_prot),
        "n_sites_by_group": {("loop" if g == LOOP_GROUP else f"bin{g}"): n
                             for g, n in n_by_grp.items()},
        "sd_of_z_all_sites": sd_z,
        "raw_curve": {
            f"bin{b}": {"mean_z": float(point[col[f'raw_bin{b}']]),
                        "ci95": ci(draws[:, col[f"raw_bin{b}"]]),
                        "n_sites": n_by_grp[b]}
            for b in range(N_BINS)},
        "raw_loop": {"mean_z": float(point[col["raw_loop"]]),
                     "ci95": ci(draws[:, col["raw_loop"]]),
                     "n_sites": n_by_grp[LOOP_GROUP]},
        "adjusted_curve_vs_interface": {
            "bin0": {"delta": 0.0, "ci95": [0.0, 0.0],
                     "note": "reference arm, zero by construction"},
            **{f"bin{b}": {"delta": float(point[col[f'adj_bin{b}_vs_bin0']]),
                           "ci95": ci(draws[:, col[f"adj_bin{b}_vs_bin0"]]),
                           "delta_standardized_within_this_scorer":
                               float(point[col[f"adj_bin{b}_vs_bin0"]]) / sd_z,
                           **meta[f"adj_bin{b}_vs_bin0"]}
               for b in range(1, N_BINS)}},
        "endpoint_gap": {
            "raw_bin4_minus_bin0": float(point[col["raw_bin4"]]
                                         - point[col["raw_bin0"]]),
            "raw_intervals_disjoint": bool(
                ci(draws[:, col["raw_bin4"]])[0] > ci(draws[:, col["raw_bin0"]])[1]),
            "paired_stratified_delta": float(point[col["adj_bin4_vs_bin0"]]),
            "paired_stratified_ci95": ci(draws[:, col["adj_bin4_vs_bin0"]]),
            **meta["adj_bin4_vs_bin0"]},
        "monotonicity": {
            "raw_curve": {
                "spearman_rho_point": float(spearman_vs_index(
                    point[raw_cols][None, :])[0]),
                "spearman_rho_bootstrap_median": float(np.median(rho_raw)),
                "spearman_rho_ci95": ci(rho_raw),
                "p_strictly_increasing_across_all_5_bins": mono_raw},
            "adjusted_curve": {
                "spearman_rho_point": float(spearman_vs_index(
                    adj_curve_point[None, :])[0]),
                "spearman_rho_bootstrap_median": float(np.median(rho_adj)),
                "spearman_rho_ci95": ci(rho_adj),
                "p_strictly_increasing_across_all_5_bins": mono_adj}},
        "loop_anchor": {
            "centre_vs_loop": {"delta": float(point[col["centre_vs_loop"]]),
                               "ci95": ci(draws[:, col["centre_vs_loop"]]),
                               **meta["centre_vs_loop"]},
            "interface_vs_loop": {"delta": float(point[col["interface_vs_loop"]]),
                                  "ci95": ci(draws[:, col["interface_vs_loop"]]),
                                  **meta["interface_vs_loop"]}},
    }
    return out


def loop_side_test(prot, nat, side, z, n_prot):
    """S3: K/R introduction, cytosolic versus non-cytosolic loops."""
    grp = np.where(side == 1, 0, 1).astype(np.int8)     # Census side 0 is cytosolic and side 1 is non-cytosolic.
    N, M, prot_of, n_prot = stratum_table(prot, nat, grp, z, n_prot, 2)
    num, den, ns, n1, n0 = contrast_arm(N, M, prot_of, n_prot, 0, 1)
    if ns == 0:
        return None
    ones = np.ones(n_prot)
    point = float((num @ ones) / (den @ ones))
    rng = np.random.default_rng(SEED)
    draws = np.empty(N_BOOT)
    for b in range(N_BOOT):
        c = np.bincount(rng.integers(0, n_prot, size=n_prot), minlength=n_prot)
        draws[b] = (c @ num) / (c @ den)
    return {"delta_cytosolic_minus_noncytosolic": point, "ci95": ci(draws),
            "n_strata": ns, "n_cytosolic_sites": n1, "n_noncytosolic_sites": n0}



def directional(estimate: float, ci95: list, direction: int) -> str:
    """B8's verdict convention: a wrong-signed significant result says so."""
    lo, hi = ci95
    if direction > 0:
        if lo > 0:
            return "MET"
        return "NOT MET - significant in the WRONG direction" if hi < 0 \
            else "NOT MET - consistent with zero"
    if hi < 0:
        return "MET"
    return "NOT MET - significant in the WRONG direction" if lo > 0 \
        else "NOT MET - consistent with zero"


def verdicts(res: dict) -> dict:
    a = res["z"]
    adj = a["adjusted_curve_vs_interface"]
    gap = a["endpoint_gap"]
    mono = a["monotonicity"]["adjusted_curve"]

    curve = [0.0] + [adj[f"bin{b}"]["delta"] for b in range(1, N_BINS)]
    strictly_increasing = all(curve[i + 1] > curve[i] for i in range(N_BINS - 1))
    raw = [a["raw_curve"][f"bin{b}"]["mean_z"] for b in range(N_BINS)]
    raw_increasing = all(raw[i + 1] > raw[i] for i in range(N_BINS - 1))

    p1 = ("MET" if strictly_increasing else "NOT MET")
    p2_raw = gap["raw_intervals_disjoint"] and raw[4] > raw[0]
    p2_paired = gap["paired_stratified_ci95"][0] > 0
    p2 = "MET" if (p2_raw and p2_paired) else "NOT MET"

    s1 = a["loop_anchor"]["centre_vs_loop"]
    s1_ok = s1["delta"] > 0 and s1["ci95"][0] > 0
    s4 = res.get("z_scaled")
    if s4 is not None:
        s4_curve = [0.0] + [s4["adjusted_curve_vs_interface"][f"bin{b}"]["delta"]
                            for b in range(1, N_BINS)]
        s4_ok = (all(s4_curve[i + 1] > s4_curve[i] for i in range(N_BINS - 1))
                 and s4["endpoint_gap"]["paired_stratified_ci95"][0] > 0)
    else:
        s4_ok = None

    kr, de = res.get("z_KR"), res.get("z_DE")

    def grad_ok(r):
        if r is None:
            return None
        c = [0.0] + [r["adjusted_curve_vs_interface"][f"bin{b}"]["delta"]
                     for b in range(1, N_BINS)]
        return (all(c[i + 1] > c[i] for i in range(N_BINS - 1))
                and r["endpoint_gap"]["paired_stratified_ci95"][0] > 0)

    def grad_verdict(r):
        """MET, qualified when the bootstrap rarely reproduces the strict rise.

        A curve can be monotone at the point estimate and still be flat across
        its upper bins. An unqualified MET would hide that, so the probability
        is carried into the verdict string itself rather than left in a detail
        block a reader may skim past.
        """
        if r is None:
            return "not computed"
        if not grad_ok(r):
            return "NOT MET"
        p = r["monotonicity"]["adjusted_curve"][
            "p_strictly_increasing_across_all_5_bins"]
        if p < 0.5:
            return (f"MET BUT WEAK - monotone at the point estimate and the "
                    f"endpoint gap is solid, but only {p:.0%} of bootstrap draws "
                    f"reproduce a strict rise across all five bins; the upper "
                    f"bins are close to flat")
        return "MET"

    def grad_detail(r):
        """Point curve, endpoint interval and how firm the monotonicity is.

        A split can be monotone at the point estimate while the bootstrap
        rarely reproduces a strict rise across all five bins. Reporting the
        probability keeps a marginal pass from hiding behind a binary MET.
        """
        if r is None:
            return None
        return {
            "adjusted_curve": [0.0] + [
                r["adjusted_curve_vs_interface"][f"bin{b}"]["delta"]
                for b in range(1, N_BINS)],
            "endpoint_paired_delta": r["endpoint_gap"]["paired_stratified_delta"],
            "endpoint_paired_ci95": r["endpoint_gap"]["paired_stratified_ci95"],
            "bootstrap_p_strictly_increasing": r["monotonicity"]["adjusted_curve"][
                "p_strictly_increasing_across_all_5_bins"]}

    # Positive-inside predicts a negative cytosolic-minus-outer K/R introduction contrast.
    s3 = res.get("S3_loop_positive_inside_zKR")
    s3_verdict = directional(s3["delta_cytosolic_minus_noncytosolic"],
                             s3["ci95"], -1) if s3 else "not computed"

    return {
        "P1_monotonic_rise_adjusted_curve": p1,
        "P1_monotonic_rise_raw_curve": "MET" if raw_increasing else "NOT MET",
        "P1_detail": {
            "adjusted_curve": curve,
            "raw_curve": raw,
            "bootstrap_p_strictly_increasing": mono[
                "p_strictly_increasing_across_all_5_bins"],
            "spearman_rho_ci95": mono["spearman_rho_ci95"]},
        "P2_centre_exceeds_interface": p2,
        "P2_detail": {
            "raw_intervals_disjoint": bool(p2_raw),
            "paired_stratified_delta": gap["paired_stratified_delta"],
            "paired_stratified_ci95": gap["paired_stratified_ci95"],
            "paired_ci_excludes_zero_upward": bool(p2_paired)},
        "S1_centre_far_above_loop": "MET" if s1_ok else "NOT MET",
        "S1_detail": s1,
        "S2_KR_split_shows_gradient": grad_verdict(kr),
        "S2_DE_split_shows_gradient": grad_verdict(de),
        "S2_detail": {"z_KR": grad_detail(kr), "z_DE": grad_detail(de)},
        "S3_loop_KR_cheaper_on_cytosolic_side": s3_verdict,
        "S3_detail": {
            "declared_direction": "negative: z_KR(cytosolic) < z_KR(non-cytosolic)",
            **(s3 or {})},
        "S4_gradient_survives_scale_normalisation": (
            ("MET" if s4_ok else "NOT MET") if s4_ok is not None else "not computed"),
        "S4_detail": grad_detail(s4),
    }


DECLARED = {
    "rule": "Charge burial depth. Introducing D, E, K or R into a transmembrane "
            "helix costs free energy that RISES toward the bilayer centre, where "
            "there is no water to solvate the charge. A scorer that has learned "
            "membrane physics should penalise charge-introducing substitutions "
            "more at deep positions than at shallow ones.",
    "role": "POSITIVE CONTROL. Both scorers are expected to pass. The point is to "
            "show the instrument can detect compliance, not only violation.",
    "P1": "PRIMARY. Mean z rises monotonically across the five depth bins, from "
          "the lipid-water interface (bin 0) to the helix centre (bin 4).",
    "P2": "PRIMARY. z at the centre bin exceeds z at the interface bin, with "
          "non-overlapping 95% intervals; and the paired stratified "
          "Delta(bin 4 - bin 0) > 0 with a 95% interval excluding zero.",
    "P3": "POSITIVE CONTROL. BOTH AlphaMissense and ESM1b satisfy P1 and P2. "
          "Divergence between the two scorers is a reportable finding.",
    "S1": "SECONDARY. z at the helix centre greatly exceeds z at loop positions.",
    "S2": "SECONDARY. Both the K/R split and the D/E split show the gradient. No "
          "prediction is made about which is larger.",
    "S3": "SECONDARY. Among loop positions, introducing K or R is penalised LESS "
          "on the cytosolic side than on the non-cytosolic side (positive-inside).",
    "S4": "SECONDARY. The gradient survives dividing each position's z by that "
          "position's own spread across all 19 substitutions, which bounds the "
          "multiplicative-conservation confound.",
    "failure_policy": "Any declared expectation that is not met is reported as a "
                      "failure in the verdicts block and is not reinterpreted.",
}


def run_scorer(tag: str, S: np.ndarray, sites: dict) -> dict:
    log(f"  computing within-position contrasts for {tag}")
    c = contrasts(S, sites["nat"].astype(np.int64))
    n_prot = len(sites["prot_names"])
    u = c["usable"]
    log(f"    usable sites: {int(u.sum()):,} of {len(u):,} "
        f"({100.0 * u.sum() / len(u):.1f}%)")

    res: dict = {}
    for name in ("z", "z_KR", "z_DE"):
        m = u & ~np.isnan(c[name])
        res[name] = analyse(sites["prot"][m], sites["nat"][m].astype(np.int64),
                            sites["grp"][m], c[name][m], n_prot, name)
    m = c["usable_scaled"] & ~np.isnan(c["z_scaled"])
    res["z_scaled"] = analyse(sites["prot"][m], sites["nat"][m].astype(np.int64),
                              sites["grp"][m], c["z_scaled"][m], n_prot, "z_scaled")

    lm = u & ~np.isnan(c["z_KR"]) & (sites["grp"] == LOOP_GROUP) & (sites["side"] > 0)
    res["S3_loop_positive_inside_zKR"] = loop_side_test(
        sites["prot"][lm], sites["nat"][lm].astype(np.int64),
        sites["side"][lm], c["z_KR"][lm], n_prot)

    a = res["z"]
    log(f"    {tag} raw depth curve (mean z, interface -> centre):")
    for b in range(N_BINS):
        r = a["raw_curve"][f"bin{b}"]
        log(f"      bin {b}: {r['mean_z']:+.6f} "
            f"[{r['ci95'][0]:+.6f}, {r['ci95'][1]:+.6f}]  n={r['n_sites']:,}")
    log(f"      loop : {a['raw_loop']['mean_z']:+.6f} "
        f"[{a['raw_loop']['ci95'][0]:+.6f}, {a['raw_loop']['ci95'][1]:+.6f}]  "
        f"n={a['raw_loop']['n_sites']:,}")
    log(f"    {tag} composition-adjusted curve (Delta vs interface bin):")
    for b in range(1, N_BINS):
        r = a["adjusted_curve_vs_interface"][f"bin{b}"]
        log(f"      bin {b}: {r['delta']:+.6f} "
            f"[{r['ci95'][0]:+.6f}, {r['ci95'][1]:+.6f}]  "
            f"std {r['delta_standardized_within_this_scorer']:+.4f}  "
            f"strata={r['n_strata']:,}")
    return res


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t_start = time.time()
    log("B9 charge burial depth -- label-free positive-control rule audit")
    log("declared predictions are fixed in DECLARED and written to the JSON verbatim")

    log("loading topology census")
    census = read_topology(ROOT / "data" / "raw" / "topology_viewer.html")
    log(f"  {len(census):,} proteins")

    log("enumerating eligible sites")
    sites = enumerate_sites(census)
    index = build_index(sites)

    results: dict = {}
    composition = {}
    for b in range(N_BINS):
        m = sites["grp"] == b
        cnt = np.bincount(sites["nat"][m].astype(np.int64), minlength=len(UNCHARGED))
        composition[f"bin{b}"] = {UNCHARGED[i]: int(cnt[i]) for i in range(len(UNCHARGED))
                                  if cnt[i]}

    log("AlphaMissense")
    S_am = load_or_extract("am", sites, index, extract_am)
    results["AlphaMissense"] = run_scorer("AlphaMissense", S_am, sites)
    del S_am

    partial = {
        "rule": "charge burial depth",
        "declared_before_running": DECLARED,
        "status": "AlphaMissense complete; ESM1b pending",
        "results": results,
        "verdicts": {"AlphaMissense": verdicts(results["AlphaMissense"])},
    }
    (OUT / "charge_burial_result.json").write_text(
        json.dumps(partial, indent=2, default=str), encoding="utf-8")
    log("  wrote interim JSON (AlphaMissense only)")

    log("ESM1b")
    S_esm = load_or_extract("esm", sites, index, extract_esm)
    results["ESM1b"] = run_scorer("ESM1b", S_esm, sites)
    del S_esm

    v = {k: verdicts(r) for k, r in results.items()}
    agree = (v["AlphaMissense"]["P1_monotonic_rise_adjusted_curve"]
             == v["ESM1b"]["P1_monotonic_rise_adjusted_curve"]
             and v["AlphaMissense"]["P2_centre_exceeds_interface"]
             == v["ESM1b"]["P2_centre_exceeds_interface"])
    p3 = "MET" if all(v[s]["P1_monotonic_rise_adjusted_curve"] == "MET"
                      and v[s]["P2_centre_exceeds_interface"] == "MET"
                      for s in ("AlphaMissense", "ESM1b")) else "NOT MET"
    disagree = sorted(k for k in v["AlphaMissense"]
                      if not k.endswith("_detail")
                      and v["AlphaMissense"][k] != v["ESM1b"][k])

    not_met = {s: {k: val for k, val in v[s].items()
                   if not k.endswith("_detail") and str(val).startswith(("NOT MET",
                                                                         "MET BUT"))}
               for s in v}
    not_met = {s: d for s, d in not_met.items() if d}

    summary = {
        "rule": "charge burial depth",
        "script": "scripts/beat_am/b9_charge_burial.py",
        "run_finished": time.strftime("%Y-%m-%d %H:%M:%S"),
        "declared_before_running": DECLARED,
        "label_free": "topology census (sequence, TM boundaries, loop sides) plus "
                      "the scorers' own outputs. No ClinVar, humsavar, gnomAD, HGMD "
                      "or deep mutational scanning data was read.",
        "design": {
            "eligible_sites": f"positions inside a TM segment of length >= {MIN_TM_LEN} "
                              f"whose native residue is one of {UNCHARGED} (15 "
                              "uncharged residues; H excluded as ambiguous)",
            "within_position_contrast": "z = mean(score to D,E,K,R) - mean(score to "
                                        "the other uncharged residues at the same "
                                        "position). Cancels additive per-position "
                                        "effects; does NOT cancel multiplicative ones.",
            "depth": "d = 1 - |(p-s) - (L-1)/2| / ((L-1)/2); 0 at either interface, "
                     "1 at the helix centre; five equal-width bins",
            "stratification": "protein x native-residue strata, stratum weight "
                              "w = n1*n0/(n1+n0)",
            "bootstrap": f"{N_BOOT} draws over WHOLE PROTEINS, seed {SEED}; all "
                         "quantities in a draw share one resample",
            "loop_anchor": f"outside EVERY annotated TM segment regardless of that "
                           f"segment's length, known side, capped at "
                           f"{MAX_LOOP_PER_STRATUM} per protein x native residue",
        },
        "counts": {
            "proteins_in_census": len(census),
            "eligible_tm_segments": sites["n_tm_eligible_segments"],
            "positions_scored": int(len(sites["pos"])),
            "proteins_contributing": len(sites["prot_names"]),
            "native_residue_composition_by_bin": composition,
        },
        "cross_scorer_comparison_caveat":
            "AlphaMissense pathogenicity and ESM1b damage (= -LLR) have different "
            "scales and different spreads. Standardised (SD-unit) values are NOT "
            "comparable between the two scorers and no such comparison is made here. "
            "Standardised values are reported only to make effect sizes readable "
            "within one scorer. The cross-scorer claims in this file concern "
            "DIRECTION and MONOTONICITY only, which are scale-free.",
        "limitations": [
            "The within-position contrast removes additive per-position effects but "
            "not multiplicative ones. If helix centres are simply more conserved and "
            "the scorer stretches its whole range at conserved positions, z would "
            "rise with depth for reasons unrelated to charge. S4 bounds this.",
            "Fractional depth is a proxy for burial depth, not a measured coordinate. "
            "Tilted or kinked helices are mis-assigned; that is noise, which flattens "
            "a gradient rather than creating one.",
            "The loop anchor compares two different structural environments at once "
            "and is reported as secondary for that reason.",
            "Loop positions were capped at 10 per protein x native residue. The cap "
            "is unbiased for a stratum mean but reduces precision on the loop arm.",
        ],
        "results": results,
        "verdicts": v,
        "P3_both_scorers_pass": p3,
        "scorers_agree_on_P1_and_P2": bool(agree),
        "declared_items_where_the_two_scorers_differ": disagree,
        "declared_expectations_not_met": not_met or "none - every declared "
                                                    "expectation was met",
        "note": "This script measures only. No correction parameter was fitted and "
                "nothing was trained.",
    }
    (OUT / "charge_burial_result.json").write_text(
        json.dumps(summary, indent=2, default=str), encoding="utf-8")

    log("")
    log("VERDICTS")
    for scorer in ("AlphaMissense", "ESM1b"):
        log(f"  {scorer}")
        for k, val in v[scorer].items():
            if k.endswith("_detail"):
                continue
            log(f"    {k:46s} {val}")
    log(f"  P3 both scorers pass: {p3}")
    log(f"  scorers agree on P1 and P2: {agree}")
    log(f"  declared items where the scorers differ: {disagree or 'none'}")
    log(f"wrote {OUT / 'charge_burial_result.json'}")
    log(f"total runtime {time.time() - t_start:.0f}s")
    (OUT / "run.log").write_text("\n".join(_LOG_LINES) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

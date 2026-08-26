"""Fit a standalone variant scorer and evaluate a positive-inside correction."""

from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
from sklearn.decomposition import PCA
from sklearn.model_selection import GroupKFold

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "scripts" / "beat_am"))

FEAT = ROOT / "results" / "beat_am_20261016" / "features"
EMB = (ROOT / "results" / "beat_am_20261016" / "_from_instance_gsq37aoc"
       / "results" / "beat_am_20261016" / "features" / "embeddings_prot_t5_all")
PANEL_EMB = ROOT / "results" / "route_b_20261012" / "site_embeddings_fp16.npy"
FROZEN = (ROOT / "results" / "beat_am_20261016" / "inventory_frozen"
          / "group_manifest_frozen.csv")
OUT = ROOT / "results" / "beat_am_20261016" / "predictor_free"

from a17_general_scorer import GENERAL, _backend, per_protein  # noqa: E402
from a9_l1_features import (GR_C, GR_P, GR_V, HYDRO, CHARGE, HELIX,  # noqa: E402
                            STRAND, grantham, blosum62)


DROPPED = ("am", "esm1b")
PREDICTOR_FREE = [f for f in GENERAL if f not in DROPPED]
KEEP_IDX = [i for i, f in enumerate(GENERAL) if f not in DROPPED]

SEED = 20261016
N_PC = 64
N_OUTER, N_INNER = 5, 4
GRID = [
    {"max_depth": 3, "learning_rate": 0.10, "n_estimators": 300},
    {"max_depth": 6, "learning_rate": 0.05, "n_estimators": 500},
]
MARGIN = 0.03
MARGIN_FIT = 0.045     # Fit the 1.5-times margin on training groups only, as in the frozen A21 rule.
DECAY_LENGTH = 8.0
PROBE = 0.05
CHUNK_SITES = 40_000
CONFIRM_ROLES = {"confirmation_membrane", "confirmation_allhuman"}
NEUTRAL = "ACFGILMNPQSTVWY"   # Use the frozen audit's 15-residue uncharged comparator.


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)



def load_dev():
    t = np.load(FEAT / "l1_table.npz", allow_pickle=True)
    d = {k: t[k] for k in t.files}
    X_all = np.column_stack([d[f].astype(np.float32) for f in GENERAL])
    X = X_all[:, KEEP_IDX]
    acc, pos = d["accession"].astype(str), d["position"].astype(int)

    E = np.zeros((len(pos), 1024), dtype=np.float32)
    missing = 0
    for a in np.unique(acc):
        f = EMB / f"{a}.npz"
        if not f.exists():
            missing += int((acc == a).sum())
            continue
        z = np.load(f, allow_pickle=True)
        tab = np.zeros((int(z["seq_len"]) + 1, 1024), dtype=np.float32)
        tab[z["position"].astype(int)] = z["embedding"].astype(np.float32)
        sel = np.flatnonzero(acc == a)
        E[sel] = tab[np.clip(pos[sel], 0, tab.shape[0] - 1)]
    if missing:
        log(f"  WARNING {missing:,} rows without embeddings")
    return {"X": X, "E": E, "y": d["target"].astype(np.float32),
            "am": d["am"].astype(np.float32), "groups": d["group_id"].astype(str),
            "accession": acc, "gene": d["gene"].astype(str),
            "dms_id": d["dms_id"].astype(str),
            "is_membrane": d["is_membrane"].astype(int),
            "distance": d["distance"].astype(np.float32),
            "side": d["side"].astype(int),
            "native": d["native"].astype(str), "alt": d["alt"].astype(str)}


def development_cv(data, make_model):
    """Out-of-fold predictions under the identical nested grouped protocol."""
    X, E, y, groups = data["X"], data["E"], data["y"], data["groups"]
    pred = np.full(len(y), np.nan, dtype=np.float32)
    chosen = []
    for fold, (tr, te) in enumerate(
            GroupKFold(n_splits=N_OUTER).split(X, y, groups), start=1):
        pca = PCA(n_components=N_PC, random_state=SEED).fit(E[tr])
        Xtr = np.hstack([X[tr], pca.transform(E[tr])]).astype(np.float32)
        Xte = np.hstack([X[te], pca.transform(E[te])]).astype(np.float32)
        g_tr = groups[tr]
        inner = GroupKFold(n_splits=min(N_INNER, len(np.unique(g_tr))))
        best, best_cfg = -np.inf, GRID[0]
        for cfg in GRID:
            ip = np.full(len(tr), np.nan)
            for itr, iva in inner.split(Xtr, y[tr], g_tr):
                m = make_model(cfg, SEED + fold)
                m.fit(Xtr[itr], y[tr][itr])
                ip[iva] = m.predict(Xtr[iva])
            full = np.full(len(y), np.nan)
            full[tr] = ip
            mask = np.zeros(len(y), bool)
            mask[tr] = True
            dd = np.array([r["difference"] for r in per_protein(full, data, mask)])
            if dd.size and dd.mean() > best:
                best, best_cfg = dd.mean(), cfg
        m = make_model(best_cfg, SEED + fold)
        m.fit(Xtr, y[tr])
        pred[te] = m.predict(Xte)
        chosen.append({"fold": fold, "inner_mean_vs_am": round(float(best), 6),
                       **best_cfg})
        log(f"  fold {fold}: inner {best:+.4f} cfg={best_cfg}")
    return pred, chosen


def summarise(rows):
    d = np.array([r["difference"] for r in rows])
    return {"n": len(d), "mean_vs_am": round(float(d.mean()), 6),
            "median_vs_am": round(float(np.median(d)), 6),
            "fraction_improved": round(float((d > 0).mean()), 4),
            "absolute_model": round(float(np.mean([r["model"] for r in rows])), 6),
            "absolute_am": round(float(np.mean([r["am"] for r in rows])), 6)}



def fit_full(data, make_model):
    """PCA + model on the entire fitting pool, for scoring the census panel."""
    pca = PCA(n_components=N_PC, random_state=SEED).fit(data["E"])
    Xf = np.hstack([data["X"], pca.transform(data["E"])]).astype(np.float32)
    model = make_model({"max_depth": 6, "learning_rate": 0.05,
                        "n_estimators": 500}, SEED)
    model.fit(Xf, data["y"])
    log(f"  full-pool scorer fitted on {len(Xf):,} rows x {Xf.shape[1]} features")
    return model, pca


def score_panel(lib, panel, model, pca, census, accs):
    """Same construction as a18, with the two predictor columns omitted."""
    n_sites = panel["scores"].shape[0]
    alphabet = list(lib.ALPHABET)

    n_tmds = np.array([len((census.get(a) or {}).get("tmds") or []) for a in accs],
                      dtype=np.float32)
    loop_len_of = {}
    for a in accs:
        rec = census.get(a)
        if not rec:
            continue
        m = {}
        for _lid, loop in (rec.get("loops") or {}).items():
            s, e = int(loop["start"]), int(loop["end"])
            for p in range(s, e + 1):
                m[p] = e - s + 1
        loop_len_of[a] = m

    emb = np.load(PANEL_EMB, mmap_mode="r")
    assert emb.shape[0] == n_sites, (emb.shape, n_sites)

    def tab(dct):
        return np.array([dct[a] for a in alphabet], dtype=np.float32)
    T_C, T_P, T_V = tab(GR_C), tab(GR_P), tab(GR_V)
    T_H, T_Q, T_HE, T_ST = tab(HYDRO), tab(CHARGE), tab(HELIX), tab(STRAND)
    BL = np.array([[blosum62().get((a, b), 0.0) for b in alphabet]
                   for a in alphabet], dtype=np.float32)
    GRM = np.array([[grantham(a, b) for b in alphabet] for a in alphabet],
                   dtype=np.float32)
    KR = np.array([1.0 if a in "KR" else 0.0 for a in alphabet], dtype=np.float32)

    scores = np.zeros((n_sites, len(alphabet)), dtype=np.float32)
    for start in range(0, n_sites, CHUNK_SITES):
        end = min(start + CHUNK_SITES, n_sites)
        ns = end - start
        prot = panel["protein"][start:end]
        pos = panel["position"][start:end]
        ref = panel["ref_aa"][start:end]
        side = panel["side"][start:end].astype(np.float32)
        dist = panel["distance"][start:end].astype(np.float32)
        tmi = panel["tmd_index"][start:end].astype(np.float32)
        pc = pca.transform(np.asarray(emb[start:end], dtype=np.float32))
        ll = np.array([loop_len_of.get(accs[p], {}).get(int(q), 0)
                       for p, q in zip(prot, pos)], dtype=np.float32)
        nt = n_tmds[prot]

        rows = []
        for a_i in range(len(alphabet)):
            f = np.column_stack([
                BL[ref, a_i], GRM[ref, a_i],
                T_C[a_i] - T_C[ref], T_P[a_i] - T_P[ref], T_V[a_i] - T_V[ref],
                T_H[a_i] - T_H[ref], T_Q[a_i] - T_Q[ref],
                T_HE[a_i] - T_HE[ref], T_ST[a_i] - T_ST[ref],
                KR[ref], np.full(ns, KR[a_i], dtype=np.float32),
                side, dist, np.zeros(ns, dtype=np.float32), tmi, ll, nt,
                np.ones(ns, dtype=np.float32),
            ]).astype(np.float32)
            rows.append(np.hstack([f, pc]).astype(np.float32))
        block = model.predict(np.vstack(rows)).astype(np.float32)
        scores[start:end] = block.reshape(len(alphabet), ns).T
        if (start // CHUNK_SITES) % 5 == 0:
            log(f"    scored {end:,}/{n_sites:,}")
    return scores


def census_labels(lib, census, kind):
    split = lib.cluster_split(census, kind=kind)
    label = dict(split["label"])
    man = {}
    with FROZEN.open(encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh):
            man[r["accession"]] = r
    cg = {r["frozen_group"] for r in man.values() if r["role"] in CONFIRM_ROLES}
    barred = {a for a, r in man.items() if r["frozen_group"] in cg}
    removed = sum(1 for a in label if a in barred and label[a] == "train")
    for a in list(label):
        if a in barred and label[a] == "train":
            label[a] = "excluded_confirmation_relative"
    return label, removed


def window_of(lib, sub):
    win = lib.window_panel(sub, window=10)
    fk = (sub["protein"].astype(np.int64) << 32) | sub["position"].astype(np.int64)
    wk = (win["protein"].astype(np.int64) << 32) | win["position"].astype(np.int64)
    o = np.argsort(fk, kind="stable")
    loc = np.searchsorted(fk[o], wk)
    assert (fk[o][loc] == wk).all()
    return win, o[loc]



def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    make_model, backend = _backend()
    log(f"backend: {backend}")
    log(f"dropped inputs: {DROPPED}")
    log(f"features: {len(PREDICTOR_FREE)} tabular + {N_PC} ProtT5 PCs "
        f"= {len(PREDICTOR_FREE) + N_PC}")

    data = load_dev()
    log(f"development rows {len(data['y']):,}  "
        f"proteins {len(np.unique(data['accession']))}")


    log("development CV (nested, homology-grouped)")
    pred, chosen = development_cv(data, make_model)
    rows = per_protein(pred, data)
    memb = [r for r in rows if r["membrane"] == 1]
    non = [r for r in rows if r["membrane"] == 0]
    flank = ((data["is_membrane"] == 1) & (data["distance"] >= 1)
             & (data["distance"] <= 60))
    frows = per_protein(pred, data, flank)
    dev = {"all": summarise(rows), "membrane": summarise(memb),
           "non_membrane": summarise(non),
           "flank_le60": summarise(frows) if frows else None}
    with (OUT / "development_per_protein.csv").open("w", encoding="utf-8",
                                                    newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    for k, v in dev.items():
        if v:
            log(f"  {k:14s} n={v['n']:3d} vs AM {v['mean_vs_am']:+.4f}  "
                f"abs {v['absolute_model']:.4f} (AM {v['absolute_am']:.4f})")


    import form_ab_lib as lib
    from run_gate1 import read_topology

    census = read_topology(ROOT / "data" / "raw" / "topology_viewer.html")
    panel = lib.P.load_panel(lib.AM_ARCHIVE, 60)
    accs = [str(a) for a in panel["accessions"]]
    log(f"panel: {panel['scores'].shape[0]:,} sites x "
        f"{panel['scores'].shape[1]} alternatives")

    cache = OUT / "panel_scores.npz"
    if cache.exists():
        log("  reusing cached panel scores")
        scores = np.load(cache)["scores"]
    else:
        model, pca = fit_full(data, make_model)
        scores = score_panel(lib, panel, model, pca, census, accs)
        np.savez_compressed(cache, scores=scores)
        log("  panel scored and cached")


    fitted = {}
    for kind in ("sequence", "symbol"):
        label, removed = census_labels(lib, census, kind)
        tr = lib.split_mask(panel, label, "train")
        sub = lib.subset_panel(panel, tr)
        win, sel = window_of(lib, sub)
        s = scores[tr][sel]
        a0 = lib.audit_matrix(win, s, draws=200)[0]
        delta0, sd = float(a0["estimate_raw"]), float(a0["penalty_sd"])
        probed = lib.apply_correction(win, PROBE, DECAY_LENGTH)
        dp = float(lib.audit_matrix(win, probed, draws=200)[0]["estimate_raw"])
        k = (delta0 - dp) / PROBE
        beta, target = lib.fit_beta(delta0, k, sd, MARGIN_FIT)
        fitted[kind] = {"train_removed": removed, "delta0_raw": delta0,
                        "delta0_standardized": delta0 / sd, "penalty_sd": sd,
                        "k": k, "raw_target": target, "beta": beta}
        log(f"  {kind:9s} train: delta0 {delta0:+.6f} ({delta0/sd:+.5f} SD) "
            f"sd {sd:.6f} k {k:.4f} -> beta {beta:.6f}")

    beta_joint = max(f["beta"] for f in fitted.values())
    log(f"joint beta = max = {beta_joint:.6f}")
    corrected = lib.apply_correction({**panel, "scores": scores},
                                     beta_joint, DECAY_LENGTH)


    held = {}
    log(f"{'split':10s} {'model':22s} {'raw':>10s} {'std':>9s} {'ci_upper':>10s}  pass")
    for kind in ("sequence", "symbol"):
        label, _ = census_labels(lib, census, kind)
        te = lib.split_mask(panel, label, "test")
        sub = lib.subset_panel(panel, te)
        win, sel = window_of(lib, sub)
        entry = {}
        for name, mat in (("alphamissense", panel["scores"][te]),
                          ("predictor_free_pre", scores[te]),
                          ("predictor_free_post", corrected[te])):
            a = lib.audit_matrix(win, mat[sel], draws=200)[0]
            raw, sd = float(a["estimate_raw"]), float(a["penalty_sd"])
            lo, hi = (float(x) for x in a["ci95_raw"])
            ok = bool(raw / sd <= -MARGIN and hi < 0)
            entry[name] = {"contrast_raw": raw, "penalty_sd": sd,
                           "contrast_standardized": raw / sd,
                           "ci95_raw": [lo, hi],
                           "ci95_standardized": [lo / sd, hi / sd],
                           "ci95_upper_below_zero": hi < 0,
                           "passes_full_condition": ok,
                           "n_sites": int(a["n_sites"]),
                           "n_proteins": int(a["n_proteins"])}
            log(f"{kind:10s} {name:22s} {raw:+10.6f} {raw/sd:+9.5f} {hi:+10.6f}  {ok}")
        held[kind] = entry

    both = all(held[k]["predictor_free_post"]["passes_full_condition"]
               for k in ("sequence", "symbol"))
    fid = lib.fidelity(scores, corrected)

    # The feature head encodes cytosol as 1, whereas the census panel encodes it as 0.
    kr_gain = np.array([(n in NEUTRAL) and (a in "KR")
                        for n, a in zip(data["native"], data["alt"])])
    sign = np.where(data["side"] == 1, -1.0, np.where(data["side"] == 2, 1.0, 0.0))
    shift = np.where(kr_gain & (data["distance"] >= 1),
                     sign * (beta_joint / 2)
                     * lib.decay(np.maximum(data["distance"], 0), DECAY_LENGTH),
                     0.0).astype(np.float32)
    post = pred + shift
    log(f"  beta {beta_joint:.6f}; dev rows shifted "
        f"{float((shift != 0).mean()):.4%}; max |shift| {float(np.abs(shift).max()):.6f}")

    kr_flank = kr_gain & (data["is_membrane"] == 1) & (data["distance"] >= 1) \
        & (data["distance"] <= 60)
    cost = {}
    for name, extra in (("all", None), ("membrane", data["is_membrane"] == 1),
                        ("non_membrane", data["is_membrane"] == 0),
                        ("flank_le60", flank), ("kr_gain_flank", kr_flank)):
        pre_rows = per_protein(pred, data, extra)
        post_rows = per_protein(post, data, extra)
        if not pre_rows:
            continue
        pre_m = float(np.mean([r["difference"] for r in pre_rows]))
        post_m = float(np.mean([r["difference"] for r in post_rows]))
        worse = sum(1 for a, b in zip(pre_rows, post_rows)
                    if b["model"] < a["model"] - 1e-12)
        cost[name] = {"n": len(pre_rows), "pre_mean_vs_am": round(pre_m, 6),
                      "post_mean_vs_am": round(post_m, 6),
                      "rule_cost": round(post_m - pre_m, 6),
                      "absolute_pre": round(float(np.mean(
                          [r["model"] for r in pre_rows])), 6),
                      "absolute_post": round(float(np.mean(
                          [r["model"] for r in post_rows])), 6),
                      "n_proteins_worsened": worse}
        log(f"  cost {name:14s} n={len(pre_rows):3d} {pre_m:+.5f} -> {post_m:+.5f} "
            f"({post_m - pre_m:+.6f})")


    out = {
        "what_this_is": "a scorer containing no audited predictor, and "
                        "enforcement of the positive-inside rule on it",
        "dropped_inputs": list(DROPPED),
        "why_both_dropped": "ESM1b already satisfies the rule (-0.073 SD), so "
                            "leaving it in would let compliance arrive through "
                            "the frontend rather than the constraint",
        "features": PREDICTOR_FREE + [f"protT5_pc{i}" for i in range(N_PC)],
        "n_features": len(PREDICTOR_FREE) + N_PC,
        "backend": backend,
        "held_identical_to_a17_a18_a21": [
            "rows, homology groups, outer/inner CV protocol, seed",
            "hyperparameter grid and selection rule",
            "beta fitted on each split's own train partition, larger taken",
            f"MARGIN_FIT={MARGIN_FIT} (1.5x), verdict target {MARGIN}",
            f"decay length {DECAY_LENGTH}",
        ],
        "scope_note": "development-set accuracy only; no reserved protein was "
                      "read, no confirmation outcome opened",
        "development": dev,
        "configs_chosen": chosen,
        "beta_joint": beta_joint,
        "per_split_fit": fitted,
        "held_out": held,
        "passes_both_splits": both,
        "fidelity_full_panel": {k: (float(v) if np.isscalar(v) else str(v))
                                for k, v in fid.items()},
        "accuracy_cost_development": cost,
    }
    (OUT / "predictor_free_result.json").write_text(
        json.dumps(out, indent=2, default=str), encoding="utf-8")
    np.savez_compressed(OUT / "panel_scores_corrected.npz", scores=corrected)

    log("")
    log(f"passes BOTH splits: {both}")
    log(f"fidelity vs uncorrected: {fid}")
    log(f"wrote {OUT / 'predictor_free_result.json'}")


if __name__ == "__main__":
    main()

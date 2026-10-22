"""A17: a GENERAL scorer, applicable to any missense variant with no assay metadata.

Why this exists. The best development model consumes two features a census panel
site cannot supply:

  * `endpoint_type` one-hot -- a panel site belongs to no assay, so it is undefined
  * the local-window embedding mean -- needs neighbours, and the panel stops at
    distance 60, so the window is truncated at the edges

The audit needs scores for all 726,444 panel sites, and the deliverable requires
`score_variants.py` to score *any* missense variant list. Both demand a scorer
with no assay-dependent and no neighbour-dependent inputs. This trains that
scorer and measures what the generality costs, so the price is visible rather
than buried.

Feature set: AM, ESM1b (damage-oriented), the substitution descriptors, the
per-site topology, and a PCA projection of the ProtT5 vector AT the position.
No endpoint, no window.

PCA is fit on outer-training rows only, per locked decision 5.

Scope: fitting pool only. Confirmation proteins are absent from l1_table.npz.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.model_selection import GroupKFold

ROOT = Path(__file__).resolve().parents[2]
FEAT = ROOT / "results" / "beat_am_20261016" / "features"
EMB = (ROOT / "results" / "beat_am_20261016" / "_from_instance_gsq37aoc"
       / "results" / "beat_am_20261016" / "features" / "embeddings_prot_t5_all")
OUT = ROOT / "results" / "beat_am_20261016" / "evaluation"

N_OUTER, N_INNER = 5, 4
SEED = 20261016
N_PC = 64

# the general feature set: no endpoint one-hot, no window term
GENERAL = [
    "am", "esm1b", "blosum62", "grantham",
    "d_composition", "d_polarity", "d_volume", "d_hydro", "d_charge",
    "d_helix", "d_strand", "native_is_kr", "alt_is_kr",
    "side", "distance", "in_tm", "tmd_index", "loop_len", "n_tmds",
    "is_membrane",
]
GRID = [
    {"max_depth": 3, "learning_rate": 0.10, "n_estimators": 300},
    {"max_depth": 6, "learning_rate": 0.05, "n_estimators": 500},
]


def load():
    t = np.load(FEAT / "l1_table.npz", allow_pickle=True)
    d = {k: t[k] for k in t.files}
    X = np.column_stack([d[f].astype(np.float32) for f in GENERAL])
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
        print(f"  WARNING {missing:,} rows without embeddings")
    return {"X": X, "E": E, "y": d["target"].astype(np.float32),
            "am": d["am"].astype(np.float32), "groups": d["group_id"].astype(str),
            "accession": acc, "gene": d["gene"].astype(str),
            "dms_id": d["dms_id"].astype(str),
            "is_membrane": d["is_membrane"].astype(int),
            "distance": d["distance"].astype(np.float32)}


def per_protein(pred, data, extra=None):
    rows = []
    for a in np.unique(data["accession"]):
        sel = data["accession"] == a
        if extra is not None:
            sel = sel & extra
        m, s = [], []
        for dm in np.unique(data["dms_id"][sel]):
            k = sel & (data["dms_id"] == dm)
            if k.sum() < 20 or np.unique(data["y"][k]).size < 3:
                continue
            m.append(spearmanr(pred[k], data["y"][k]).statistic)
            s.append(spearmanr(data["am"][k], data["y"][k]).statistic)
        if m:
            rows.append({"accession": a, "gene": data["gene"][sel][0],
                         "membrane": int(data["is_membrane"][sel][0]),
                         "model": float(np.mean(m)), "am": float(np.mean(s)),
                         "difference": float(np.mean(m) - np.mean(s))})
    return rows


def _backend():
    """lightgbm if its native library loads, else sklearn's equivalent.

    The two agree to 0.0002 on the L1 headline (verified on the instance), and
    lightgbm's DLL does not load on this Windows host. Whichever is used is
    recorded in the summary.
    """
    try:
        import lightgbm as lgb

        def make(cfg, seed):
            return lgb.LGBMRegressor(random_state=seed, n_jobs=-1, verbose=-1, **cfg)
        return make, "lightgbm"
    except Exception:
        from sklearn.ensemble import HistGradientBoostingRegressor

        def make(cfg, seed):
            return HistGradientBoostingRegressor(
                max_depth=cfg["max_depth"], learning_rate=cfg["learning_rate"],
                max_iter=cfg["n_estimators"], random_state=seed,
                early_stopping=False)
        return make, "sklearn HistGradientBoostingRegressor"


def main() -> None:
    make_model, backend = _backend()
    print(f"backend: {backend}", flush=True)

    OUT.mkdir(parents=True, exist_ok=True)
    data = load()
    X, E, y, groups = data["X"], data["E"], data["y"], data["groups"]
    print(f"rows {len(y):,}  general features {X.shape[1]}  PCA {N_PC}", flush=True)

    pred = np.full(len(y), np.nan, dtype=np.float32)
    for fold, (tr, te) in enumerate(
            GroupKFold(n_splits=N_OUTER).split(X, y, groups), start=1):
        pca = PCA(n_components=N_PC, random_state=SEED)
        pca.fit(E[tr])
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
        print(f"  fold {fold}: inner {best:+.4f} cfg={best_cfg}", flush=True)

    rows = per_protein(pred, data)
    memb = [r for r in rows if r["membrane"] == 1]
    non = [r for r in rows if r["membrane"] == 0]
    flank = ((data["is_membrane"] == 1) & (data["distance"] >= 1)
             & (data["distance"] <= 60))
    frows = per_protein(pred, data, flank)

    def st(rs):
        d = np.array([r["difference"] for r in rs])
        return {"n": len(d), "mean": round(float(d.mean()), 6),
                "median": round(float(np.median(d)), 6),
                "fraction_improved": round(float((d > 0).mean()), 4),
                "worst": round(float(d.min()), 6),
                "absolute_model": round(float(np.mean([r["model"] for r in rs])), 6),
                "absolute_am": round(float(np.mean([r["am"] for r in rs])), 6)}

    with (OUT / "general_scorer_per_protein.csv").open("w", encoding="utf-8",
                                                       newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    FULL = {"all": 0.0350, "membrane": 0.0106, "non_membrane": 0.0472,
            "flank_le60": -0.0328}
    summary = {
        "purpose": "a scorer with no assay-dependent and no neighbour-dependent "
                   "features, so it can be applied to every census panel site and "
                   "to any arbitrary missense variant list",
        "dropped_vs_best_dev_model": ["endpoint_type one-hot (undefined for a "
                                      "panel site)",
                                      "local-window embedding mean (truncated at "
                                      "the panel edge)"],
        "features": GENERAL + [f"protT5_pc{i}" for i in range(N_PC)],
        "all": st(rows), "membrane": st(memb), "non_membrane": st(non),
        "flank_le60": st(frows) if frows else None,
        "backend": backend,
        "reference_full_feature_model_pc64": FULL,
        "cost_of_generality": {
            k: round(st(v)["mean"] - FULL[k], 6)
            for k, v in (("all", rows), ("membrane", memb),
                         ("non_membrane", non), ("flank_le60", frows)) if v
        },
    }
    (OUT / "general_scorer_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    print(f"\n{'scope':14s} {'n':>3s} {'mean vs AM':>11s} {'abs model':>10s} "
          f"{'abs AM':>8s}   full-feature was   cost")
    for name, rs in (("all", rows), ("membrane", memb),
                     ("non-membrane", non), ("flank<=60", frows)):
        if not rs:
            continue
        s = st(rs)
        key = {"all": "all", "membrane": "membrane",
               "non-membrane": "non_membrane", "flank<=60": "flank_le60"}[name]
        print(f"{name:14s} {s['n']:3d} {s['mean']:+11.4f} {s['absolute_model']:+10.4f} "
              f"{s['absolute_am']:+8.4f}   {FULL[key]:+.4f}   "
              f"{s['mean'] - FULL[key]:+.4f}")


if __name__ == "__main__":
    main()

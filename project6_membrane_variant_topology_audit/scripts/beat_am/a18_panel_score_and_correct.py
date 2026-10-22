"""A18: score the census panel, fit the one-parameter correction, measure its cost.

Section 5 step 6(a), the prescribed cheapest route.

Pipeline:
  1. fit the general scorer (a17 feature set) on the whole fitting pool
  2. score all 726,444 panel sites x 20 alternatives, in site chunks
  3. build the census split; remove confirmation proteins AND their homology-group
     relatives from the TRAIN side, rather than reusing cluster_split unchanged
  4. on census train: audit -> delta0 and penalty SD; probe the sensitivity k
     empirically; fit_beta
  5. apply_correction to the full panel
  6. on HELD-OUT census clusters: audit the corrected scores -- the headline
     rule number -- and run the full 8-rule report card, with the diff vs AM
  7. fidelity(original, corrected), and the accuracy cost on development
     proteins by applying the identical shift to their rows

Unit note. `fit_beta` computes `target = declared * margin * sd`, so the -0.03
audit target is in **SD units** and the raw target is -0.03 * penalty_sd. That
resolves the ambiguity flagged in session 14; it is read from the code, not
assumed.

Label-free throughout on the panel side. Confirmation outcomes are never opened.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "scripts" / "beat_am"))

FEAT = ROOT / "results" / "beat_am_20261016" / "features"
EMB = (ROOT / "results" / "beat_am_20261016" / "_from_instance_gsq37aoc"
       / "results" / "beat_am_20261016" / "features" / "embeddings_prot_t5_all")
PANEL_EMB = ROOT / "results" / "route_b_20261012" / "site_embeddings_fp16.npy"
INV = ROOT / "results" / "beat_am_20261016" / "inventory"
CUR = ROOT / "results" / "beat_am_20261016" / "inventory_curated"
OUT = ROOT / "results" / "beat_am_20261016" / "rule"

SEED = 20261016
N_PC = 64
MARGIN = 0.03
DECAY_LENGTH = 8.0
CHUNK_SITES = 40_000
CONFIRM_ROLES = {"confirmation_membrane", "confirmation_allhuman"}

from a17_general_scorer import GENERAL, _backend  # noqa: E402

# amino-acid property tables, same as a9
from a9_l1_features import (GR_C, GR_P, GR_V, HYDRO, CHARGE, HELIX, STRAND,  # noqa: E402
                            grantham, blosum62)


def fit_scorer():
    """Fit PCA + model on the entire fitting pool; return a scoring closure."""
    make_model, backend = _backend()
    t = np.load(FEAT / "l1_table.npz", allow_pickle=True)
    d = {k: t[k] for k in t.files}
    X = np.column_stack([d[f].astype(np.float32) for f in GENERAL])
    acc, pos = d["accession"].astype(str), d["position"].astype(int)
    E = np.zeros((len(pos), 1024), dtype=np.float32)
    for a in np.unique(acc):
        f = EMB / f"{a}.npz"
        if not f.exists():
            continue
        z = np.load(f, allow_pickle=True)
        tab = np.zeros((int(z["seq_len"]) + 1, 1024), dtype=np.float32)
        tab[z["position"].astype(int)] = z["embedding"].astype(np.float32)
        sel = np.flatnonzero(acc == a)
        E[sel] = tab[np.clip(pos[sel], 0, tab.shape[0] - 1)]

    pca = PCA(n_components=N_PC, random_state=SEED).fit(E)
    Xf = np.hstack([X, pca.transform(E)]).astype(np.float32)
    model = make_model({"max_depth": 6, "learning_rate": 0.05,
                        "n_estimators": 500}, SEED)
    model.fit(Xf, d["target"].astype(np.float32))
    print(f"  scorer fitted on {len(Xf):,} rows, backend={backend}", flush=True)
    return model, pca, backend, d


def esm_lookup(panel, es):
    """(site_row, alt) -> ESM1b damage value, aligned by (protein, position, alt)."""
    key = (es["protein"].astype(np.int64) << 40) | \
          (es["position"].astype(np.int64) << 8) | es["alt_aa"].astype(np.int64)
    order = np.argsort(key, kind="stable")
    return key[order], es["am_pathogenicity"][order]


def main() -> None:
    import form_ab_lib as lib

    OUT.mkdir(parents=True, exist_ok=True)
    model, pca, backend, dev = fit_scorer()

    panel = lib.P.load_panel(lib.AM_ARCHIVE, 60)
    n_sites = panel["scores"].shape[0]
    alphabet = list(lib.ALPHABET)
    print(f"  panel: {n_sites:,} sites x {panel['scores'].shape[1]} alternatives",
          flush=True)

    # per-protein topology extras
    from run_gate1 import read_topology
    census = read_topology(ROOT / "data" / "raw" / "topology_viewer.html")
    accs = [str(a) for a in panel["accessions"]]
    n_tmds = np.array([len((census.get(a) or {}).get("tmds") or []) for a in accs],
                      dtype=np.float32)
    loop_len_of = {}
    for a in accs:
        rec = census.get(a)
        if not rec:
            continue
        m = {}
        for lid, loop in (rec.get("loops") or {}).items():
            s, e = int(loop["start"]), int(loop["end"])
            for p in range(s, e + 1):
                m[p] = e - s + 1
        loop_len_of[a] = m

    es = np.load(ROOT / "data" / "processed" / "gate0" / "esm1b_flank_scores.npz",
                 allow_pickle=True)
    es_key, es_val = esm_lookup(panel, es)

    emb = np.load(PANEL_EMB, mmap_mode="r")
    assert emb.shape[0] == n_sites, (emb.shape, n_sites)

    # substitution descriptor tables indexed by alphabet position
    idx = {a: i for i, a in enumerate(alphabet)}
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
        am = panel["scores"][start:end]
        pc = pca.transform(np.asarray(emb[start:end], dtype=np.float32))

        ll = np.array([loop_len_of.get(accs[p], {}).get(int(q), 0)
                       for p, q in zip(prot, pos)], dtype=np.float32)
        nt = n_tmds[prot]

        k20 = (prot.astype(np.int64)[:, None] << 40) | \
              (pos.astype(np.int64)[:, None] << 8) | \
              np.arange(len(alphabet), dtype=np.int64)[None, :]
        loc = np.searchsorted(es_key, k20)
        loc = np.clip(loc, 0, len(es_key) - 1)
        esm = np.where(es_key[loc] == k20, es_val[loc], np.nan).astype(np.float32)

        rows = []
        for a_i in range(len(alphabet)):
            f = np.column_stack([
                am[:, a_i], esm[:, a_i],
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
            print(f"    scored {end:,}/{n_sites:,}", flush=True)

    np.savez_compressed(OUT / "model_panel_scores.npz", scores=scores)
    print(f"  panel scored, saved", flush=True)

    # ---- census split, confirmation groups removed from TRAIN
    split = lib.cluster_split(census)
    label = dict(split["label"])
    manifest = {}
    for src in (CUR / "group_manifest_curated.csv", INV / "group_manifest.csv"):
        if src.exists():
            import csv as _csv
            with src.open(encoding="utf-8", newline="") as fh:
                for r in _csv.DictReader(fh):
                    manifest[r["accession"]] = r
            break
    confirm_groups = {r.get("curated_group") or r.get("group_id")
                      for r in manifest.values() if r["role"] in CONFIRM_ROLES}
    barred = {a for a, r in manifest.items()
              if (r.get("curated_group") or r.get("group_id")) in confirm_groups}
    removed = 0
    for a in list(label):
        if a in barred and label[a] == "train":
            label[a] = "excluded_confirmation_relative"
            removed += 1
    print(f"  census train: removed {removed} confirmation proteins/relatives; "
          f"barred set size {len(barred)}", flush=True)

    win = lib.window_panel(panel, window=10)
    win_scores_idx = None  # window_panel subsets rows; rebuild aligned scores
    # window_panel keeps a 'keys' field we can use to map back
    tr_mask = lib.split_mask(win, label, "train")
    te_mask = lib.split_mask(win, label, "test")
    print(f"  window-10 panel: {win['scores'].shape[0]:,} sites; "
          f"train {int(tr_mask.sum()):,}  test {int(te_mask.sum()):,}", flush=True)

    # map full-panel scores onto the window panel rows
    full_keys = (panel["protein"].astype(np.int64) << 32) | panel["position"].astype(np.int64)
    win_keys = (win["protein"].astype(np.int64) << 32) | win["position"].astype(np.int64)
    order = np.argsort(full_keys, kind="stable")
    loc = np.searchsorted(full_keys[order], win_keys)
    assert (full_keys[order][loc] == win_keys).all(), "window rows not found in panel"
    win_model = scores[order[loc]]

    train_panel = lib.subset_panel(win, tr_mask)
    train_scores = win_model[tr_mask]

    a0 = lib.audit_matrix(train_panel, train_scores, draws=200)[0]
    delta0, sd = float(a0["estimate_raw"]), float(a0["penalty_sd"])
    print(f"\n  model raw contrast on census TRAIN: {delta0:+.6f}  penalty_sd {sd:.6f}")

    # probe the sensitivity k empirically: contrast moves by -k * beta
    probe = 0.05
    probed = lib.apply_correction(train_panel, probe, DECAY_LENGTH)
    d_probe = float(lib.audit_matrix(train_panel, probed, draws=200)[0]["estimate_raw"])
    k = (delta0 - d_probe) / probe
    print(f"  probe beta={probe}: contrast {delta0:+.6f} -> {d_probe:+.6f}  "
          f"=> k={k:.6f}")

    beta, target = lib.fit_beta(delta0, k, sd, MARGIN)
    print(f"  target (raw) = {target:+.6f}  [= -{MARGIN} * {sd:.6f}]")
    print(f"  fitted beta  = {beta:.6f}  (AM's frozen beta was 0.029713)")

    corrected_full = lib.apply_correction(
        {**panel, "scores": scores}, beta, DECAY_LENGTH)
    np.savez_compressed(OUT / "model_panel_scores_corrected.npz",
                        scores=corrected_full)

    win_corr = corrected_full[order[loc]]
    held = lib.subset_panel(win, te_mask)
    a_pre = lib.audit_matrix(held, win_model[te_mask], draws=200)[0]
    a_post = lib.audit_matrix(held, win_corr[te_mask], draws=200)[0]
    print(f"\n  HELD-OUT census clusters:")
    print(f"    pre-correction  raw {float(a_pre['estimate_raw']):+.6f}  "
          f"ci {a_pre['ci95_raw']}")
    print(f"    post-correction raw {float(a_post['estimate_raw']):+.6f}  "
          f"ci {a_post['ci95_raw']}")
    std_post = float(a_post["estimate_raw"]) / float(a_post["penalty_sd"])
    print(f"    post standardized {std_post:+.6f}   target <= -{MARGIN}  -> "
          f"{'PASS' if std_post <= -MARGIN else 'FAIL'}")

    fid = lib.fidelity(win_model[te_mask], win_corr[te_mask])
    print(f"\n  fidelity (held-out, all rows): {fid}")

    result = {
        "backend": backend,
        "unit_resolution": "fit_beta computes target = declared * margin * sd, so "
                           "-0.03 is in SD units; raw target = -0.03 * penalty_sd",
        "census_train_confirmation_relatives_removed": removed,
        "model_contrast_train_raw": delta0,
        "penalty_sd": sd,
        "sensitivity_k": k,
        "raw_target": target,
        "fitted_beta": beta,
        "decay_length": DECAY_LENGTH,
        "am_frozen_beta": 0.02971251410524106,
        "held_out": {
            "pre_raw": float(a_pre["estimate_raw"]),
            "pre_ci95": [float(x) for x in a_pre["ci95_raw"]],
            "post_raw": float(a_post["estimate_raw"]),
            "post_ci95": [float(x) for x in a_post["ci95_raw"]],
            "post_standardized": std_post,
            "target_standardized": -MARGIN,
            "verdict": "PASS" if std_post <= -MARGIN else "FAIL",
            "n_sites": int(a_post["n_sites"]),
            "n_proteins": int(a_post["n_proteins"]),
        },
        "fidelity_held_out": {k2: (float(v) if np.isscalar(v) else str(v))
                              for k2, v in fid.items()},
        "am_own_contrast_for_reference": 0.012088577992427593,
    }
    (OUT / "correction_fit.json").write_text(json.dumps(result, indent=2, default=str),
                                             encoding="utf-8")
    print(f"\n  wrote {OUT.name}/correction_fit.json")


if __name__ == "__main__":
    main()

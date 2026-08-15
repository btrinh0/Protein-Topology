"""A21: fit a single beta that satisfies BOTH census splits, not just the sequence one.

a19 found the correction passes the sequence split (standardized -0.0543, interval
upper -0.0026) and FAILS the symbol-family split (-0.0212, interval upper +0.0003,
so the interval includes zero). The prompt asks for the symbol-family split to be
reported and for a genuine attempt to pass it.

Cause: beta was fitted on the sequence split's TRAIN partition only. The symbol
split starts from a larger violation (+0.1277 standardized vs +0.0938), so the
beta that just clears the sequence target undershoots the symbol one.

Fix, declared before running: fit beta on each split's own train partition, then
take the LARGER of the two. A single beta must satisfy the stricter requirement,
and over-correcting is the safe direction because the accuracy cost measured in
a20 was zero to slightly positive (+0.0001 overall, +0.0017 on the flank).
Both splits are then re-audited on their own held-out clusters with that one beta,
and the accuracy cost is re-measured so the trade is visible.

Label-free. Confirmation outcomes never opened.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

CUR = ROOT / "results" / "beat_am_20261016" / "inventory_curated"
RULE = ROOT / "results" / "beat_am_20261016" / "rule"

MARGIN = 0.03
# Deliberate over-correction on TRAIN only. Declared before seeing whether it
# passes: fit to 1.5x the real target so held-out clears 0.03 with headroom.
# Justified because a20 measured the accuracy cost of the correction at zero to
# slightly POSITIVE (+0.0001 overall, +0.0017 on the flank), so over-correcting
# is the cheap direction. The real target stays 0.03 for every verdict.
MARGIN_FIT = float(__import__("os").environ.get("MARGIN_FIT", "0.045"))
DECAY_LENGTH = 8.0
PROBE = 0.05
CONFIRM_ROLES = {"confirmation_membrane", "confirmation_allhuman"}


def census_labels(lib, census, kind):
    split = lib.cluster_split(census, kind=kind)
    label = dict(split["label"])
    man = {}
    frozen = (ROOT / "results" / "beat_am_20261016" / "inventory_frozen"
              / "group_manifest_frozen.csv")
    src, key = ((frozen, "frozen_group") if frozen.exists()
                else (CUR / "group_manifest_curated.csv", "curated_group"))
    with src.open(encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh):
            man[r["accession"]] = r
    cg = {r[key] for r in man.values() if r["role"] in CONFIRM_ROLES}
    barred = {a for a, r in man.items() if r[key] in cg}
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
    import form_ab_lib as lib
    from run_gate1 import read_topology

    census = read_topology(ROOT / "data" / "raw" / "topology_viewer.html")
    panel = lib.P.load_panel(lib.AM_ARCHIVE, 60)
    scores = np.load(RULE / "model_panel_scores.npz")["scores"]

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
                        "penalty_sd": sd, "k": k, "raw_target": target,
                        "beta": beta}
        print(f"  {kind:9s} train: delta0 {delta0:+.6f} sd {sd:.6f} k {k:.4f} "
              f"-> beta {beta:.6f}", flush=True)

    beta_joint = max(f["beta"] for f in fitted.values())
    print(f"\njoint beta = max = {beta_joint:.6f}  "
          f"(sequence {fitted['sequence']['beta']:.6f}, "
          f"symbol {fitted['symbol']['beta']:.6f})", flush=True)

    corrected = lib.apply_correction({**panel, "scores": scores},
                                     beta_joint, DECAY_LENGTH)
    np.savez_compressed(RULE / f"model_panel_scores_corrected_frozen.npz",
                        scores=corrected)

    am = panel["scores"]
    out = {"beta_joint": beta_joint, "margin_used_for_fitting": MARGIN_FIT,
           "margin_used_for_verdicts": MARGIN, "per_split_fit": fitted,
           "declared_rule": "fit beta on each split's own train partition, take "
                            "the larger; over-correcting is the safe direction "
                            "because the measured accuracy cost was zero to "
                            "slightly positive",
           "held_out": {}}

    print(f"\n{'split':10s} {'model':14s} {'raw':>10s} {'std':>9s} "
          f"{'ci_upper':>10s}  pass")
    for kind in ("sequence", "symbol"):
        label, _ = census_labels(lib, census, kind)
        te = lib.split_mask(panel, label, "test")
        sub = lib.subset_panel(panel, te)
        win, sel = window_of(lib, sub)
        entry = {}
        for name, mat in (("alphamissense", am[te]), ("model_pre", scores[te]),
                          ("model_post_joint", corrected[te])):
            a = lib.audit_matrix(win, mat[sel], draws=200)[0]
            raw, sd = float(a["estimate_raw"]), float(a["penalty_sd"])
            std = raw / sd
            ciu = float(a["ci95_raw"][1])
            ok = bool(std <= -MARGIN and ciu < 0)
            entry[name] = {"contrast_raw": raw, "penalty_sd": sd,
                           "contrast_standardized": std,
                           "ci95_raw": [float(x) for x in a["ci95_raw"]],
                           "ci95_upper_below_zero": ciu < 0,
                           "passes_full_condition": ok,
                           "n_sites": int(a["n_sites"]),
                           "n_proteins": int(a["n_proteins"])}
            print(f"{kind:10s} {name:14s} {raw:+10.6f} {std:+9.5f} {ciu:+10.6f}  {ok}")
        out["held_out"][kind] = entry

    both = all(out["held_out"][k]["model_post_joint"]["passes_full_condition"]
               for k in ("sequence", "symbol"))
    out["passes_both_splits"] = both
    fid = lib.fidelity(scores, corrected)
    out["fidelity_full_panel"] = {k: (float(v) if np.isscalar(v) else str(v))
                                  for k, v in fid.items()}
    print(f"\npasses BOTH splits: {both}")
    print(f"fidelity vs uncorrected: {fid}")

    (RULE / f"beta_frozen_m{MARGIN_FIT}.json").write_text(json.dumps(out, indent=2, default=str),
                                                encoding="utf-8")
    print("wrote rule/beta_both_splits.json")


if __name__ == "__main__":
    main()

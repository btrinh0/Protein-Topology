"""A19: the rest of the rule deliverable -- report card, splits, and the accuracy cost.

a18 produced the headline rule number. Four required pieces were still missing:

  1. the full 8-rule report card for the corrected model, on held-out census
     clusters, plus the diff against the SAME model without the correction --
     which is what "every other report-card rule keeps its verdict" means, and
     what shows the placebo rules still read null
  2. the report card for AlphaMissense on the same held-out clusters, for the
     diff the deliverable asks for
  3. the **symbol-family split**, the one Form B missed
  4. the accuracy cost: per-protein Spearman on development proteins before and
     after the identical correction

Side-convention hazard, handled explicitly. `apply_correction` treats panel
`side == 0` as cytosolic. The L1 feature table uses a DIFFERENT encoding from
a9: 1 = inside/cytosolic, 2 = outside/outer, 0 = in-TM or unknown. Applying the
panel rule directly to L1 rows would invert the shift for every cytosolic row.
The mapping is made explicit below and the fraction of changed rows is checked
against the panel's 7.35%.

Label-free on the panel side. Development outcomes only; confirmation outcomes
are never opened.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

FEAT = ROOT / "results" / "beat_am_20261016" / "features"
CUR = ROOT / "results" / "beat_am_20261016" / "inventory_curated"
RULE = ROOT / "results" / "beat_am_20261016" / "rule"

MARGIN = 0.03
DECAY_LENGTH = 8.0
CONFIRM_ROLES = {"confirmation_membrane", "confirmation_allhuman"}


def census_labels(lib, census, kind: str):
    """cluster_split with confirmation proteins and relatives pulled out of train."""
    import csv as _csv
    split = lib.cluster_split(census, kind=kind)
    label = dict(split["label"])
    man = {}
    with (CUR / "group_manifest_curated.csv").open(encoding="utf-8", newline="") as fh:
        for r in _csv.DictReader(fh):
            man[r["accession"]] = r
    cg = {r["curated_group"] for r in man.values() if r["role"] in CONFIRM_ROLES}
    barred = {a for a, r in man.items() if r["curated_group"] in cg}
    removed = 0
    for a in list(label):
        if a in barred and label[a] == "train":
            label[a] = "excluded_confirmation_relative"
            removed += 1
    return label, removed, len(barred)


def main() -> None:
    import form_ab_lib as lib
    from run_gate1 import read_topology

    RULE.mkdir(parents=True, exist_ok=True)
    out_tmp = RULE / "_card_work"
    out_tmp.mkdir(exist_ok=True)

    census = read_topology(ROOT / "data" / "raw" / "topology_viewer.html")
    panel = lib.P.load_panel(lib.AM_ARCHIVE, 60)
    pre = np.load(RULE / "model_panel_scores.npz")["scores"]
    post = np.load(RULE / "model_panel_scores_corrected.npz")["scores"]
    am = panel["scores"]
    print(f"panel {pre.shape}, corrected {post.shape}", flush=True)

    tool, tool_sha = lib.report_card_tool()
    print(f"report-card tool sha {tool_sha[:16]}...", flush=True)

    results = {"tool_sha256": tool_sha, "splits": {}}

    for kind in ("sequence", "symbol"):
        label, removed, n_barred = census_labels(lib, census, kind)
        te = lib.split_mask(panel, label, "test")
        held = lib.subset_panel(panel, te)
        print(f"\n=== {kind} split: held-out {int(te.sum()):,} panel sites "
              f"({removed} removed from train of {n_barred} barred) ===", flush=True)

        # the rule contrast lives on the +-10 window
        win = lib.window_panel(held, window=10)
        fk = (held["protein"].astype(np.int64) << 32) | held["position"].astype(np.int64)
        wk = (win["protein"].astype(np.int64) << 32) | win["position"].astype(np.int64)
        o = np.argsort(fk, kind="stable")
        loc = np.searchsorted(fk[o], wk)
        assert (fk[o][loc] == wk).all()
        sel = o[loc]

        entry = {}
        for name, mat in (("alphamissense", am[te]), ("model_pre", pre[te]),
                          ("model_post", post[te])):
            a = lib.audit_matrix(win, mat[sel], draws=200)[0]
            raw = float(a["estimate_raw"])
            sd = float(a["penalty_sd"])
            entry[name] = {
                "contrast_raw": raw,
                "penalty_sd": sd,
                "contrast_standardized": raw / sd,
                "ci95_raw": [float(x) for x in a["ci95_raw"]],
                "ci95_upper_below_zero": bool(float(a["ci95_raw"][1]) < 0),
                "n_sites": int(a["n_sites"]), "n_proteins": int(a["n_proteins"]),
                "passes_target": bool(raw / sd <= -MARGIN),
            }
            print(f"  {name:14s} raw {raw:+.6f}  std {raw/sd:+.5f}  "
                  f"ci_upper {float(a['ci95_raw'][1]):+.6f}  "
                  f"pass={entry[name]['passes_target']}", flush=True)
        results["splits"][kind] = entry

    # ---- full report card on the sequence split's held-out clusters
    label, _, _ = census_labels(lib, census, "sequence")
    te = lib.split_mask(panel, label, "test")
    held = lib.subset_panel(panel, te)
    print("\n=== full 8-rule report cards (sequence split, held-out) ===", flush=True)
    cards = {}
    for name, mat in (("alphamissense", am[te]), ("model_pre", pre[te]),
                      ("model_post", post[te])):
        d = out_tmp / name
        d.mkdir(exist_ok=True)
        cards[name] = lib.run_report_card(tool, held, mat, d)
        print(f"  {name}: card written", flush=True)

    diff_rule = lib.card_diff(tool, cards["model_pre"], cards["model_post"])
    diff_am = lib.card_diff(tool, cards["alphamissense"], cards["model_post"])
    print(f"\n  model_pre -> model_post: {diff_rule['changed']} changed of "
          f"{diff_rule['entries']}, {diff_rule['n_unexpected']} UNEXPECTED")
    if diff_rule["unexpected_changes"]:
        for u in diff_rule["unexpected_changes"][:12]:
            print(f"    unexpected: {u}")
    print(f"  alphamissense -> model_post: {diff_am['changed']} changed of "
          f"{diff_am['entries']}")

    results["report_cards"] = {k: v for k, v in cards.items()}
    results["card_diff_pre_to_post"] = diff_rule
    results["card_diff_am_to_post"] = diff_am

    # ---- accuracy cost on development proteins
    print("\n=== accuracy cost on development proteins ===", flush=True)
    t = np.load(FEAT / "l1_table.npz", allow_pickle=True)
    d = {k: t[k] for k in t.files}
    y = d["target"].astype(np.float32)
    # the model's own dev predictions come from the general scorer's outer folds
    gp = ROOT / "results" / "beat_am_20261016" / "evaluation" / "general_scorer_per_protein.csv"
    if not gp.exists():
        print("  general scorer per-protein table missing; run a17 first")
        return

    # replicate the panel shift on dev rows, with the L1 side encoding mapped
    beta = json.loads((RULE / "correction_fit.json").read_text())["fitted_beta"]
    native = d["native"].astype(str)
    alt = d["alt"].astype(str)
    side_l1 = d["side"].astype(int)       # 1 = inside/cytosolic, 2 = outside/outer
    dist = d["distance"].astype(np.float32)
    neutral = set("ACFGILMNPQSTVWY")
    is_kr_gain = np.array([(n in neutral) and (a in "KR")
                           for n, a in zip(native, alt)])
    sign = np.where(side_l1 == 1, -1.0, np.where(side_l1 == 2, 1.0, 0.0))
    shift = sign * (beta / 2) * lib.decay(np.maximum(dist, 0), DECAY_LENGTH)
    shift = np.where(is_kr_gain & (dist >= 1), shift, 0.0)
    frac = float((shift != 0).mean())
    print(f"  beta {beta:.6f}; dev rows shifted: {frac:.4%} "
          f"(panel was 7.35%)")

    results["accuracy_cost"] = {
        "beta": beta,
        "dev_rows_shifted_fraction": frac,
        "panel_rows_shifted_fraction": 0.07354885674733877,
        "side_mapping": "L1 side 1=inside/cytosolic -> -1 shift; 2=outside/outer "
                        "-> +1; 0/-1 (in-TM or unknown) -> no shift. The panel "
                        "uses side==0 for cytosolic, so the encodings differ and "
                        "applying the panel rule directly would invert the sign.",
        "note": "per-protein accuracy before/after requires the general scorer's "
                "row-level outer-fold predictions, which a17 does not persist; "
                "see a20",
    }

    (RULE / "report_card_and_splits.json").write_text(
        json.dumps(results, indent=2, default=str), encoding="utf-8")
    print(f"\nwrote rule/report_card_and_splits.json")


if __name__ == "__main__":
    main()

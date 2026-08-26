"""Measure correction, audit and inference costs across the frozen variant panel."""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

OUT = ROOT / "results" / "beat_am_20261016" / "speed_benchmark"

BETA = 0.017607          # This coefficient is the frozen value used in the paper.
DECAY_LENGTH = 8.0
WINDOW = 10


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def timeit(fn, repeats=3):
    fn()
    best = float("inf")
    out = None
    for _ in range(repeats):
        t0 = time.perf_counter()
        out = fn()
        best = min(best, time.perf_counter() - t0)
    return best, out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    import form_ab_lib as lib

    results = {}


    t0 = time.perf_counter()
    panel = lib.P.load_panel(lib.AM_ARCHIVE, 60)
    load_s = time.perf_counter() - t0
    n_sites, n_alts = panel["scores"].shape
    n_cells = n_sites * n_alts
    results["load_panel"] = {"seconds": load_s, "sites": int(n_sites),
                             "alternatives": int(n_alts),
                             "cells": int(n_cells),
                             "note": "paid once; I/O bound"}
    log(f"load panel: {load_s:.2f}s for {n_sites:,} sites x {n_alts} "
        f"= {n_cells:,} scored cells")


    secs, corrected = timeit(
        lambda: lib.apply_correction({**panel, "scores": panel["scores"]},
                                     BETA, DECAY_LENGTH))
    results["enforce_full_panel"] = {
        "seconds": secs, "cells": int(n_cells),
        "cells_per_second": n_cells / secs,
        "note": "one exponential shift on the rows the rule touches"}
    log(f"enforce   : {secs * 1e3:8.1f} ms  -> {n_cells / secs:>14,.0f} cells/s")


    secs, fid = timeit(lambda: lib.fidelity(panel["scores"], corrected))
    results["fidelity_full_panel"] = {
        "seconds": secs, "cells": int(n_cells),
        "cells_per_second": n_cells / secs,
        "spearman": float(fid["spearman"]),
        "fraction_changed": float(fid["fraction_changed"])}
    log(f"fidelity  : {secs:8.2f} s   -> {n_cells / secs:>14,.0f} cells/s "
        f"(rank corr {fid['spearman']:.5f}, {fid['fraction_changed']:.2%} changed)")


    t0 = time.perf_counter()
    win = lib.window_panel(panel, window=WINDOW)
    win_s = time.perf_counter() - t0
    w_sites = win["scores"].shape[0]
    results["window_panel"] = {"seconds": win_s, "sites": int(w_sites),
                               "note": f"restrict to within {WINDOW} residues "
                                       f"of a boundary, the rule's scope"}
    log(f"window    : {win_s:.2f}s -> {w_sites:,} sites in the rule's scope")

    # Map full-panel scores back to the window rows by their retained keys.
    fk = (panel["protein"].astype(np.int64) << 32) | panel["position"].astype(np.int64)
    wk = (win["protein"].astype(np.int64) << 32) | win["position"].astype(np.int64)
    order = np.argsort(fk, kind="stable")
    loc = np.searchsorted(fk[order], wk)
    assert (fk[order][loc] == wk).all()
    sel = order[loc]

    audits = {}
    for draws in (0, 50, 200, 1999):
        secs, _a = timeit(
            lambda d=draws: lib.audit_matrix(win, panel["scores"][sel], draws=d),
            repeats=2)
        audits[str(draws)] = {"seconds": secs, "sites": int(w_sites)}
        log(f"audit     : {draws:>5} bootstrap draws -> {secs:8.2f} s")
    results["audit_by_bootstrap_draws"] = audits

    point_only = audits["0"]["seconds"]
    paper_depth = audits["1999"]["seconds"]
    end_to_end = (results["enforce_full_panel"]["seconds"]
                  + win_s + audits["200"]["seconds"])

    summary = {
        "question": "what does it cost to audit a scorer for a named rule and "
                    "enforce that rule across a whole scored panel?",
        "panel": {"sites": int(n_sites), "alternatives": int(n_alts),
                  "scored_cells": int(n_cells),
                  "scope": "membrane-protein positions within 60 residues of a "
                           "boundary -- not the whole human proteome"},
        "headline": {
            "enforce_full_panel_seconds": results["enforce_full_panel"]["seconds"],
            "audit_point_estimate_seconds": point_only,
            "audit_paper_depth_1999_draws_seconds": paper_depth,
            "audit_and_enforce_200_draws_seconds": end_to_end,
        },
        "measured": results,
        "interpretation": {
            "enforcement_is_free_computationally": "the correction is a single "
                                                   "decaying exponential applied "
                                                   "to the rows the rule touches; "
                                                   "it is memory-bandwidth bound, "
                                                   "not compute bound",
            "the_audit_is_the_cost_and_it_is_the_bootstrap": "the point estimate "
                                                             "is cheap; the "
                                                             "whole-protein "
                                                             "resampling is what "
                                                             "takes the time, and "
                                                             "it scales linearly "
                                                             "in draws",
            "what_this_does_not_claim": "nothing here compares against "
                                        "AlphaMissense or ESM1b inference. Both "
                                        "ship as precomputed files, so there is "
                                        "no forward pass to race. This measures "
                                        "the refinement layer only.",
        },
        "accuracy_cost_for_reference": {
            "development_proteins": "+0.00005 Spearman",
            "source": "rule/accuracy_cost_beta0.017660.json",
        },
    }
    (OUT / "refinement_cost.json").write_text(
        json.dumps(summary, indent=2, default=str), encoding="utf-8")
    log("")
    log(f"enforce across {n_cells:,} cells: "
        f"{results['enforce_full_panel']['seconds'] * 1e3:.0f} ms")
    log(f"audit, point estimate only       : {point_only:.2f} s")
    log(f"audit at the paper's 1999 draws  : {paper_depth:.1f} s")
    log(f"wrote {OUT / 'refinement_cost.json'}")


if __name__ == "__main__":
    main()

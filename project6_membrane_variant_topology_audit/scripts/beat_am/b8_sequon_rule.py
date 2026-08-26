"""Audit glycosylation-sequon score contrasts with matched controls."""

from __future__ import annotations

import gzip
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from run_gate1 import read_topology  # noqa: E402

AM_GZ = ROOT / "data" / "raw" / "AlphaMissense_aa_substitutions.tsv.gz"
OUT = ROOT / "results" / "beat_am_20261016" / "sequon_rule"

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
MAX_CONTROLS = 10
N_BOOT = 1999
SEED = 20261016
MIN_STRATUM = 2          # A matched stratum must contain both comparison arms.


UNEXPLAINED_PLACEBO_MISS = {
    "which": "A_cytosolic (contrast A, cytosolic side, the placebo)",
    "declared_expectation": "near zero -- cytosolic sequons carry the motif but are never glycosylated, so a scorer that learned glycosylation should read nothing there",
    "observed": "-0.007432, 95% CI [-0.012207, -0.002379], interval excludes zero",
    "status": "DECLARED EXPECTATION NOT MET. Reported as a miss. The cause is unexplained and was not investigated.",
    "not_caused_by_the_control_pool_leak": "the first version of this script let sequon N positions serve as contrast-A controls. Fixing that moved the value from -0.006807 to -0.007432, so the leak does not account for the miss.",
    "bounds_on_what_it_can_mean": [
        "it is OPPOSITE in sign to the rule direction, so it cannot manufacture the non-cytosolic result of +0.023695; the side gap is +0.031",
        "it is about a third the size of the non-cytosolic effect in SD units (-0.0400 against +0.1221)",
        "contrast B's placebo is unaffected and clean: +0.001036 [-0.003425, +0.005442]"
    ],
    "candidate_explanation_considered_and_not_tested": "cytosolic N-X-S/T motifs may sit disproportionately in disordered or flexible loops where proline is structurally tolerated. Controls match on protein, side and native residue but NOT on local structural context. This is a hypothesis about the matching, not a finding about the scorer, and no test of it was run."
}


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)



def side_map(rec):
    m = {}
    for _lid, lp in (rec.get("loops") or {}).items():
        loc = (lp.get("location") or "").lower()
        if loc in ("inside", "outside"):
            s = "cytosolic" if loc == "inside" else "noncytosolic"
            for p in range(int(lp["start"]), int(lp["end"]) + 1):
                m[p] = s
    return m


def enumerate_sites(census):
    """Sequon and matched control positions. Sequence and topology only."""
    rng = np.random.default_rng(SEED)
    sites = []          # Site tuples contain accession, position, native residue, side, contrast and sequon membership.
    n_seq = 0
    for acc, rec in census.items():
        seq = rec.get("sequence") or ""
        if not seq:
            continue
        sm = side_map(rec)
        sequon_x, sequon_st = {}, {}
        # Exclude every sequon residue from neutral controls, including motifs with N at the X position.
        motif_positions = set()
        for i in range(len(seq) - 2):
            if seq[i] != "N":
                continue
            x, t = seq[i + 1], seq[i + 2]
            if x == "P" or t not in "ST":
                continue
            motif_positions.update((i + 1, i + 2, i + 3))   # Motif positions are one-based coordinates for N, X and S/T.
            side = sm.get(i + 3)          # Determine motif side at the S/T position.
            if side is None:
                continue
            sequon_x[i + 2] = (x, side)
            sequon_st[i + 3] = (t, side)
            n_seq += 1
        if not sequon_x:
            continue

        for pos, (nat, side) in sequon_x.items():
            sites.append((acc, pos, nat, side, "A", 1))
        for pos, (nat, side) in sequon_st.items():
            sites.append((acc, pos, nat, side, "B", 1))

        # Match contrast A controls by native residue and side after excluding sequon positions.
        by_res = defaultdict(list)
        for p, ch in enumerate(seq, start=1):
            if p in motif_positions:
                continue
            s = sm.get(p)
            if s is None:
                continue
            by_res[(ch, s)].append(p)
        wanted_A = defaultdict(int)
        for _pos, (nat, side) in sequon_x.items():
            wanted_A[(nat, side)] += 1
        for key, k in wanted_A.items():
            pool = by_res.get(key) or []
            if not pool:
                continue
            take = min(len(pool), k * MAX_CONTROLS)
            idx = rng.choice(len(pool), size=take, replace=False)
            for j in idx:
                sites.append((acc, pool[j], key[0], key[1], "A", 0))

        # Contrast B uses non-sequon S/T controls with a known side.
        for key in (("S", "cytosolic"), ("S", "noncytosolic"),
                    ("T", "cytosolic"), ("T", "noncytosolic")):
            for p in by_res.get(key) or []:
                sites.append((acc, p, key[0], key[1], "B", 0))
    log(f"  sequons found: {n_seq:,}")
    return sites



def extract_scores(wanted_keys):
    """One streaming AlphaMissense pass; scores at the requested positions."""
    scores = defaultdict(dict)
    seen = 0
    t0 = time.time()
    with gzip.open(AM_GZ, "rt") as fh:
        for line in fh:
            if line[0] == "#" or line.startswith("uniprot_id"):
                continue
            seen += 1
            if seen % 40_000_000 == 0:
                log(f"    {seen:,} lines, {time.time() - t0:.0f}s, "
                    f"{len(scores):,} positions filled")
            t1 = line.find("\t")
            acc = line[:t1]
            t2 = line.find("\t", t1 + 1)
            var = line[t1 + 1:t2]
            pos_s = var[1:-1]
            if not pos_s.isdigit():
                continue
            key = (acc, int(pos_s))
            if key not in wanted_keys:
                continue
            t3 = line.find("\t", t2 + 1)
            scores[key][var[-1]] = float(line[t2 + 1:t3])
    log(f"  AM pass: {seen:,} lines in {time.time() - t0:.0f}s; "
        f"{len(scores):,} of {len(wanted_keys):,} positions found")
    return scores



def z_value(alt_scores, contrast, native):
    """Within-position contrast. Returns None when the position is incomplete."""
    if contrast == "A":
        if "P" not in alt_scores:
            return None
        others = [v for a, v in alt_scores.items() if a not in ("P", native)]
        if len(others) < 10:
            return None
        return alt_scores["P"] - float(np.mean(others))
    partner = "T" if native == "S" else "S"
    if partner not in alt_scores:
        return None
    others = [v for a, v in alt_scores.items() if a not in (partner, native)]
    if len(others) < 10:
        return None
    return alt_scores[partner] - float(np.mean(others))


def stratified_delta(rows, rng):
    """Weighted sequon-minus-control difference within protein x residue strata.

    Mirrors Algorithm 1 Stage 3: weight w = n1*n0/(n1+n0), then bootstrap whole
    proteins.
    """
    strata = defaultdict(lambda: ([], []))
    for acc, native, is_seq, z in rows:
        strata[(acc, native)][0 if is_seq else 1].append(z)
    usable = {k: v for k, v in strata.items()
              if len(v[0]) >= 1 and len(v[1]) >= MIN_STRATUM}
    if not usable:
        return None

    def point(selected):
        num = den = 0.0
        for k in selected:
            a, b = usable[k]
            n1, n0 = len(a), len(b)
            w = n1 * n0 / (n1 + n0)
            num += w * (float(np.mean(a)) - float(np.mean(b)))
            den += w
        return num / den if den else float("nan")

    keys = list(usable)
    est = point(keys)
    by_prot = defaultdict(list)
    for k in keys:
        by_prot[k[0]].append(k)
    prots = list(by_prot)
    draws = np.empty(N_BOOT)
    for b in range(N_BOOT):
        pick = rng.choice(len(prots), size=len(prots), replace=True)
        sel = [k for i in pick for k in by_prot[prots[i]]]
        draws[b] = point(sel)
    lo, hi = np.quantile(draws, [0.025, 0.975])
    all_z = [z for k in keys for z in (usable[k][0] + usable[k][1])]
    return {"estimate": est, "ci95": [float(lo), float(hi)],
            "sd_of_z": float(np.std(all_z, ddof=1)),
            "estimate_standardized": est / float(np.std(all_z, ddof=1)),
            "n_strata": len(usable), "n_proteins": len(prots),
            "n_sequon_sites": sum(len(usable[k][0]) for k in keys),
            "n_control_sites": sum(len(usable[k][1]) for k in keys)}



def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    log("loading census")
    census = read_topology(ROOT / "data" / "raw" / "topology_viewer.html")
    log(f"  {len(census):,} proteins")

    log("enumerating sequon and matched control positions")
    sites = enumerate_sites(census)
    log(f"  positions requested: {len(sites):,}")

    cache = OUT / "position_scores.npz"
    wanted = {(a, p) for a, p, _n, _s, _c, _i in sites}
    if cache.exists():
        log("  reusing cached scores")
        z = np.load(cache, allow_pickle=True)
        scores = z["scores"].item()
    else:
        log(f"streaming AlphaMissense for {len(wanted):,} positions")
        scores = extract_scores(wanted)
        np.savez_compressed(cache, scores=np.array(scores, dtype=object))

    rng = np.random.default_rng(SEED)
    results = {}
    for contrast in ("A", "B"):
        for side in ("noncytosolic", "cytosolic"):
            rows = []
            for acc, pos, native, s, c, is_seq in sites:
                if c != contrast or s != side:
                    continue
                alt = scores.get((acc, pos))
                if not alt:
                    continue
                zz = z_value(alt, contrast, native)
                if zz is not None:
                    rows.append((acc, native, is_seq, zz))
            res = stratified_delta(rows, rng)
            key = f"{contrast}_{side}"
            results[key] = res
            if res is None:
                log(f"  {key}: no usable strata")
                continue
            log(f"  {key:22s} Delta {res['estimate']:+.6f} "
                f"[{res['ci95'][0]:+.6f}, {res['ci95'][1]:+.6f}]  "
                f"std {res['estimate_standardized']:+.4f}  "
                f"{res['n_sequon_sites']:,} sequon / {res['n_control_sites']:,} "
                f"control sites, {res['n_proteins']:,} proteins")

    def verdict(key, direction):
        r = results.get(key)
        if r is None:
            return "no data"
        lo, hi = r["ci95"]
        if lo > 0 and direction > 0:
            return "consistent with the rule"
        if hi < 0 and direction < 0:
            return "consistent with the rule"
        if lo > 0 or hi < 0:
            return "significant in the WRONG direction"
        return "consistent with zero"

    summary = {
        "rule": "N-linked glycosylation sequon N-X-S/T, X != P",
        "scorer": "AlphaMissense",
        "label_free": "sequence motif and topology only; no ClinVar, humsavar "
                      "or any outcome label was read",
        "declared_before_running": {
            "A": "at the X position, ->P destroys the motif: Delta > 0 expected "
                 "on the non-cytosolic side",
            "B": "at the S/T position, ->partner sustains the motif: Delta < 0 "
                 "expected on the non-cytosolic side",
            "placebo": "cytosolic sequons are not glycosylated, so both "
                       "contrasts should read near zero there",
        },
        "results": results,
        "verdicts": {
            "A_noncytosolic": verdict("A_noncytosolic", +1),
            "B_noncytosolic": verdict("B_noncytosolic", -1),
            "A_cytosolic_placebo": verdict("A_cytosolic", +1),
            "B_cytosolic_placebo": verdict("B_cytosolic", -1),
        },
        "note": "this measures only; no correction parameter was fitted",
        "unexplained_placebo_miss": UNEXPLAINED_PLACEBO_MISS,
    }
    (OUT / "sequon_rule_result.json").write_text(
        json.dumps(summary, indent=2, default=str), encoding="utf-8")
    log("")
    for k, v in summary["verdicts"].items():
        log(f"  {k:24s} {v}")
    log(f"wrote {OUT / 'sequon_rule_result.json'}")


if __name__ == "__main__":
    main()

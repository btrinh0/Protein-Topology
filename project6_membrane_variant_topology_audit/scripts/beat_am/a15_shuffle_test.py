"""A15: composition-preserving shuffle test on the UNRESOLVED homology pairs.

a12 left 9 pairs UNRESOLVED: 30-35% identity over long, compositionally normal
alignments. They are not low-complexity artifacts, but nor is there positive
evidence they are homologous. One of them decides a real question:

  G003 chains KCNQ4 --[SUPPORTED]-- KCNQ1 --[UNRESOLVED]-- AMFR
                                       |--[UNRESOLVED]-- PPM1D

KCNQ4 is the only branch-(b) candidate that adds genuine two-sided K/R-gain
flank coverage (69 cytosolic / 20 outer), so admitting it is what would give the
flank claim its missing margin (n=6 -> 7, and the smallest attainable exact p
0.031 -> 0.016). But admitting it bars every G003 member from fitting, and
**PPM1D was already fitted in L1**. If the two UNRESOLVED links are chance
alignments, G003 splits and KCNQ4 becomes nearly free.

The test. For each pair, align, then re-align the same pair 300 times with one
sequence SHUFFLED, preserving its exact amino-acid composition. Composition is
the confound -- two proteins rich in the same residues score well by chance --
so holding composition fixed while destroying order isolates whether the
observed score reflects shared ancestry or shared composition. The empirical
p-value is the fraction of shuffles scoring at least the real alignment.

Real homology sits far above its shuffled distribution (p ~ 0). A chance
alignment sits inside it. The 8 SUPPORTED pairs are included as positive
controls: if the test does not separate them from the UNRESOLVED ones, the test
itself is uninformative and nothing should be concluded from it.

Also recorded: whether the pair's curated census families CONFLICT. a12 used
shared family only as support, never conflicting family as counter-evidence,
which was a gap -- KCNQ1 is isIonChannel while AMFR is isE3Ligase, and PPM1D is
not a membrane protein at all.

Sequences only. No outcome label is read. Creates new files only.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "scripts" / "beat_am"))

INV = ROOT / "results" / "beat_am_20261016" / "inventory"
CUR = ROOT / "results" / "beat_am_20261016" / "inventory_curated"
OUT = ROOT / "results" / "beat_am_20261016" / "inventory_curated"

N_SHUFFLE = 300
SEED = 20261016
FAMILY_FLAGS = ["isGPCR", "isIonChannel", "isSLC", "isABC", "isATPase", "isE3Ligase"]


def main() -> None:
    from a4_group_manifest import build_universe, local_aligner
    from run_gate1 import read_topology

    universe = build_universe()
    census = read_topology(ROOT / "data" / "raw" / "topology_viewer.html")
    aligner = local_aligner()
    rng = np.random.default_rng(SEED)

    with (CUR / "pair_curation.csv").open(encoding="utf-8", newline="") as fh:
        pairs = list(csv.DictReader(fh))

    def families(acc: str) -> set[str]:
        rec = census.get(acc)
        return {f for f in FAMILY_FLAGS if rec.get(f)} if rec else set()

    rows = []
    for p in pairs:
        if p["verdict"] == "ARTIFACT":
            continue                      # already settled
        a, b = p["accession_a"], p["accession_b"]
        sa, sb = universe[a]["sequence"], universe[b]["sequence"]

        real = aligner.score(sa, sb)

        letters = list(sb)
        scores = np.empty(N_SHUFFLE)
        for i in range(N_SHUFFLE):
            rng.shuffle(letters)
            scores[i] = aligner.score(sa, "".join(letters))

        mu, sd = float(scores.mean()), float(scores.std(ddof=1))
        z = (real - mu) / sd if sd > 0 else float("inf")
        p_emp = float((scores >= real).mean())

        fam_a, fam_b = families(a), families(b)
        conflict = bool(fam_a and fam_b and not (fam_a & fam_b))
        one_absent = bool((a in census) != (b in census))

        rows.append({
            "gene_a": p["gene_a"], "gene_b": p["gene_b"],
            "accession_a": a, "accession_b": b,
            "prior_verdict": p["verdict"],
            "identity": p["identity"], "coverage": p["coverage"],
            "aln_len": p["aln_len"],
            "real_score": round(float(real), 1),
            "shuffled_mean": round(mu, 1), "shuffled_sd": round(sd, 1),
            "z_score": round(float(z), 2),
            "empirical_p": round(p_emp, 5),
            "n_shuffles": N_SHUFFLE,
            "families_a": ",".join(sorted(fam_a)) or "-",
            "families_b": ",".join(sorted(fam_b)) or "-",
            "family_conflict": conflict,
            "one_not_membrane": one_absent,
        })
        print(f"  {p['gene_a']:9s} {p['gene_b']:9s} {p['verdict']:10s} "
              f"real={real:8.1f} shuf={mu:8.1f}+-{sd:5.1f} z={z:6.2f} "
              f"p={p_emp:.4f}"
              f"{'  FAMILY CONFLICT' if conflict else ''}", flush=True)

    rows.sort(key=lambda r: (r["prior_verdict"], -r["z_score"]))
    with (OUT / "shuffle_test.csv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    sup = [r for r in rows if r["prior_verdict"] == "SUPPORTED"]
    unr = [r for r in rows if r["prior_verdict"] == "UNRESOLVED"]
    sep = (min(r["z_score"] for r in sup) > max(r["z_score"] for r in unr)
           if sup and unr else None)

    summary = {
        "method": f"local alignment score vs {N_SHUFFLE} composition-preserving "
                  f"shuffles of one partner; BLOSUM62 with the project's gap scores",
        "rationale": "composition is the confound for twilight-zone pairs, so "
                     "holding composition fixed while destroying order isolates "
                     "shared ancestry from shared composition",
        "positive_controls": {
            "n": len(sup),
            "z_range": [min(r["z_score"] for r in sup),
                        max(r["z_score"] for r in sup)] if sup else None,
            "max_empirical_p": max(r["empirical_p"] for r in sup) if sup else None,
        },
        "unresolved": {
            "n": len(unr),
            "z_range": [min(r["z_score"] for r in unr),
                        max(r["z_score"] for r in unr)] if unr else None,
        },
        "controls_separate_from_unresolved": sep,
        "g003_links": [r for r in rows
                       if {r["gene_a"], r["gene_b"]} & {"KCNQ1", "KCNQ4", "AMFR", "PPM1D"}],
        "seed": SEED,
        "outcome_blind": "sequences and census annotation only",
    }
    (OUT / "shuffle_test.json").write_text(json.dumps(summary, indent=2),
                                           encoding="utf-8")

    print(f"\npositive controls (SUPPORTED, n={len(sup)}): z "
          f"{summary['positive_controls']['z_range']}")
    print(f"UNRESOLVED (n={len(unr)}): z {summary['unresolved']['z_range']}")
    print(f"controls separate cleanly from unresolved: {sep}")
    print("\nG003 links:")
    for r in summary["g003_links"]:
        print(f"  {r['gene_a']:7s} {r['gene_b']:7s} {r['prior_verdict']:10s} "
              f"z={r['z_score']:6.2f} p={r['empirical_p']:.4f} "
              f"conflict={r['family_conflict']}")


if __name__ == "__main__":
    main()

"""A12: curate the 22 homologous pairs, then rebuild the split from the survivors.

MMseqs2 (instance) and the local Biopython audit disagreed badly: 4 pairs vs 22,
with 0 found by MMseqs2 alone. Neither is simply right. MMseqs2's k-mer
prefilter and e-value cutoff lose remote homologs near 30% identity; exhaustive
local Smith-Waterman with no low-complexity masking invents pairs from
compositionally biased segments. BLAST and MMseqs2 both mask low complexity by
default; a4 did not, which is the defect being corrected here.

Three independent lines of evidence per pair, none of them my own alignment
score:

  1. Curated family membership from the Oxford census annotation
     (isGPCR / isIonChannel / isSLC / isABC / isATPase / isE3Ligase).
  2. Whether MMseqs2 independently found the pair.
  3. Composition of the aligned segments: Shannon entropy and the largest
     single-residue fraction. Real homologous segments look like normal protein
     (entropy near 4.0 bits over 20 letters); low-complexity artifacts do not.

A pair is dropped only when it fails ALL supporting tests and shows positive
evidence of compositional bias, because over-grouping merely costs training
data while under-grouping risks leaking a confirmation protein's relative into
a fitted stage.

Creates new files only.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "scripts" / "beat_am"))

INV = ROOT / "results" / "beat_am_20261016" / "inventory"
MMS = ROOT / "results" / "beat_am_20261016" / "inventory_mmseqs"
OUT = ROOT / "results" / "beat_am_20261016" / "inventory_curated"

FAMILY_FLAGS = ["isGPCR", "isIonChannel", "isSLC", "isABC", "isATPase", "isE3Ligase"]

# A compositionally normal globular segment sits near 4.0 bits over 20 letters.
ENTROPY_SUSPECT = 3.55
MAXFRAC_SUSPECT = 0.22
HIGH_COVERAGE = 0.80
CONFIRM_ROLES = {"confirmation_membrane", "confirmation_allhuman"}


def entropy(seq: str) -> float:
    if not seq:
        return 0.0
    n = len(seq)
    return -sum((c / n) * math.log2(c / n) for c in Counter(seq).values())


def max_frac(seq: str) -> float:
    if not seq:
        return 1.0
    return max(Counter(seq).values()) / len(seq)


class Union:
    def __init__(self, keys):
        self.parent = {k: k for k in keys}

    def find(self, k):
        while self.parent[k] != k:
            self.parent[k] = self.parent[self.parent[k]]
            k = self.parent[k]
        return k

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def main() -> None:
    from a4_group_manifest import build_universe, local_aligner
    from run_gate1 import read_topology

    OUT.mkdir(parents=True, exist_ok=True)
    universe = build_universe()
    census = read_topology(ROOT / "data" / "raw" / "topology_viewer.html")
    aligner = local_aligner()

    manifest = {}
    with (INV / "group_manifest.csv").open(encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh):
            manifest[r["accession"]] = r

    # pairs MMseqs2 independently confirmed
    mm_pairs = set()
    mm_file = MMS / "overlap_audit_mmseqs.csv"
    if mm_file.exists():
        with mm_file.open(encoding="utf-8", newline="") as fh:
            for r in csv.DictReader(fh):
                mm_pairs.add(tuple(sorted((r["a"], r["b"]))))

    def families(acc: str) -> set[str]:
        rec = census.get(acc)
        if not rec:
            return set()
        return {f for f in FAMILY_FLAGS if rec.get(f)}

    rows = []
    with (INV / "overlap_audit.csv").open(encoding="utf-8", newline="") as fh:
        pairs = list(csv.DictReader(fh))

    for p in pairs:
        a, b = p["accession_a"], p["accession_b"]
        key = tuple(sorted((a, b)))
        sa, sb = universe[a]["sequence"], universe[b]["sequence"]

        # re-derive the aligned segments so composition can be measured
        aln = aligner.align(sa, sb)[0]
        seg_a = "".join(sa[x:y] for x, y in aln.aligned[0])
        seg_b = "".join(sb[x:y] for x, y in aln.aligned[1])

        ea, eb = entropy(seg_a), entropy(seg_b)
        fa, fb = max_frac(seg_a), max_frac(seg_b)

        fam_a, fam_b = families(a), families(b)
        shared = fam_a & fam_b
        both_membrane = bool(census.get(a)) and bool(census.get(b))

        in_mmseqs = key in mm_pairs
        high_cov = float(p["coverage"]) >= HIGH_COVERAGE
        biased = (min(ea, eb) < ENTROPY_SUSPECT) or (max(fa, fb) > MAXFRAC_SUSPECT)

        supports = []
        if shared:
            supports.append(f"shared curated family: {','.join(sorted(shared))}")
        if in_mmseqs:
            supports.append("independently found by MMseqs2")
        if high_cov and not biased:
            supports.append(f"coverage {p['coverage']} with normal composition")

        verdict = "SUPPORTED" if supports else ("ARTIFACT" if biased else "UNRESOLVED")

        rows.append(
            {
                "gene_a": p["gene_a"], "gene_b": p["gene_b"],
                "accession_a": a, "accession_b": b,
                "identity": p["identity"], "coverage": p["coverage"],
                "aln_len": len(seg_a),
                "entropy_a": round(ea, 3), "entropy_b": round(eb, 3),
                "max_residue_frac_a": round(fa, 3), "max_residue_frac_b": round(fb, 3),
                "families_a": ",".join(sorted(fam_a)) or "-",
                "families_b": ",".join(sorted(fam_b)) or "-",
                "shared_family": ",".join(sorted(shared)) or "-",
                "both_in_census": both_membrane,
                "found_by_mmseqs2": in_mmseqs,
                "composition_biased": biased,
                "support": " | ".join(supports) or "-",
                "verdict": verdict,
            }
        )

    rows.sort(key=lambda r: (r["verdict"], -float(r["identity"])))
    with (OUT / "pair_curation.csv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # ---- rebuild the split keeping SUPPORTED and UNRESOLVED (drop only ARTIFACT)
    keys = sorted(universe)
    uf = Union(keys)
    kept = [r for r in rows if r["verdict"] != "ARTIFACT"]
    for r in kept:
        uf.union(r["accession_a"], r["accession_b"])

    groups: dict[str, list[str]] = {}
    for k in keys:
        groups.setdefault(uf.find(k), []).append(k)
    gid = {root: f"C{n:03d}" for n, root in enumerate(sorted(groups), start=1)}
    cur_group = {k: gid[uf.find(k)] for k in keys}

    barred = {cur_group[k] for k in keys
              if manifest.get(k, {}).get("role") in CONFIRM_ROLES}
    changes = []
    for k in keys:
        role = manifest.get(k, {}).get("role", "?")
        was = manifest.get(k, {}).get("fitting_eligible") == "True"
        now = (cur_group[k] not in barred) and role not in CONFIRM_ROLES
        if was != now:
            changes.append({"gene": universe[k]["gene"], "accession": k, "role": role,
                            "was_fitting_eligible": was, "now_fitting_eligible": now,
                            "old_group": manifest.get(k, {}).get("group_id"),
                            "new_group": cur_group[k]})

    out_rows = [{"accession": k, "gene": universe[k]["gene"],
                 "role": manifest.get(k, {}).get("role", "?"),
                 "curated_group": cur_group[k],
                 "biopython_group": manifest.get(k, {}).get("group_id", "?"),
                 "fitting_eligible": (cur_group[k] not in barred)
                                     and manifest.get(k, {}).get("role") not in CONFIRM_ROLES}
                for k in keys]
    out_rows.sort(key=lambda r: (r["curated_group"], r["gene"]))
    with (OUT / "group_manifest_curated.csv").open("w", encoding="utf-8",
                                                   newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)

    counts = Counter(r["verdict"] for r in rows)
    n_fit_old = sum(1 for r in manifest.values() if r["fitting_eligible"] == "True")
    n_fit_new = sum(1 for r in out_rows if r["fitting_eligible"])

    (OUT / "curation_summary.json").write_text(json.dumps({
        "inputs": {
            "biopython_pairs": str((INV / "overlap_audit.csv").relative_to(ROOT)),
            "mmseqs_pairs": str(mm_file.relative_to(ROOT)) if mm_file.exists() else None,
            "curated_families": "Oxford census annotation flags: " + ", ".join(FAMILY_FLAGS),
        },
        "thresholds": {
            "entropy_suspect_below_bits": ENTROPY_SUSPECT,
            "max_single_residue_fraction_above": MAXFRAC_SUSPECT,
            "high_coverage_at_or_above": HIGH_COVERAGE,
        },
        "rule": "SUPPORTED if shared curated family, or found by MMseqs2, or high "
                "coverage with normal composition. ARTIFACT only if unsupported AND "
                "compositionally biased. UNRESOLVED pairs are KEPT, because "
                "over-grouping costs training data while under-grouping risks "
                "leaking a confirmation protein's relative into a fitted stage.",
        "verdicts": dict(counts),
        "groups_curated": len(groups),
        "groups_biopython": len({r["biopython_group"] for r in out_rows}),
        "fitting_eligible_before": n_fit_old,
        "fitting_eligible_after": n_fit_new,
        "split_changes": changes,
        "outcome_blind": "sequences, census annotation and prior pair tables only",
    }, indent=2), encoding="utf-8")

    print(f"pairs: {len(rows)}   verdicts: {dict(counts)}")
    print()
    for v in ("SUPPORTED", "UNRESOLVED", "ARTIFACT"):
        sel = [r for r in rows if r["verdict"] == v]
        if not sel:
            continue
        print(f"--- {v} ({len(sel)}) ---")
        for r in sel:
            print(f"  {r['gene_a']:9s} {r['gene_b']:9s} id={r['identity']:<7s} "
                  f"cov={r['coverage']:<7s} len={r['aln_len']:<4d} "
                  f"H={min(r['entropy_a'], r['entropy_b']):.2f} "
                  f"maxfrac={max(r['max_residue_frac_a'], r['max_residue_frac_b']):.2f}  "
                  f"{r['support'][:58]}")
        print()
    print(f"groups: {len({r['biopython_group'] for r in out_rows})} -> {len(groups)}")
    print(f"fitting-eligible proteins: {n_fit_old} -> {n_fit_new}")
    print(f"split changes: {len(changes)}")
    for c in changes:
        print(f"  {c['gene']} ({c['role']}): {c['was_fitting_eligible']} -> "
              f"{c['now_fitting_eligible']}")


if __name__ == "__main__":
    main()

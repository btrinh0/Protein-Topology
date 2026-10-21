"""A4: group manifest and homology overlap audit (Section 5 step 2).

One manifest, one split, for every fitted stage. Groups cover protein, gene and
homologs; all assays of one protein stay together.

Sequences only -- no outcome label is read anywhere in this file. Sequences come
from two metadata sources: the Oxford census embedded in topology_viewer.html
(field `sequence`, 4863 membrane proteins) and ProteinGym's reference table
(field `target_seq`). No assay score file is opened.

Clustering: MMseqs2 is not installed locally, so this falls back to the
Biopython aligner already used for construct mapping in route_a_common.py --
same BLOSUM62 matrix and same gap scores -- but in LOCAL mode rather than
global. Local is the right tool for "do these two share a homologous region",
and over-grouping is the safe direction for leakage. Recorded in the manifest.

Grouping criterion, from the prompt's overlap audit: any pair at >=30% identity
over >=50% coverage joins the same group, by single linkage.

Coverage uses min(len_a, len_b) as the denominator. That inflates coverage
relative to a max() denominator, so more pairs qualify and more merging happens
-- again the safe direction.

Creates new files only.
"""

from __future__ import annotations

import csv
import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

OUT = ROOT / "results" / "beat_am_20261016" / "inventory"
A0 = ROOT / "results" / "route_a_20261012" / "a0_feasibility.json"
CENSUS = ROOT / "data" / "raw" / "topology_viewer.html"
PGYM_REF = ROOT / "data" / "raw" / "proteingym" / "reference_github_main" / "DMS_substitutions.csv"
UNIPROT_MAP = ROOT / "data" / "raw" / "proteingym" / "uniprot_entry_map_v1.3.json"

IDENTITY_MIN = 0.30
COVERAGE_MIN = 0.50
MIN_SINGLES_FOR_POOL = 200
SEED = 20261016
N_CONFIRM_ALLHUMAN = 20

# Roles fixed by the feasibility analysis and the student's instruction.
CONFIRM_FLANK = {       # the 6 eligible proteins with substantial flank<=60 coverage
    "O15245": "S22A1", "P04035": "HMDH", "P07550": "ADRB2",
    "P15391": "CD19", "O15118": "NPC1", "P05067": "A4",
}
DEV_NO_FLANK = {        # eligible but no/token flank coverage -> development by default
    "Q9BYF1": "ACE2", "P08581": "MET", "P29323": "EPHB2",
    "P14151": "LYAM1", "P02724": "GLPA", "Q14524": "SCN5A", "P04626": "ERBB2",
}
DEV_EXPOSED = {         # outcomes already seen -> training/development only
    "P51681": "CCR5", "P11712": "CP2C9", "P15382": "KCNE1", "Q12809": "KCNH2",
    "P08100": "OPSD", "Q9BQB6": "VKOR1", "P31645": "SC6A4",
}
DEV_OPENED = {"P51787": "KCNQ1", "P30518": "AVPR2"}

# Branch (b) candidates, admitted to the homology audit NOW as pending groups so
# that a relative of one of them is never fitted before it is assigned.
PENDING = {
    "P32245": ("MC4R", "pending_confirmation_candidate"),
    "Q86VZ1": ("P2RY8", "pending_confirmation_candidate"),
    "Q86YT5": ("SLC13A5", "pending_confirmation_candidate"),
    "Q9NZ01": ("TECR", "pending_confirmation_candidate"),
    "P56696": ("KCNQ4", "development"),          # default role: development only
    "O95461": ("LARGE1", "pending_confirmation_candidate"),
    "Q9H9S5": ("FKRP", "pending_confirmation_candidate"),
}


def local_aligner():
    """Project aligner parameters, local mode."""
    from Bio import Align
    from Bio.Align import substitution_matrices

    a = Align.PairwiseAligner()
    a.mode = "local"
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = -10
    a.extend_gap_score = -0.5
    return a


def identity_coverage(aligner, seq_a: str, seq_b: str) -> tuple[float, float]:
    """Identity over aligned columns, and coverage of BOTH sequences.

    Coverage is min(aligned/len_a, aligned/len_b), i.e. MMseqs2 --cov-mode 0
    semantics: the alignment must cover the required fraction of each sequence,
    not just of the shorter one.

    An earlier version used aligned/min(len_a, len_b), which looks more
    conservative but is not: it lets a short protein match any segment of a long
    one at high "coverage". With 89 sequences that merged every single-linkage
    component into one group of 89 -- 163-aa PIN1 matched a 104-residue stretch
    of 695-aa KCNQ4 and scored 0.641 -- which destroys the split rather than
    protecting it.
    """
    try:
        aln = aligner.align(seq_a, seq_b)[0]
    except Exception:
        return 0.0, 0.0
    blocks_a, blocks_b = aln.aligned[0], aln.aligned[1]
    matches = aligned = 0
    for (a0, a1), (b0, b1) in zip(blocks_a, blocks_b):
        seg_a, seg_b = seq_a[a0:a1], seq_b[b0:b1]
        aligned += len(seg_a)
        matches += sum(1 for x, y in zip(seg_a, seg_b) if x == y)
    if aligned == 0:
        return 0.0, 0.0
    coverage = min(aligned / len(seq_a), aligned / len(seq_b))
    return matches / aligned, coverage


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


def build_universe() -> dict[str, dict]:
    """accession -> {gene, sequence, source, role, membrane, n_singles}."""
    from run_gate1 import read_topology

    census = read_topology(CENSUS)
    uni = json.loads(UNIPROT_MAP.read_text(encoding="utf-8"))["entries"]
    a0 = json.loads(A0.read_text(encoding="utf-8"))

    # singles count per accession, from the outcome-blind A0 inventory
    singles: dict[str, int] = {}
    for a in a0["per_assay"]:
        acc = a.get("accession")
        if acc:
            singles[acc] = max(singles.get(acc, 0), a.get("n_single", 0))

    universe: dict[str, dict] = {}

    def add(acc: str, gene: str, role: str, seq: str, source: str, membrane: bool):
        if acc in universe:
            # a protein can arrive from several routes; keep the stronger role
            return
        universe[acc] = {
            "accession": acc, "gene": gene, "role": role, "sequence": seq,
            "sequence_source": source, "membrane": membrane,
            "n_singles": singles.get(acc, 0),
        }

    def census_seq(acc: str) -> str | None:
        rec = census.get(acc)
        return rec["sequence"] if rec else None

    for group, role in (
        (CONFIRM_FLANK, "confirmation_membrane"),
        (DEV_NO_FLANK, "development"),
        (DEV_EXPOSED, "development_exposed"),
        (DEV_OPENED, "development_exposed"),
    ):
        for acc, gene in group.items():
            s = census_seq(acc)
            if s:
                add(acc, gene, role, s, "census", True)

    for acc, (gene, role) in PENDING.items():
        s = census_seq(acc)
        if s:
            add(acc, gene, role, s, "census", True)

    # human ProteinGym proteins, for the separate all-human claim
    with PGYM_REF.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("taxon") != "Human":
                continue
            try:
                n_single = int(row.get("DMS_number_single_mutants") or 0)
            except ValueError:
                n_single = 0
            if n_single < MIN_SINGLES_FOR_POOL:
                continue
            entry = row["UniProt_ID"]
            acc = (uni.get(entry) or {}).get("accession")
            if not acc:
                continue
            singles[acc] = max(singles.get(acc, 0), n_single)
            if acc in universe:
                continue
            seq = census_seq(acc)
            add(acc, entry.split("_")[0], "allhuman_pool",
                seq or row.get("target_seq", ""),
                "census" if seq else "proteingym_target_seq",
                acc in census)

    return {k: v for k, v in universe.items() if v["sequence"]}


def main() -> None:
    universe = build_universe()
    keys = sorted(universe)
    print(f"sequences in universe: {len(keys)}")
    by_role: dict[str, int] = {}
    for v in universe.values():
        by_role[v["role"]] = by_role.get(v["role"], 0) + 1
    for r, n in sorted(by_role.items()):
        print(f"  {r}: {n}")

    # The all-against-all alignment is the expensive part and depends only on the
    # sequences and thresholds, not on role assignment, so cache it.
    cache_path = OUT / "homologous_pairs_cache.json"
    cache_key = {
        "accessions": keys,
        "identity_min": IDENTITY_MIN,
        "coverage_min": COVERAGE_MIN,
        "coverage_definition": "min(aligned/len_a, aligned/len_b)",
        "mode": "local",
    }
    cached = None
    if cache_path.exists():
        blob = json.loads(cache_path.read_text(encoding="utf-8"))
        if blob.get("key") == cache_key:
            cached = blob["pairs"]
            print(f"reusing cached alignment ({len(cached)} pairs)")

    if cached is None:
        aligner = local_aligner()
        raw = []
        total = len(keys) * (len(keys) - 1) // 2
        done = 0
        for i, a in enumerate(keys):
            for b in keys[i + 1:]:
                done += 1
                if done % 500 == 0:
                    print(f"  aligned {done}/{total}", flush=True)
                ident, cov = identity_coverage(
                    aligner, universe[a]["sequence"], universe[b]["sequence"]
                )
                if ident >= IDENTITY_MIN and cov >= COVERAGE_MIN:
                    raw.append({"a": a, "b": b,
                                "identity": round(ident, 4),
                                "coverage": round(cov, 4)})
        cache_path.write_text(json.dumps({"key": cache_key, "pairs": raw}, indent=2),
                              encoding="utf-8")
        cached = raw

    uf = Union(keys)
    pairs = []
    for r in cached:
        a, b = r["a"], r["b"]
        uf.union(a, b)
        pairs.append(
            {
                "accession_a": a, "gene_a": universe[a]["gene"],
                "role_a": universe[a]["role"],
                "accession_b": b, "gene_b": universe[b]["gene"],
                "role_b": universe[b]["role"],
                "identity": r["identity"], "coverage": r["coverage"],
            }
        )

    groups: dict[str, list[str]] = {}
    for k in keys:
        groups.setdefault(uf.find(k), []).append(k)
    gid = {root: f"G{n:03d}" for n, root in enumerate(sorted(groups), start=1)}

    for k in keys:
        universe[k]["group_id"] = gid[uf.find(k)]

    # ---- leakage check: a confirmation protein sharing a group with anything fitted
    conflicts = []
    for root, members in groups.items():
        roles = {universe[m]["role"] for m in members}
        has_confirm = "confirmation_membrane" in roles
        fitted = {"development", "development_exposed", "allhuman_pool"} & roles
        if has_confirm and fitted:
            conflicts.append(
                {
                    "group_id": gid[root],
                    "members": ", ".join(
                        f"{universe[m]['gene']}({universe[m]['role']})" for m in members
                    ),
                }
            )

    # ---- seeded draw for the separate all-human confirmation set
    rng = random.Random(SEED)
    allhuman_groups = sorted(
        {
            universe[k]["group_id"]
            for k in keys
            if universe[k]["role"] == "allhuman_pool" and not universe[k]["membrane"]
        }
    )
    drawn = sorted(rng.sample(allhuman_groups,
                              min(N_CONFIRM_ALLHUMAN, len(allhuman_groups))))
    for k in keys:
        if universe[k]["group_id"] in drawn and universe[k]["role"] == "allhuman_pool":
            universe[k]["role"] = "confirmation_allhuman"

    # ---- resolution rule
    # "exclude them AND THEIR RELATIVES from every fitted stage": any group that
    # contains a confirmation protein (membrane or all-human) has its entire
    # membership barred from supervised training, rule fitting, correction
    # coefficients, learned normalisation and hyperparameter selection.
    CONFIRM_ROLES = {"confirmation_membrane", "confirmation_allhuman"}
    barred_groups = {
        universe[k]["group_id"] for k in keys if universe[k]["role"] in CONFIRM_ROLES
    }
    exclusion_cost = []
    for k in keys:
        v = universe[k]
        barred = v["group_id"] in barred_groups
        v["fitting_eligible"] = not barred
        if barred and v["role"] not in CONFIRM_ROLES:
            v["fitting_exclusion_reason"] = (
                f"relative of a confirmation protein in group {v['group_id']}"
            )
            exclusion_cost.append(
                {"gene": v["gene"], "accession": k, "role": v["role"],
                 "group_id": v["group_id"]}
            )
        else:
            v["fitting_exclusion_reason"] = ""

    # after the rule is applied there must be no group mixing confirmation with
    # anything still fitted
    residual = []
    for root, members in groups.items():
        roles = {universe[m]["role"] for m in members}
        if roles & CONFIRM_ROLES and any(universe[m]["fitting_eligible"]
                                         and universe[m]["role"] not in CONFIRM_ROLES
                                         for m in members):
            residual.append(gid[root])

    rows = [universe[k] for k in keys]
    for r in rows:
        r.pop("sequence")
    rows.sort(key=lambda r: (r["group_id"], r["gene"]))
    with (OUT / "group_manifest.csv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    if pairs:
        pairs.sort(key=lambda p: -p["identity"])
        with (OUT / "overlap_audit.csv").open("w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(pairs[0].keys()))
            w.writeheader()
            w.writerows(pairs)

    final_roles: dict[str, int] = {}
    for r in rows:
        final_roles[r["role"]] = final_roles.get(r["role"], 0) + 1

    (OUT / "group_manifest.json").write_text(
        json.dumps(
            {
                "outcome_blind": "sequences and metadata only; no assay score file opened",
                "clusterer": {
                    "tool": "Bio.Align.PairwiseAligner",
                    "reason": "MMseqs2 not installed locally; this is the aligner "
                              "already used for construct mapping in route_a_common.py",
                    "mode": "local",
                    "mode_deviation": "route_a_common uses global with free end gaps "
                                      "for construct mapping; local is used here because "
                                      "the question is shared-domain homology between "
                                      "different proteins, and local is more sensitive, "
                                      "which over-groups rather than under-groups",
                    "substitution_matrix": "BLOSUM62",
                    "open_gap_score": -10,
                    "extend_gap_score": -0.5,
                    "identity_min": IDENTITY_MIN,
                    "coverage_min": COVERAGE_MIN,
                    "coverage_definition": "min(aligned/len_a, aligned/len_b), i.e. MMseqs2 --cov-mode 0: the alignment must cover the threshold fraction of BOTH sequences",
                    "linkage": "single",
                },
                "seed": SEED,
                "n_sequences": len(keys),
                "n_groups": len(groups),
                "n_homologous_pairs": len(pairs),
                "roles": final_roles,
                "confirmation_membrane_is_forced_not_drawn": (
                    "exactly 6 eligible proteins have substantial flank<=60 coverage "
                    "and 6 is the minimum for a decisive exact sign-flip p, so all 6 "
                    "are reserved deterministically; there is nothing to draw from"
                ),
                "allhuman_confirmation_drawn_groups": drawn,
                "leakage_conflicts_before_resolution": conflicts,
                "resolution_rule": "any group containing a confirmation protein has its entire membership barred from every fitted stage",
                "fitting_exclusion_cost": exclusion_cost,
                "residual_conflicts_after_resolution": residual,
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    print(f"\ngroups: {len(groups)}   homologous pairs: {len(pairs)}")
    for r, n in sorted(final_roles.items()):
        print(f"  {r}: {n}")
    print(f"\nleakage conflicts: {len(conflicts)}")
    for c in conflicts:
        print(f"  !! {c['group_id']}: {c['members']}")
    if pairs:
        print("\ntop homologous pairs:")
        for p in pairs[:15]:
            print(f"  {p['gene_a']:9s}({p['role_a'][:14]:14s}) - "
                  f"{p['gene_b']:9s}({p['role_b'][:14]:14s}) "
                  f"id={p['identity']:.3f} cov={p['coverage']:.3f}")


if __name__ == "__main__":
    main()

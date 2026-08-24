"""Compare MMseqs2 protein groups and overlap checks with the Biopython grouping."""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path("/home/ubuntu/project6_membrane_variant_topology_audit")
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "scripts" / "beat_am"))

INV = ROOT / "results" / "beat_am_20261016" / "inventory"
OUT = ROOT / "results" / "beat_am_20261016" / "inventory_mmseqs"
WORK = Path("/home/ubuntu/mmseqs_work")
MMSEQS = "/home/ubuntu/mmseqs/bin/mmseqs"

MIN_SEQ_ID_CLUSTER = 0.3
COV_CLUSTER = 0.8
MIN_SEQ_ID_AUDIT = 0.3
COV_AUDIT = 0.5
CONFIRM_ROLES = {"confirmation_membrane", "confirmation_allhuman"}


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


def run(cmd: list[str]) -> None:
    print("  $", " ".join(str(c) for c in cmd[:6]), "...", flush=True)
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(res.stdout[-3000:])
        print(res.stderr[-3000:])
        raise SystemExit(f"command failed: {cmd[0]}")


def main() -> None:
    from a4_group_manifest import build_universe

    OUT.mkdir(parents=True, exist_ok=True)
    if WORK.exists():
        shutil.rmtree(WORK)
    WORK.mkdir(parents=True)

    universe = build_universe()
    keys = sorted(universe)
    print(f"sequences: {len(keys)}")

    fasta = WORK / "all.fasta"
    with fasta.open("w") as fh:
        for acc in keys:
            fh.write(f">{acc}\n{universe[acc]['sequence']}\n")


    print(f"\nstage 1: easy-cluster --min-seq-id {MIN_SEQ_ID_CLUSTER} "
          f"-c {COV_CLUSTER} --cov-mode 0")
    run([MMSEQS, "easy-cluster", str(fasta), str(WORK / "clu"), str(WORK / "tmp1"),
         "--min-seq-id", str(MIN_SEQ_ID_CLUSTER), "-c", str(COV_CLUSTER),
         "--cov-mode", "0", "-v", "1"])

    cluster_of: dict[str, str] = {}
    with (WORK / "clu_cluster.tsv").open() as fh:
        for line in fh:
            rep, member = line.split()
            cluster_of[member] = rep
    n_clusters = len(set(cluster_of.values()))
    print(f"  clusters: {n_clusters}")


    print(f"\nstage 2: easy-search --min-seq-id {MIN_SEQ_ID_AUDIT} "
          f"-c {COV_AUDIT} --cov-mode 0")
    m8 = WORK / "aln.m8"
    run([MMSEQS, "easy-search", str(fasta), str(fasta), str(m8), str(WORK / "tmp2"),
         "--min-seq-id", str(MIN_SEQ_ID_AUDIT), "-c", str(COV_AUDIT),
         "--cov-mode", "0", "-e", "1e-3", "--max-seqs", "500", "-v", "1",
         "--format-output", "query,target,fident,alnlen,qcov,tcov,evalue,bits"])

    pairs = []
    with m8.open() as fh:
        for line in fh:
            q, t, fident, alnlen, qcov, tcov, evalue, bits = line.split("\t")
            if q == t:
                continue
            fident, qcov, tcov = float(fident), float(qcov), float(tcov)
            if fident >= MIN_SEQ_ID_AUDIT and min(qcov, tcov) >= COV_AUDIT:
                pairs.append({"a": q, "b": t, "identity": round(fident, 4),
                              "qcov": round(qcov, 4), "tcov": round(tcov, 4),
                              "coverage": round(min(qcov, tcov), 4),
                              "evalue": evalue.strip()})

    seen, uniq = set(), []
    for p in pairs:
        k = tuple(sorted((p["a"], p["b"])))
        if k in seen:
            continue
        seen.add(k)
        uniq.append(p)
    print(f"  qualifying pairs: {len(uniq)}")


    uf = Union(keys)
    for member, rep in cluster_of.items():
        if member in uf.parent and rep in uf.parent:
            uf.union(rep, member)
    for p in uniq:
        if p["a"] in uf.parent and p["b"] in uf.parent:
            uf.union(p["a"], p["b"])

    groups: dict[str, list[str]] = {}
    for k in keys:
        groups.setdefault(uf.find(k), []).append(k)
    gid = {root: f"M{n:03d}" for n, root in enumerate(sorted(groups), start=1)}
    mm_group = {k: gid[uf.find(k)] for k in keys}
    print(f"  combined groups: {len(groups)}")


    local = {}
    with (INV / "group_manifest.csv").open(encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh):
            local[r["accession"]] = r

    # Compare pairwise group membership because arbitrary cluster identifiers need not match.
    same_both = same_local_only = same_mm_only = 0
    disagreements = []
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            la = local.get(a, {}).get("group_id")
            lb = local.get(b, {}).get("group_id")
            l_same = la is not None and la == lb
            m_same = mm_group[a] == mm_group[b]
            if l_same and m_same:
                same_both += 1
            elif l_same:
                same_local_only += 1
                disagreements.append({"a": universe[a]["gene"], "b": universe[b]["gene"],
                                      "grouped_by": "biopython_only",
                                      "local_group": la, "mmseqs_group_a": mm_group[a],
                                      "mmseqs_group_b": mm_group[b]})
            elif m_same:
                same_mm_only += 1
                disagreements.append({"a": universe[a]["gene"], "b": universe[b]["gene"],
                                      "grouped_by": "mmseqs_only",
                                      "local_group_a": la, "local_group_b": lb,
                                      "mmseqs_group": mm_group[a]})


    barred_mm = {mm_group[k] for k in keys
                 if local.get(k, {}).get("role") in CONFIRM_ROLES}
    changed = []
    for k in keys:
        role = local.get(k, {}).get("role", "?")
        was_eligible = local.get(k, {}).get("fitting_eligible") == "True"
        now_eligible = (mm_group[k] not in barred_mm) and role not in CONFIRM_ROLES
        if was_eligible != now_eligible:
            changed.append({"gene": universe[k]["gene"], "accession": k, "role": role,
                            "biopython_fitting_eligible": was_eligible,
                            "mmseqs_fitting_eligible": now_eligible,
                            "local_group": local.get(k, {}).get("group_id"),
                            "mmseqs_group": mm_group[k]})

    rows = [{"accession": k, "gene": universe[k]["gene"],
             "role": local.get(k, {}).get("role", "?"),
             "mmseqs_group": mm_group[k],
             "mmseqs_cluster_rep": cluster_of.get(k, k),
             "biopython_group": local.get(k, {}).get("group_id", "?")}
            for k in keys]
    rows.sort(key=lambda r: (r["mmseqs_group"], r["gene"]))
    with (OUT / "group_manifest_mmseqs.csv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    if uniq:
        with (OUT / "overlap_audit_mmseqs.csv").open("w", encoding="utf-8",
                                                     newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(uniq[0].keys()))
            w.writeheader()
            w.writerows(sorted(uniq, key=lambda p: -p["identity"]))

    summary = {
        "sequences": len(keys),
        "mmseqs_version": subprocess.run([MMSEQS, "version"], capture_output=True,
                                         text=True).stdout.strip(),
        "stage1_clustering": {
            "command": f"easy-cluster --min-seq-id {MIN_SEQ_ID_CLUSTER} "
                       f"-c {COV_CLUSTER} --cov-mode 0",
            "n_clusters": n_clusters,
        },
        "stage2_overlap_audit": {
            "command": f"easy-search --min-seq-id {MIN_SEQ_ID_AUDIT} "
                       f"-c {COV_AUDIT} --cov-mode 0 -e 1e-3",
            "criterion": "identity >= 0.30 and min(qcov, tcov) >= 0.50",
            "qualifying_pairs": len(uniq),
        },
        "combined_groups": len(groups),
        "biopython_groups": len({r["biopython_group"] for r in rows}),
        "pair_level_agreement": {
            "co_grouped_in_both": same_both,
            "co_grouped_biopython_only": same_local_only,
            "co_grouped_mmseqs_only": same_mm_only,
            "note": "counted over all 89*88/2 accession pairs",
        },
        "split_changes": changed,
        "split_unchanged": not changed,
        "disagreements": disagreements[:60],
        "outcome_blind": "sequences and metadata only; no assay score read",
    }
    (OUT / "mmseqs_vs_biopython.json").write_text(json.dumps(summary, indent=2),
                                                  encoding="utf-8")

    print(f"\n--- comparison ---")
    print(f"  groups: mmseqs {len(groups)}  biopython "
          f"{len({r['biopython_group'] for r in rows})}")
    print(f"  pairs co-grouped in both        : {same_both}")
    print(f"  pairs co-grouped biopython only : {same_local_only}")
    print(f"  pairs co-grouped mmseqs only    : {same_mm_only}")
    print(f"\n  SPLIT CHANGES: {len(changed)}")
    for c in changed:
        print(f"    {c['gene']} ({c['role']}): fitting_eligible "
              f"{c['biopython_fitting_eligible']} -> {c['mmseqs_fitting_eligible']}")
    if not changed:
        print("    none - the split is unchanged, no rebuild needed")
    if disagreements:
        print(f"\n  example disagreements ({len(disagreements)} total):")
        for d in disagreements[:12]:
            print(f"    {d['a']} / {d['b']}: grouped by {d['grouped_by']}")


if __name__ == "__main__":
    main()

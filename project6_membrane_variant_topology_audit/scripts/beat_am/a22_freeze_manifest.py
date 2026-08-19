"""Freeze the curated protein split and verify exclusion of confirmation relatives."""

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
OUT = ROOT / "results" / "beat_am_20261016" / "inventory_frozen"
FEAT = ROOT / "results" / "beat_am_20261016" / "features"

SHUFFLE_P_MAX = 0.05
KCNQ4 = "P56696"
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


def main() -> None:
    from a4_group_manifest import build_universe

    OUT.mkdir(parents=True, exist_ok=True)
    universe = build_universe()
    keys = sorted(universe)

    prior = {}
    with (CUR / "group_manifest_curated.csv").open(encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh):
            prior[r["accession"]] = r

    shuffle = {}
    with (CUR / "shuffle_test.csv").open(encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh):
            shuffle[tuple(sorted((r["accession_a"], r["accession_b"])))] = r

    kept, dropped = [], []
    for pair, r in shuffle.items():
        p = float(r["empirical_p"])
        (kept if p <= SHUFFLE_P_MAX else dropped).append((pair, r, p))
    print(f"shuffle-test filter at p <= {SHUFFLE_P_MAX}: "
          f"kept {len(kept)}, dropped {len(dropped)}")
    for pair, r, p in sorted(dropped, key=lambda x: -x[2]):
        print(f"  dropped  {r['gene_a']:8s} {r['gene_b']:8s} z={float(r['z_score']):6.2f} "
              f"p={p:.4f}  ({r['prior_verdict']})")
    print("  (the 5 ARTIFACT pairs were already excluded before the shuffle test)")

    uf = Union(keys)
    for pair, r, p in kept:
        a, b = pair
        if a in uf.parent and b in uf.parent:
            uf.union(a, b)

    groups: dict[str, list[str]] = {}
    for k in keys:
        groups.setdefault(uf.find(k), []).append(k)
    gid = {root: f"F{n:03d}" for n, root in enumerate(sorted(groups), start=1)}
    gof = {k: gid[uf.find(k)] for k in keys}


    role = {k: prior.get(k, {}).get("role", "unassigned") for k in keys}
    before = role.get(KCNQ4)
    role[KCNQ4] = "confirmation_membrane"
    print(f"\nKCNQ4 ({KCNQ4}): role {before} -> confirmation_membrane")
    print(f"  its group is now {gof[KCNQ4]} = "
          f"{sorted(universe[m]['gene'] for m in groups[uf.find(KCNQ4)])}")

    confirm_groups = {gof[k] for k in keys if role[k] in CONFIRM_ROLES}
    barred = {k for k in keys if gof[k] in confirm_groups}


    t = np.load(FEAT / "l1_table.npz", allow_pickle=True)
    fitted = set(t["accession"].astype(str))
    violations = sorted(barred & fitted)
    print(f"\n=== GATE: barred-from-fitting proteins that were actually fitted ===")
    print(f"  barred set: {len(barred)}   fitted set: {len(fitted)}")
    if violations:
        for v in violations:
            print(f"  !! VIOLATION {universe[v]['gene']} ({v}) role={role[v]} "
                  f"group={gof[v]}")
    else:
        print("  NONE -- no confirmation protein's relative was ever fitted")

    rows = []
    for k in keys:
        rows.append({"accession": k, "gene": universe[k]["gene"],
                     "role": role[k], "frozen_group": gof[k],
                     "prior_curated_group": prior.get(k, {}).get("curated_group", ""),
                     "fitting_eligible": (gof[k] not in confirm_groups)
                                         and role[k] not in CONFIRM_ROLES,
                     "was_fitted_in_l1": k in fitted})
    rows.sort(key=lambda r: (r["frozen_group"], r["gene"]))
    with (OUT / "group_manifest_frozen.csv").open("w", encoding="utf-8",
                                                  newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    conf_memb = sorted(universe[k]["gene"] for k in keys
                       if role[k] == "confirmation_membrane")
    conf_all = sorted(universe[k]["gene"] for k in keys
                      if role[k] == "confirmation_allhuman")
    n_fit = sum(1 for r in rows if r["fitting_eligible"])

    summary = {
        "shuffle_filter": {"p_max": SHUFFLE_P_MAX, "pairs_kept": len(kept),
                           "pairs_dropped": len(dropped),
                           "dropped": [{"a": r["gene_a"], "b": r["gene_b"],
                                        "z": float(r["z_score"]),
                                        "p": float(r["empirical_p"])}
                                       for _, r, _ in dropped],
                           "rationale": "applied uniformly to all 17 non-artifact "
                                        "pairs; measured in session 11 before this "
                                        "decision existed"},
        "kcnq4_admitted": {"accession": KCNQ4, "role_before": before,
                           "group": gof[KCNQ4],
                           "group_members": sorted(universe[m]["gene"] for m in
                                                   groups[uf.find(KCNQ4)])},
        "groups": len(groups),
        "confirmation_membrane": conf_memb,
        "n_confirmation_membrane": len(conf_memb),
        "confirmation_allhuman_n": len(conf_all),
        "barred_from_fitting": len(barred),
        "fitting_eligible": n_fit,
        "gate_violations": [{"gene": universe[v]["gene"], "accession": v,
                             "role": role[v]} for v in violations],
        "gate_passed": not violations,
        "beta_must_be_refit": True,
        "beta_refit_reason": "AMFR (shuffle z=4.70, kept) shares KCNQ4's group and "
                             "sits in the census sequence-TRAIN partition that beta "
                             "was fitted on",
    }
    (OUT / "freeze_summary.json").write_text(json.dumps(summary, indent=2),
                                             encoding="utf-8")

    print(f"\ngroups {len(groups)}   fitting-eligible {n_fit}   barred {len(barred)}")
    print(f"confirmation_membrane now {len(conf_memb)}: {conf_memb}")
    print(f"confirmation_allhuman: {len(conf_all)} proteins")
    print(f"\nGATE PASSED: {not violations}")


if __name__ == "__main__":
    main()

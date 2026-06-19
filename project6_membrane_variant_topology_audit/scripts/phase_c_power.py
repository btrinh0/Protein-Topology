"""Phase C feasibility, computed from Phase B counts only.

Reads the DMS inventory, which records how many variants fall in each flank cell.
It does not read any measured effect, so running this does not break the Phase C
blind on DMS outcomes.

The binding constraint in Phase C is not the variant count, it is the number of
independent proteins. Three GPCRs assayed by surface expression contribute most of
the data, so this script reports the leave-one-family-out picture alongside the
full one.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INVENTORY = ROOT / "results" / "phaseB_20260924" / "dms_inventory.csv"
DEFAULT_OUT = ROOT / "results" / "prereg_phase_c"

GPCR_FAMILY = {"CCR5", "CXCR4", "GPR68", "RHO"}
CELLS = {
    "positive_charge_gain": ("flank_10_cytosolic_positive_charge_gain",
                             "flank_10_outer_positive_charge_gain"),
    "positive_charge_loss": ("flank_10_cytosolic_positive_charge_loss",
                             "flank_10_outer_positive_charge_loss"),
}
ALPHA, POWER = 0.05, 0.80
Z_ALPHA, Z_POWER = 1.959964, 0.8416212


def best_per_gene(rows: list[dict]) -> list[dict]:
    """One score set per gene: the one covering the most flank charge changes.

    A gene can appear several times in the inventory, once per deposited score set,
    and those score sets measure overlapping variants. Summing them would both
    double-count variants and inflate the apparent number of independent clusters,
    which is exactly the quantity Phase C is constrained by.
    """
    chosen: dict[str, dict] = {}
    for row in rows:
        gene = row["gene"]
        coverage = int(row["flank_10_charge_change"])
        if gene not in chosen or coverage > int(chosen[gene]["flank_10_charge_change"]):
            chosen[gene] = row
    return list(chosen.values())


def cell_feasibility(rows: list[dict], cytosolic_field: str, outer_field: str) -> dict:
    per_protein = []
    for row in best_per_gene(rows):
        cytosolic, outer = int(row[cytosolic_field]), int(row[outer_field])
        if cytosolic + outer > 0:
            per_protein.append({
                "gene": row["gene"],
                "cytosolic": cytosolic,
                "outer": outer,
                "both_sides": cytosolic > 0 and outer > 0,
            })
    n1 = sum(p["cytosolic"] for p in per_protein)
    n2 = sum(p["outer"] for p in per_protein)
    both = [p for p in per_protein if p["both_sides"]]
    total = n1 + n2
    mean_cluster = total / len(per_protein) if per_protein else 0.0
    largest = max((p["cytosolic"] + p["outer"] for p in per_protein), default=0)

    detectable = {}
    for icc in (0.0, 0.02, 0.05, 0.10, 0.20):
        design_effect = 1.0 + (mean_cluster - 1.0) * icc
        e1, e2 = n1 / design_effect, n2 / design_effect
        detectable[f"icc_{icc:g}"] = {
            "design_effect": round(design_effect, 2),
            "effective_n": [round(e1, 1), round(e2, 1)],
            "detectable_standardised_difference": (
                round((Z_ALPHA + Z_POWER) * math.sqrt(1 / e1 + 1 / e2), 3)
                if e1 > 0 and e2 > 0 else None),
        }
    return {
        "n_proteins": len(per_protein),
        "n_proteins_sampling_both_sides": len(both),
        "n_cytosolic": n1,
        "n_outer": n2,
        "mean_measurements_per_protein": round(mean_cluster, 1),
        "largest_protein_share": round(largest / total, 3) if total else None,
        "per_protein": sorted(per_protein, key=lambda p: -(p["cytosolic"] + p["outer"])),
        "detectable_at_80pct_power": detectable,
    }


def run(args: argparse.Namespace) -> dict:
    with Path(args.inventory).open(encoding="utf-8", newline="") as handle:
        rows = [row for row in csv.DictReader(handle)
                if row.get("expression_readout", "").strip().lower() == "true"]

    scopes = {
        "all_proteins": rows,
        "leave_one_family_out_no_gpcr": [r for r in rows if r["gene"] not in GPCR_FAMILY],
    }
    results = {
        scope: {name: cell_feasibility(subset, *fields) for name, fields in CELLS.items()}
        for scope, subset in scopes.items()
    }

    summary = {
        "analysis": "Phase C feasibility from Phase B counts (no DMS effects read)",
        "window_residues": 10,
        "alpha": ALPHA, "target_power": POWER,
        "caveat": "Normal-approximation detectable effects are optimistic. With about a dozen "
                  "protein clusters, cluster-robust standard errors are downward biased; "
                  "Phase C uses a wild cluster bootstrap-t and should expect wider intervals.",
        "results": results,
    }
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "phase_c_feasibility.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    for scope, block in results.items():
        print(f"\n=== {scope} ===")
        for name, cell in block.items():
            print(f"  {name}: n={cell['n_cytosolic']}/{cell['n_outer']}  "
                  f"proteins={cell['n_proteins']} (both sides {cell['n_proteins_sampling_both_sides']})  "
                  f"largest share {cell['largest_protein_share']}")
            for icc, values in cell["detectable_at_80pct_power"].items():
                print(f"      {icc:<9} DE={values['design_effect']:<6} "
                      f"detectable d={values['detectable_standardised_difference']}")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", default=str(DEFAULT_INVENTORY))
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    run(parser.parse_args())


if __name__ == "__main__":
    main()

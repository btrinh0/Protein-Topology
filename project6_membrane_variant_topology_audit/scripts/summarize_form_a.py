"""Summarize Form A estimates and evaluate the declared pass/fail criteria."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

import form_ab_lib as lib

PLANTED = ("P_neg05", "P_neg10", "P_pos10", "P_decay")
NULLS = ("N0", "S_main")
GATES = {"s0_exact": 1e-9, "s2_identity": 1e-9, "coverage_min": 0.90, "width_ratio": (1.7, 2.3),
         "bias_se_multiple": 4.0, "slope": (0.7, 1.1), "sign_min": 0.95, "null_within_min": 0.90}


def read(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        row["replicate"] = int(row["replicate"])
        for name in ("oracle_raw", "independent_raw", "estimate_raw", "ci_lo", "ci_hi", "penalty_sd"):
            row[name] = float(row[name])
    return rows


def table(rows: list[dict], sd_z0: float) -> list[dict]:
    out = []
    for scenario in lib.SCENARIOS:
        for student in ("S0", "S0n", "S1", "S2", "S3"):
            cell = [r for r in rows if r["scenario"] == scenario and r["student"] == student]
            if not cell:
                continue
            est = np.array([r["estimate_raw"] for r in cell])
            oracle = cell[0]["oracle_raw"]
            lo, hi = np.array([r["ci_lo"] for r in cell]), np.array([r["ci_hi"] for r in cell])
            out.append({
                "scenario": scenario, "student": student, "n": len(cell), "oracle": oracle,
                "mean_estimate": float(est.mean()), "sd_across_replicates": float(est.std(ddof=1)) if len(est) > 1 else None,
                "bias_to_oracle": float(est.mean() - oracle),
                "recovery_ratio": float(est.mean() / oracle) if oracle else None,
                "ci_coverage_of_oracle": float(np.mean((lo <= oracle) & (oracle <= hi))),
                "mean_ci_half_width": float(np.mean((hi - lo) / 2)),
                "fraction_ci_excludes_zero": float(np.mean((lo > 0) | (hi < 0))),
                "fraction_correct_sign": float(np.mean(np.sign(est) == np.sign(oracle))) if oracle else None,
                "fraction_within_margin": float(np.mean(np.abs(est) <= lib.EQUIVALENCE * sd_z0)),
                "max_abs_estimate": float(np.abs(est).max()),
                "max_abs_error_to_oracle": float(np.abs(est - oracle).max())})
    return out


def criteria(rows: list[dict], tab: list[dict], sd_z0: float) -> dict:
    by = {(t["scenario"], t["student"]): t for t in tab}
    g = GATES
    out = {"gates": {k: list(v) if isinstance(v, tuple) else v for k, v in g.items()}}
    s0 = [r for r in rows if r["student"] == "S0"]
    out["S0_exact"] = {"max_abs_error": max(abs(r["estimate_raw"] - r["oracle_raw"]) for r in s0),
                       "passed": all(abs(r["estimate_raw"] - r["oracle_raw"]) <= g["s0_exact"] for r in s0)}
    s0n = [t for t in tab if t["student"] == "S0n"]
    s0n_rows = [r for r in rows if r["student"] == "S0n"]
    cover = float(np.mean([r["ci_lo"] <= r["oracle_raw"] <= r["ci_hi"] for r in s0n_rows]))
    pooled_sd = float(np.sqrt(np.mean([t["sd_across_replicates"] ** 2 for t in s0n])))
    ratio = float(np.mean([(r["ci_hi"] - r["ci_lo"]) / 2 for r in s0n_rows]) / pooled_sd)
    bias_ok = all(abs(t["bias_to_oracle"]) <= g["bias_se_multiple"] * t["sd_across_replicates"] / np.sqrt(t["n"])
                  for t in s0n)
    out["S0n_calibration"] = {"pooled_coverage": cover, "half_width_over_pooled_sd": ratio,
                              "pooled_sd": pooled_sd, "bias_within_4_se": bias_ok,
                              "passed": bool(cover >= g["coverage_min"] and g["width_ratio"][0] <= ratio <= g["width_ratio"][1]
                                             and bias_ok)}
    s1 = {}
    est = np.array([by[(s, "S1")]["mean_estimate"] for s in PLANTED if (s, "S1") in by])
    ora = np.array([by[(s, "S1")]["oracle"] for s in PLANTED if (s, "S1") in by])
    slope = float((est * ora).sum() / (ora ** 2).sum())
    s1["slope_through_origin"] = slope
    s1["sign_fraction"] = {s: by[(s, "S1")]["fraction_correct_sign"] for s in PLANTED if (s, "S1") in by}
    s1["null_within_margin"] = {s: by[(s, "S1")]["fraction_within_margin"] for s in NULLS if (s, "S1") in by}
    s1["recovery_ratio"] = {s: by[(s, "S1")]["recovery_ratio"] for s in PLANTED if (s, "S1") in by}
    s1["passed"] = bool(g["slope"][0] <= slope <= g["slope"][1]
                        and all(v >= g["sign_min"] for v in s1["sign_fraction"].values())
                        and all(v >= g["null_within_min"] for v in s1["null_within_margin"].values()))
    out["S1"] = s1
    s2 = [r for r in rows if r["student"] == "S2"]
    out["S2_identity"] = {"max_abs_estimate": max(abs(r["estimate_raw"]) for r in s2),
                          "passed": all(abs(r["estimate_raw"]) <= g["s2_identity"] for r in s2)}
    out["independent_arithmetic"] = {
        "max_abs_difference": max(abs(r["estimate_raw"] - r["independent_raw"]) for r in rows),
        "n_runs": len(rows), "passed": all(abs(r["estimate_raw"] - r["independent_raw"]) <= g["s0_exact"]
                                           for r in rows)}
    out["S3_reported"] = {s: {"mean_estimate": by[(s, "S3")]["mean_estimate"], "oracle": by[(s, "S3")]["oracle"],
                              "recovery_ratio": by[(s, "S3")]["recovery_ratio"]}
                          for s in lib.SCENARIOS if (s, "S3") in by}
    out["form_b_may_start"] = bool(out["S0_exact"]["passed"] and out["S0n_calibration"]["passed"]
                                   and out["S1"]["passed"] and out["S2_identity"]["passed"]
                                   and out["independent_arithmetic"]["passed"])
    return out


def naive_table(path: Path, sd_z0: float) -> list[dict]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    out = []
    for scenario in NULLS:
        for student in ("S1", "S2"):
            cell = [r for r in rows if r["scenario"] == scenario and r["student"] == student]
            if cell:
                out.append({"scenario": scenario, "student": student, "n": len(cell),
                            **{name: float(np.mean([float(r[name]) for r in cell])) for name in
                               ("matched", "charged_arm_only", "pooled_paired")},
                            "margin": lib.EQUIVALENCE * sd_z0, "sd_z0": sd_z0})
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads((args.results / "run_manifest.json").read_text(encoding="utf-8"))
    sd = manifest["sd_z0"]
    rows = read(args.results / "estimates.csv")
    tab = table(rows, sd)
    verdict = criteria(rows, tab, sd)
    naive = naive_table(args.results / "naive_comparators.csv", sd)
    with (args.results / "summary_tables.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(tab[0]))
        writer.writeheader()
        writer.writerows(tab)
    (args.results / "criteria.json").write_text(json.dumps({**verdict, "naive_comparators": naive,
                                                            "stage": manifest["stage"]}, indent=2) + "\n",
                                                encoding="utf-8")
    print(json.dumps({k: (v if k != "gates" else "...") for k, v in verdict.items()
                      if k not in ("S3_reported",)}, indent=2))


if __name__ == "__main__":
    main()

"""Design quantities for the preregistration, computed without unblinding outcomes.

Reads the blinded flank-event table. The only outcome information it uses is the
pooled P/LP fraction across all flank events, which Gate 1 already reported and
which is not split by side or charge class.

Reports, per preregistered cell: cluster structure (events per gene), the number of
conditional strata that carry information, and the difference a two-sided test can
detect at 80% power once gene clustering is taken into account.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EVENTS = ROOT / "data" / "processed" / "gate0" / "blinded_flank_events.csv"
DEFAULT_OUT = ROOT / "results" / "prereg_v1"

# Gate 1, pooled across all flank events and both charge classes. Not a per-cell rate.
POOLED_PATHOGENIC = 935
POOLED_LABELED = 1724

ALPHA = 0.05
TARGET_POWER = 0.80

POSITIVE = set("KR")
NEGATIVE = set("DE")
AMBIGUOUS = set("H")


def substitution_class(ref: str, alt: str) -> str | None:
    """Same class definitions as run_gate0.py, so the two gates line up."""
    if ref in POSITIVE and alt not in POSITIVE:
        return "positive_charge_loss"
    if ref not in POSITIVE and alt in POSITIVE:
        return "positive_charge_gain"
    if ref in NEGATIVE and alt not in NEGATIVE:
        return "negative_charge_loss"
    if ref not in NEGATIVE and alt in NEGATIVE:
        return "negative_charge_gain"
    if not ({ref, alt} & (POSITIVE | NEGATIVE | AMBIGUOUS)):
        return "charge_neutral_placebo"
    return None


def normal_cdf(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def normal_quantile(p: float) -> float:
    low, high = -12.0, 12.0
    for _ in range(200):
        middle = (low + high) / 2.0
        if normal_cdf(middle) < p:
            low = middle
        else:
            high = middle
    return (low + high) / 2.0


def power_for_difference(delta: float, n1: float, n2: float, baseline: float) -> float:
    if n1 <= 0 or n2 <= 0:
        return 0.0
    weight1, weight2 = n1 / (n1 + n2), n2 / (n1 + n2)
    p1 = baseline + delta * weight2
    p2 = baseline - delta * weight1
    if not (0.0 < p1 < 1.0 and 0.0 < p2 < 1.0):
        return float("nan")
    pooled = weight1 * p1 + weight2 * p2
    se_null = math.sqrt(pooled * (1 - pooled) * (1 / n1 + 1 / n2))
    se_alt = math.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
    critical = normal_quantile(1 - ALPHA / 2) * se_null
    return normal_cdf((delta - critical) / se_alt) + normal_cdf((-delta - critical) / se_alt)


def detectable_difference(n1: float, n2: float, baseline: float) -> float | None:
    low, high = 0.0, min(baseline, 1 - baseline) * 2.0
    if power_for_difference(high, n1, n2, baseline) < TARGET_POWER:
        return None
    for _ in range(200):
        middle = (low + high) / 2.0
        if power_for_difference(middle, n1, n2, baseline) < TARGET_POWER:
            low = middle
        else:
            high = middle
    return (low + high) / 2.0


def odds_ratio(delta: float, baseline: float, n1: float, n2: float) -> float:
    weight1, weight2 = n1 / (n1 + n2), n2 / (n1 + n2)
    p1, p2 = baseline + delta * weight2, baseline - delta * weight1
    return (p1 / (1 - p1)) / (p2 / (1 - p2))


def load_events(path: Path) -> list[dict]:
    import csv
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def cell_report(events: list[dict], window: int, class_name: str, baseline: float) -> dict:
    subset = [event for event in events
              if int(event["distance_to_tmd"]) <= window
              and substitution_class(event["ref_aa"], event["alt_aa"]) == class_name]
    by_side = defaultdict(list)
    for event in subset:
        by_side[event["side"]].append(event)
    counts = {side: len(rows) for side, rows in by_side.items()}
    n1, n2 = counts.get("cytosolic", 0), counts.get("outer", 0)

    genes = {side: Counter(event["gene"] for event in rows) for side, rows in by_side.items()}
    all_genes = Counter()
    for counter in genes.values():
        all_genes.update(counter)
    both_sides = set(genes.get("cytosolic", {})) & set(genes.get("outer", {}))

    strata = defaultdict(lambda: {"cytosolic": 0, "outer": 0})
    for event in subset:
        strata[(event["gene"], f"{event['ref_aa']}>{event['alt_aa']}")][event["side"]] += 1
    informative = {key: value for key, value in strata.items()
                   if value["cytosolic"] > 0 and value["outer"] > 0}

    mean_cluster = (sum(all_genes.values()) / len(all_genes)) if all_genes else 0.0
    report = {
        "n_cytosolic": n1,
        "n_outer": n2,
        "n_genes": len(all_genes),
        "n_genes_with_both_sides": len(both_sides),
        "events_in_genes_with_both_sides": sum(
            genes["cytosolic"][gene] + genes["outer"][gene] for gene in both_sides),
        "mean_events_per_gene": round(mean_cluster, 3),
        "largest_gene_share_of_cell": round(max(all_genes.values()) / (n1 + n2), 4) if all_genes else None,
        "conditional_strata_gene_x_substitution_both_sides": len(informative),
        "events_in_conditional_strata": sum(
            value["cytosolic"] + value["outer"] for value in informative.values()),
        "detectable_difference_80pct_power": {},
    }
    for icc in (0.0, 0.02, 0.05, 0.10):
        design_effect = 1.0 + (mean_cluster - 1.0) * icc
        delta = detectable_difference(n1 / design_effect, n2 / design_effect, baseline)
        report["detectable_difference_80pct_power"][f"icc_{icc:g}"] = {
            "design_effect": round(design_effect, 3),
            "effective_n_cytosolic": round(n1 / design_effect, 1),
            "effective_n_outer": round(n2 / design_effect, 1),
            "percentage_points": round(delta * 100, 1) if delta else None,
            "equivalent_odds_ratio": round(odds_ratio(delta, baseline, n1, n2), 2) if delta else None,
        }
    return report


def did_report(cells: dict, window: int, left: str, right: str, baseline: float) -> dict:
    """Power for a difference of two independent risk differences.

    Gate 0 showed a raw side comparison is confounded by a side baseline that also
    affects substitutions changing no charge, so the primary clinical contrast is
    (side difference in one class) minus (side difference in a reference class).
    Under the null both risk differences are zero, so their standard errors add in
    quadrature and the detectable effects combine the same way.
    """
    def standard_error(cell: dict) -> float | None:
        n1, n2 = cell["n_cytosolic"], cell["n_outer"]
        if n1 == 0 or n2 == 0:
            return None
        return math.sqrt(baseline * (1 - baseline) * (1 / n1 + 1 / n2))

    left_cell, right_cell = cells[f"window{window}_{left}"], cells[f"window{window}_{right}"]
    se_left, se_right = standard_error(left_cell), standard_error(right_cell)
    if se_left is None or se_right is None:
        return {"detectable_difference_in_differences_pp": None}
    multiplier = normal_quantile(1 - ALPHA / 2) + normal_quantile(TARGET_POWER)
    return {
        "left": left, "right": right,
        "left_n": [left_cell["n_cytosolic"], left_cell["n_outer"]],
        "right_n": [right_cell["n_cytosolic"], right_cell["n_outer"]],
        "detectable_difference_in_differences_pp": round(
            multiplier * math.sqrt(se_left ** 2 + se_right ** 2) * 100, 1),
    }


def run(args: argparse.Namespace) -> dict:
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    events = load_events(Path(args.events))
    baseline = POOLED_PATHOGENIC / POOLED_LABELED

    classes = ["positive_charge_loss", "positive_charge_gain", "negative_charge_loss",
               "negative_charge_gain", "charge_neutral_placebo"]
    cells = {}
    for window in (5, 10, 15):
        for class_name in classes:
            cells[f"window{window}_{class_name}"] = cell_report(
                events, window, class_name, baseline)

    contrasts = {}
    for window in (5, 10, 15):
        for left, right in (("positive_charge_loss", "charge_neutral_placebo"),
                            ("positive_charge_loss", "positive_charge_gain"),
                            ("negative_charge_loss", "charge_neutral_placebo")):
            contrasts[f"window{window}_{left}_minus_{right}"] = did_report(
                cells, window, left, right, baseline)

    summary = {
        "analysis": "Preregistration design quantities (outcome-blind)",
        "difference_in_differences_power": contrasts,
        "baseline_pathogenic_fraction": round(baseline, 4),
        "baseline_scope": "pooled across all flank events, not split by side or charge class",
        "alpha": ALPHA,
        "target_power": TARGET_POWER,
        "blinded_events_rows": len(events),
        "cells": cells,
    }
    (out_dir / "power_and_clustering.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    for name, report in cells.items():
        print(f"\n{name}: n={report['n_cytosolic']}/{report['n_outer']} "
              f"genes={report['n_genes']} (both sides {report['n_genes_with_both_sides']}) "
              f"mean events/gene={report['mean_events_per_gene']} "
              f"strata={report['conditional_strata_gene_x_substitution_both_sides']}")
        for icc, values in report["detectable_difference_80pct_power"].items():
            print(f"    {icc:<9} DE={values['design_effect']:<6} "
                  f"detectable {values['percentage_points']} pp  (OR {values['equivalent_odds_ratio']})")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--events", default=str(DEFAULT_EVENTS))
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    run(parser.parse_args())


if __name__ == "__main__":
    main()

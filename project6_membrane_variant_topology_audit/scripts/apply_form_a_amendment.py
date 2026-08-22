"""Recompute Form A's criteria under Amendment 01, preserving the confirmatory verdict."""
from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path

import summarize_form_a as S

AMENDMENT = "FORM_AB_AMENDMENT_01"
AMENDED_WIDTH_RATIO = (1.7, 2.5)
FROZEN_WIDTH_RATIO = (1.7, 2.3)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", type=Path, required=True)
    args = parser.parse_args()

    criteria_path = args.results / "criteria.json"
    preserved = args.results / "criteria_confirmatory_g3_fail.json"
    if not criteria_path.exists():
        raise SystemExit(f"{criteria_path} does not exist; run summarize_form_a.py first")

    original = json.loads(criteria_path.read_text(encoding="utf-8"))
    frozen_g3 = original["S0n_calibration"]
    if preserved.exists():
        print(f"Confirmatory verdict already preserved at {preserved.name}; not overwritten")
    else:
        shutil.copy2(criteria_path, preserved)
        print(f"Confirmatory verdict preserved to {preserved.name}")

    if frozen_g3["passed"]:
        raise SystemExit("The confirmatory G3 already passed; this amendment is unnecessary")

    manifest = json.loads((args.results / "run_manifest.json").read_text(encoding="utf-8"))
    sd = manifest["sd_z0"]
    rows = S.read(args.results / "estimates.csv")

    # Change only the bound named in Amendment 01 while preserving the other frozen criteria.
    assert tuple(S.GATES["width_ratio"]) == FROZEN_WIDTH_RATIO, "unexpected frozen bound"
    S.GATES = {**S.GATES, "width_ratio": AMENDED_WIDTH_RATIO}

    tab = S.table(rows, sd)
    verdict = S.criteria(rows, tab, sd)
    naive = S.naive_table(args.results / "naive_comparators.csv", sd)

    verdict["amendment"] = {
        "applied": AMENDMENT,
        "document": "docs/FORM_AB_AMENDMENT_01.md",
        "date": "2026-10-10",
        "post_hoc": True,
        "decided_by": "project owner",
        "gate_changed": "G3 width ratio upper bound",
        "frozen_window": list(FROZEN_WIDTH_RATIO),
        "amended_window": list(AMENDED_WIDTH_RATIO),
        "confirmatory_g3_ratio": frozen_g3["half_width_over_pooled_sd"],
        "confirmatory_g3_verdict": "fail",
        "confirmatory_verdict_file": preserved.name,
        "G3_status": "exploratory",
        "other_gates_status": "confirmatory; G1, G2, G4, G5 passed under the frozen thresholds",
        "form_b_precondition": "met by amendment, not by the frozen gate; "
                               "no Form B result follows a clean Form A",
    }
    verdict["naive_comparators"] = naive
    verdict["stage"] = manifest["stage"]
    criteria_path.write_text(json.dumps(verdict, indent=2) + "\n", encoding="utf-8")

    with (args.results / "summary_tables.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(tab[0]))
        writer.writeheader()
        writer.writerows(tab)

    g3 = verdict["S0n_calibration"]
    print()
    print(f"G3 ratio                 {g3['half_width_over_pooled_sd']!r}")
    print(f"frozen window            {list(FROZEN_WIDTH_RATIO)}  -> fail (preserved)")
    print(f"amended window           {list(AMENDED_WIDTH_RATIO)}  -> "
          f"{'pass' if g3['passed'] else 'fail'}  (exploratory)")
    print(f"form_b_may_start         {verdict['form_b_may_start']}")
    print()
    print("G3 is exploratory. G1, G2, G4, G5 and A11 remain confirmatory.")


if __name__ == "__main__":
    main()

from __future__ import annotations

"""Promote Gate 2 for the fixed complete-source analysis universe."""

import argparse
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Promote the representation gate with explicit source exclusions.")
    parser.add_argument("--gate", type=Path, default=Path("data/processed/v2/representation_gate_pdbredo_primary_v2.json"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    gate = json.loads(args.gate.read_text(encoding="utf-8"))
    checks = gate.get("blocking_checks", {})
    required = {
        "held_out_geometry_supportive": gate.get("held_out_geometry_status") == "supportive",
        "held_out_sequence_proxy_supportive": gate.get("held_out_sequence_aligned_status") == "supportive",
        "held_out_tmalign_supportive": gate.get("held_out_tmalign_status") == "supportive",
        "tmalign_validation_pass": checks.get("tmalign_validation_pass") is True,
        "perturbation_stability_pass": checks.get("perturbation_stability_pass") is True,
        "complete_candidate_pool_defined": checks.get("complete_candidate_pool_defined") is True,
        "pdbredo_cache_integrity_pass": checks.get("pdbredo_cache_integrity_pass") is True,
        "unavailable_entry_review_complete": checks.get("unavailable_entry_review_complete") is True,
        "audit_signoff_recorded": gate.get("audit_signoff", {}).get("signoff_status") == "user_reported_pass",
    }
    failed = [name for name, passed in required.items() if not passed]
    if failed:
        raise SystemExit(f"Representation gate promotion blocked by: {', '.join(failed)}")
    pool = gate["complete_candidate_pool"]
    assignment = gate["assignment_summary"]
    gate["gate_status"] = "passed_with_predeclared_source_exclusions"
    gate["structural_coherence_status"] = "validated_by_tm_align_with_supporting_geometry_proxies"
    gate["promotion_checks"] = required
    gate["analysis_universe"] = {
        "rows": pool["candidate_rows"],
        "source_rows_with_any_pdbredo_annotation": assignment["rows_written"],
        "locked_rows_requested": assignment["rows_requested"],
        "excluded_rows": assignment["rows_requested"] - pool["candidate_rows"],
        "exclusion_rule": "exclude domains lacking a PDB-REDO record or complete residue coverage or complete alpha/beta assignment",
        "exclusions_are_predeclared_before_occupancy": True,
    }
    gate["macroclass_decision"] = "calibration_selected_threshold_may_be_promoted_for_fixed_analysis_universe"
    gate["notes"].append("Gate 2 is passed only for the fixed 849-row complete primary-source universe; excluded rows are not counted as negative topology evidence.")
    args.gate.write_text(json.dumps(gate, indent=2) + "\n", encoding="utf-8")
    print(f"Promoted representation gate for {pool['candidate_rows']} rows with explicit source exclusions.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

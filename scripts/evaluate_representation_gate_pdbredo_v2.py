from __future__ import annotations

"""Summarize primary-source SSE evidence for the V2 representation gate."""

import argparse
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate PDB-REDO DSSP representation evidence without assigning final classes.")
    parser.add_argument("--assignment", type=Path, default=Path("data/processed/v2/domain_topology_pdbredo_dssp_summary_v2.json"))
    parser.add_argument("--calibration", type=Path, default=Path("data/processed/v2/topology_signature_calibration_pdbredo_primary_summary_v2.json"))
    parser.add_argument("--geometry", type=Path, default=Path("data/processed/v2/structural_coherence_pdbredo_summary_v2.json"))
    parser.add_argument("--sequence", type=Path, default=Path("data/processed/v2/sequence_aligned_structure_pdbredo_summary_v2.json"))
    parser.add_argument("--audit-signoff", type=Path, default=Path("data/processed/v2/domain_mapping_audit_user_signoff_v2.json"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/representation_gate_pdbredo_v2.json"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    assignment = json.loads(args.assignment.read_text(encoding="utf-8"))
    calibration = json.loads(args.calibration.read_text(encoding="utf-8"))
    geometry = json.loads(args.geometry.read_text(encoding="utf-8"))
    sequence = json.loads(args.sequence.read_text(encoding="utf-8"))
    audit = json.loads(args.audit_signoff.read_text(encoding="utf-8"))
    held_out_geometry = [row for row in geometry["results"] if row["split"] == "held_out"]
    held_out_sequence = [
        row for row in sequence["results"]
        if row["split"] == "held_out" and row["pair_type"] == "effect_between_minus_within"
    ]
    output = {
        "assignment_summary": assignment,
        "audit_signoff": audit,
        "held_out_geometry_effects": [row["coherence_effect_between_minus_within"] for row in held_out_geometry],
        "held_out_geometry_status": "supportive" if held_out_geometry and all(row["coherence_effect_between_minus_within"] > 0 for row in held_out_geometry) else "not_supportive",
        "held_out_sequence_aligned_effects_angstrom": [row["mean_kabsch_ca_rmsd_angstrom"] for row in held_out_sequence],
        "held_out_sequence_aligned_status": "supportive" if held_out_sequence and all(row["mean_kabsch_ca_rmsd_angstrom"] > 0 for row in held_out_sequence) else "not_supportive",
        "held_out_calibration_metrics": [
            {
                "threshold": row["threshold"],
                "nmi_scope_fold": row["nmi_scope_fold"],
                "nmi_cath_topology": row["nmi_cath_topology"],
                "recurring_class_count": row["recurring_class_count"],
            }
            for row in calibration["candidate_metrics"]
            if row["split"] == "held_out"
        ],
        "gate_status": "conditional_evidence_pending_structural_similarity_and_unavailable_entry_review",
        "macroclass_decision": "do_not_assign_final_macroclass",
        "blocking_checks": {
            "complete_domain_source_coverage": assignment["rows_written"] == assignment["rows_requested"],
            "all_pdbredo_domains_have_complete_residue_coverage": assignment["complete_domain_coverage_fraction"] == 1.0,
            "sequence_aligned_structural_method_is_external_structure_proxy": True,
            "manual_audit_csv_row_fields_complete": audit.get("audit_template_pending_rows_at_recording") == 0,
        },
        "notes": [
            "PDB-REDO DSSP 4.6.1 provides the primary-source candidate assignment for available entries.",
            "Eight domains lack a PDB-REDO legacy record and remain outside primary-source calibration.",
            "Both geometry checks are proxies; TM-align or Foldseek has not been run.",
            "The user-reported audit pass is retained, while the row-level audit CSV remains pending.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(f"PDB-REDO representation gate remains {output['gate_status']}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

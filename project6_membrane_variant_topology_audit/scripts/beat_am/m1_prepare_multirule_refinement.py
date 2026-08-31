"""Prepare annotation features for standalone charge, glycosylation and disulfide
refinements.
"""
from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_gate1 import read_topology
from membrane_rule_audit.refinement import (ALPHABET, RULES, FeatureSpec,
    prepare_record, variant_features, write_json, ablation_coefficients)
from m2_standalone_refinement_scorer import build_tabular, FEATURES

OUT = ROOT/"results/standalone_multirule_20261025"
TABLE = ROOT/"results/beat_am_20261016/features/l1_table.npz"
TOPOLOGY = ROOT/"data/raw/topology_viewer.html"


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    if (OUT/"preparation_manifest.json").exists():
        raise SystemExit("Preparation is frozen; use a new dated directory for a changed design")
    census = read_topology(TOPOLOGY)
    source_dirs = [ROOT/"data/external/proline_rule_20261025/uniprot",
                   ROOT/"data/external/disulfide_rule_20261025/uniprot"]
    # Access NPZ columns by name without opening target or audited-predictor columns.
    with np.load(TABLE,allow_pickle=False) as data:
        accession = data["accession"].astype(str)
        membrane = data["is_membrane"]
        positions,natives,alternatives = data["position"],data["native"],data["alt"]
        tabular = np.column_stack([data[name] for name in FEATURES]).astype(np.float32)
    records, sources, parity, counts = {}, [], {}, Counter()
    arrays = {name:[] for name in ("accession","position","native","alt","features")}
    for acc in sorted(set(accession[membrane == 1])):
        path = next((directory/f"{acc}.json" for directory in source_dirs
                     if (directory/f"{acc}.json").exists()),None)
        if path is None:
            raise FileNotFoundError(f"No cached annotation for {acc}; no silent fallback")
        record = prepare_record(acc,census[acc],json.loads(path.read_text(encoding="utf-8")),spec=FeatureSpec())
        selected = accession == acc
        rebuilt = build_tabular(record,positions[selected],natives[selected],alternatives[selected])
        if not np.array_equal(rebuilt,tabular[selected]):
            maximum = np.max(np.abs(rebuilt-tabular[selected]),axis=0)
            raise ValueError(f"Supervised-feature parity failed for {acc}: {dict(zip(FEATURES,maximum.tolist()))}")
        parity[acc] = dict(rows=int(selected.sum()),all_18_features_exactly_match=True)
        records[acc] = record
        sources.append(dict(accession=acc,path=str(path.resolve()),sha256=sha(path)))
        position = np.repeat(np.arange(1,len(record["sequence"])+1),19)
        native = np.repeat(np.array(list(record["sequence"])),19)
        alt = np.concatenate([np.array([a for a in ALPHABET if a != n]) for n in record["sequence"]])
        feature = variant_features(record,position,native,alt)
        arrays["accession"].append(np.repeat(acc,len(position)))
        arrays["position"].append(position)
        arrays["native"].append(native)
        arrays["alt"].append(alt)
        arrays["features"].append(feature)
        counts["annotated_N_linked_sites"] += len(record["glycosites"])
        counts["annotated_disulfide_bonds"] += len(record["disulfides"])
        counts["canonical_positions"] += len(record["sequence"])
    joined = {name:np.concatenate(value) for name,value in arrays.items()}
    joined["feature_names"] = np.array(RULES)
    np.savez_compressed(OUT/"prepared_variant_features.npz",**joined)
    write_json(OUT/"annotation_bank.json",records)
    write_json(OUT/"feature_parity.json",parity)
    counts.update({f"variants_with_{rule}_feature":int(np.count_nonzero(joined["features"][:,i])) for i,rule in enumerate(RULES)})
    manifest = dict(prepared_utc=datetime.now(timezone.utc).isoformat(),
        state="implemented_and_annotation_prepared;_real_calibration_and_benchmark_deferred",
        model="standalone_DMS_trained_ProtT5_plus_18_tabular_features",
        predictor_score_inputs=[],rules=list(RULES),feature_spec=vars(FeatureSpec()),
        preparation_cohort="existing 13 fitted membrane proteins; annotation/feature integration only, never a fresh holdout",
        counts=dict(counts),variant_rows=len(joined["position"]),sources=sources,
        table_metadata_source=dict(path=str(TABLE.resolve()),sha256=sha(TABLE),
            columns_read=["accession","position","native","alt",*FEATURES],
            target_array_read=False,am_array_read=False,esm1b_array_read=False),
        artifacts={name:sha(OUT/name) for name in ("prepared_variant_features.npz","annotation_bank.json","feature_parity.json")},
        source_sha256={str(path.relative_to(ROOT)):sha(path) for path in
            [Path(__file__),ROOT/"scripts/beat_am/m2_standalone_refinement_scorer.py",ROOT/"tools/membrane_rule_audit/refinement.py",TOPOLOGY]},
        disulfide=dict(mode="control_free_annotation_feature",verified_free_controls_required=False,
            interpretation="later held-out functional prediction test; not isolated causal bond-loss estimation",
            endpoint_score_equalization=False),
        old_positive_inside_beta_reused=False,real_model_loaded=False,real_model_scored=False,
        real_coefficients_fitted=False,experimental_benchmark_run=False,
        old_model_or_paper_modified=False,
        final_ablation_names=list(ablation_coefficients({rule:0.0 for rule in RULES})),
        final_stage=["freeze common family roles and prior-exposure ledger",
            "fit baseline/preprocessing on supervised training families only, or use verified group-out-of-fold predictions",
            "fit correction coefficients on calibration families; passing rules can be locked at zero",
            "evaluate every declared ablation once on group-disjoint evaluation families",
            "report functional accuracy, all rule audits, affected regions, per-protein harm and collateral changes"],
        caveats=["directional features do not guarantee total-score rule compliance",
                 "annotation counts are not performance results",
                 "no probability recalibration or speedup claim",
                 "the old predictor-free enforcement verdict remains invalidated by its review"])
    write_json(OUT/"preparation_manifest.json",manifest)
    print(json.dumps({k:manifest[k] for k in ("state","counts","variant_rows","real_model_scored","real_coefficients_fitted","experimental_benchmark_run")},indent=2))


if __name__ == "__main__":
    main()

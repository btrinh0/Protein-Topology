"""Independently reconstruct the glycosylation repair results and collateral checks."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import m33_glycosylation_constraint_repair as repair
from m2_standalone_refinement_scorer import refine_position_scores
from membrane_rule_audit import panel


def reconstruct(matrix, sites, index, bank, groups, eligible):
    records = []
    alphabet = np.array(list(repair.audit.ALPHABET))
    occupied = {a: {g["position"] for g in r["glycosites"]} for a, r in bank.items()}
    for a, p, native, side, contrast, target in sites:
        if eligible and target and p - (1 if contrast == "A" else 2) not in occupied[a]:
            continue
        score = matrix[index[(a, p)]]
        special = "P" if contrast == "A" else ("T" if native == "S" else "S")
        others = (alphabet != native) & (alphabet != special) & np.isfinite(score)
        if not np.isfinite(score[repair.audit.ALPHABET.index(special)]) or others.sum() < 10:
            continue
        z = float(score[repair.audit.ALPHABET.index(special)] - score[others].mean())
        records.append(dict(accession=a, group=groups[a], native=native, side=side,
                            contrast=contrast, target=target, z=z))
    frame = pd.DataFrame(records)
    output = {}
    for contrast in ("A", "B"):
        for side in ("noncytosolic", "cytosolic"):
            f = frame[(frame.contrast == contrast) & (frame.side == side)]
            cells = f.groupby(["accession", "native", "group", "target"]).z.agg(["count", "mean"])
            wide = cells.unstack("target")
            key = contrast + "_" + side
            if ("count", 0) not in wide or ("count", 1) not in wide:
                output[key] = None
                continue
            wide = wide[(wide[("count", 1)] >= 1) & (wide[("count", 0)] >= 2)]
            if wide.empty:
                output[key] = None
                continue
            n1, n0 = wide[("count", 1)], wide[("count", 0)]
            weight = n1 * n0 / (n1 + n0)
            delta = wide[("mean", 1)] - wide[("mean", 0)]
            parts = pd.DataFrame(dict(num=weight * delta, den=weight)).groupby(level="group").sum().sort_index()
            num, den = parts.num.to_numpy(), parts.den.to_numpy()
            counts = np.random.default_rng(20261016).multinomial(len(parts), np.full(len(parts), 1 / len(parts)), size=1999)
            interval = np.quantile((counts @ num) / (counts @ den), [.025, .975]).tolist() if len(parts) >= 2 else None
            output[key] = dict(estimate=float(num.sum() / den.sum()), ci95=interval,
                               n_groups=len(parts), n_strata=len(wide), n_target_sites=int(n1.sum()),
                               n_reference_sites=int(n0.sum()),
                               n_proteins=wide.index.get_level_values("accession").nunique())
    return output


def main(out):
    bank, census, groups, folds, cached, index, sites = repair.metadata()
    expected = repair.read(out / "audit_results.json")["guarded_group_out"]
    with np.load(out / "guarded_position_scores.npz", allow_pickle=False) as z:
        scores = {name: z[name].copy() for name in repair.MODELS}
    maximum_error, checked = 0., 0
    for model, matrix in scores.items():
        for scope, eligible in [("original_all_motifs", False), ("feature_eligible", True)]:
            independent = reconstruct(matrix, sites, index, bank, groups, eligible)
            for key, value in independent.items():
                original = expected[model][scope][key]["group_resampling"]
                if value is None:
                    if original is not None:
                        raise AssertionError("An unsupported contrast was reported")
                    continue
                for field in ("n_groups", "n_strata", "n_target_sites", "n_reference_sites", "n_proteins"):
                    if value[field] != original[field]:
                        raise AssertionError((model, scope, key, field))
                error = max(abs(value["estimate"] - original["estimate"]),
                            float(np.max(np.abs(np.array(value["ci95"]) - original["ci95"]))))
                maximum_error = max(maximum_error, error)
                if error > 1e-12:
                    raise AssertionError("Independent reconstruction disagrees")
                checked += 1

    profile = repair.read(out / "saved_full_glycosylation_profile.json")
    exported = repair.read(out / "saved_full_glycosylation_coefficients.json")
    from_raw = cached["saved_full_fit"] + profile["beta"] * cached["features"][:, :, repair.RULE_INDEX]
    api_error = 0.
    for a, record in bank.items():
        keep = cached["accession"] == a
        scored = dict(accession=cached["accession"][keep], position=cached["position"][keep],
                      native=cached["native"][keep], alphabet=cached["alphabet"],
                      scores=cached["saved_full_fit"][keep], score_state="raw")
        applied = refine_position_scores(record, scored, exported["coefficients"])
        api_error = max(api_error, float(np.nanmax(np.abs(applied["scores"] - from_raw[keep]))))
    if api_error > 1e-12:
        raise AssertionError("The exported profile does not reproduce the repair")
    # Apply the original positive-inside estimator directly to the score matrices.
    accessions = sorted(bank)
    pidx = {a: i for i, a in enumerate(accessions)}
    maps = {a: repair.audit.refinement.topology_map(r["topology"]) for a, r in bank.items()}
    geo = [maps[str(a)].get(int(p), dict(side=-1, distance=-1, in_tm=-1))
           for a, p in zip(cached["accession"], cached["position"])]
    side = np.array([0 if g["side"] == 1 else 1 if g["side"] == 2 else -1 for g in geo])
    distance = np.array([g["distance"] for g in geo])
    flank = (side >= 0) & (distance >= 1) & (distance <= 10)
    tm = np.array([g["in_tm"] == 1 for g in geo])
    positive_inside = {}
    for name in ("raw", "existing_PI_DS", "PI_DS_constraint_glyco"):
        source = dict(scores=scores[name][flank], side=side[flank],
                      protein=np.array([pidx[str(a)] for a in cached["accession"][flank]]),
                      ref_aa=np.array([repair.audit.ALPHABET.index(str(n)) for n in cached["native"][flank]]))
        contrast = panel.gain_contrast(source)
        positive_inside[name] = panel.summarize_parts(
            panel.stratum_parts(source, contrast["penalty"], contrast["eligible"]), 1999, 20261003)
    for candidate, base in [("constraint_glyco", "raw"), ("PI_DS_constraint_glyco", "existing_PI_DS")]:
        if not np.array_equal(scores[candidate][tm], scores[base][tm], equal_nan=True):
            raise AssertionError("Charge-burial TM scores were changed")
    per_protein = pd.read_csv(out / "per_protein_contrasts.csv")
    eligible_directions = {}
    for name in ("existing_PI_glyco_DS", "PI_DS_constraint_glyco"):
        subset = per_protein[(per_protein.model == name) & (per_protein.scope == "feature_eligible")]
        signed = np.where(subset.contrast == "A", 1., -1.) * subset.baseline_delta
        eligible_directions[name] = dict(
            comparisons=len(subset), expected_or_zero=int((signed >= -1e-12).sum()),
            opposite=int((signed < -1e-12).sum()), minimum_signed_delta=float(signed.min()))
    # Match the recorded audit implementation, including its documented accuracy-stage optimization.
    plan = repair.read(out / "PLAN.json")
    snapshots = []
    for name, digest in plan["input_sha256"].items():
        current = repair.ROOT / name
        if repair.audit.sha(current) == digest:
            continue
        snapshot = out / "audit_runner_snapshot.py"
        if current.name == "m33_glycosylation_constraint_repair.py" and snapshot.exists() and repair.audit.sha(snapshot) == digest:
            snapshots.append(str(snapshot.relative_to(repair.ROOT)))
        else:
            raise AssertionError(f"Changed source input: {name}")
    result = dict(
        independently_reconstructed_audit_rows=checked,
        maximum_absolute_reconstruction_error=maximum_error,
        deployment_API_maximum_error=api_error,
        all_TM_scores_exactly_unchanged=True, positive_inside=positive_inside,
        feature_eligible_protein_directions=eligible_directions,
        input_hashes_verified=True, preserved_implementation_snapshots=snapshots,
        new_model_or_PCA_fits=0, new_experimental_outcome_reads=False,
    )
    (out / "INDEPENDENT_VERIFICATION.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    main(args.out.resolve())

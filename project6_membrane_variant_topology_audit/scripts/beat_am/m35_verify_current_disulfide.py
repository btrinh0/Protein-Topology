"""Verify disulfide-loss contrasts in the current saved score stack."""
from __future__ import annotations

from collections import defaultdict
from contextlib import ExitStack, redirect_stdout
import csv
from datetime import datetime, timezone
import io
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np

import m7_audit_standalone_rules as audit

ROOT = audit.ROOT
SOURCE = ROOT / "results/glycosylation_constraint_repair_20261031/run2/guarded_position_scores.npz"
OUT = ROOT / "results/disulfide_refinement_verification_20261031"
MODELS = ("raw", "existing_PI_DS", "existing_PI_glyco_DS", "PI_DS_constraint_glyco")


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    bank = read(audit.BASE / "annotation_bank.json")
    plan = read(audit.BENCH / "BENCHMARK_PLAN.json")
    groups = {a: g for g, members in plan["dependence_members"].items() for a in members}
    census_all = audit.read_topology(ROOT / "data/raw/topology_viewer.html")
    census = {a: census_all[a] for a in sorted(bank)}
    paths = [Path(__file__), SOURCE, audit.BASE / "annotation_bank.json",
             audit.BENCH / "BENCHMARK_PLAN.json", audit.b10.UNIPROT_TSV,
             Path(audit.__file__), Path(audit.b10.__file__),
             ROOT / "data/raw/topology_viewer.html"]
    before_hashes = {str(p.relative_to(ROOT)): audit.sha(p) for p in paths}
    with np.load(SOURCE, allow_pickle=False) as z:
        acc, pos = z["accession"].copy(), z["position"].copy()
        native = z["native"].astype(str)
        if "".join(z["alphabet"]) != audit.ALPHABET:
            raise AssertionError("Unexpected score alphabet")
        matrices = {name: z[name].copy() for name in MODELS}
    index = {(str(a), int(p)): i for i, (a, p) in enumerate(zip(acc, pos))}
    with redirect_stdout(io.StringIO()):
        uni = {a: r for a, r in audit.b10.load_uniprot().items() if a in bank}
        sites = audit.b10.enumerate_sites(census, uni)
    targets = {(a, b[k]) for a, r in bank.items() for b in r["disulfides"] for k in ("start", "end")}
    target = np.array([(str(a), int(p)) in targets for a, p in zip(sites["acc"], sites["pos"])])
    keep = target | (sites["bonded"] == 0)
    strata = audit.b10.stratum_ids(sites, ("comp",))
    values = {name: audit.score_rows(matrix, index, sites["acc"], sites["pos"])
              for name, matrix in matrices.items()}
    results, protein_rows = {}, []
    with ExitStack() as stack:
        for obj, method in ((audit.PCA, "fit"), (audit.PCA, "fit_transform"),
                            (audit.HistGradientBoostingRegressor, "fit"),
                            (audit.refinement, "fit_coefficients")):
            stack.enter_context(patch.object(obj, method, side_effect=audit.forbidden_fit))
        for name, matrix in values.items():
            results[name] = {}
            for replacement in ("S", "A", "all19_mean"):
                v = np.nanmean(matrix, axis=1) if replacement == "all19_mean" else matrix[:, audit.ALPHABET.index(replacement)]
                whole = audit.b10.matched_gap(
                    sites["prot"][keep], strata[keep], target[keep].astype(int), v[keep],
                    len(sites["prot_names"]), np.random.default_rng(audit.b10.SEED))
                lo, hi = whole["ci95"]
                results[name][replacement] = dict(
                    whole_protein=whole,
                    status="pass_internal_loss_criterion" if lo > 0 else
                           "opposite_direction" if hi < 0 else "inconclusive")
                if replacement != "S":
                    continue
                grouped_rows, buckets = [], defaultdict(lambda: [[], []])
                for i in np.flatnonzero(keep & np.isfinite(v)):
                    a, comp = str(sites["acc"][i]), int(sites["comp"][i])
                    grouped_rows.append((a, "C/" + str(comp), bool(target[i]), float(v[i])))
                    buckets[(a, comp)][int(target[i])].append(float(v[i]))
                results[name][replacement]["dependence_group_sensitivity"] = audit.matched_cluster(grouped_rows, groups)
                parts = defaultdict(lambda: [0., 0., 0, 0])
                for (a, comp), (reference, endpoint) in buckets.items():
                    if len(endpoint) < 1 or len(reference) < audit.b10.MIN_REF_ARM:
                        continue
                    w = len(endpoint) * len(reference) / (len(endpoint) + len(reference))
                    parts[a][0] += w * (np.mean(endpoint) - np.mean(reference))
                    parts[a][1] += w
                    parts[a][2] += len(endpoint)
                    parts[a][3] += len(reference)
                for a, (numerator, denominator, n1, n0) in sorted(parts.items()):
                    protein_rows.append(dict(model=name, accession=a, gene=census[a]["gene"],
                        delta=float(numerator / denominator), n_endpoint_sites=n1,
                        n_background_sites=n0, weight=float(denominator),
                        expected_direction=bool(numerator > 0)))
                reconstructed = sum(p[0] for p in parts.values()) / sum(p[1] for p in parts.values())
                if abs(reconstructed - whole["estimate"]) > 1e-12:
                    raise AssertionError("Per-protein reconstruction disagrees with the original instrument")
    # Reproduce the archived result before comparing the current score stack.
    old = read(audit.BASE / "disulfide_internal_consistency_20261025/internal_consistency.json")
    historical = old["results"]["guarded_oof_PI_plus_disulfide"]["S"]["measurement"]
    measured = results["existing_PI_DS"]["S"]["whole_protein"]
    if not np.isclose(measured["estimate"], historical["estimate"], atol=1e-12, rtol=0):
        raise AssertionError("The old disulfide result was not reproduced")
    if not np.allclose(measured["ci95"], historical["ci95"], atol=1e-12, rtol=0):
        raise AssertionError("The original uncertainty protocol changed")
    si = audit.ALPHABET.index("S")
    cysteine = native == "C"
    equal = np.array_equal(matrices["existing_PI_DS"][cysteine, si],
                           matrices["PI_DS_constraint_glyco"][cysteine, si], equal_nan=True)
    if not equal:
        raise AssertionError("The glycosylation repair changed C-to-S disulfide scores")
    current_rows = [r for r in protein_rows if r["model"] == "PI_DS_constraint_glyco"]
    final_pass = results["PI_DS_constraint_glyco"]["S"]["status"] == "pass_internal_loss_criterion"
    summary = dict(
        completed_utc=datetime.now(timezone.utc).isoformat(),
        criterion="existing annotation-defined matched endpoint-minus-background C-to-S gap; original whole-protein interval must be positive",
        evaluation_kind="internal_annotation_consistency",
        same_annotations_select_correction_and_evaluation_targets=True,
        background_is_verified_free_cysteine=False,
        current_model_passes=final_pass,
        current_proteins_positive=sum(r["expected_direction"] for r in current_rows),
        current_contributing_proteins=len(current_rows),
        no_new_refinement_needed=final_pass,
        C_to_S_scores_unchanged_after_glycosylation_repair=equal,
        new_model_or_coefficient_fits=0, experimental_outcome_reads=0,
        results=results, source_sha256=before_hashes)
    if not all(audit.sha(ROOT / p) == digest for p, digest in before_hashes.items()):
        raise AssertionError("An input changed during the verification")
    (OUT / "verification.json").write_text(json.dumps(audit.clean(summary), indent=2) + "\n", encoding="utf-8")
    with (OUT / "per_protein.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(protein_rows[0]))
        writer.writeheader()
        writer.writerows(protein_rows)
    print(json.dumps(dict(pass_internal_loss_criterion=final_pass,
        current_measurement=results["PI_DS_constraint_glyco"]["S"],
        current_per_protein=current_rows, no_new_refinement_needed=final_pass), indent=2))


if __name__ == "__main__":
    main()

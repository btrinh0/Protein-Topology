"""Phase C blind-spot test: lab minus predictor on identical variants.

Required by the frozen preregistration section 7. The headline claim is about a
*difference between two measurements*, so it is tested as one: the same within-site
estimator is run on the lab scores and on each predictor, restricted to exactly the
same (protein, position, substitution) set, and the two are resampled in the SAME
bootstrap draw so their correlation is preserved.

Both sides are standardised by their own within-site penalty spread before
differencing, because the lab scores, AlphaMissense and ESM1b are on three
different scales.

ESM1b is the positive control, with a prediction of no difference from the lab.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from run_gate0 import load_frame
from run_phase_c import (CHARGED, NEGATIVE, POSITIVE, PRIMARY_WINDOW, SCORE_SETS,
                         WEBB, assemble, site_penalties)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "results" / "phase_c_20260928"
ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


def predictor_sites(npz: Path, wanted: dict[str, set[tuple[int, str]]],
                    accession_of: dict[str, str], kind: str, group: set[str]) -> list[dict]:
    """Within-site penalties from a predictor, restricted to the lab's own sites."""
    frame, _ = load_frame(npz)
    blob = np.load(npz, allow_pickle=False)
    accessions = list(blob["accessions"])
    index_of = {a: i for i, a in enumerate(accessions)}
    gene_of = {index_of[a]: g for g, a in accession_of.items() if a in index_of}

    keep = frame[frame["protein"].isin(gene_of)]
    keep = keep[keep["distance"] <= PRIMARY_WINDOW]
    alphabet = np.array(list(ALPHABET))
    refs = alphabet[keep["ref_aa"].to_numpy()]
    alts = alphabet[keep["alt_aa"].to_numpy()]

    grouped: dict = defaultdict(lambda: defaultdict(dict))
    for protein, position, ref, alt, score, side, distance in zip(
            keep["protein"].to_numpy(), keep["position"].to_numpy(), refs, alts,
            keep["score"].to_numpy(), keep["side"].to_numpy(), keep["distance"].to_numpy()):
        gene = gene_of[protein]
        if (int(position), str(ref)) not in wanted.get(gene, ()):
            continue
        grouped[(gene, int(side), int(distance))][(int(position), str(ref))][str(alt)] = float(score)

    rows = []
    for (gene, side, distance), variants in grouped.items():
        for site in site_penalties(variants, kind, group):
            rows.append({"gene": gene, "ref": site["ref"], "position": site["position"],
                         "penalty": site["penalty"], "cytosolic": side == 0,
                         "distance": distance})
    return rows


def stratum_stats(rows: list[dict]):
    strata: dict = defaultdict(lambda: {"cyto": [], "outer": []})
    for row in rows:
        strata[(row["gene"], row["ref"])]["cyto" if row["cytosolic"] else "outer"].append(row["penalty"])
    keys, weight, diff, gene = [], [], [], []
    for (g, ref), arms in sorted(strata.items()):
        n1, n0 = len(arms["cyto"]), len(arms["outer"])
        if n1 == 0 or n0 == 0:
            continue
        keys.append((g, ref))
        weight.append(n1 * n0 / (n1 + n0))
        diff.append(float(np.mean(arms["cyto"]) - np.mean(arms["outer"])))
        gene.append(g)
    return keys, np.array(weight), np.array(diff), gene


def paired_contrast(lab: list[dict], predictor: list[dict], draws: int,
                    weights_kind: str, seed: int) -> dict:
    """Standardised lab minus predictor, resampled in one shared bootstrap draw."""
    lab_keys, lab_w, lab_d, lab_genes = stratum_stats(lab)
    pred_keys, pred_w, pred_d, _ = stratum_stats(predictor)
    shared = sorted(set(lab_keys) & set(pred_keys))
    if len(shared) < 2:
        return {"n_shared_strata": len(shared)}
    lab_index = {k: i for i, k in enumerate(lab_keys)}
    pred_index = {k: i for i, k in enumerate(pred_keys)}
    li = np.array([lab_index[k] for k in shared])
    pi = np.array([pred_index[k] for k in shared])
    weight = lab_w[li]
    lab_diff, pred_diff = lab_d[li], pred_d[pi]
    genes = sorted({k[0] for k in shared})
    gene_index = {g: i for i, g in enumerate(genes)}
    cluster = np.array([gene_index[k[0]] for k in shared])

    lab_sd = float(np.std([r["penalty"] for r in lab], ddof=1)) or 1.0
    pred_sd = float(np.std([r["penalty"] for r in predictor], ddof=1)) or 1.0

    def contrast(lab_values, pred_values):
        total = weight.sum()
        a = float((weight * lab_values).sum() / total) / lab_sd
        b = float((weight * pred_values).sum() / total) / pred_sd
        return a, b, a - b

    lab_point, pred_point, delta = contrast(lab_diff, pred_diff)
    rng = np.random.default_rng(seed)
    samples = np.empty(draws)
    lab_centred, pred_centred = lab_diff - lab_diff.mean(), pred_diff - pred_diff.mean()
    for i in range(draws):
        flips = (rng.choice((-1.0, 1.0), size=len(genes)) if weights_kind == "rademacher"
                 else rng.choice(WEBB, size=len(genes)))
        multiplier = flips[cluster]
        _, _, samples[i] = contrast(lab_centred * multiplier, pred_centred * multiplier)
    extreme = int((np.abs(samples) >= abs(delta)).sum())
    return {
        "n_shared_strata": len(shared), "n_genes": len(genes),
        "lab_standardised": lab_point, "predictor_standardised": pred_point,
        "difference_lab_minus_predictor": delta,
        "p_value": (extreme + 1) / (draws + 1), "weights": weights_kind,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", default=str(ROOT / "data" / "raw" / "topology_viewer.html"))
    parser.add_argument("--intramem", default=str(ROOT / "data" / "external" / "uniprot" / "human_intramem.tsv"))
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    parser.add_argument("--wild-draws", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260928)
    args = parser.parse_args()

    sites, meta = assemble(args)
    lab = [s for s in sites
           if s["class"] == "positive_charge_gain" and s["distance"] <= PRIMARY_WINDOW]
    wanted: dict = defaultdict(set)
    for site in lab:
        wanted[site["gene"]].add((site["position"], site["ref"]))
    accession_of = {gene: accession for gene, (_, accession, _) in SCORE_SETS.items()
                    if meta.get(gene, {}).get("status") == "included"}

    out: dict = {}
    for name, npz in (("alphamissense", ROOT / "data/processed/gate0/am_full_flank_scores.npz"),
                      ("esm1b", ROOT / "data/processed/gate0/esm1b_flank_scores.npz")):
        predictor = predictor_sites(npz, wanted, accession_of, "gain", POSITIVE)
        out[name] = {
            "n_predictor_sites": len(predictor),
            "webb": paired_contrast(lab, predictor, args.wild_draws, "webb", args.seed),
            "rademacher": paired_contrast(lab, predictor, args.wild_draws, "rademacher", args.seed),
        }

    summary = {
        "analysis": "Phase C blind-spot test: lab minus predictor on identical variants",
        "preregistration_sha256": "754a9a713a2afc7fd2e3e43be9af55bb762805feeb3d3991be0667c207635d79",
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "n_lab_sites": len(lab),
        "decision_rule": "Blind spot confirmed only if the lab side difference is negative "
                         "with an interval excluding zero AND this difference excludes zero.",
        "results": out,
    }
    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    (Path(args.outdir) / "phase_c_blind_spot.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

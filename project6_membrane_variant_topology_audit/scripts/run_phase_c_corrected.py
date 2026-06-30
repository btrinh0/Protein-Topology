"""Phase C, rerun under the frozen correction plan.

Correction plan: docs/PHASE_C_CORRECTION_PLAN.md
SHA-256 c74e03ec741bcd08c4062e6a282db972e145d486b4e54540c4295d0a331dca41
frozen before any corrected number was computed.

Protocol unchanged: docs/PHASE_C_PREREGISTRATION.md, SHA-256 754a9a71...
This repairs the engine, not the protocol.

F1 matched substitutions · F2 weights through inference · F3 shared estimator
F4 interval by test inversion · F5 honest contributing counts
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from extract_am_flank_scores import build_flank_map
from phase_c_engine import build_strata, fit, statistic, wild_interval, wild_p_value
from run_gate0 import load_frame
from run_gate1 import read_topology
from run_phase_c import (CHARGED, CROSS_SPECIES, GPCR_FAMILY, NEGATIVE, POSITIVE,
                         PRIMARY_WINDOW, SCORE_SETS, classify, direction_for,
                         score_rows, site_penalties)
from run_within_site_controls import load_pore_spans, sequon_positions

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "results" / "phase_c_corrected_20260929"
ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
PORE_MARGIN = 10


def assemble_variants(args) -> tuple[dict, dict, dict]:
    """Lab variants kept at substitution level, so F1 can intersect on the alternate too."""
    proteins = read_topology(Path(args.topology))
    pore_spans = load_pore_spans(Path(args.intramem))
    lab: dict = {}
    context: dict = {}
    meta: dict = {}

    for gene, (urn, accession, offset) in SCORE_SETS.items():
        rows = score_rows(urn)
        if not rows:
            meta[gene] = {"status": "score file not cached"}
            continue
        synonymous, nonsense, missense = classify(rows)
        direction = direction_for(gene, synonymous, nonsense, missense)
        meta[gene] = {"urn": urn, "uniprot": accession, **direction}
        if direction["flip"] is None:
            meta[gene]["status"] = "EXCLUDED: no usable anchors"
            continue
        if abs(direction["separation"] - 0.5) < 0.05:
            meta[gene]["status"] = "EXCLUDED: anchors do not separate"
            continue
        meta[gene]["status"] = "included"

        record = proteins[accession]
        sequence = record["sequence"]
        flanks = build_flank_map(record, 60)
        sequons = sequon_positions(sequence)
        blocked: set[int] = set()
        for start, end in pore_spans.get(accession, ()):
            blocked |= set(range(start - PORE_MARGIN, end + PORE_MARGIN + 1))

        scores = [s for _, _, _, s in missense]
        centre, spread = float(np.mean(scores)), float(np.std(scores, ddof=1))
        for ref, position, alt, score in missense:
            census_position = position + offset
            if census_position < 1 or census_position > len(sequence):
                continue
            if sequence[census_position - 1] != ref:
                continue
            flank = flanks.get(census_position)
            if flank is None:
                continue
            value = (score - centre) / spread if spread > 0 else 0.0
            if direction["flip"]:
                value = -value
            key = (gene, census_position, ref)
            lab.setdefault(key, {})[alt] = value
            side, distance, tmd_index = flank
            context[key] = {
                "gene": gene, "uniprot": accession, "cytosolic": side == 0,
                "distance": distance, "tmd_index": tmd_index,
                "in_sequon": census_position in sequons,
                "near_pore": census_position in blocked,
                "ion_channel": bool(record.get("isIonChannel")),
                "gpcr": gene in GPCR_FAMILY, "cross_species": gene in CROSS_SPECIES,
            }
    return lab, context, meta


def penalties_from(variants: dict, context: dict, kind: str, group: set[str],
                   label: str) -> list[dict]:
    rows = []
    by_site = defaultdict(dict)
    for (gene, position, ref), alts in variants.items():
        by_site[(gene, position, ref)] = alts
    for (gene, position, ref), alts in by_site.items():
        computed = site_penalties({(position, ref): alts}, kind, group)
        if not computed:
            continue
        info = context[(gene, position, ref)]
        rows.append({**info, "class": label, "position": position, "ref": ref,
                     "penalty": computed[0]["penalty"]})
    return rows


def predictor_variants(npz: Path, accession_of: dict[str, str]) -> dict:
    frame, _ = load_frame(npz)
    blob = np.load(npz, allow_pickle=False)
    index_of = {a: i for i, a in enumerate(list(blob["accessions"]))}
    gene_of = {index_of[a]: g for g, a in accession_of.items() if a in index_of}
    keep = frame[frame["protein"].isin(gene_of)]
    alphabet = np.array(list(ALPHABET))
    out: dict = defaultdict(dict)
    for protein, position, ref_code, alt_code, score in zip(
            keep["protein"].to_numpy(), keep["position"].to_numpy(),
            keep["ref_aa"].to_numpy(), keep["alt_aa"].to_numpy(), keep["score"].to_numpy()):
        out[(gene_of[protein], int(position), str(alphabet[ref_code]))][
            str(alphabet[alt_code])] = float(score)
    return out


def intersect(lab: dict, predictor: dict) -> tuple[dict, dict, dict]:
    """F1: keep only substitutions measured by both, then rebuild both penalties from them."""
    lab_out, pred_out = {}, {}
    counts = {"lab_substitutions": 0, "predictor_substitutions": 0, "shared_substitutions": 0,
              "sites_before": len(lab), "sites_after": 0}
    for key, lab_alts in lab.items():
        counts["lab_substitutions"] += len(lab_alts)
        pred_alts = predictor.get(key)
        if not pred_alts:
            continue
        shared = set(lab_alts) & set(pred_alts)
        if not shared:
            continue
        counts["shared_substitutions"] += len(shared)
        lab_out[key] = {a: lab_alts[a] for a in shared}
        pred_out[key] = {a: pred_alts[a] for a in shared}
    counts["predictor_substitutions"] = sum(len(v) for v in predictor.values())
    counts["sites_after"] = len(lab_out)
    return lab_out, pred_out, counts


def standardise(rows: list[dict]) -> list[dict]:
    spread = float(np.std([r["penalty"] for r in rows], ddof=1)) if len(rows) > 1 else 1.0
    spread = spread or 1.0
    return [{**r, "penalty": r["penalty"] / spread} for r in rows]


def paired_blind_spot(lab_rows, pred_rows, draws, kind, seed) -> dict:
    """Lab minus predictor on identical substitutions, resampled in one shared draw."""
    lab_strata = build_strata(standardise(lab_rows))
    pred_strata = build_strata(standardise(pred_rows))
    if lab_strata["n_strata"] < 2 or lab_strata["n_strata"] != pred_strata["n_strata"]:
        return {"note": "strata sets differ after matching",
                "lab_strata": lab_strata["n_strata"], "pred_strata": pred_strata["n_strata"]}
    delta = {**lab_strata, "diff": lab_strata["diff"] - pred_strata["diff"]}
    beta, _ = statistic(delta)
    return {
        "lab_standardised": statistic(lab_strata)[0],
        "predictor_standardised": statistic(pred_strata)[0],
        "difference_lab_minus_predictor": beta,
        "p_value": wild_p_value(delta, 0.0, draws, kind, seed),
        "ci95_inverted": wild_interval(delta, max(draws // 10, 199), kind, seed),
        "n_strata": delta["n_strata"], "n_genes": len(delta["genes"]), "weights": kind,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", default=str(ROOT / "data/raw/topology_viewer.html"))
    parser.add_argument("--intramem", default=str(ROOT / "data/external/uniprot/human_intramem.tsv"))
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    parser.add_argument("--wild-draws", type=int, default=9999)
    parser.add_argument("--boot-draws", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260929)
    args = parser.parse_args()

    lab, context, meta = assemble_variants(args)
    included = {g for g, m in meta.items() if m.get("status") == "included"}
    accession_of = {g: SCORE_SETS[g][1] for g in included}

    def rows_for(cls: str, kind: str, group: set[str], window: int = PRIMARY_WINDOW,
                 source: dict | None = None, ctx: dict | None = None) -> list[dict]:
        rows = penalties_from(source if source is not None else lab,
                              ctx if ctx is not None else context, kind, group, cls)
        return [r for r in rows if r["distance"] <= window]

    gain = rows_for("positive_charge_gain", "gain", POSITIVE)
    results = {
        "primary": fit(gain, "K/R gain, window <=10", args.wild_draws, args.boot_draws, args.seed),
    }
    results["sensitivities"] = {
        "window_5": fit(rows_for("positive_charge_gain", "gain", POSITIVE, 5),
                        "window <=5", args.wild_draws, args.boot_draws, args.seed),
        "window_15": fit(rows_for("positive_charge_gain", "gain", POSITIVE, 15),
                         "window <=15", args.wild_draws, args.boot_draws, args.seed),
        "no_cross_species": fit([r for r in gain if not r["cross_species"]],
                                "without cross-species", args.wild_draws, args.boot_draws, args.seed),
        "no_ion_channels": fit([r for r in gain if not r["ion_channel"]],
                               "without ion channels", args.wild_draws, args.boot_draws, args.seed),
        "pore_excluded": fit([r for r in gain if not r["near_pore"]],
                             "pore regions excluded", args.wild_draws, args.boot_draws, args.seed),
        "no_sequon": fit([r for r in gain if not r["in_sequon"]],
                         "sequons excluded", args.wild_draws, args.boot_draws, args.seed),
        "gene_equal_weight": fit(gain, "gene-equal weighting", args.wild_draws,
                                 args.boot_draws, args.seed, equal_gene_weight=True),
    }
    results["robustness"] = {
        "leave_one_family_out_no_gpcr": fit([r for r in gain if not r["gpcr"]],
                                            "GPCR family removed", args.wild_draws,
                                            args.boot_draws, args.seed),
        "leave_one_protein_out": {
            g: fit([r for r in gain if r["gene"] != g], f"without {g}",
                   args.wild_draws, args.boot_draws, args.seed)["estimate"]
            for g in sorted({r["gene"] for r in gain})
        },
    }
    results["secondary_classes"] = {
        "positive_charge_loss": fit(rows_for("positive_charge_loss", "loss", POSITIVE),
                                    "K/R loss", args.wild_draws, args.boot_draws, args.seed),
        "negative_charge_gain": fit(rows_for("negative_charge_gain", "gain", NEGATIVE),
                                    "D/E gain", args.wild_draws, args.boot_draws, args.seed),
        "negative_charge_loss": fit(rows_for("negative_charge_loss", "loss", NEGATIVE),
                                    "D/E loss", args.wild_draws, args.boot_draws, args.seed),
    }

    blind: dict = {}
    for name, npz in (("alphamissense", ROOT / "data/processed/gate0/am_full_flank_scores.npz"),
                      ("esm1b", ROOT / "data/processed/gate0/esm1b_flank_scores.npz")):
        predictor = predictor_variants(npz, accession_of)
        lab_matched, pred_matched, counts = intersect(lab, predictor)
        lab_rows = rows_for("positive_charge_gain", "gain", POSITIVE,
                            source=lab_matched, ctx=context)
        pred_rows = rows_for("positive_charge_gain", "gain", POSITIVE,
                             source=pred_matched, ctx=context)
        blind[name] = {"matching": counts,
                       "n_lab_sites": len(lab_rows), "n_predictor_sites": len(pred_rows)}
        for kind in ("webb", "rademacher"):
            blind[name][kind] = paired_blind_spot(lab_rows, pred_rows,
                                                  args.wild_draws, kind, args.seed)
    results["blind_spot_matched_substitutions"] = blind

    summary = {
        "analysis": "Phase C corrected under the frozen correction plan",
        "correction_plan_sha256": "c74e03ec741bcd08c4062e6a282db972e145d486b4e54540c4295d0a331dca41",
        "preregistration_sha256": "754a9a713a2afc7fd2e3e43be9af55bb762805feeb3d3991be0667c207635d79",
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "fixes_applied": ["F1 matched substitutions", "F2 weights through inference",
                          "F3 shared estimator", "F4 inverted-test interval",
                          "F5 contributing counts"],
        "score_sets": meta,
        "results": results,
    }
    out_dir = Path(args.outdir); out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "phase_c_corrected.json").write_text(json.dumps(summary, indent=2) + "\n",
                                                    encoding="utf-8")

    p = results["primary"]
    print("=== CORRECTED PRIMARY: K/R gain, window <=10 ===")
    print(f"  estimate {p['estimate']:+.4f}   t={p['t_statistic']:+.3f}")
    print(f"  contributing sites {p['n_contributing_sites']} "
          f"({p['n_cytosolic_sites']}/{p['n_outer_sites']})  dropped single-sided "
          f"{p['n_dropped_single_sided_sites']}  strata {p['n_strata']}  genes {p['n_genes']}")
    for kind in ("webb", "rademacher"):
        w = p[f"wild_{kind}"]
        print(f"  wild {kind:<11} p={w['p_value']:.4f}  inverted CI95 "
              f"[{w['ci95_inverted'][0]:+.4f}, {w['ci95_inverted'][1]:+.4f}]")
    print(f"  percentile CI95 (superseded) {p['percentile_ci95_superseded']}")
    print("\n=== sensitivities ===")
    for name, v in results["sensitivities"].items():
        if v.get("n_strata", 0) >= 2:
            print(f"  {name:<20} {v['estimate']:+.4f}  webb p={v['wild_webb']['p_value']:.4f}  "
                  f"n={v['n_contributing_sites']}  genes={v['n_genes']}")
    lofo = results["robustness"]["leave_one_family_out_no_gpcr"]
    if lofo.get("n_strata", 0) >= 2:
        print(f"\n=== leave-one-family-out ===\n  no GPCR {lofo['estimate']:+.4f}  "
              f"webb p={lofo['wild_webb']['p_value']:.4f}  genes={lofo['n_genes']}  "
              f"n={lofo['n_contributing_sites']}")
    print("\n=== blind-spot test, identical substitutions ===")
    for name, block in blind.items():
        m = block["matching"]
        print(f"  {name}: shared substitutions {m['shared_substitutions']} "
              f"(lab had {m['lab_substitutions']}), sites {m['sites_before']} -> {m['sites_after']}")
        w = block.get("webb", {})
        if "difference_lab_minus_predictor" in w:
            print(f"     lab {w['lab_standardised']:+.4f}  predictor {w['predictor_standardised']:+.4f}"
                  f"  diff {w['difference_lab_minus_predictor']:+.4f}  p={w['p_value']:.4f}"
                  f"  CI [{w['ci95_inverted'][0]:+.3f}, {w['ci95_inverted'][1]:+.3f}]")
    print("\n=== secondary classes ===")
    for cls, v in results["secondary_classes"].items():
        if v.get("n_strata", 0) >= 2:
            print(f"  {cls:<24} {v['estimate']:+.4f}  webb p={v['wild_webb']['p_value']:.4f}  "
                  f"n={v['n_contributing_sites']}")


if __name__ == "__main__":
    main()

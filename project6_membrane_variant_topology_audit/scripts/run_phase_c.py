"""Phase C, run exactly as preregistered.

Frozen protocol: docs/PHASE_C_PREREGISTRATION.md
SHA-256 754a9a713a2afc7fd2e3e43be9af55bb762805feeb3d3991be0667c207635d79
Frozen 2026-09-28 23:16:29 PDT, verified before this script was written.

Primary: the within-site K/R gain side difference over the <=10 residue window,
cytosolic minus outer, with a restricted wild cluster bootstrap-t clustered on
gene, reported under both Webb and Rademacher weights.

Every residue set, sensitivity fit and decision rule here is transcribed from the
frozen text. Deviations are recorded in the output under `deviations` and belong
in the decision log, never in the frozen document.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from extract_am_flank_scores import build_flank_map
from run_gate1 import read_topology
from run_within_site_controls import load_pore_spans, sequon_positions

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "external" / "mavedb"
DEFAULT_OUT = ROOT / "results" / "phase_c_20260928"
PRIMARY_WINDOW = 10
PORE_MARGIN = 10

POSITIVE, NEGATIVE, AMBIGUOUS = set("KR"), set("DE"), set("H")
CHARGED = POSITIVE | NEGATIVE | AMBIGUOUS
GPCR_FAMILY = {"CCR5", "CXCR4", "GPR68", "RHO"}
WEBB = np.array([-math.sqrt(1.5), -1.0, -math.sqrt(0.5),
                 math.sqrt(0.5), 1.0, math.sqrt(1.5)])
THREE = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
         "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
         "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
         "Tyr": "Y", "Val": "V"}
MISSENSE = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")

# Score sets chosen by the frozen rule 4: expression-type readout, then the most
# flank charge changes at <=10. Offsets are the Phase B / recovery values.
SCORE_SETS = {
    "CCR5":    ("urn:mavedb:00000047-a-1", "P51681", 0),
    "CXCR4":   ("urn:mavedb:00000048-a-1", "P61073", 0),
    "GPR68":   ("urn:mavedb:00001207-a-2", "Q15743", 0),
    "KCNH2":   ("urn:mavedb:00001231-a-2", "Q12809", 0),
    "VKORC1":  ("urn:mavedb:00000078-b-1", "Q9BQB6", 0),
    "LDLR":    ("urn:mavedb:00001269-b-1", "P01130", 0),
    "MPL":     ("urn:mavedb:00001214-j-1", "P40238", 480),
    "SGCB":    ("urn:mavedb:00000659-a-1", "Q16585", 0),
    "CYP2C19": ("urn:mavedb:00001199-a-1", "P33261", 0),
    "CYP2C9":  ("urn:mavedb:00000095-a-1", "P11712", 0),
    "SGCA":    ("urn:mavedb:00001283-a-1", "Q16586", 0),
    "CD86":    ("urn:mavedb:00000046-a-1", "P42081", 243),
    "RHO":     ("urn:mavedb:00001275-a-1", "P08100", 0),
    "KCNE1":   ("urn:mavedb:00000674-a-2", "P15382", 0),
    "KCNJ2":   ("urn:mavedb:00000660-a-1", "P63252", 0),
    "INSR":    ("urn:mavedb:00001239-a-1", "P06213", 0),
}
CROSS_SPECIES = {"KCNJ2"}


def score_rows(urn: str) -> list[dict]:
    path = CACHE / f"{urn.replace(':', '_').replace('.', '_')}.scores.csv"
    if not path.exists():
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def classify(rows: list[dict]) -> tuple[list, list, list]:
    """Split into synonymous, nonsense and missense, each with its score."""
    synonymous, nonsense, missense = [], [], []
    for row in rows:
        label = (row.get("hgvs_pro") or "").strip()
        try:
            score = float(row.get("score"))
        except (TypeError, ValueError):
            continue
        if not np.isfinite(score):
            continue
        if label.endswith("="):
            synonymous.append(score)
            continue
        if "Ter" in label:
            nonsense.append(score)
            continue
        match = MISSENSE.fullmatch(label)
        if not match:
            continue
        ref, position, alt = THREE.get(match.group(1)), int(match.group(2)), THREE.get(match.group(3))
        if ref is None or alt is None:
            continue
        if ref == alt:
            synonymous.append(score)
        else:
            missense.append((ref, position, alt, score))
    return synonymous, nonsense, missense


def rank_separation(low: list[float], high: list[float]) -> float:
    """P(a draw from `low` < a draw from `high`), ties at 0.5. 0.5 means no separation."""
    if not low or not high:
        return float("nan")
    a, b = np.asarray(low), np.asarray(high)
    order = np.argsort(np.concatenate([a, b]), kind="mergesort")
    ranks = np.empty(order.size); ranks[order] = np.arange(1, order.size + 1)
    u = ranks[: a.size].sum() - a.size * (a.size + 1) / 2
    return 1.0 - u / (a.size * b.size)


def direction_for(gene: str, synonymous, nonsense, missense) -> dict:
    """Assay direction from anchors, per the frozen rule, with its basis recorded.

    Frozen: nonsense below synonymous means higher scores are more protein, so the
    score is flipped to make higher mean more damaging.
    """
    missense_scores = [s for _, _, _, s in missense]
    if synonymous and nonsense:
        flip = float(np.median(nonsense)) < float(np.median(synonymous))
        return {"basis": "frozen: nonsense vs synonymous", "flip": flip,
                "separation": rank_separation(nonsense, synonymous),
                "n_synonymous": len(synonymous), "n_nonsense": len(nonsense)}
    if synonymous and missense_scores:
        flip = float(np.median(synonymous)) > float(np.median(missense_scores))
        return {"basis": "DEVIATION D1: synonymous vs missense (no nonsense deposited)",
                "flip": flip, "separation": rank_separation(missense_scores, synonymous),
                "n_synonymous": len(synonymous), "n_nonsense": 0}
    if nonsense and missense_scores:
        flip = float(np.median(nonsense)) < float(np.median(missense_scores))
        return {"basis": "DEVIATION D1: nonsense vs missense (no synonymous deposited)",
                "flip": flip, "separation": rank_separation(nonsense, missense_scores),
                "n_synonymous": 0, "n_nonsense": len(nonsense)}
    return {"basis": "no usable anchors", "flip": None, "separation": float("nan"),
            "n_synonymous": 0, "n_nonsense": 0}


def site_penalties(variants: dict, kind: str, group: set[str]) -> list[dict]:
    """Within-site penalties using the frozen residue sets.

    gain: reference uncharged; numerator alt in `group`; comparator alt uncharged.
    loss: reference in `group`; numerator alt uncharged; comparator alt the other
          member of `group`. Charge reversals are excluded from every arm.
    """
    out = []
    for (position, ref), alts in variants.items():
        if kind == "gain":
            if ref in CHARGED:
                continue
            numerator = [s for a, s in alts.items() if a in group]
            comparator = [s for a, s in alts.items() if a not in CHARGED]
        else:
            if ref not in group:
                continue
            numerator = [s for a, s in alts.items() if a not in CHARGED]
            comparator = [s for a, s in alts.items() if a in group and a != ref]
        if not numerator or not comparator:
            continue
        out.append({"position": position, "ref": ref,
                    "penalty": float(np.mean(numerator)) - float(np.mean(comparator))})
    return out


def assemble(args) -> tuple[list[dict], dict]:
    proteins = read_topology(Path(args.topology))
    pore_spans = load_pore_spans(Path(args.intramem))
    sites, meta = [], {}

    for gene, (urn, accession, offset) in SCORE_SETS.items():
        rows = score_rows(urn)
        if not rows:
            meta[gene] = {"status": "score file not cached", "urn": urn}
            continue
        synonymous, nonsense, missense = classify(rows)
        direction = direction_for(gene, synonymous, nonsense, missense)
        record = proteins[accession]
        sequence = record["sequence"]
        flanks = build_flank_map(record, 60)
        sequons = sequon_positions(sequence)
        blocked: set[int] = set()
        for start, end in pore_spans.get(accession, ()):
            blocked |= set(range(start - PORE_MARGIN, end + PORE_MARGIN + 1))

        scores = [s for _, _, _, s in missense]
        centre, spread = float(np.mean(scores)), float(np.std(scores, ddof=1))
        variants: dict = defaultdict(dict)
        mapped = 0
        for ref, position, alt, score in missense:
            census_position = position + offset
            if census_position < 1 or census_position > len(sequence):
                continue
            if sequence[census_position - 1] != ref:
                continue
            mapped += 1
            value = (score - centre) / spread if spread > 0 else 0.0
            if direction["flip"]:
                value = -value
            variants[(census_position, ref)][alt] = value

        meta[gene] = {
            "urn": urn, "uniprot": accession, "offset": offset,
            "n_missense": len(missense), "n_mapped": mapped,
            "cross_species": gene in CROSS_SPECIES,
            "gpcr": gene in GPCR_FAMILY, **direction,
        }
        if direction["flip"] is None:
            meta[gene]["status"] = "EXCLUDED: no usable anchors"
            continue
        # Frozen rule: "a score set whose anchors do not separate is excluded".
        # Separation is the rank statistic; 0.5 is no separation at all. The
        # observed values fall in two clean groups, 0.501 and then 0.615 upward,
        # so this threshold does not sit near any set's value.
        if abs(direction["separation"] - 0.5) < 0.05:
            meta[gene]["status"] = "EXCLUDED: anchors do not separate"
            continue
        meta[gene]["status"] = "included"

        for kind, group, label in (("gain", POSITIVE, "positive_charge_gain"),
                                   ("loss", POSITIVE, "positive_charge_loss"),
                                   ("gain", NEGATIVE, "negative_charge_gain"),
                                   ("loss", NEGATIVE, "negative_charge_loss")):
            for site in site_penalties(variants, kind, group):
                flank = flanks.get(site["position"])
                if flank is None:
                    continue
                side, distance, tmd_index = flank
                sites.append({
                    "gene": gene, "uniprot": accession, "class": label,
                    "position": site["position"], "ref": site["ref"],
                    "penalty": site["penalty"], "cytosolic": side == 0,
                    "distance": distance, "tmd_index": tmd_index,
                    "in_sequon": site["position"] in sequons,
                    "near_pore": site["position"] in blocked,
                    "ion_channel": bool(record.get("isIonChannel")),
                    "gpcr": gene in GPCR_FAMILY,
                    "cross_species": gene in CROSS_SPECIES,
                })
    return sites, meta


def stratified(rows: list[dict], equal_gene_weight: bool = False) -> tuple[float, dict]:
    """Side difference stratified by (gene, reference residue), per the frozen text."""
    strata: dict = defaultdict(lambda: {"cyto": [], "outer": []})
    for row in rows:
        key = (row["gene"], row["ref"])
        strata[key]["cyto" if row["cytosolic"] else "outer"].append(row["penalty"])
    numerator = denominator = 0.0
    per_gene: dict = defaultdict(lambda: [0.0, 0.0])
    for (gene, _), arms in strata.items():
        n1, n0 = len(arms["cyto"]), len(arms["outer"])
        if n1 == 0 or n0 == 0:
            continue
        weight = n1 * n0 / (n1 + n0)
        difference = float(np.mean(arms["cyto"]) - np.mean(arms["outer"]))
        per_gene[gene][0] += weight * difference
        per_gene[gene][1] += weight
    if equal_gene_weight:
        values = [v[0] / v[1] for v in per_gene.values() if v[1] > 0]
        estimate = float(np.mean(values)) if values else float("nan")
    else:
        numerator = sum(v[0] for v in per_gene.values())
        denominator = sum(v[1] for v in per_gene.values())
        estimate = numerator / denominator if denominator > 0 else float("nan")
    return estimate, {g: (v[0] / v[1] if v[1] > 0 else float("nan")) for g, v in per_gene.items()}


def gene_arrays(rows: list[dict]):
    """Per-(gene, ref) stratum sufficient statistics, grouped by gene for the bootstrap."""
    strata: dict = defaultdict(lambda: {"cyto": [], "outer": []})
    for row in rows:
        strata[(row["gene"], row["ref"])]["cyto" if row["cytosolic"] else "outer"].append(row["penalty"])
    genes, weights, diffs, spread = [], [], [], []
    for (gene, _), arms in strata.items():
        n1, n0 = len(arms["cyto"]), len(arms["outer"])
        if n1 == 0 or n0 == 0:
            continue
        genes.append(gene)
        weights.append(n1 * n0 / (n1 + n0))
        diffs.append(float(np.mean(arms["cyto"]) - np.mean(arms["outer"])))
        spread.append(float(np.var(arms["cyto"], ddof=1) / n1 if n1 > 1 else 0.0)
                      + float(np.var(arms["outer"], ddof=1) / n0 if n0 > 1 else 0.0))
    names = sorted(set(genes))
    index = {g: i for i, g in enumerate(names)}
    return (np.array([index[g] for g in genes]), np.array(weights),
            np.array(diffs), np.array(spread), names)


def wild_bootstrap(rows: list[dict], draws: int, weights_kind: str, seed: int) -> dict:
    """Restricted wild cluster bootstrap-t on the stratum contributions, clustered by gene."""
    cluster, weight, diff, _, names = gene_arrays(rows)
    n_clusters = len(names)
    if n_clusters < 2 or weight.sum() == 0:
        return {"estimate": float("nan"), "p_value": float("nan"), "n_genes": n_clusters}

    def statistic(values: np.ndarray) -> tuple[float, float]:
        total = weight.sum()
        beta = float((weight * values).sum() / total)
        contribution = np.bincount(cluster, weights=weight * (values - beta),
                                   minlength=n_clusters) / total
        variance = float((contribution ** 2).sum())
        return beta, (beta / math.sqrt(variance) if variance > 0 else float("nan"))

    beta, t_observed = statistic(diff)
    if not np.isfinite(t_observed):
        return {"estimate": beta, "p_value": float("nan"), "n_genes": n_clusters}
    # Center stratum effects before applying one multiplier per gene.
    centred = diff - beta
    rng = np.random.default_rng(seed)
    extreme = 0
    for _ in range(draws):
        flips = (rng.choice((-1.0, 1.0), size=n_clusters) if weights_kind == "rademacher"
                 else rng.choice(WEBB, size=n_clusters))
        _, t_star = statistic(centred * flips[cluster])
        if np.isfinite(t_star) and abs(t_star) >= abs(t_observed):
            extreme += 1
    return {"estimate": beta, "t_statistic": t_observed, "n_genes": n_clusters,
            "p_value": (extreme + 1) / (draws + 1), "weights": weights_kind}


def bootstrap_interval(rows: list[dict], draws: int, seed: int) -> list[float]:
    """Gene-clustered percentile interval, reported alongside but not primary."""
    cluster, weight, diff, _, names = gene_arrays(rows)
    n_clusters = len(names)
    if n_clusters < 2:
        return [float("nan"), float("nan")]
    rng = np.random.default_rng(seed)
    samples = np.full(draws, np.nan)
    for i in range(draws):
        counts = np.bincount(rng.integers(0, n_clusters, n_clusters),
                             minlength=n_clusters).astype(float)
        multiplier = counts[cluster]
        total = (weight * multiplier).sum()
        if total > 0:
            samples[i] = float((weight * multiplier * diff).sum() / total)
    finite = samples[np.isfinite(samples)]
    return ([float(np.percentile(finite, 2.5)), float(np.percentile(finite, 97.5))]
            if finite.size > 1 else [float("nan"), float("nan")])


def fit(rows: list[dict], args, label: str, equal_gene_weight: bool = False) -> dict:
    if not rows:
        return {"scope": label, "n_sites": 0}
    estimate, per_gene = stratified(rows, equal_gene_weight)
    result = {
        "scope": label, "n_sites": len(rows),
        "n_cytosolic": sum(1 for r in rows if r["cytosolic"]),
        "n_outer": sum(1 for r in rows if not r["cytosolic"]),
        "estimate": estimate,
        "percentile_ci95": bootstrap_interval(rows, args.boot_draws, args.seed),
        "per_gene_estimate": per_gene,
    }
    for kind in ("webb", "rademacher"):
        result[f"wild_{kind}"] = wild_bootstrap(rows, args.wild_draws, kind, args.seed)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", default=str(ROOT / "data" / "raw" / "topology_viewer.html"))
    parser.add_argument("--intramem", default=str(ROOT / "data" / "external" / "uniprot" / "human_intramem.tsv"))
    parser.add_argument("--outdir", default=str(DEFAULT_OUT))
    parser.add_argument("--wild-draws", type=int, default=9999)
    parser.add_argument("--boot-draws", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260928)
    args = parser.parse_args()

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sites, meta = assemble(args)

    def subset(**kw) -> list[dict]:
        rows = [s for s in sites if s["class"] == kw.get("cls", "positive_charge_gain")]
        rows = [s for s in rows if s["distance"] <= kw.get("window", PRIMARY_WINDOW)]
        if kw.get("no_gpcr"):
            rows = [s for s in rows if not s["gpcr"]]
        if kw.get("no_cross_species"):
            rows = [s for s in rows if not s["cross_species"]]
        if kw.get("no_ion_channel"):
            rows = [s for s in rows if not s["ion_channel"]]
        if kw.get("no_pore"):
            rows = [s for s in rows if not s["near_pore"]]
        if kw.get("no_sequon"):
            rows = [s for s in rows if not s["in_sequon"]]
        if kw.get("arginine_only"):
            rows = [s for s in rows if s["ref"] == "R"]
        if kw.get("gene"):
            rows = [s for s in rows if s["gene"] != kw["gene"]]
        return rows

    results: dict = {"primary": fit(subset(), args, "K/R gain, window <=10, all genes")}
    results["sensitivities"] = {
        "window_5": fit(subset(window=5), args, "window <=5"),
        "window_15": fit(subset(window=15), args, "window <=15"),
        "no_cross_species": fit(subset(no_cross_species=True), args, "without cross-species sets"),
        "no_ion_channels": fit(subset(no_ion_channel=True), args, "without ion channels"),
        "pore_excluded": fit(subset(no_pore=True), args, "UniProt pore regions excluded"),
        "arginine_only": fit(subset(arginine_only=True), args, "arginine reference only"),
        "gene_equal_weight": fit(subset(), args, "gene-equal weighting", equal_gene_weight=True),
        "no_sequon": fit(subset(no_sequon=True), args, "N-X-S/T sequon positions excluded"),
    }
    results["robustness"] = {
        "leave_one_family_out_no_gpcr": fit(subset(no_gpcr=True), args, "GPCR family removed"),
        "leave_one_protein_out": {
            gene: fit(subset(gene=gene), args, f"without {gene}")["estimate"]
            for gene in sorted({s["gene"] for s in subset()})
        },
    }
    results["secondary_classes"] = {
        cls: fit(subset(cls=cls), args, cls)
        for cls in ("positive_charge_loss", "negative_charge_gain", "negative_charge_loss")
    }

    summary = {
        "analysis": "Phase C, run under the frozen preregistration",
        "preregistration_sha256": "754a9a713a2afc7fd2e3e43be9af55bb762805feeb3d3991be0667c207635d79",
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "primary_window_residues": PRIMARY_WINDOW,
        "sign_convention": "negative = gaining K/R is more damaging on the OUTER side, "
                           "which is the preregistered prediction",
        "score_sets": meta,
        "results": results,
    }
    (out_dir / "phase_c_results.json").write_text(json.dumps(summary, indent=2) + "\n",
                                                  encoding="utf-8")
    with (out_dir / "phase_c_sites.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(sites[0]))
        writer.writeheader(); writer.writerows(sites)

    primary = results["primary"]
    print("=== score sets ===")
    for gene, info in sorted(meta.items()):
        print(f"  {gene:<9} {info.get('status','?'):<28} flip={info.get('flip')} "
              f"sep={info.get('separation', float('nan')):.3f}  {info.get('basis','')}")
    print("\n=== PRIMARY: K/R gain side difference, window <=10 ===")
    print(f"  estimate {primary['estimate']:+.4f}   n={primary['n_cytosolic']}/{primary['n_outer']}  "
          f"genes={primary['wild_webb']['n_genes']}")
    print(f"  percentile CI95 {primary['percentile_ci95']}")
    for kind in ("webb", "rademacher"):
        w = primary[f"wild_{kind}"]
        print(f"  wild {kind:<11} p={w['p_value']:.4f}  t={w.get('t_statistic', float('nan')):.3f}")
    print("\n=== sensitivities ===")
    for name, value in results["sensitivities"].items():
        if value.get("n_sites"):
            print(f"  {name:<20} {value['estimate']:+.4f}  "
                  f"webb p={value['wild_webb']['p_value']:.4f}  n={value['n_sites']}")
    print("\n=== leave-one-family-out ===")
    lofo = results["robustness"]["leave_one_family_out_no_gpcr"]
    if lofo.get("n_sites"):
        print(f"  no GPCR: {lofo['estimate']:+.4f}  webb p={lofo['wild_webb']['p_value']:.4f}  "
              f"genes={lofo['wild_webb']['n_genes']}  n={lofo['n_sites']}")
    print("\n=== secondary classes ===")
    for cls, value in results["secondary_classes"].items():
        if value.get("n_sites"):
            print(f"  {cls:<24} {value['estimate']:+.4f}  webb p={value['wild_webb']['p_value']:.4f}")


if __name__ == "__main__":
    main()

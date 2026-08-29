"""Audit disulfide-loss scores with matched cysteine controls."""

from __future__ import annotations

import csv
import gzip
import json
import re
import sys
import time
import warnings
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from run_gate1 import read_topology  # noqa: E402

AM_GZ = ROOT / "data" / "raw" / "AlphaMissense_aa_substitutions.tsv.gz"
ESM_ZIP = ROOT / "data" / "raw" / "ESM1b_ALL_hum_isoforms_LLR.zip"
ESM_PREFIX = "content/ALL_hum_isoforms_ESM1b_LLR/"
UNIPROT_TSV = (ROOT / "data" / "external" / "disulfide_rule_20261025"
               / "human_reviewed_features.tsv")
UNIPROT_MANIFEST = (ROOT / "data" / "external" / "disulfide_rule_20261025"
                    / "source_manifest.json")
OUT = ROOT / "results" / "beat_am_20261016" / "disulfide_rule"

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
AIDX = {a: i for i, a in enumerate(ALPHABET)}

PROBES = ("S", "A")             # Keep primary C-to-S and sensitivity C-to-A probes separate.
MIN_REF_ARM = 2                 # Require at least two free controls as in the frozen B8 stratum definition.
N_BOOT = 1999
N_PERM = 1999
SEED = 20261016
HYDRO_WINDOW = 7
N_HYDRO_BINS = 4
MIN_BONDS_FOR_PAIRING = 2       # At least two complete bonds are needed to permute partners within a protein.

# Use Kyte-Doolittle hydropathy as the S3 burial proxy.
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
      "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
      "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
      "Y": -1.3, "V": 4.2}

COMPARTMENTS = ("noncytosolic", "cytosolic", "tm", "unknown")
CIDX = {c: i for i, c in enumerate(COMPARTMENTS)}

DECLARED = {
    "note_on_ordering":
        "The eligibility census (UniProt feature annotation only: how many "
        "bonded and free cysteines exist, their evidence codes, their "
        "compartments) was computed BEFORE these predictions were written, "
        "because a prediction has to be answerable by the sample that exists. "
        "No AlphaMissense or ESM1b value was read until after these were fixed "
        "in the script. S5 is declared unavailable on the basis of that census.",
    "P1_primary_magnitude":
        "score(C->S) is higher -- more damaging -- at bonded cysteines than at "
        "matched free cysteines (same protein, same compartment), with a 95% "
        "interval excluding zero.",
    "P2_primary_pairing":
        "the mean within-bond |difference| in score(C->S) between the two "
        "endpoints of the same bond is LOWER than the mean under within-protein "
        "permuted re-pairing of the same endpoints.",
    "S1_secondary_C_to_A":
        "C->A shows the same direction as P1. Reported separately; C->S and "
        "C->A are never averaged.",
    "S2_secondary_direct_evidence":
        "P1 survives restricting the bonded set to endpoints with direct "
        "experimental evidence (ECO:0000269, or ECO:0007744 with a PDB source).",
    "S3_secondary_burial_proxy_match":
        "P1 survives adding a local-hydrophobicity quartile to the stratum, so "
        "bonded and free cysteines are matched on a sequence proxy for burial "
        "as well as on protein and compartment. If P1 disappears under this "
        "match, the primary result is substantially a burial artefact.",
    "S4_secondary_free_set_quality":
        "P1 survives restricting free controls to proteins carrying at least "
        "two annotated disulfide bonds, whose disulfide pattern a curator has "
        "clearly characterised.",
    "S5_control_compartment_placebo":
        "DECLARED UNAVAILABLE BEFORE RUNNING. A cytosolic placebo arm would "
        "test whether the effect appears where disulfides cannot form. The "
        "eligibility census finds 16,517 non-cytosolic bonded cysteines against "
        "16 cytosolic ones among the census proteins, so the arm does not "
        "exist. Reported as NOT AVAILABLE, never as a pass.",
    "X1_cross_scorer":
        "AlphaMissense and ESM1b are reported separately and compared on "
        "DIRECTION and pass/fail pattern only. No standardised (SD-unit) "
        "quantity is compared between scorers.",
}

LIMITATIONS = [
    "BURIAL AND CONSERVATION CONFOUND, the main one. Bonded cysteines are on "
    "average more buried and more conserved than free ones, so a positive P1 is "
    "consistent both with the scorer having learned bond loss and with it having "
    "learned only that these positions matter. Matching within protein and "
    "compartment helps and is not sufficient. Using annotation alone this "
    "matched-site design CANNOT fully separate bond loss from burial. S3 bounds "
    "the burial part with a weak sequence proxy. P2, the pairing test, is the "
    "one result here that this confound cannot produce.",
    "THE FREE SET IS AN UNKNOWN SET. Absence of a disulfide annotation is not "
    "evidence that a cysteine is free. The annotation source supplies no usable "
    "set of positively evidenced free thiols. A genuinely bonded cysteine "
    "sitting in the free arm biases P1 toward zero, so this cost works against "
    "the declared direction rather than manufacturing it.",
    "THE WITHIN-POSITION ESTIMATOR OF B8 AND B9 IS NOT AVAILABLE FOR THIS RULE. "
    "Every ordinary non-Cys substitution at a bonded cysteine removes the bond, "
    "so there is no motif-preserving comparator at the same position and no "
    "exact cancellation of per-position effects. This audit is a weaker "
    "instrument than B8 and B9 for that reason, not for a reason of "
    "implementation.",
    "NO STRUCTURE WAS READ. No solvent accessibility, no PDB coordinates, no "
    "predicted structure. Burial enters only through the S3 hydrophobicity "
    "proxy, which is weak.",
    "EVIDENCE MIXTURE. Most annotated endpoints are inferred from family rules "
    "rather than observed directly. The breakdown is reported and S2 is the "
    "experimental-only arm.",
    "NO CYTOSOLIC PLACEBO. Unlike B8's sequon audit, this rule has no built-in "
    "arm where the chemistry is impossible but the annotation is present. See S5.",
]

POST_HOC = {
    "status": "Observations made AFTER the declared results were computed. "
              "They are not predictions, nothing was redefined to accommodate "
              "them, and the declared primary is reported unchanged.",
    "the_pooled_primary_is_inflated_by_the_TM_compartment": {
        "observation":
            "The declared primary pools all compartments as strata. The "
            "per-compartment breakdown shows the non-cytosolic arm -- where "
            "structural disulfides actually form, and where both arms are "
            "large -- gives a SMALLER gap than the pooled figure, while the "
            "transmembrane arm gives a much larger one and pulls the pooled "
            "estimate up. See results[scorer]['by_compartment'].",
        "why_the_TM_arm_is_not_credible":
            "A structural disulfide inside the hydrocarbon core is "
            "biophysically implausible. Those ~600 bonded sites are most "
            "likely boundary disagreements between the census TM limits and "
            "the UniProt disulfide annotation: a cysteine just outside the "
            "real membrane boundary assigned to the TM compartment. Their "
            "free comparators, by contrast, are genuine core cysteines, for "
            "which a substitution is unremarkable. That makes the TM arm a "
            "comparison of two different environments, not a bond-loss test.",
        "how_to_read_the_result":
            "The non-cytosolic gap is the interpretable primary number. The "
            "pooled figure should not be quoted on its own.",
    },
    "the_cytosolic_arm_argues_against_the_bond_loss_reading": {
        "observation":
            "S5 was declared unavailable before running, and it stays "
            "unavailable: 16 cytosolic bonded cysteines over 5 strata cannot "
            "support a placebo test. But the gap those few sites give is "
            "POSITIVE and its interval excludes zero, for both scorers.",
        "why_this_matters":
            "Disulfides cannot form in the reducing cytosol, so a genuine "
            "bond-loss effect should be absent there. An effect that appears "
            "anyway is what the burial-and-conservation confound predicts. "
            "With 16 sites this is far too thin to conclude anything, and it "
            "is reported here precisely because it points AGAINST the "
            "comfortable reading rather than for it. It is the one piece of "
            "evidence in this file that bears on the placebo question, and it "
            "is not reassuring.",
    },
}

_LOG_LINES: list[str] = []


def log(m: str) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {m}"
    _LOG_LINES.append(line)
    print(line, flush=True)



def features(text: str, kind: str):
    """Yield (start, end, detail) for one UniProt feature type.

    The detail string runs to the next feature of the same kind, so evidence
    codes stay attached to the feature they belong to.
    """
    pat = (rf'{re.escape(kind)} (\d+)(?:\.\.(\d+))?;?'
           rf'(.*?)(?={re.escape(kind)} |\Z)')
    for m in re.finditer(pat, text or ""):
        yield int(m[1]), int(m[2] or m[1]), m[3]


def direct_evidence(detail: str) -> bool:
    """Direct experimental support, as c2_disulfide_inventory.py defines it."""
    return ("ECO:0000269" in detail
            or bool(re.search(r"ECO:0007744\|PDB:", detail)))


def evidence_code(detail: str) -> str:
    m = re.search(r"(ECO:\d{7})", detail or "")
    return m[1] if m else "none"


def load_uniprot() -> dict:
    csv.field_size_limit(10_000_000)
    with UNIPROT_TSV.open(encoding="utf-8", newline="") as fh:
        return {r["Entry"]: r for r in csv.DictReader(fh, delimiter="\t")}



def compartment_map(rec: dict) -> dict:
    m = {}
    for _lid, lp in (rec.get("loops") or {}).items():
        loc = (lp.get("location") or "").lower()
        if loc in ("inside", "outside"):
            c = "cytosolic" if loc == "inside" else "noncytosolic"
            for p in range(int(lp["start"]), int(lp["end"]) + 1):
                m[p] = c
    for tm in (rec.get("tmds") or []):
        for p in range(int(tm["start"]), int(tm["end"]) + 1):
            m[p] = "tm"
    return m


def window_hydropathy(seq: str, pos: int) -> float:
    lo = max(0, pos - 1 - HYDRO_WINDOW)
    hi = min(len(seq), pos + HYDRO_WINDOW)
    vals = [KD.get(c, 0.0) for c in seq[lo:hi]]
    return float(np.mean(vals)) if vals else 0.0


def enumerate_sites(census: dict, uni: dict) -> dict:
    """Eligible cysteine sites plus the bond graph. Annotation only, no scores."""
    counts = defaultdict(int)
    ev_breakdown = defaultdict(int)
    acc_l, pos_l, state_l, comp_l, hyd_l, direct_l, nbonds_l = [], [], [], [], [], [], []
    partner_l = []
    bonds_by_prot: dict[str, list[tuple[int, int, bool]]] = {}

    shared = sorted(set(census) & set(uni))
    counts["census_proteins"] = len(census)
    counts["census_in_reviewed_uniprot"] = len(shared)

    for acc in shared:
        rec, row = census[acc], uni[acc]
        seq = rec.get("sequence") or ""
        if not seq:
            counts["protein_no_census_sequence"] += 1
            continue
        if seq != row["Sequence"]:
            # Require annotation and score coordinates to refer to the same sequence.
            counts["protein_sequence_mismatch_excluded"] += 1
            continue
        cysset = {p for p, ch in enumerate(seq, 1) if ch == "C"}
        if not cysset:
            continue
        counts["proteins_considered"] += 1
        counts["cysteines_total"] += len(cysset)


        excl: dict[int, str] = {}
        for kind, col in (("BINDING", "Binding site"),
                          ("ACT_SITE", "Active site"),
                          ("MOD_RES", "Modified residue"),
                          ("SIGNAL", "Signal peptide")):
            for s, e, _d in features(row.get(col, ""), kind):
                for p in range(s, e + 1):
                    if p in cysset and p not in excl:
                        excl[p] = col
                        counts[f"excluded_{col.replace(' ', '_')}"] += 1


        partner: dict[int, int] = {}
        ev_direct: dict[int, bool] = {}
        uncertain: set[int] = set()
        pairs: list[tuple[int, int, bool]] = []
        for s, e, detail in features(row.get("Disulfide bond", ""), "DISULFID"):
            counts["disulfide_features"] += 1
            if s == e:
                counts["feature_singleton_no_partner"] += 1
                uncertain.add(s)
                continue
            if s < 1 or e > len(seq) or seq[s - 1] != "C" or seq[e - 1] != "C":
                # Count unexpected residue mismatches even after the full-sequence equality check.
                counts["feature_position_not_cysteine"] += 1
                continue
            if "Interchain" in detail:
                counts["feature_interchain_excluded"] += 1
                uncertain.update((s, e))
                continue
            d = direct_evidence(detail)
            conflict = False
            for p, q in ((s, e), (e, s)):
                if p in partner and partner[p] != q:
                    uncertain.add(p)
                    counts["feature_alternative_partner"] += 1
                    conflict = True
            if conflict:
                continue
            partner[s], partner[e] = e, s
            ev_direct[s] = ev_direct[e] = d
            ev_breakdown["direct_experimental" if d else "inferred"] += 2
            ev_breakdown[evidence_code(detail)] += 2
            pairs.append((s, e, d))

        n_bonds_here = len([1 for s, e, _d in pairs
                            if s not in uncertain and e not in uncertain])

        cmap = compartment_map(rec)
        complete: list[tuple[int, int, bool]] = []
        for p in sorted(cysset):
            if p in uncertain:
                counts["site_uncertain_excluded"] += 1
                continue
            if p in excl:
                continue
            bonded = p in partner
            acc_l.append(acc)
            pos_l.append(p)
            state_l.append(1 if bonded else 0)
            comp_l.append(CIDX[cmap.get(p, "unknown")])
            hyd_l.append(window_hydropathy(seq, p))
            direct_l.append(1 if (bonded and ev_direct.get(p)) else 0)
            nbonds_l.append(n_bonds_here)
            partner_l.append(partner.get(p, 0) if bonded else 0)
            counts["bonded_sites" if bonded else "free_sites"] += 1

        for s, e, d in pairs:
            if s in uncertain or e in uncertain:
                continue
            if (s in excl) != (e in excl):
                counts["bond_orphaned_one_endpoint_excluded"] += 1
                continue
            if s in excl and e in excl:
                counts["bond_both_endpoints_excluded"] += 1
                continue
            complete.append((s, e, d))
        counts["bonds_complete_both_endpoints_eligible"] += len(complete)
        if complete:
            bonds_by_prot[acc] = complete

    prot_names = sorted({a for a in acc_l})
    pidx = {a: i for i, a in enumerate(prot_names)}
    sites = {
        "acc": np.array(acc_l),
        "pos": np.array(pos_l, dtype=np.int32),
        "prot": np.array([pidx[a] for a in acc_l], dtype=np.int32),
        "bonded": np.array(state_l, dtype=np.int8),
        "comp": np.array(comp_l, dtype=np.int8),
        "hyd": np.array(hyd_l, dtype=np.float32),
        "direct": np.array(direct_l, dtype=np.int8),
        "nbonds": np.array(nbonds_l, dtype=np.int16),
        "partner": np.array(partner_l, dtype=np.int32),
        "prot_names": prot_names,
        "pidx": pidx,
        "counts": dict(counts),
        "evidence_breakdown": dict(ev_breakdown),
        "bonds_by_prot": bonds_by_prot,
    }
    log(f"  proteins considered: {counts['proteins_considered']:,}  "
        f"sites: {len(pos_l):,}  bonded: {counts['bonded_sites']:,}  "
        f"free: {counts['free_sites']:,}")
    log(f"  complete bonds (both endpoints eligible): "
        f"{counts['bonds_complete_both_endpoints_eligible']:,}")

    # Compute S3 hydropathy quartiles over the eligible control pool.
    q = np.quantile(sites["hyd"], [0.25, 0.5, 0.75])
    sites["hyd_bin"] = np.searchsorted(q, sites["hyd"]).astype(np.int8)
    sites["hyd_quartile_edges"] = [float(x) for x in q]
    return sites


# Stream AlphaMissense scores while indexing ESM1b matrices by position.
def build_index(sites: dict) -> dict:
    return {(a, p): i for i, (a, p) in
            enumerate(zip(sites["acc"].tolist(), sites["pos"].tolist()))}


def extract_am(index: dict, n_rows: int) -> np.ndarray:
    S = np.full((n_rows, 20), np.nan, dtype=np.float32)
    seen = filled = 0
    t0 = time.time()
    with gzip.open(AM_GZ, "rt") as fh:
        for line in fh:
            if line[0] == "#" or line.startswith("uniprot_id"):
                continue
            seen += 1
            if seen % 40_000_000 == 0:
                log(f"    {seen:,} lines, {time.time() - t0:.0f}s, {filled:,} cells")
            t1 = line.find("\t")
            t2 = line.find("\t", t1 + 1)
            var = line[t1 + 1:t2]
            pos_s = var[1:-1]
            if not pos_s.isdigit():
                continue
            row = index.get((line[:t1], int(pos_s)))
            if row is None:
                continue
            col = AIDX.get(var[-1])
            if col is None:
                continue
            t3 = line.find("\t", t2 + 1)
            S[row, col] = float(line[t2 + 1:t3])
            filled += 1
    log(f"  AM pass: {seen:,} lines in {time.time() - t0:.0f}s, {filled:,} scores")
    return S


def extract_esm(sites: dict, index: dict, n_rows: int) -> np.ndarray:
    """Per-accession CSVs out of the zip. Writes damage = -LLR (negated)."""
    S = np.full((n_rows, 20), np.nan, dtype=np.float32)
    wanted: dict[str, set[int]] = defaultdict(set)
    for (a, p) in index:
        wanted[a].add(p)
    found = missing = filled = 0
    t0 = time.time()
    with zipfile.ZipFile(ESM_ZIP) as z:
        available = set(z.namelist())
        for n_done, acc in enumerate(sorted(wanted), start=1):
            if n_done % 500 == 0:
                log(f"    {n_done:,} accessions, {time.time() - t0:.0f}s, "
                    f"{filled:,} cells")
            name = f"{ESM_PREFIX}{acc}_LLR.csv"
            if name not in available:
                missing += 1
                continue
            with z.open(name) as fh:
                rows = fh.read().decode("utf-8").splitlines()
            if len(rows) < 2:
                missing += 1
                continue
            found += 1
            want = wanted[acc]
            header = rows[0].split(",")
            cols: list[tuple[int, int]] = []
            for c, cell in enumerate(header[1:], start=1):
                parts = cell.strip().split()
                if len(parts) != 2:
                    continue
                p = int(parts[1])
                if p not in want:
                    continue
                r = index.get((acc, p))
                if r is not None:
                    cols.append((c, r))
            if not cols:
                continue
            for line in rows[1:]:
                cells = line.split(",")
                col = AIDX.get(cells[0].strip())
                if col is None:
                    continue
                for c, r in cols:
                    v = cells[c]
                    if not v or v == "nan":
                        continue
                    try:
                        S[r, col] = -float(v)      # Negate fitness-oriented LLR values to obtain damage-oriented scores.
                    except ValueError:
                        continue
                    filled += 1
    log(f"  ESM1b: {found:,} accessions read, {missing:,} absent from the zip, "
        f"{filled:,} scores in {time.time() - t0:.0f}s")
    return S


def load_or_extract(tag: str, sites: dict, index: dict, fn) -> np.ndarray:
    cache = OUT / f"position_scores_{tag}.npz"
    n = len(sites["pos"])
    if cache.exists():
        z = np.load(cache, allow_pickle=False)
        if (z["pos"].shape == sites["pos"].shape
                and bool((z["pos"] == sites["pos"]).all())
                and bool((z["acc"] == sites["acc"]).all())):
            log(f"  reusing cached {tag} scores from {cache.name}")
            return z["scores"]
        log(f"  cached {tag} scores do not match the site table; re-extracting")
    S = fn(index, n) if tag == "am" else fn(sites, index, n)
    np.savez_compressed(cache, scores=S, acc=sites["acc"], pos=sites["pos"])
    log(f"  cached {tag} scores -> {cache.name}")
    return S



def ci(draws: np.ndarray) -> list[float]:
    lo, hi = np.nanquantile(draws, [0.025, 0.975])
    return [float(lo), float(hi)]


def matched_gap(prot, stratum, bonded, v, n_prot, rng) -> dict | None:
    """Weighted bonded-minus-free gap in v within strata; whole-protein bootstrap.

    `stratum` is an integer stratum id that already encodes the protein, so a
    protein resample moves all of its strata together -- and therefore moves
    both endpoints of every intrachain bond together.
    """
    ok = np.isfinite(v)
    if ok.sum() == 0:
        return None
    prot, stratum, bonded, v = prot[ok], stratum[ok], bonded[ok], v[ok]

    key = stratum.astype(np.int64) * 2 + bonded.astype(np.int64)
    uniq, inv = np.unique(key, return_inverse=True)
    cnt = np.bincount(inv).astype(np.float64)
    tot = np.bincount(inv, weights=v.astype(np.float64))
    sid_of_key = uniq // 2
    arm_of_key = (uniq % 2).astype(np.int64)
    sids, sinv = np.unique(sid_of_key, return_inverse=True)
    N = np.zeros((len(sids), 2))
    M = np.zeros((len(sids), 2))
    N[sinv, arm_of_key] = cnt
    M[sinv, arm_of_key] = np.divide(tot, cnt, out=np.zeros_like(tot), where=cnt > 0)

    # Derive each stratum's protein index from its contributing sites.
    prot_of_stratum = np.zeros(len(sids), dtype=np.int64)
    order = np.argsort(stratum, kind="stable")
    s_sorted = stratum[order]
    first = np.searchsorted(s_sorted, sids)
    prot_of_stratum[:] = prot[order][first]

    usable = (N[:, 1] >= 1) & (N[:, 0] >= MIN_REF_ARM)
    if not usable.any():
        return None
    n1, n0 = N[usable, 1], N[usable, 0]
    w = n1 * n0 / (n1 + n0)
    d = M[usable, 1] - M[usable, 0]
    pof = prot_of_stratum[usable]
    num = np.bincount(pof, weights=w * d, minlength=n_prot)
    den = np.bincount(pof, weights=w, minlength=n_prot)

    est = float(num.sum() / den.sum())
    live = np.flatnonzero(den > 0)
    draws = np.empty(N_BOOT)
    for b in range(N_BOOT):
        pick = live[rng.integers(0, len(live), size=len(live))]
        dd = den[pick].sum()
        draws[b] = num[pick].sum() / dd if dd else np.nan
    sd = float(np.std(v, ddof=1)) if len(v) > 1 else float("nan")
    return {
        "estimate": est,
        "ci95": ci(draws),
        "sd_of_values": sd,
        "estimate_standardized_within_this_scorer": (est / sd) if sd else None,
        "n_strata": int(usable.sum()),
        "n_proteins": int(len(live)),
        "n_bonded_sites": int(n1.sum()),
        "n_free_sites": int(n0.sum()),
        "mean_bonded": float(np.average(M[usable, 1], weights=w)),
        "mean_free": float(np.average(M[usable, 0], weights=w)),
    }


def stratum_ids(sites: dict, keys: tuple[str, ...]) -> np.ndarray:
    """Integer stratum id from the protein plus the named site attributes."""
    out = sites["prot"].astype(np.int64)
    for k in keys:
        width = {"comp": len(COMPARTMENTS), "hyd_bin": N_HYDRO_BINS}[k]
        out = out * width + sites[k].astype(np.int64)
    return out



def pairing_test(sites: dict, v: np.ndarray, rng) -> dict:
    """Within-bond |difference| against a within-protein re-pairing null."""
    row_of = {(a, p): i for i, (a, p) in
              enumerate(zip(sites["acc"].tolist(), sites["pos"].tolist()))}
    per_prot: list[np.ndarray] = []
    n_bonds = n_prot_used = 0
    dropped_unscored = 0
    bonds_per_prot = defaultdict(int)
    for acc, bonds in sites["bonds_by_prot"].items():
        vals = []
        for s, e, _d in bonds:
            rs, re_ = row_of.get((acc, s)), row_of.get((acc, e))
            if rs is None or re_ is None:
                dropped_unscored += 1
                continue
            a, b = v[rs], v[re_]
            if not (np.isfinite(a) and np.isfinite(b)):
                dropped_unscored += 1
                continue
            vals.append((a, b))
        if len(vals) < MIN_BONDS_FOR_PAIRING:
            continue
        per_prot.append(np.asarray(vals, dtype=np.float64))
        n_bonds += len(vals)
        n_prot_used += 1
        bonds_per_prot[min(len(vals), 10)] += 1

    if not per_prot:
        return {"status": "no protein had two complete scored bonds"}

    # Use the same bond-weighted denominator for observed and permuted statistics.
    observed = float(
        np.sum([np.abs(a[:, 0] - a[:, 1]).sum() for a in per_prot]) / n_bonds)

    flat = [a.reshape(-1).copy() for a in per_prot]
    null = np.empty(N_PERM)
    for d in range(N_PERM):
        total = 0.0
        for arr in flat:
            rng.shuffle(arr)
            total += np.abs(arr[0::2] - arr[1::2]).sum()
        null[d] = total / n_bonds
    n_le = int((null <= observed).sum())
    return {
        "observed_mean_within_bond_abs_diff": observed,
        "null_mean": float(null.mean()),
        "null_ci95": ci(null),
        "difference_observed_minus_null_mean": observed - float(null.mean()),
        "permutation_p_one_sided_observed_lower": (1 + n_le) / (1 + N_PERM),
        "n_permutations": N_PERM,
        "n_bonds": n_bonds,
        "n_proteins": n_prot_used,
        "bonds_per_protein_histogram_capped_at_10":
            {str(k): int(bonds_per_prot[k]) for k in sorted(bonds_per_prot)},
        "bond_endpoints_dropped_unscored": dropped_unscored,
        "what_the_null_holds_fixed":
            "the null re-pairs the SAME endpoint scores within the same "
            "protein, so protein level, compartment, burial and conservation "
            "are all held fixed and only the pairing varies",
        "conservatism":
            "a protein with exactly two bonds has three possible matchings, "
            "one of them the true one, so the null contains the observed "
            "pairing with positive probability; this makes the test "
            "conservative rather than liberal",
    }


def proximity_control(sites: dict, v: np.ndarray, rng) -> dict:
    """POST-HOC control for P2: is the within-bond similarity just proximity?

    Added after P2 was computed, for a confound P2's own null does not remove.
    The two endpoints of a disulfide are often close together in sequence, and
    neighbouring residues tend to receive similar scores for reasons that have
    nothing to do with the bond. P2's null re-pairs endpoints at random, which
    pairs cysteines that are typically FURTHER apart than real partners, so some
    of P2 could be sequence proximity rather than bond representation.

    This control matches on separation. For each true bond, a control pair is
    drawn from the same protein's own bonded endpoints -- never a true bond --
    with the closest available sequence separation. Ties are broken by position
    order, so the choice is deterministic, and each control pair is used once
    where possible.

    DIRECTION DECLARED BEFORE COMPUTING: if P2 reflects the bond rather than
    proximity, the true-bond mean |difference| stays BELOW the separation-
    matched control mean |difference|.
    """
    row_of = {(a, p): i for i, (a, p) in
              enumerate(zip(sites["acc"].tolist(), sites["pos"].tolist()))}
    rows = []
    for acc, bonds in sites["bonds_by_prot"].items():
        pts, pairs = {}, []
        for s, e, _d in bonds:
            rs, re_ = row_of.get((acc, s)), row_of.get((acc, e))
            if rs is None or re_ is None:
                continue
            if not (np.isfinite(v[rs]) and np.isfinite(v[re_])):
                continue
            pairs.append((s, e, rs, re_))
            pts[s], pts[e] = rs, re_
        if len(pairs) < MIN_BONDS_FOR_PAIRING:
            continue
        positions = sorted(pts)
        truth = {frozenset((s, e)) for s, e, _r, _r2 in pairs}
        cand = [(abs(a - b), a, b)
                for i, a in enumerate(positions)
                for b in positions[i + 1:]
                if frozenset((a, b)) not in truth]
        if not cand:
            continue
        cand_sep = np.asarray([c[0] for c in cand], dtype=np.float64)
        used: set[int] = set()
        s_true = s_ctrl = sep_true = sep_ctrl = 0.0
        n = 0
        for s, e, rs, re_ in pairs:
            sep = abs(s - e)
            order = np.argsort(np.abs(cand_sep - sep), kind="stable")
            pick = next((int(i) for i in order if int(i) not in used),
                        int(order[0]))
            used.add(pick)
            csep, a, b = cand[pick]
            s_true += abs(v[rs] - v[re_])
            s_ctrl += abs(v[pts[a]] - v[pts[b]])
            sep_true += sep
            sep_ctrl += csep
            n += 1
        rows.append((s_true, s_ctrl, float(n), sep_true, sep_ctrl))

    if not rows:
        return {"status": "no protein supplied a separation-matched control pair"}
    A = np.asarray(rows, dtype=np.float64)
    tot_n = A[:, 2].sum()
    est_true = A[:, 0].sum() / tot_n
    est_ctrl = A[:, 1].sum() / tot_n
    draws = np.empty(N_BOOT)
    for b in range(N_BOOT):
        pick = rng.integers(0, len(A), size=len(A))
        sub = A[pick]
        d = sub[:, 2].sum()
        draws[b] = (sub[:, 0].sum() - sub[:, 1].sum()) / d if d else np.nan
    return {
        "declared_direction_before_computing":
            "true-bond mean |difference| BELOW the separation-matched control "
            "mean |difference| (a negative difference)",
        "post_hoc":
            "this control was added after P2 was computed, to address a "
            "confound P2's own permutation null does not remove. It is "
            "reported as post hoc, with its direction declared before it "
            "was computed.",
        "true_bond_mean_abs_diff": float(est_true),
        "separation_matched_control_mean_abs_diff": float(est_ctrl),
        "difference_true_minus_control": float(est_true - est_ctrl),
        "ci95_of_difference": ci(draws),
        "mean_sequence_separation_true_bonds": float(A[:, 3].sum() / tot_n),
        "mean_sequence_separation_control_pairs": float(A[:, 4].sum() / tot_n),
        "n_bonds": int(tot_n),
        "n_proteins": int(len(A)),
    }



def directional(res: dict | None, direction: int) -> str:
    """B8/B9 verdict convention: a wrong-signed significant result says so."""
    if res is None:
        return "NO DATA"
    lo, hi = res["ci95"]
    if direction > 0:
        if lo > 0:
            return "MET"
        return ("NOT MET - significant in the WRONG direction" if hi < 0
                else "NOT MET - consistent with zero")
    if hi < 0:
        return "MET"
    return ("NOT MET - significant in the WRONG direction" if lo > 0
            else "NOT MET - consistent with zero")


def run_scorer(tag: str, S: np.ndarray, sites: dict) -> dict:
    rng = np.random.default_rng(SEED)
    n_prot = len(sites["prot_names"])
    res: dict = {"coverage": {}}

    for probe in PROBES:
        v = S[:, AIDX[probe]].astype(np.float64)
        res["coverage"][f"C_to_{probe}"] = {
            "sites_with_a_score": int(np.isfinite(v).sum()),
            "sites_total": int(len(v)),
            "bonded_with_a_score": int(np.isfinite(v)[sites["bonded"] == 1].sum()),
            "free_with_a_score": int(np.isfinite(v)[sites["bonded"] == 0].sum()),
        }


    base = stratum_ids(sites, ("comp",))
    for probe in PROBES:
        v = S[:, AIDX[probe]].astype(np.float64)
        res[f"primary_C_to_{probe}"] = matched_gap(
            sites["prot"], base, sites["bonded"], v, n_prot, rng)

    v_s = S[:, AIDX["S"]].astype(np.float64)


    keep = (sites["bonded"] == 0) | (sites["direct"] == 1)
    res["S2_direct_evidence_only"] = matched_gap(
        sites["prot"][keep], base[keep], sites["bonded"][keep],
        v_s[keep], n_prot, rng)


    hyd = stratum_ids(sites, ("comp", "hyd_bin"))
    res["S3_burial_proxy_matched"] = matched_gap(
        sites["prot"], hyd, sites["bonded"], v_s, n_prot, rng)


    keep4 = (sites["bonded"] == 1) | (sites["nbonds"] >= 2)
    res["S4_free_set_well_characterised_proteins"] = matched_gap(
        sites["prot"][keep4], base[keep4], sites["bonded"][keep4],
        v_s[keep4], n_prot, rng)


    res["by_compartment"] = {}
    for cname, cid in CIDX.items():
        m = sites["comp"] == cid
        if m.sum() == 0:
            continue
        res["by_compartment"][cname] = {
            "n_bonded_sites": int((sites["bonded"][m] == 1).sum()),
            "n_free_sites": int((sites["bonded"][m] == 0).sum()),
            "gap": matched_gap(sites["prot"][m], base[m], sites["bonded"][m],
                               v_s[m], n_prot, rng),
        }


    res["P2_pairing_C_to_S"] = pairing_test(sites, v_s, rng)
    res["P2_proximity_control_post_hoc"] = proximity_control(sites, v_s, rng)


    with warnings.catch_warnings():
        # An unscored alternative set legitimately produces an empty row.
        warnings.simplefilter("ignore", category=RuntimeWarning)
        mean19 = np.nanmean(np.where(np.arange(20) == AIDX["C"], np.nan, S), axis=1)
    res["diagnostic_all19_mean_gap"] = {
        "gap": matched_gap(sites["prot"], base, sites["bonded"],
                           mean19.astype(np.float64), n_prot, rng),
        "how_to_read_this":
            "the same matched gap computed on the mean damage over all 19 "
            "substitutions instead of C->S alone. It is NOT a specificity "
            "test and no prediction was declared for it. At a bonded cysteine "
            "EVERY non-Cys substitution removes the bond, so a real bond-loss "
            "effect inflates this quantity too. A gap here of similar size to "
            "the C->S gap is therefore consistent BOTH with generic positional "
            "importance AND with uniform bond loss across substitutions, and "
            "cannot distinguish them. It is reported because it tells the "
            "reader how much of the C->S gap is specific to that substitution "
            "rather than present across the board.",
    }
    return res


def verdicts(res: dict) -> dict:
    p1 = directional(res.get("primary_C_to_S"), +1)
    pair = res.get("P2_pairing_C_to_S") or {}
    if "observed_mean_within_bond_abs_diff" not in pair:
        p2 = "NO DATA"
    else:
        lower = pair["difference_observed_minus_null_mean"] < 0
        sig = pair["permutation_p_one_sided_observed_lower"] <= 0.05
        if lower and sig:
            p2 = "MET"
        elif lower:
            p2 = "NOT MET - lower than the null but not significantly so"
        else:
            p2 = ("NOT MET - HIGHER than the null, the opposite of the "
                  "declared direction")
    prox = res.get("P2_proximity_control_post_hoc") or {}
    if "difference_true_minus_control" not in prox:
        p2p = "NO DATA"
    else:
        lo, hi = prox["ci95_of_difference"]
        if hi < 0:
            p2p = "MET"
        elif lo > 0:
            p2p = ("NOT MET - HIGHER than the separation-matched control, the "
                   "opposite of the declared direction")
        else:
            p2p = "NOT MET - consistent with zero"
    return {
        "P1_primary_magnitude": p1,
        "P2_primary_pairing": p2,
        "P2_proximity_control_post_hoc": p2p,
        "S1_secondary_C_to_A": directional(res.get("primary_C_to_A"), +1),
        "S2_secondary_direct_evidence": directional(
            res.get("S2_direct_evidence_only"), +1),
        "S3_secondary_burial_proxy_match": directional(
            res.get("S3_burial_proxy_matched"), +1),
        "S4_secondary_free_set_quality": directional(
            res.get("S4_free_set_well_characterised_proteins"), +1),
        "S5_control_compartment_placebo":
            "NOT AVAILABLE - declared unavailable before running; too few "
            "cytosolic bonded cysteines exist to form an arm",
    }



def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t_start = time.time()
    log("B10 disulfide bond loss -- label-free matched-site rule audit")
    log("declared predictions are fixed in DECLARED and written to JSON verbatim")

    log("loading topology census")
    census = read_topology(ROOT / "data" / "raw" / "topology_viewer.html")
    log(f"  {len(census):,} proteins")

    log("loading cached UniProt feature table (no network fetch)")
    uni = load_uniprot()
    log(f"  {len(uni):,} reviewed human entries")

    log("enumerating eligible cysteine sites")
    sites = enumerate_sites(census, uni)
    index = build_index(sites)

    free_positive_evidence = {
        "sites_with_positive_free_thiol_evidence": 0,
        "sites_merely_unannotated": int(sites["counts"].get("free_sites", 0)),
        "explanation":
            "A cysteine that carries positive evidence of a free thiol -- a "
            "modified-residue, lipidation or active-site annotation -- has a "
            "role that is not 'ordinary free structural cysteine', so those "
            "sites are EXCLUDED from both arms by the eligibility rule rather "
            "than used as free controls. That leaves no positively evidenced "
            "free control set at all: every site in the free arm is merely "
            "unannotated. An independent earlier pass over the same source "
            "reached the same conclusion (results/disulfide_rule_20261025/"
            "inventory_summary.json, 'verified_free_controls_in_import': 0).",
        "counts_of_the_excluded_positively_evidenced_cysteines": {
            k: v for k, v in sites["counts"].items() if k.startswith("excluded_")
        },
    }

    results: dict = {}
    log("AlphaMissense")
    S_am = load_or_extract("am", sites, index, extract_am)
    results["AlphaMissense"] = run_scorer("AlphaMissense", S_am, sites)
    del S_am

    def assemble(status: str) -> dict:
        v = {k: verdicts(r) for k, r in results.items()}
        not_met = {s: {k: val for k, val in d.items()
                       if str(val).startswith(("NOT MET", "NO DATA"))}
                   for s, d in v.items()}
        return {
            "rule": "disulfide bond loss",
            "script": "scripts/beat_am/b10_disulfide_rule.py",
            "run_finished": time.strftime("%Y-%m-%d %H:%M:%S"),
            "status": status,
            "declared_before_running": DECLARED,
            "label_free":
                "Oxford topology census (sequence, TM boundaries, loop sides), "
                "UniProt sequence and feature annotation, and the scorers' own "
                "outputs. No ClinVar, humsavar, gnomAD, HGMD or deep mutational "
                "scanning data was read. UniProt feature annotation is sequence "
                "and structure annotation, not clinical data.",
            "why_not_the_within_position_estimator":
                "B8 and B9 subtract, at one position, one class of substitution "
                "from another, which cancels per-position effects exactly. That "
                "is impossible for this rule: at a bonded cysteine EVERY "
                "ordinary non-Cys replacement removes the linkage, so there is "
                "no motif-preserving comparator to subtract against and the "
                "within-position mean already contains the bond-loss penalty in "
                "every term. Reusing that estimator would subtract the quantity "
                "being measured and report approximately zero for a scorer that "
                "had learned the chemistry perfectly. No comparator was invented "
                "to avoid this. The matched-site design used instead is weaker, "
                "and the cost is the burial/conservation confound named in "
                "limitations.",
            "design": {
                "primary_statistic":
                    "the RAW damage score for C->S at each eligible cysteine, "
                    "not a contrast. C->A is run separately as S1 and the two "
                    "are never averaged.",
                "primary_comparison":
                    "bonded versus free cysteines within protein x compartment "
                    "strata, stratum weight w = n1*n0/(n1+n0), bootstrapped "
                    f"over WHOLE PROTEINS, {N_BOOT} draws, seed {SEED}",
                "bond_partner_dependence":
                    "every retained bond is intrachain, so both endpoints "
                    "belong to the same protein and a whole-protein resample "
                    "carries them together automatically -- either both are in "
                    "a draw or neither is. Interchain bonds are excluded for "
                    "exactly this reason: no protein-level resample could "
                    "preserve the dependence of endpoints in different "
                    "molecules. Bonds with one endpoint excluded on role "
                    "grounds are counted as orphaned endpoints and the eligible "
                    "endpoint is kept in the bonded arm.",
                "compartment":
                    "from the census topology: loop location 'inside' is "
                    "cytosolic, 'outside' is non-cytosolic, positions inside a "
                    "TM segment form their own compartment. Always a stratum "
                    "variable, so a bonded cysteine is only ever compared "
                    "against free cysteines on the same side.",
                "bonded_definition":
                    "carries a DISULFID feature with a defined intrachain "
                    "partner; interchain, singleton and conflicting-partner "
                    "features excluded as uncertain",
                "free_definition":
                    "carries NO disulfide annotation. THIS IS AN UNKNOWN, NOT "
                    "EVIDENCE OF A FREE THIOL. See free_set_assumption.",
                "excluded_from_every_arm":
                    "cysteines with a binding-site, active-site or "
                    "modified-residue annotation, and cysteines inside a signal "
                    "peptide. Removed from EVERY arm, not only from the arm "
                    "being scored.",
                "numbering_reconciliation":
                    "only accessions whose census sequence is byte-identical to "
                    "the UniProt sequence are used; every annotated disulfide "
                    "position is checked to be a cysteine and mismatches are "
                    "counted, not dropped silently",
                "burial_proxy_for_S3":
                    f"mean Kyte-Doolittle hydropathy over a +/-{HYDRO_WINDOW} "
                    f"residue window, binned into {N_HYDRO_BINS} quartiles of "
                    "the eligible pool",
                "pairing_test":
                    "mean over bonds of |score(C->S, p) - score(C->S, q)| "
                    "against a null that shuffles the same endpoint scores "
                    "within each protein and re-pairs them consecutively, "
                    f"{N_PERM} draws; proteins need >= "
                    f"{MIN_BONDS_FOR_PAIRING} complete scored bonds",
            },
            "counts": sites["counts"],
            "bond_evidence_breakdown_endpoints": sites["evidence_breakdown"],
            "free_set_assumption": free_positive_evidence,
            "hydrophobicity_quartile_edges": sites["hyd_quartile_edges"],
            "uniprot_source":
                "data/external/disulfide_rule_20261025/human_reviewed_features.tsv, "
                "fetched from the UniProt REST API by c1_disulfide_sources.py "
                "with URL and SHA-256 recorded in that directory's "
                "source_manifest.json. Not re-fetched on a rerun.",
            "cross_scorer_comparison_caveat":
                "AlphaMissense pathogenicity and ESM1b damage (= -LLR) have "
                "different scales and different spreads. Standardised (SD-unit) "
                "values are NOT comparable between the two scorers and no such "
                "comparison is made here. Standardised values appear only to "
                "make effect sizes readable within one scorer. Cross-scorer "
                "claims concern DIRECTION and pass/fail pattern only.",
            "limitations": LIMITATIONS,
            "post_hoc_observations": POST_HOC,
            "results": results,
            "verdicts": v,
            "declared_expectations_not_met":
                {s: d for s, d in not_met.items() if d} or
                "none - every declared expectation was met",
            "note": "This script measures only. No correction parameter was "
                    "fitted and nothing was trained.",
        }

    (OUT / "disulfide_rule_result.json").write_text(
        json.dumps(assemble("AlphaMissense complete; ESM1b pending"),
                   indent=2, default=str), encoding="utf-8")
    log("  wrote interim JSON (AlphaMissense only)")

    log("ESM1b")
    try:
        S_esm = load_or_extract("esm", sites, index, extract_esm)
        results["ESM1b"] = run_scorer("ESM1b", S_esm, sites)
        del S_esm
        status = "complete"
    except Exception as exc:                      # noqa: BLE001
        log(f"  ESM1b FAILED: {exc!r}")
        status = f"AlphaMissense complete; ESM1b failed: {exc!r}"

    summary = assemble(status)
    (OUT / "disulfide_rule_result.json").write_text(
        json.dumps(summary, indent=2, default=str), encoding="utf-8")

    log("")
    log("VERDICTS")
    for scorer, d in summary["verdicts"].items():
        log(f"  {scorer}")
        for k, val in d.items():
            log(f"    {k:34s} {val}")
    log(f"wrote {OUT / 'disulfide_rule_result.json'}")
    log(f"total runtime {time.time() - t_start:.0f}s")
    (OUT / "run.log").write_text("\n".join(_LOG_LINES) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

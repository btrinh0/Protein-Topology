"""Audit disulfide pairing with domain-matched cysteine controls."""

from __future__ import annotations

import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "scripts" / "beat_am"))

import b10_disulfide_rule as b10                   # noqa: E402
from run_gate1 import read_topology                # noqa: E402

OUT = ROOT / "results" / "beat_am_20261016" / "conservation_confound"
DR = ROOT / "results" / "beat_am_20261016" / "disulfide_rule"

AIDX = b10.AIDX
SEED = b10.SEED
N_BOOT = b10.N_BOOT
N_PERM = b10.N_PERM
MIN_BONDS = b10.MIN_BONDS_FOR_PAIRING

PRIOR_BOUNDS = {                       # The E3 baseline comes from B10's frozen summary and pairing controls.
    "AlphaMissense": {"constraint_matched_share": 0.21,
                      "separation_matched_share": 0.82},
    "ESM1b": {"constraint_matched_share": 0.02,
              "separation_matched_share": 0.60},
}

DECLARED = {
    "rule_under_scrutiny": "B10, disulfide bond loss -- the P2 pairing test",
    "why_this_control": "B10C's constraint proxy is derived from the audited scorer's "
                        "own outputs and still contains bond-loss signal, so it "
                        "over-corrects and gives only a lower bound. UniProt DOMAIN "
                        "annotation is independent of both scorers and carries no "
                        "bond-loss signal.",
    "E1": "true bond partners are far more likely to share an annotated domain than "
          "randomly re-paired bonded cysteines of the same protein. "
          "Confound-friendly; expected to hold; established from annotation alone.",
    "E2": "among intra-domain pairs only, true bonds still show a lower mean "
          "|difference| than intra-domain non-partner controls, 95% interval "
          "excluding zero",
    "E3": "the share of the subset P2 gap surviving domain matching falls between "
          "B10C's constraint-matched 21% and B10's separation-matched 82% for "
          "AlphaMissense",
    "E4": "no prediction declared for ESM1b; reported either way",
    "cross_domain_arm": "declared NOT AVAILABLE before running: only 26 of 7688 true "
                        "bonds are inter-domain against 5500 intra-domain",
    "failure_policy": "any declared item that fails is reported as a failure and is "
                      "not reinterpreted",
}

_LOG: list[str] = []


def log(m: str) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {m}"
    _LOG.append(line)
    print(line, flush=True)


def ci(draws: np.ndarray) -> list[float]:
    lo, hi = np.nanquantile(draws, [0.025, 0.975])
    return [float(lo), float(hi)]



def domain_intervals(uni: dict) -> dict:
    """Annotated DOMAIN intervals per accession, using B10's own parser."""
    out = {}
    for acc, rec in uni.items():
        iv = [(s, e) for s, e, _d in
              b10.features(rec.get("Domain [FT]") or "", "DOMAIN")]
        if iv:
            out[acc] = iv
    return out


def domain_of(dom: dict, acc: str, pos: int):
    """Index of the annotated domain containing pos, or None."""
    for i, (s, e) in enumerate(dom.get(acc, ())):
        if s <= pos <= e:
            return i
    return None


def collect(sites: dict, dom: dict, v: np.ndarray):
    """Per protein: intra-domain true bonds and intra-domain candidate controls."""
    row_of = {(a, p): i for i, (a, p) in
              enumerate(zip(sites["acc"].tolist(), sites["pos"].tolist()))}
    per_prot = []
    counts = defaultdict(int)
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
        if len(pairs) < MIN_BONDS:
            continue

        truth = {frozenset((s, e)) for s, e, _r, _r2 in pairs}
        true_intra = []
        for s, e, rs, re_ in pairs:
            ds, de = domain_of(dom, acc, s), domain_of(dom, acc, e)
            if ds is None or de is None:
                counts["true_outside_domain"] += 1
            elif ds == de:
                counts["true_intra_domain"] += 1
                true_intra.append((s, e, rs, re_, ds))
            else:
                counts["true_inter_domain"] += 1

        positions = sorted(pts)
        cand_intra = []
        for i, a in enumerate(positions):
            for b in positions[i + 1:]:
                if frozenset((a, b)) in truth:
                    continue
                da, db = domain_of(dom, acc, a), domain_of(dom, acc, b)
                if da is None or db is None:
                    continue
                if da == db:
                    counts["cand_intra_domain"] += 1
                    cand_intra.append((a, b, da))
                else:
                    counts["cand_inter_domain"] += 1

        cand_all = []
        for i, a in enumerate(positions):
            for b in positions[i + 1:]:
                if frozenset((a, b)) in truth:
                    continue
                da, db = domain_of(dom, acc, a), domain_of(dom, acc, b)
                cand_all.append((a, b, da if (da is not None and da == db) else None))

        if true_intra and cand_intra:
            per_prot.append((acc, true_intra, cand_intra, pts, cand_all))
    return per_prot, dict(counts)


def matched_on_key(per_prot, v, k, key: str) -> dict:
    """Same intra-domain true bonds, same FULL candidate pool, different key.

    Added so that the separation-, constraint- and domain-matched controls can be
    compared on one denominator and one set of bonds. B10's separation control and
    B10C's constraint control were each computed over ALL bonds, whose P2 gap is
    0.01445, while the domain control runs on the intra-domain subset, whose gap is
    0.01208. Quoting the three shares side by side without this is wrong, so the
    comparison table in the write-up uses these figures instead.
    """
    rows = []
    for _acc, true_intra, _cand_intra, pts, cand_all in per_prot:
        c_sep = np.array([abs(a - b) for a, b, _d in cand_all], dtype=np.float64)
        c_k = np.array([abs(k[pts[a]] - k[pts[b]]) for a, b, _d in cand_all],
                       dtype=np.float64)
        used: set[int] = set()
        s_true = s_ctrl = 0.0
        n = 0
        for s, e, rs, re_, _d in true_intra:
            target = abs(s - e) if key == "separation" else abs(k[rs] - k[re_])
            dist = np.abs((c_sep if key == "separation" else c_k) - target)
            order = np.argsort(dist, kind="stable")
            pick = next((int(i) for i in order if int(i) not in used),
                        int(order[0]))
            used.add(pick)
            a, b, _dd = cand_all[pick]
            s_true += abs(v[rs] - v[re_])
            s_ctrl += abs(v[pts[a]] - v[pts[b]])
            n += 1
        rows.append((s_true, s_ctrl, float(n)))
    A = np.asarray(rows, dtype=np.float64)
    tot = A[:, 2].sum()
    est_true, est_ctrl = A[:, 0].sum() / tot, A[:, 1].sum() / tot
    rng = np.random.default_rng(SEED)
    draws = np.empty(N_BOOT)
    for b in range(N_BOOT):
        sub = A[rng.integers(0, len(A), size=len(A))]
        dd = sub[:, 2].sum()
        draws[b] = (sub[:, 0].sum() - sub[:, 1].sum()) / dd if dd else np.nan
    return {"matching_key": key,
            "candidate_pool": "all non-bond pairs of bonded cysteines in the protein",
            "true_bond_mean_abs_diff": float(est_true),
            "matched_control_mean_abs_diff": float(est_ctrl),
            "difference_true_minus_control": float(est_true - est_ctrl),
            "ci95_of_difference": ci(draws),
            "n_bonds": int(tot), "n_proteins": int(len(A))}


def premise(per_prot_all) -> dict:
    """E1: intra-domain rate for true bonds vs within-protein random re-pairing."""
    obs_intra = obs_n = 0
    null_tot = np.zeros(N_PERM)
    null_n = 0
    rng = np.random.default_rng(SEED)
    per_prot_arrays = []
    for acc, dom_of_pos, bond_list in per_prot_all:
        for a, b in bond_list:
            da, db = dom_of_pos.get(a), dom_of_pos.get(b)
            obs_n += 1
            if da is not None and db is not None and da == db:
                obs_intra += 1
        per_prot_arrays.append(
            np.array([(p, dom_of_pos.get(p) if dom_of_pos.get(p) is not None else -1)
                      for p in sorted(dom_of_pos)], dtype=np.int64))
        null_n += len(bond_list)
    for d in range(N_PERM):
        tot = 0
        for arr in per_prot_arrays:
            idx = rng.permutation(len(arr))
            dd = arr[idx, 1]
            a, b = dd[0::2], dd[1::2]
            tot += int(np.sum((a == b) & (a >= 0)))
        null_tot[d] = tot
    obs_rate = obs_intra / obs_n
    null_rate = null_tot / null_n
    n_ge = int((null_rate >= obs_rate).sum())
    return {
        "observed_intra_domain_rate_true_bonds": float(obs_rate),
        "null_intra_domain_rate_random_repairing": float(null_rate.mean()),
        "null_ci95": ci(null_rate),
        "ratio_observed_over_null": float(obs_rate / null_rate.mean())
        if null_rate.mean() else float("nan"),
        "permutation_p_one_sided_observed_higher": (1 + n_ge) / (1 + N_PERM),
        "n_true_bonds_scored": int(obs_n),
        "reading": "a large ratio means partners are co-located, so P2's permutation "
                   "null -- which re-pairs across domains -- does NOT hold the "
                   "structural environment fixed",
    }


def subset_p2(per_prot, v) -> dict:
    """P2 recomputed on the intra-domain true bonds only, for a matching denominator."""
    obs_sum = n = 0.0
    flat = []
    for _acc, true_intra, _cand, _pts, _all in per_prot:
        vals = []
        for _s, _e, rs, re_, _d in true_intra:
            obs_sum += abs(v[rs] - v[re_])
            n += 1
            vals.extend((v[rs], v[re_]))
        flat.append(np.array(vals, dtype=np.float64))
    observed = obs_sum / n
    rng = np.random.default_rng(SEED)
    null = np.empty(N_PERM)
    for d in range(N_PERM):
        tot = 0.0
        for arr in flat:
            rng.shuffle(arr)
            tot += np.abs(arr[0::2] - arr[1::2]).sum()
        null[d] = tot / n
    n_le = int((null <= observed).sum())
    return {"observed_mean_within_bond_abs_diff": float(observed),
            "null_mean": float(null.mean()), "null_ci95": ci(null),
            "gap_null_minus_observed": float(null.mean() - observed),
            "permutation_p_one_sided_observed_lower": (1 + n_le) / (1 + N_PERM),
            "n_bonds": int(n), "n_proteins": len(per_prot)}


def domain_matched(per_prot, v, prefer_separation: bool) -> dict:
    """Intra-domain true bonds against intra-domain non-partner control pairs."""
    rows = []
    n_same_instance = n_total = 0
    for _acc, true_intra, cand, pts, _all in per_prot:
        c_sep = np.array([abs(a - b) for a, b, _d in cand], dtype=np.float64)
        c_dom = np.array([d for _a, _b, d in cand], dtype=np.int64)
        used: set[int] = set()
        s_true = s_ctrl = sep_true = sep_ctrl = 0.0
        n = 0
        for s, e, rs, re_, d in true_intra:
            same_inst = np.flatnonzero(c_dom == d)
            pool = same_inst if same_inst.size else np.arange(len(cand))
            if prefer_separation:
                order = pool[np.argsort(np.abs(c_sep[pool] - abs(s - e)),
                                        kind="stable")]
            else:
                order = pool                    # The residues are already ordered by canonical position.
            pick = next((int(i) for i in order if int(i) not in used),
                        int(order[0]))
            used.add(pick)
            a, b, dd = cand[pick]
            n_total += 1
            if dd == d:
                n_same_instance += 1
            s_true += abs(v[rs] - v[re_])
            s_ctrl += abs(v[pts[a]] - v[pts[b]])
            sep_true += abs(s - e)
            sep_ctrl += abs(a - b)
            n += 1
        rows.append((s_true, s_ctrl, float(n), sep_true, sep_ctrl))

    A = np.asarray(rows, dtype=np.float64)
    tot = A[:, 2].sum()
    est_true, est_ctrl = A[:, 0].sum() / tot, A[:, 1].sum() / tot
    rng = np.random.default_rng(SEED)
    draws = np.empty(N_BOOT)
    for b in range(N_BOOT):
        sub = A[rng.integers(0, len(A), size=len(A))]
        dd = sub[:, 2].sum()
        draws[b] = (sub[:, 0].sum() - sub[:, 1].sum()) / dd if dd else np.nan
    return {
        "control_definition": ("intra-domain non-bond pair, same domain instance "
                               "preferred, closest sequence separation"
                               if prefer_separation else
                               "intra-domain non-bond pair, same domain instance "
                               "preferred, first by position order"),
        "declared_direction_before_computing":
            "true-bond mean |difference| BELOW the domain-matched control",
        "true_bond_mean_abs_diff": float(est_true),
        "domain_matched_control_mean_abs_diff": float(est_ctrl),
        "difference_true_minus_control": float(est_true - est_ctrl),
        "ci95_of_difference": ci(draws),
        "mean_sequence_separation_true_bonds": float(A[:, 3].sum() / tot),
        "mean_sequence_separation_control_pairs": float(A[:, 4].sum() / tot),
        "share_of_controls_in_the_same_domain_instance":
            n_same_instance / n_total if n_total else float("nan"),
        "n_bonds": int(tot),
        "n_proteins": int(len(A)),
    }



def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    log("B10D -- domain-matched control for B10's pairing test")
    log("declared predictions are fixed in DECLARED and written verbatim")

    log("loading topology census and B10's cached UniProt features")
    census = read_topology(ROOT / "data" / "raw" / "topology_viewer.html")
    uni = b10.load_uniprot()
    dom = domain_intervals(uni)
    log(f"  {len(census):,} census proteins, {len(uni):,} UniProt entries, "
        f"{len(dom):,} with DOMAIN annotation")
    log("re-enumerating B10's sites (same seed, same code path)")
    sites = b10.enumerate_sites(census, uni)

    results: dict[str, dict] = {}
    for tag, fname in (("AlphaMissense", "position_scores_am.npz"),
                       ("ESM1b", "position_scores_esm.npz")):
        path = DR / fname
        if not path.exists():
            log(f"{tag}: cache {fname} missing -- skipped")
            continue
        log("")
        log(f"{tag}: loading cached scores from {fname} (no new pass)")
        z = np.load(path, allow_pickle=False)
        if not (np.array_equal(z["acc"], sites["acc"])
                and np.array_equal(z["pos"], sites["pos"])):
            raise SystemExit(f"{tag} cache rows do not match re-enumerated sites")
        S = z["scores"]
        v = S[:, AIDX["S"]].astype(np.float64)
        kcols = [AIDX[a] for a in b10.ALPHABET if a not in ("S", "C")]
        ksub = S[:, kcols].astype(np.float64)
        with np.errstate(invalid="ignore"):
            k = np.where(np.sum(np.isfinite(ksub), axis=1) >= 15,
                         np.nanmean(np.where(np.isfinite(ksub), ksub, np.nan),
                                    axis=1), np.nan)

        per_prot, counts = collect(sites, dom, v)
        log(f"  true bonds: intra-domain {counts.get('true_intra_domain', 0):,}, "
            f"inter-domain {counts.get('true_inter_domain', 0):,}, "
            f"outside any domain {counts.get('true_outside_domain', 0):,}")
        log(f"  candidate non-bond pairs: intra-domain "
            f"{counts.get('cand_intra_domain', 0):,}, inter-domain "
            f"{counts.get('cand_inter_domain', 0):,}")
        log(f"  proteins usable for the intra-domain arm: {len(per_prot):,}")

        R: dict[str, object] = {"counts": counts}

        # Define E1 eligibility from annotations alone.
        prem_input = []
        for acc, bonds in sites["bonds_by_prot"].items():
            pts = {}
            bl = []
            for s, e, _d in bonds:
                pts[s] = domain_of(dom, acc, s)
                pts[e] = domain_of(dom, acc, e)
                bl.append((s, e))
            if len(bl) >= MIN_BONDS:
                prem_input.append((acc, pts, bl))
        pr = premise(prem_input)
        R["E1_premise_partners_share_a_domain"] = pr
        log(f"  E1 intra-domain rate: true bonds "
            f"{pr['observed_intra_domain_rate_true_bonds']:.4f} vs random re-pairing "
            f"{pr['null_intra_domain_rate_random_repairing']:.4f} "
            f"({pr['ratio_observed_over_null']:.1f}x), "
            f"p = {pr['permutation_p_one_sided_observed_higher']:.4f}")

        sp2 = subset_p2(per_prot, v)
        R["subset_P2_on_intra_domain_bonds"] = sp2
        log(f"  subset P2 (intra-domain bonds only): observed "
            f"{sp2['observed_mean_within_bond_abs_diff']:.5f} vs null "
            f"{sp2['null_mean']:.5f}, gap {sp2['gap_null_minus_observed']:.5f}, "
            f"{sp2['n_bonds']:,} bonds")

        for key, pref in (("primary_separation_tiebreak", True),
                          ("sensitivity_position_order", False)):
            mc = domain_matched(per_prot, v, pref)
            gap = sp2["gap_null_minus_observed"]
            mc["share_of_subset_P2_gap_surviving"] = (
                float(abs(mc["difference_true_minus_control"]) / gap)
                if gap else float("nan"))
            R[key] = mc
            log(f"  {key:30s} true {mc['true_bond_mean_abs_diff']:.5f} vs control "
                f"{mc['domain_matched_control_mean_abs_diff']:.5f}, diff "
                f"{mc['difference_true_minus_control']:+.5f} "
                f"[{mc['ci95_of_difference'][0]:+.5f}, "
                f"{mc['ci95_of_difference'][1]:+.5f}]  "
                f"({mc['share_of_subset_P2_gap_surviving'] * 100:.0f}% of subset gap)  "
                f"sep true {mc['mean_sequence_separation_true_bonds']:.1f} vs ctrl "
                f"{mc['mean_sequence_separation_control_pairs']:.1f}, "
                f"same-instance {mc['share_of_controls_in_the_same_domain_instance']:.2f}")

        # Evaluate all three controls on the same bonds and denominator.
        log("  comparable arms on the intra-domain subset "
            "(same bonds, same candidate pool, different key):")
        comp = {}
        for key in ("separation", "constraint"):
            mk = matched_on_key(per_prot, v, k, key)
            gap = sp2["gap_null_minus_observed"]
            mk["share_of_subset_P2_gap_surviving"] = (
                float(abs(mk["difference_true_minus_control"]) / gap)
                if gap else float("nan"))
            comp[key] = mk
            log(f"    {key:12s} true {mk['true_bond_mean_abs_diff']:.5f} vs control "
                f"{mk['matched_control_mean_abs_diff']:.5f}, diff "
                f"{mk['difference_true_minus_control']:+.5f} "
                f"[{mk['ci95_of_difference'][0]:+.5f}, "
                f"{mk['ci95_of_difference'][1]:+.5f}]  "
                f"({mk['share_of_subset_P2_gap_surviving'] * 100:.0f}% of subset gap)")
        dm = R["primary_separation_tiebreak"]
        comp["domain"] = {
            "matching_key": "domain",
            "difference_true_minus_control": dm["difference_true_minus_control"],
            "ci95_of_difference": dm["ci95_of_difference"],
            "share_of_subset_P2_gap_surviving":
                dm["share_of_subset_P2_gap_surviving"],
            "note": "the primary arm above; its candidate pool is restricted to "
                    "intra-domain pairs, which IS the matching",
        }
        log(f"    {'domain':12s} "
            f"({comp['domain']['share_of_subset_P2_gap_surviving'] * 100:.0f}% "
            f"of subset gap, from the primary arm)")
        R["comparable_controls_on_intra_domain_subset"] = {
            "why": "B10's separation control and B10C's constraint control were each "
                   "computed over ALL bonds (P2 gap 0.01445) while the domain control "
                   "runs on the intra-domain subset (gap 0.01208). These arms put all "
                   "three on the same bonds and the same denominator so the shares "
                   "are comparable.",
            "subset_P2_gap": sp2["gap_null_minus_observed"],
            "arms": comp,
        }

        R["cross_domain_arm"] = {
            "status": "NOT AVAILABLE -- declared so before running",
            "n_true_inter_domain_bonds": counts.get("true_inter_domain", 0),
            "n_true_intra_domain_bonds": counts.get("true_intra_domain", 0),
            "why": "too few inter-domain true bonds to carry a comparison; no number "
                   "was computed from this arm",
        }
        results[tag] = R


    verdicts: dict[str, dict] = {}
    for tag, R in results.items():
        pr = R["E1_premise_partners_share_a_domain"]
        e1 = ("MET" if (pr["ratio_observed_over_null"] > 1
                        and pr["permutation_p_one_sided_observed_higher"] < 0.05)
              else "NOT MET")
        mc = R["primary_separation_tiebreak"]
        lo, hi = mc["ci95_of_difference"]
        e2 = ("MET" if (mc["difference_true_minus_control"] < 0 and hi < 0)
              else "NOT MET")
        share = mc["share_of_subset_P2_gap_surviving"]
        b = PRIOR_BOUNDS[tag]
        e3 = ("MET" if (b["constraint_matched_share"] < share
                        < b["separation_matched_share"]) else "NOT MET")
        verdicts[tag] = {
            "E1_partners_share_an_annotated_domain": e1,
            "E1_detail": {
                "true_bond_rate": pr["observed_intra_domain_rate_true_bonds"],
                "random_repairing_rate":
                    pr["null_intra_domain_rate_random_repairing"],
                "ratio": pr["ratio_observed_over_null"],
                "consequence": "P2's permutation null re-pairs across domains, so it "
                               "does NOT hold the structural environment fixed",
            },
            "E2_P2_survives_domain_matching": e2,
            "E2_detail": {"difference": mc["difference_true_minus_control"],
                          "ci95": [lo, hi]},
            "E3_share_between_constraint_and_separation_bounds": e3
            if tag == "AlphaMissense" else "no prediction declared (E4)",
            "E3_detail": {
                "share_surviving_domain_matching": share,
                "share_surviving_constraint_matching_b10c":
                    b["constraint_matched_share"],
                "share_surviving_separation_matching_b10":
                    b["separation_matched_share"],
            },
            "sensitivity_agrees_with_primary": bool(
                (R["sensitivity_position_order"]["difference_true_minus_control"] < 0)
                == (mc["difference_true_minus_control"] < 0)),
        }

    summary = {
        "script": "scripts/beat_am/b10d_pairing_domain_control.py",
        "rule_under_scrutiny": "B10 disulfide bond loss -- the P2 pairing test",
        "run_finished": time.strftime("%Y-%m-%d %H:%M:%S"),
        "declared_before_running": DECLARED,
        "why_annotation_beats_the_score_derived_proxy":
            "UniProt DOMAIN annotation is independent of both audited scorers and "
            "carries no bond-loss signal, so unlike B10C's constraint proxy it cannot "
            "over-correct. It is an estimate rather than a bound.",
        "domain_source": "Domain [FT] column of "
                         "data/external/disulfide_rule_20261025/"
                         "human_reviewed_features.tsv, parsed with b10.features(); "
                         "the response B10 had already cached. No network request.",
        "label_free": "topology census, cached UniProt features, and the scorers' own "
                      "cached outputs. No ClinVar, humsavar, gnomAD, HGMD or DMS read.",
        "results": results,
        "verdicts": verdicts,
        "note": "measurement only; nothing was fitted, trained or tuned",
    }
    (OUT / "b10d_pairing_domain_result.json").write_text(
        json.dumps(summary, indent=2, default=float), encoding="utf-8")

    log("")
    for tag, vv in verdicts.items():
        log(f"  {tag}")
        for k, val in vv.items():
            if isinstance(val, str):
                log(f"    {k:56s} {val}")
    (OUT / "b10d_run.log").write_text("\n".join(_LOG) + "\n", encoding="utf-8")
    log(f"wrote {OUT / 'b10d_pairing_domain_result.json'}")


if __name__ == "__main__":
    main()

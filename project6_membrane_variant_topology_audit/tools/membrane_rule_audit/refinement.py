"""Annotation-driven additions to a predictor that consumes no predictor scores.

These are directional correction features, not clinical probabilities or a claim
that every mutation following a chemical pattern has the same functional effect.
Real coefficients must be calibrated and the complete predictions evaluated
later. In particular, the control-free disulfide arm tests predictive utility;
it does not identify a causal disulfide-specific penalty.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
NEUTRAL = "ACFGILMNPQSTVWY"
RULES = ("positive_inside", "charge_burial", "glycosylation_loss", "disulfide_loss")


@dataclass(frozen=True)
class FeatureSpec:
    decay_length: float = 8.0
    maximum_flank_distance: int = 60
    minimum_tm_length: int = 17
    annotation_evidence: str = "direct"
    disulfide_mode: str = "control_free_annotation_feature"

    def __post_init__(self):
        if not np.isfinite(self.decay_length) or self.decay_length <= 0:
            raise ValueError("A positive finite decay length is required")
        if self.maximum_flank_distance < 1 or self.minimum_tm_length < 3:
            raise ValueError("Invalid membrane window")
        if self.annotation_evidence not in {"direct", "all_annotated"}:
            raise ValueError("Declare direct or all_annotated evidence")
        if self.disulfide_mode != "control_free_annotation_feature":
            raise ValueError("The controlled disulfide audit is a separate interface")


def sequence_digest(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def topology_map(record):
    """Use the original supervised-table conventions, never panel 0/1 sides.

    Cytosolic=1, outer=2, TM=0, unknown=-1. Unknown is not cytosolic.
    Geometry remains an annotation proxy for membrane depth.
    """
    sequence = record["sequence"]
    out = {}
    tms = record.get("tmds") or []
    tm_positions = set()
    for index, segment in enumerate(tms, 1):
        start, end = int(segment["start"]), int(segment["end"])
        if not 1 <= start <= end <= len(sequence):
            raise ValueError("TM coordinates outside canonical sequence")
        for position in range(start, end + 1):
            if position in tm_positions:
                raise ValueError("Overlapping TM segments")
            tm_positions.add(position)
            depth = 1 - abs(position - (start + end) / 2) / ((end - start) / 2) if end > start else 0
            out[position] = dict(side=0, distance=0, in_tm=1, tmd_index=index,
                                 loop_len=0, tm_length=end-start+1, depth=float(depth))
    for key, loop in (record.get("loops") or {}).items():
        index = int(key)
        start, end = int(loop["start"]), int(loop["end"])
        if not 1 <= start <= end <= len(sequence):
            raise ValueError("Loop coordinates outside canonical sequence")
        side = {"inside": 1, "outside": 2}.get(str(loop.get("location", "")).strip().lower(), -1)
        bounds = []
        if index > 1 and index - 2 < len(tms):
            bounds.append(("left", int(tms[index - 2]["end"])))
        if index <= len(tms):
            bounds.append(("right", int(tms[index - 1]["start"])))
        for position in range(start, end + 1):
            if position in tm_positions or position in out:
                raise ValueError("Overlapping loop/TM annotations")
            distances = [position-b if direction == "left" else b-position for direction,b in bounds]
            distances = [x for x in distances if x > 0]
            if not distances:
                continue
            out[position] = dict(side=side, distance=min(distances), in_tm=0,
                                 tmd_index=index, loop_len=end-start+1, tm_length=0, depth=0.0)
    return out


def _exact_location(feature):
    start = feature.get("location", {}).get("start", {})
    end = feature.get("location", {}).get("end", {})
    if start.get("modifier") != "EXACT" or end.get("modifier") != "EXACT":
        return None
    try:
        return int(start["value"]), int(end["value"])
    except (KeyError, TypeError, ValueError):
        return None


def _direct(feature):
    return any(e.get("evidenceCode") == "ECO:0000269" or
               (e.get("evidenceCode") == "ECO:0007744" and e.get("source") == "PDB")
               for e in feature.get("evidences", []))


def prepare_record(accession, topology, annotation, *, spec=FeatureSpec()):
    """Parse source annotations without reading any variant-effect score.

    No annotation is interpreted as verified absence of a bond. The alternative
    disulfide design intentionally needs no free-cysteine bank.
    """
    sequence = topology["sequence"]
    if not sequence or any(a not in ALPHABET for a in sequence):
        raise ValueError("A canonical standard-amino-acid sequence is required")
    if annotation.get("primaryAccession") != accession:
        raise ValueError("Annotation accession mismatch")
    if annotation.get("sequence", {}).get("value") != sequence:
        raise ValueError("Annotation and model-reference sequences must match exactly")
    geometry = topology_map(topology)
    exclusions = Counter()
    banned = set()
    features = annotation.get("features", [])
    for feature in features:
        location = _exact_location(feature)
        if location and feature.get("type") in {"Signal", "Active site", "Binding site", "Modified residue"}:
            banned.update(range(location[0], location[1]+1))
    glycosites, bonds = [], []
    for feature in features:
        kind = feature.get("type")
        if kind not in {"Glycosylation", "Disulfide bond"}:
            continue
        location = _exact_location(feature)
        if location is None:
            exclusions["inexact_or_unknown_feature_location"] += 1
            continue
        start, end = location
        if not 1 <= start <= end <= len(sequence):
            exclusions["out_of_range_feature"] += 1
            continue
        direct = _direct(feature)
        if spec.annotation_evidence == "direct" and not direct:
            exclusions["not_direct_evidence"] += 1
            continue
        evidence = dict(direct_evidence=direct, evidences=feature.get("evidences", []),
                        description=feature.get("description", ""))
        if kind == "Glycosylation":
            if start != end or not feature.get("description", "").lower().startswith("n-linked"):
                exclusions["not_single_N_linked_site"] += 1
                continue
            motif = sequence[start-1:start+2]
            if len(motif) != 3 or motif[0] != "N" or motif[1] == "P" or motif[2] not in "ST":
                exclusions["noncanonical_or_absent_sequon"] += 1
                continue
            positions = range(start, start+3)
            if any(p in banned for p in positions) or any(geometry.get(p, {}).get("side") != 2 for p in positions):
                exclusions["glycosite_not_resolved_outer_mature_loop"] += 1
                continue
            glycosites.append(dict(position=start, motif=motif, **evidence))
        else:
            description = feature.get("description", "").lower()
            if start == end or sequence[start-1] != "C" or sequence[end-1] != "C":
                exclusions["not_two_canonical_cysteine_endpoints"] += 1
                continue
            if start in banned or end in banned or any(term in description for term in ("redox", "transient", "interchain", "inter-chain")):
                exclusions["competing_function_signal_or_nonstructural_bond"] += 1
                continue
            bonds.append(dict(start=start, end=end, role="annotated_bond_not_independently_functionally_validated", **evidence))
    # Duplicate source records are harmless; competing partners are not.
    glycosites = list({g["position"]: g for g in glycosites}.values())
    bonds = list({(b["start"],b["end"]): b for b in bonds}.values())
    partners = defaultdict(set)
    for bond in bonds:
        partners[bond["start"]].add(bond["end"])
        partners[bond["end"]].add(bond["start"])
    ambiguous = {p for p, q in partners.items() if len(q)>1}
    kept = [b for b in bonds if b["start"] not in ambiguous and b["end"] not in ambiguous]
    exclusions["ambiguous_partner_bonds"] += len(bonds)-len(kept)
    return dict(accession=accession, sequence=sequence, sequence_sha256=sequence_digest(sequence),
                topology={k:topology.get(k) for k in ("sequence", "tmds", "loops")},
                glycosites=glycosites, disulfides=kept, exclusions=dict(exclusions),
                feature_spec=asdict(spec), free_cysteine_controls_required=False,
                interpretation="annotation-driven refinement candidate; no effect-size or pathogenicity conclusion")


def variant_features(record, positions, natives, alternatives):
    """Return four named features, rejecting coordinate/native/identity mistakes."""
    positions = np.asarray(positions)
    natives, alternatives = np.asarray(natives).astype(str), np.asarray(alternatives).astype(str)
    if positions.ndim != 1 or natives.shape != positions.shape or alternatives.shape != positions.shape:
        raise ValueError("Variant columns must have equal one-dimensional shape")
    sequence = record["sequence"]
    if sequence_digest(sequence) != record["sequence_sha256"]:
        raise ValueError("Sequence hash mismatch")
    if not np.isfinite(positions).all() or np.any(positions != positions.astype(int)):
        raise ValueError("Positions must be finite integers")
    pos = positions.astype(int)
    if np.any(pos < 1) or np.any(pos > len(sequence)):
        raise ValueError("Variant outside canonical sequence")
    reference = np.array(list(sequence))[pos-1]
    if np.any(natives != reference) or not np.isin(alternatives, list(ALPHABET)).all():
        raise ValueError("Native mismatch or nonstandard alternative")
    if np.any(natives == alternatives):
        raise ValueError("Only missense substitutions are accepted")
    spec = FeatureSpec(**record["feature_spec"])
    geometry = topology_map(record["topology"])
    matrix = np.zeros((len(pos), len(RULES)), dtype=float)
    bonded = {b[k] for b in record["disulfides"] for k in ("start", "end")}
    glyco = defaultdict(list)
    for site in record["glycosites"]:
        for offset in range(3):
            glyco[site["position"]+offset].append(offset)
    for i, (position, native, alt) in enumerate(zip(pos,natives,alternatives)):
        context = geometry.get(int(position), {})
        neutral_gain = native in NEUTRAL
        distance = context.get("distance", -1)
        if neutral_gain and alt in "KR" and context.get("side") in (1,2) and 1 <= distance <= spec.maximum_flank_distance:
            matrix[i,0] = (-0.5 if context["side"] == 1 else 0.5) * np.exp(-(distance-1)/spec.decay_length)
        if neutral_gain and alt in "DEKR" and context.get("in_tm") == 1 and context["tm_length"] >= spec.minimum_tm_length:
            matrix[i,1] = context["depth"]
        # Binary union: a single mutation is not counted twice for overlapping motifs.
        matrix[i,2] = float(any((offset == 0 and alt != "N") or (offset == 1 and alt == "P") or
                               (offset == 2 and alt not in "ST") for offset in glyco.get(int(position), [])))
        matrix[i,3] = float(int(position) in bonded)
    return matrix


def apply_refinement(baseline, features, coefficients, *, already_refined=False):
    """One additive application, with explicit direction, scope and no clipping."""
    if already_refined:
        raise ValueError("Refinement must be applied once to the raw baseline")
    if set(coefficients) != set(RULES):
        raise ValueError("Every named coefficient must be explicitly supplied")
    beta = np.array([coefficients[r] for r in RULES], dtype=float)
    x, baseline = np.asarray(features,dtype=float), np.asarray(baseline,dtype=float)
    if x.shape != (len(baseline),len(RULES)) or baseline.ndim != 1:
        raise ValueError("Mismatched score/feature shapes")
    if not np.isfinite(x).all() or not np.isfinite(beta).all() or np.any(beta<0):
        raise ValueError("Finite features and nonnegative directional coefficients required")
    if np.isinf(baseline).any():
        raise ValueError("Infinite baseline score")
    contribution = x*beta
    result = baseline + contribution.sum(axis=1)
    return dict(scores=result, contributions=contribution, rules=list(RULES),
                coefficients=dict(coefficients), baseline_nan_preserved=True,
                total_rule_compliance_guaranteed=False, score_scale="unchanged_ranking_score")


def ablation_coefficients(coefficients):
    """Raw baseline plus PI and all eight subsets of the three additions."""
    if set(coefficients) != set(RULES):
        raise ValueError("Every named coefficient must be explicitly supplied")
    result = {"raw": {r:0.0 for r in RULES}}
    for n in range(4):
        for selected in itertools.combinations(RULES[1:], n):
            active = {RULES[0], *selected}
            key = "+".join(r for r in RULES if r in active)
            result[key] = {r:float(coefficients[r]) if r in active else 0.0 for r in RULES}
    return result


def fit_coefficients(features, baseline, target, *, accession, assay, family, role,
                     target_scale, label_provenance, ridge=0.001, minimum_families=3,
                     locked_zero_rules=(), tolerance=1e-10, max_iterations=10000):
    """Nonnegative residual calibration on declared calibration families only.

    This API is prepared for the final stage. It is not called by preparation.
    Equal family/protein/assay weights and within-assay centering remove assay
    offsets; they do not replace a chemically matched control experiment.
    """
    x = np.asarray(features,dtype=float)
    baseline,target = np.asarray(baseline,dtype=float),np.asarray(target,dtype=float)
    n = len(baseline)
    columns = [np.asarray(c).astype(str) for c in (accession,assay,family,role)]
    if x.shape != (n,4) or target.shape != (n,) or any(c.shape != (n,) for c in columns) or not n:
        raise ValueError("Calibration dimensions do not agree")
    if not np.isfinite(x).all() or not np.isfinite(baseline).all() or not np.isfinite(target).all():
        raise ValueError("Calibration values must be finite; do not impute missing outcomes")
    acc, assays, families, roles = columns
    if any(any(not value.strip() for value in column) for column in (acc,assays,families)):
        raise ValueError("Calibration accession, assay and family identifiers are required")
    if set(roles) != {"calibration"}:
        raise ValueError("Only calibration rows may fit coefficients; evaluation roles are forbidden")
    if not target_scale or not label_provenance:
        raise ValueError("Experimental target scale and label provenance must be declared")
    if ridge <= 0 or not np.isfinite(ridge) or minimum_families < 1:
        raise ValueError("Declare a positive ridge penalty and family threshold")
    if len(set(families)) < minimum_families:
        raise ValueError("Too few calibration families")
    if any(len(set(families[acc == a])) != 1 for a in set(acc)):
        raise ValueError("A protein cannot belong to multiple calibration families")
    locked = set(locked_zero_rules)
    if not locked <= set(RULES):
        raise ValueError("Unknown locked rule")
    xc, residual, weight = x.copy(),target-baseline,np.zeros(n)
    for fam in sorted(set(families)):
        proteins = sorted(set(acc[families == fam]))
        for protein in proteins:
            ap = acc == protein
            names = sorted(set(assays[ap]))
            for name in names:
                selected = ap & (assays == name)
                weight[selected] = 1/(len(set(families))*len(proteins)*len(names)*selected.sum())
                xc[selected] -= xc[selected].mean(axis=0)
                residual[selected] -= residual[selected].mean()
    gram = xc.T @ (weight[:,None]*xc)
    response = xc.T @ (weight*residual)
    unsupported = {RULES[j] for j in range(4) if gram[j,j] <= np.finfo(float).eps}
    gram += ridge*np.eye(4)
    beta = np.zeros(4)
    converged = False
    for iteration in range(max_iterations):
        previous = beta.copy()
        for j,rule in enumerate(RULES):
            if rule in locked or rule in unsupported:
                continue
            beta[j] = max(0.0,(response[j] - gram[j]@beta + gram[j,j]*beta[j])/gram[j,j])
        if np.max(np.abs(beta-previous)) <= tolerance:
            converged = True
            break
    if not converged:
        raise RuntimeError("Calibration did not converge; no coefficients released")
    return dict(coefficients=dict(zip(RULES,beta.tolist())), converged=True,
                iterations=iteration+1, ridge=ridge, target_scale=target_scale,
                label_provenance=label_provenance, calibration_families=sorted(set(families)),
                locked_zero_rules=sorted(locked), unsupported_rules=sorted(unsupported),
                disulfide_design="control_free_predictive_refinement_not_causal_effect_estimation",
                total_rule_compliance_guaranteed=False,
                evaluation_required="group-disjoint accuracy, rule audits, ablations and collateral changes")


def write_json(path, value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+"\n",encoding="utf-8")


def main(argv=None):
    """Explicit final-stage calibration/application; never runs on import."""
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command",required=True)
    fit = commands.add_parser("fit",help="Fit only a supplied calibration table")
    fit.add_argument("--calibration",type=Path,required=True)
    fit.add_argument("--target-scale",required=True)
    fit.add_argument("--label-provenance",required=True)
    fit.add_argument("--minimum-families",type=int,default=3)
    fit.add_argument("--ridge",type=float,default=0.001)
    fit.add_argument("--lock-zero",choices=RULES,action="append",default=[])
    fit.add_argument("--out",type=Path,required=True)
    apply = commands.add_parser("apply",help="Apply to supplied raw ranking scores")
    apply.add_argument("--input",type=Path,required=True)
    apply.add_argument("--coefficients",type=Path,required=True)
    apply.add_argument("--score-direction",choices=["higher-is-damaging"],required=True)
    apply.add_argument("--out",type=Path,required=True)
    args = parser.parse_args(argv)
    source = args.calibration if args.command == "fit" else args.input
    protected = [source] + ([args.coefficients] if args.command == "apply" else [])
    if args.out.resolve() in {p.resolve() for p in protected}:
        raise ValueError("Outputs must not replace inputs")
    with np.load(source,allow_pickle=False) as data:
        if tuple(data["feature_names"].astype(str)) != RULES:
            raise ValueError("Unexpected feature names or order")
        if args.command == "fit":
            result = fit_coefficients(data["features"],data["baseline"],data["target"],
                **{name:data[name] for name in ("accession","assay","family","role")},
                target_scale=args.target_scale,label_provenance=args.label_provenance,
                ridge=args.ridge,minimum_families=args.minimum_families,locked_zero_rules=args.lock_zero)
            args.out.parent.mkdir(parents=True,exist_ok=True)
            write_json(args.out,result)
        else:
            if str(data["score_state"]) != "raw":
                raise ValueError("Input must explicitly declare score_state=raw")
            identifiers = {name:data[name] for name in ("accession","position","native","alt")}
            n = len(data["baseline"])
            if any(values.shape != (n,) for values in identifiers.values()):
                raise ValueError("Variant identifiers must align with raw scores")
            keys = list(zip(*(identifiers[name].tolist() for name in ("accession","position","native","alt"))))
            if len(set(keys)) != n:
                raise ValueError("Duplicate scored variants")
            result = apply_refinement(data["baseline"],data["features"],
                json.loads(args.coefficients.read_text(encoding="utf-8"))["coefficients"])
            args.out.parent.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(args.out,scores=result["scores"],contributions=result["contributions"],
                feature_names=np.array(RULES),score_state=np.array("refined"),**identifiers)
    print(json.dumps(dict(command=args.command,out=str(args.out),
        benchmark_run=False,clinical_probability_claim=False)))


if __name__ == "__main__":
    main()

"""Construct variant features and score the saved predictor-free head."""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(Path(__file__).resolve().parent))
from membrane_rule_audit.refinement import ALPHABET, topology_map, sequence_digest
from a9_l1_features import GR_C, GR_P, GR_V, HYDRO, CHARGE, HELIX, STRAND, grantham, blosum62

FEATURES = ["blosum62", "grantham", "d_composition", "d_polarity", "d_volume",
            "d_hydro", "d_charge", "d_helix", "d_strand", "native_is_kr",
            "alt_is_kr", "side", "distance", "in_tm", "tmd_index", "loop_len",
            "n_tmds", "is_membrane"]


def build_tabular(record, positions, natives, alternatives):
    """Construct the 18 original supervised features with the correct codes."""
    pos = np.asarray(positions)
    native, alt = np.asarray(natives).astype(str),np.asarray(alternatives).astype(str)
    if pos.ndim != 1 or native.shape != pos.shape or alt.shape != pos.shape:
        raise ValueError("Mismatched variant columns")
    seq = record["sequence"]
    if sequence_digest(seq) != record["sequence_sha256"]:
        raise ValueError("Canonical sequence hash mismatch")
    if not np.isfinite(pos).all() or np.any(pos != pos.astype(int)) or np.any(pos<1) or np.any(pos>len(seq)):
        raise ValueError("Invalid canonical position")
    if np.any(native != np.array(list(seq))[pos.astype(int)-1]) or not np.isin(alt,list(ALPHABET)).all():
        raise ValueError("Native mismatch or invalid alternative")
    if np.any(native == alt):
        raise ValueError("Unchanged residues must not enter the missense comparator")
    bl = blosum62()
    context = topology_map(record["topology"])
    tms = record["topology"].get("tmds") or []
    rows = []
    for p,n,a in zip(pos.astype(int),native,alt):
        geo = context.get(int(p), dict(side=-1,distance=-1,in_tm=-1,tmd_index=-1,loop_len=-1))
        rows.append([bl[n,a],grantham(n,a),GR_C[a]-GR_C[n],GR_P[a]-GR_P[n],
                     GR_V[a]-GR_V[n],HYDRO[a]-HYDRO[n],CHARGE[a]-CHARGE[n],
                     HELIX[a]-HELIX[n],STRAND[a]-STRAND[n],float(n in "KR"),float(a in "KR"),
                     geo["side"],geo["distance"],geo["in_tm"],geo["tmd_index"],geo["loop_len"],
                     len(tms),float(bool(tms))])
    return np.asarray(rows,dtype=np.float32).reshape(-1,len(FEATURES))


def load_bundle(path):
    """Load a locally trusted serialized model only when explicitly invoked."""
    import joblib
    bundle = joblib.load(path)
    if bundle.get("features") != FEATURES or bundle.get("n_pc") != 64:
        raise ValueError("Bundle is not the declared predictor-free 18+64 feature head")
    if bundle["model"].n_features_in_ != 82 or bundle["pca"].n_features_in_ != 1024:
        raise ValueError("Bundle feature dimensions do not agree")
    return bundle


def score_positions(bundle, record, positions, embedding_positions, embeddings):
    """Optional final-stage inference. Does not access AM/ESM scores or geometry.

    Returned native cells are NaN. Missing embeddings fail rather than becoming
    zero vectors. Callers must verify embedding reference provenance separately.
    """
    if bundle.get("features") != FEATURES:
        raise ValueError("Predictor inputs must be the declared standalone features")
    positions = np.asarray(positions)
    embedding_positions = np.asarray(embedding_positions)
    embeddings = np.asarray(embeddings,dtype=np.float32)
    if positions.ndim != 1 or not len(positions) or len(set(positions)) != len(positions):
        raise ValueError("Unique canonical positions required")
    if not np.isfinite(positions).all() or np.any(positions != positions.astype(int)):
        raise ValueError("Canonical positions must be finite integers")
    if embedding_positions.ndim != 1 or not np.isfinite(embedding_positions).all() or np.any(embedding_positions != embedding_positions.astype(int)):
        raise ValueError("Embedding positions must be finite integers")
    positions,embedding_positions = positions.astype(int),embedding_positions.astype(int)
    if embeddings.shape != (len(embedding_positions),1024) or len(set(embedding_positions)) != len(embedding_positions):
        raise ValueError("Invalid embedding shape or duplicate embedding position")
    if not np.isfinite(embeddings).all():
        raise ValueError("Nonfinite embedding")
    lookup = {int(p):i for i,p in enumerate(embedding_positions)}
    if any(int(p) not in lookup for p in positions):
        raise ValueError("Missing embedding for a requested position")
    sequence = record["sequence"]
    if np.any(positions<1) or np.any(positions>len(sequence)):
        raise ValueError("Position outside canonical sequence")
    native = np.array(list(sequence))[positions-1]
    vp = np.repeat(positions,19)
    vn = np.repeat(native,19)
    va = np.concatenate([np.array([a for a in ALPHABET if a != n]) for n in native])
    tabular = build_tabular(record,vp,vn,va)
    projected = bundle["pca"].transform(embeddings[[lookup[int(p)] for p in positions]])
    X = np.hstack([tabular,np.repeat(projected,19,axis=0)]).astype(np.float32)
    predicted = bundle["model"].predict(X)
    scores = np.full((len(positions),20),np.nan)
    for i,n in enumerate(native):
        scores[i,[j for j,a in enumerate(ALPHABET) if a != n]] = predicted[i*19:(i+1)*19]
    return dict(accession=np.repeat(record["accession"],len(positions)),position=positions,
                native=native,alphabet=np.array(list(ALPHABET)),scores=scores)


def same_model_probe(baseline, features, rule, amount):
    """A calibration probe always starts from the supplied standalone scores."""
    from membrane_rule_audit.refinement import RULES, apply_refinement
    if rule not in RULES:
        raise ValueError("Unknown correction rule")
    coefficients = {r:amount if r == rule else 0.0 for r in RULES}
    return apply_refinement(baseline,features,coefficients)["scores"]


def refine_position_scores(record, scored, coefficients):
    """Connect the standalone score matrix to all four additive components.

    This is a final-stage API. Existing missing alternatives stay missing and
    native identity cells are never substituted for real missense scores.
    """
    from membrane_rule_audit.refinement import variant_features,apply_refinement,RULES
    if "score_state" in scored and scored["score_state"] != "raw":
        raise ValueError("Only a raw standalone score matrix may enter refinement")
    if "".join(scored["alphabet"]) != ALPHABET:
        raise ValueError("Unexpected score alphabet")
    positions = np.asarray(scored["position"])
    scores = np.asarray(scored["scores"],dtype=float)
    native = np.asarray(scored["native"]).astype(str)
    if scores.shape != (len(positions),20) or native.shape != positions.shape:
        raise ValueError("Malformed standalone score matrix")
    if positions.ndim != 1 or not np.isfinite(positions).all() or len(set(positions)) != len(positions):
        raise ValueError("Unique finite canonical positions required")
    accession = np.asarray(scored["accession"]).astype(str)
    if accession.shape != positions.shape or np.any(accession != record["accession"]):
        raise ValueError("Score accession mismatch")
    if np.isinf(scores).any():
        raise ValueError("Infinite standalone scores")
    for i,(p,n) in enumerate(zip(positions,native)):
        if not 1 <= p <= len(record["sequence"]) or int(p) != p or record["sequence"][int(p)-1] != n:
            raise ValueError("Native/canonical position mismatch")
        if np.isfinite(scores[i,ALPHABET.index(n)]):
            raise ValueError("Unchanged residues must be masked before refinement")
    row,col = np.where(np.isfinite(scores))
    features = variant_features(record,positions[row],native[row],np.array(list(ALPHABET))[col])
    result = apply_refinement(scores[row,col],features,coefficients)
    updated = scores.copy()
    updated[row,col] = result["scores"]
    contributions = np.full((len(positions),20,len(RULES)),np.nan)
    contributions[row,col] = result["contributions"]
    return {**scored,"scores":updated,"score_state":"refined","contributions":contributions,
            "coefficients":dict(coefficients),"total_rule_compliance_guaranteed":False}

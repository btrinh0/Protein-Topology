"""Shared implementation for the Form A and Form B studies.

Form A reads site geometry without loading AlphaMissense scores. Estimation uses
the installed validated panel module. The working-copy audit module supplies the
report card and is pinned by SHA-256.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


def _validated_panel():
    """The installed validated panel.py, loaded under its own module name.

    The working-copy tools package has the same import name. Importing the package here would make the
    project's own tests (which put tools/ first on sys.path) pick up whichever copy was imported first.
    """
    import importlib.metadata as metadata
    location = metadata.distribution("membrane-rule-audit").locate_file("membrane_rule_audit/panel.py")
    spec = importlib.util.spec_from_file_location("membrane_rule_audit_validated_panel", location)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


P = _validated_panel()

ALPHABET, NEUTRAL = P.ALPHABET, P.NEUTRAL
NEUTRAL_INDEX = np.array([ALPHABET.index(a) for a in NEUTRAL])
KR_INDEX = np.array([ALPHABET.index("K"), ALPHABET.index("R")])
# Row i lists the 19 alternatives of native residue i.
ALT_TABLE = np.array([[a for a in range(20) if a != r] for r in range(20)], dtype=np.int64)

AM_ARCHIVE = ROOT / "data/processed/gate0/am_full_flank_scores.npz"
ESM_ARCHIVE = ROOT / "data/processed/gate0/esm1b_flank_scores.npz"
TOPOLOGY = ROOT / "data/raw/topology_viewer.html"
WHEEL = ROOT / "dist/membrane-rule-audit-v0.1.0-validated-20261004/wheel/membrane_rule_audit-0.1.0-py3-none-any.whl"
MANIFEST = ROOT / "dist/membrane-rule-audit-v0.1.0-validated-20261004/RELEASE_MANIFEST.json"
WORKING_AUDIT = ROOT / "tools/membrane_rule_audit/audit.py"
GEOMETRY = ("protein", "position", "side", "distance", "tmd_index", "ref_aa")

MASTER_SEED, PILOT_SEED = 20261008, 20261001
# Keep clustering aligned with the composition-control analysis.
MINHASH_SEED = 20260929
PANEL_SEED, PANEL_DRAWS = 20261003, 1999
REPORT_CARD = {"max_distance": 60, "folds": 5, "draws": 999, "seed": 20261006}
SPLIT_PERSON = b"p6-ab-split-v1"
CONTEXT_RADIUS = 5
CONTEXT_CLASSES = ("KR", "DE", "ACFILMVW")
NOISE_SD, SITE_SD = 0.5, 0.5
# Express the equivalence margin as a fraction of SD_z0.
EQUIVALENCE = 0.02

SCENARIOS = {
    "N0": dict(rho=0.0, lam=None, side_main=0.0),
    "S_main": dict(rho=0.0, lam=None, side_main=0.30),
    "P_neg05": dict(rho=-0.05, lam=None, side_main=0.0),
    "P_neg10": dict(rho=-0.10, lam=None, side_main=0.0),
    "P_pos10": dict(rho=0.10, lam=None, side_main=0.0),
    "P_decay": dict(rho=-0.20, lam=8.0, side_main=0.0),
}
STUDENT_FEATURES = {
    "S1": ("ref", "alt", "side", "distance"),
    "S2": ("ref", "alt"),
    "S3": ("ref", "alt", "kr", "de", "hyd"),
}
STREAMS = {"train": 1, "validation": 2, "test": 3, "model": 4, "subsample": 5}


def validated_tool() -> dict:
    """Confirm the installed package is the validated release, byte for byte."""
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))["files"]
    installed = Path(P.__file__).resolve().parent
    for name in ("panel.py", "audit.py", "__init__.py", "__main__.py"):
        if P.sha256(installed / name) != manifest[f"source/membrane_rule_audit/{name}"]["sha256"]:
            raise RuntimeError(f"Installed {name} differs from the validated release")
    return {"wheel_sha256": P.sha256(WHEEL), "panel_sha256": P.sha256(installed / "panel.py"),
            "validated_audit_sha256": P.sha256(installed / "audit.py")}


def report_card_tool():
    """The working-copy audit.py (placebo rules, declared directions, verdicts), loaded by path."""
    spec = importlib.util.spec_from_file_location("mra_report_card", WORKING_AUDIT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module, P.sha256(WORKING_AUDIT)


def package_versions() -> dict:
    import importlib.metadata as md
    return {name: md.version(name) for name in
            ("numpy", "pandas", "scikit-learn", "scipy", "torch", "biopython", "membrane-rule-audit")}


def load_sites(path: Path = AM_ARCHIVE) -> tuple[dict, np.ndarray]:
    """Unique sites with geometry only. The score array is never read."""
    with np.load(path, allow_pickle=False) as blob:
        if "".join(blob["aa_alphabet"].tolist()) != ALPHABET:
            raise ValueError("Unexpected archive alphabet")
        arrays = {name: blob[name] for name in GEOMETRY}
        accessions = blob["accessions"].copy()
    code = (arrays["protein"].astype(np.int64) << 32) + arrays["position"].astype(np.int64)
    _, first = np.unique(code, return_index=True)
    return {name: arrays[name][first].astype(np.int64) for name in GEOMETRY}, accessions


def take(sites: dict, index) -> dict:
    return {name: values[index] for name, values in sites.items()}


def context_counts(sites: dict, accessions, census: dict, radius: int = CONTEXT_RADIUS) -> np.ndarray:
    """Counts of K/R, D/E and hydrophobic residues within +/-radius of each site, site excluded."""
    protein = sites["protein"]
    out = np.zeros((len(protein), len(CONTEXT_CLASSES)), dtype=np.int64)
    bounds = np.flatnonzero(np.r_[True, protein[1:] != protein[:-1], True])
    for start, stop in zip(bounds[:-1], bounds[1:]):
        sequence = census[str(accessions[protein[start]])]["sequence"]
        codes = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
        index = sites["position"][start:stop] - 1
        if index.max() >= len(sequence):
            raise ValueError("Site beyond the census sequence")
        lo, hi = np.maximum(index - radius, 0), np.minimum(index + radius, len(sequence) - 1)
        for column, letters in enumerate(CONTEXT_CLASSES):
            member = np.isin(codes, np.frombuffer(letters.encode("ascii"), dtype=np.uint8)).astype(np.int64)
            running = np.concatenate([[0], np.cumsum(member)])
            out[start:stop, column] = running[hi + 1] - running[lo] - member[index]
    return out


def split_bucket(key: str) -> int:
    digest = hashlib.blake2b(key.encode("ascii"), digest_size=8, person=SPLIT_PERSON).digest()
    return int.from_bytes(digest, "little") % 10


def bucket_label(bucket: int) -> str:
    return "train" if bucket <= 5 else "validation" if bucket == 6 else "test"


def cluster_split(census: dict, seed: int = MINHASH_SEED, kind: str = "sequence") -> dict:
    """Assign whole clusters to train (buckets 0-5), validation (6) or test (7-9).

    A cluster is keyed by its lexicographically smallest accession (sequence clusters) or its
    gene-symbol family root (symbol families), then hashed with blake2b.
    """
    from run_step2_composition_control import minhash_clusters, symbol_families
    if kind == "sequence":
        ordered = sorted(census)
        clusters = minhash_clusters(census, seed)
        key_of = {accession: ordered[root] for accession, root in clusters.items()}
    elif kind == "symbol":
        key_of = {accession: root.encode("ascii", "replace").decode("ascii")
                  for accession, root in symbol_families(census).items()}
    else:
        raise ValueError(kind)
    label = {accession: bucket_label(split_bucket(key)) for accession, key in key_of.items()}
    members = defaultdict(set)
    for accession, key in key_of.items():
        members[key].add(accession)
    counts = {name: {"clusters": sum(1 for key in members if label[next(iter(members[key]))] == name),
                     "proteins": sum(1 for v in label.values() if v == name)}
              for name in ("train", "validation", "test")}
    lines = "\n".join(f"{a}\t{label[a]}" for a in sorted(label))
    return {"label": label, "key_of": key_of, "kind": kind, "seed": seed, "counts": counts,
            "assignment_sha256": hashlib.sha256(lines.encode("ascii")).hexdigest()}


def site_pools(sites: dict, accessions, label: dict) -> dict:
    """Site indices for each split; a protein sits in exactly one."""
    per_protein = np.array([label[str(a)] for a in accessions])
    per_site = per_protein[sites["protein"]]
    return {name: np.flatnonzero(per_site == name) for name in ("train", "validation", "test")}


def blosum_base() -> np.ndarray:
    """Higher means more damaging: -BLOSUM62, z-scored over the 380 off-diagonal pairs (population SD)."""
    from Bio.Align import substitution_matrices
    matrix = substitution_matrices.load("BLOSUM62")
    raw = np.array([[matrix[a, b] for b in ALPHABET] for a in ALPHABET], dtype=float)
    off = ~np.eye(20, dtype=bool)
    base = -raw
    return (base - base[off].mean()) / base[off].std()


def arm_difference(base: np.ndarray) -> np.ndarray:
    """K/R-arm mean minus uncharged-arm mean of the baseline, by native residue (NaN if not eligible)."""
    out = np.full(20, np.nan)
    for ref in NEUTRAL_INDEX:
        comparator = [a for a in NEUTRAL_INDEX if a != ref]
        out[ref] = base[ref, KR_INDEX].mean() - base[ref, comparator].mean()
    return out


def sd_z0(sites: dict, base: np.ndarray) -> float:
    """SD of unplanted K/R-gain site penalties on this geometry: baseline spread plus arm noise."""
    z = arm_difference(base)[sites["ref_aa"]]
    z = z[np.isfinite(z)]
    return float(np.sqrt(np.var(z, ddof=1) + NOISE_SD ** 2 * (1 / 2 + 1 / 14)))


def rng_for(master: int, scenario: str, replicate: int, stream: str) -> np.random.Generator:
    index = list(SCENARIOS).index(scenario)
    return np.random.default_rng(np.random.SeedSequence([master, index, replicate, STREAMS[stream]]))


def planted_h(sites: dict, amp: float, lam) -> np.ndarray:
    """Planted contribution to the K/R-gain penalty at each site (zero unless the native residue is uncharged)."""
    sign = np.where(sites["side"] == 0, 1.0, -1.0)
    decay = np.ones(len(sign)) if lam is None else np.exp(-(sites["distance"] - 1) / lam)
    return np.where(np.isin(sites["ref_aa"], NEUTRAL_INDEX), (amp / 2) * sign * decay, 0.0)


def draw_target(sites: dict, base: np.ndarray, amp: float, lam, side_main: float,
                rng: np.random.Generator) -> dict:
    """clean = baseline + side main effect + planted rule; site is shared across a site's alternatives."""
    n = len(sites["ref_aa"])
    ref, sign = sites["ref_aa"], np.where(sites["side"] == 0, 1.0, -1.0)
    clean = base[ref] + side_main * sign[:, None]
    h = planted_h(sites, amp, lam)
    clean[:, KR_INDEX[0]] += h
    clean[:, KR_INDEX[1]] += h
    site = rng.normal(0.0, SITE_SD, n)[:, None]
    noise = rng.normal(0.0, NOISE_SD, (n, 20))
    return {"clean": clean, "site": site, "noise": noise, "h": h}


def mask_native(matrix: np.ndarray, ref: np.ndarray) -> np.ndarray:
    out = matrix.copy()
    out[np.arange(len(ref)), ref] = np.nan
    return out


def oracle_delta(sites: dict, h: np.ndarray) -> tuple[float, float]:
    """The estimator's stratified weighting applied to the planted term alone, in plain Python."""
    eligible = np.isin(sites["ref_aa"], NEUTRAL_INDEX)
    cyto, outer = defaultdict(list), defaultdict(list)
    for protein, ref, side, value in zip(sites["protein"][eligible], sites["ref_aa"][eligible],
                                         sites["side"][eligible], h[eligible]):
        (cyto if side == 0 else outer)[(protein, ref)].append(value)
    numerator = denominator = 0.0
    for key in cyto.keys() & outer.keys():
        n1, n0 = len(cyto[key]), len(outer[key])
        weight = n1 * n0 / (n1 + n0)
        numerator += weight * (np.mean(cyto[key]) - np.mean(outer[key]))
        denominator += weight
    return numerator / denominator, denominator


def independent_delta(sites: dict, scores: np.ndarray) -> float:
    """The estimator from arm means in plain Python, for any model's output, sharing no code with the tool."""
    penalty = np.nanmean(scores[:, KR_INDEX], axis=1) - np.nanmean(scores[:, NEUTRAL_INDEX], axis=1)
    return oracle_delta(sites, np.where(np.isfinite(penalty), penalty, 0.0))[0]


def audit_matrix(sites: dict, scores: np.ndarray, draws: int = PANEL_DRAWS, seed: int = PANEL_SEED):
    """The packaged estimator on a (sites x 20) matrix, the same code run_panel_audit uses after loading."""
    panel = {"protein": sites["protein"], "ref_aa": sites["ref_aa"], "side": sites["side"], "scores": scores}
    contrast = P.gain_contrast(panel)
    parts = P.stratum_parts(panel, contrast["penalty"], contrast["eligible"])
    return P.summarize_parts(parts, draws, seed), parts, contrast


def write_scores_npz(path: Path, sites: dict, scores: np.ndarray, accessions) -> None:
    """Mapped archive schema with a `score` array and one row per finite alternative."""
    rows, alts = np.nonzero(np.isfinite(scores))
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, protein=sites["protein"][rows].astype(np.int32),
             position=sites["position"][rows].astype(np.int32),
             side=sites["side"][rows].astype(np.int8), distance=sites["distance"][rows].astype(np.int16),
             tmd_index=sites["tmd_index"][rows].astype(np.int16), ref_aa=sites["ref_aa"][rows].astype(np.int8),
             alt_aa=alts.astype(np.int8), score=scores[rows, alts], accessions=accessions,
             aa_alphabet=np.array(list(ALPHABET)))


def naive_comparators(sites: dict, scores: np.ndarray) -> dict:
    """(i) charged-arm-only side difference; (ii) same-position subtraction pooled without strata."""
    eligible = np.isin(sites["ref_aa"], NEUTRAL_INDEX)
    cyto = eligible & (sites["side"] == 0)
    outer = eligible & (sites["side"] == 1)
    arm = np.nanmean(scores[:, KR_INDEX], axis=1)
    penalty = arm - np.nanmean(scores[:, NEUTRAL_INDEX], axis=1)
    return {"charged_arm_only": float(arm[cyto].mean() - arm[outer].mean()),
            "pooled_paired": float(penalty[cyto].mean() - penalty[outer].mean())}


def site_columns(sites: dict, context: np.ndarray) -> dict:
    return {"ref": sites["ref_aa"], "side": sites["side"], "distance": sites["distance"],
            "kr": context[:, 0], "de": context[:, 1], "hyd": context[:, 2]}


def design_matrix(columns: dict, index, names) -> np.ndarray:
    """One row per (site, alternative); alternatives follow ALT_TABLE order, site-major."""
    ref = columns["ref"][index]
    built = {"ref": np.repeat(ref, 19), "alt": ALT_TABLE[ref].ravel()}
    for name in names:
        if name not in built:
            built[name] = np.repeat(columns[name][index], 19)
    return np.column_stack([built[name] for name in names]).astype(np.float32)


def row_targets(matrix: np.ndarray, ref: np.ndarray) -> np.ndarray:
    return np.take_along_axis(matrix, ALT_TABLE[ref], axis=1).ravel()


def predictions_to_matrix(pred: np.ndarray, ref: np.ndarray) -> np.ndarray:
    matrix = np.full((len(ref), 20), np.nan)
    np.put_along_axis(matrix, ALT_TABLE[ref], pred.reshape(-1, 19), axis=1)
    return matrix


def fit_student(X, y, Xv, yv, seed: int, max_iter: int = 300):
    """Boosted trees with fixed settings; the iteration count is chosen on validation clusters only."""
    from sklearn.ensemble import HistGradientBoostingRegressor
    model = HistGradientBoostingRegressor(max_iter=max_iter, learning_rate=0.1, categorical_features=[0, 1],
                                          early_stopping=False, random_state=seed)
    model.fit(X, y)
    curve = np.array([np.mean((p - yv) ** 2) for p in model.staged_predict(Xv)])
    return model, int(curve.argmin()) + 1, curve


def predict_at(model, X, n_iter: int) -> np.ndarray:
    for count, pred in enumerate(model.staged_predict(X), start=1):
        if count == n_iter:
            return pred
    raise ValueError("n_iter exceeds the fitted iterations")


def wilson(k: int, n: int) -> list[float]:
    z = 1.959963984540054
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [float(centre - half), float(centre + half)]


# Fixed ten-residue K/R-gain penalty SD per predictor (verified_values.json), the paper's unit.
SD_FIXED = {"alphamissense": 0.15656056894897286, "esm1b": 2.0195231488036884}
CORRECTION = {"margin": 0.05, "decay_length": 8.0}


def subset_panel(panel: dict, mask: np.ndarray) -> dict:
    """Subset a panel and widen geometry to prevent int8 cell-index overflow."""
    out = {name: panel[name][mask].astype(np.int64) for name in GEOMETRY}
    out["scores"] = panel["scores"][mask]
    out["accessions"] = panel["accessions"]
    return out


def split_mask(panel: dict, label: dict, name: str) -> np.ndarray:
    per_protein = np.array([label[str(a)] for a in panel["accessions"]])
    return per_protein[panel["protein"]] == name


def window_panel(panel: dict, window: int = 10) -> dict:
    return subset_panel(panel, panel["distance"] <= window)


def decay(distance: np.ndarray, length) -> np.ndarray:
    return np.ones(len(distance)) if length is None else np.exp(-(distance - 1) / length)


def _strata(panel: dict, values: np.ndarray, keep: np.ndarray):
    """Protein x native-residue strata with both sides present: weights and per-arm sums of `values`."""
    group = panel["protein"][keep] * 20 + panel["ref_aa"][keep]
    cyto = panel["side"][keep] == 0
    y = values[keep]
    unique, inverse = np.unique(group, return_inverse=True)
    size = len(unique)
    n1, n0 = np.bincount(inverse[cyto], minlength=size), np.bincount(inverse[~cyto], minlength=size)
    s1 = np.bincount(inverse[cyto], weights=y[cyto], minlength=size)
    s0 = np.bincount(inverse[~cyto], weights=y[~cyto], minlength=size)
    shared = (n1 > 0) & (n0 > 0)
    weight = n1[shared] * n0[shared] / (n1[shared] + n0[shared])
    return unique[shared], weight, s1[shared] / n1[shared], s0[shared] / n0[shared]


def correction_weight(window: dict, length) -> float:
    """K in Delta'(beta) = Delta0 - beta*K: the weighted mean of the two arms' average decay."""
    contrast = P.gain_contrast(window)
    keep = contrast["eligible"] & np.isfinite(contrast["penalty"])
    _, weight, cyto, outer = _strata(window, decay(window["distance"], length), keep)
    return float((weight * (cyto + outer) / 2).sum() / weight.sum())


def fit_beta(delta0: float, k: float, sd: float, margin: float, declared: int = -1) -> tuple[float, float]:
    """One-sided closed form: a scorer already at or beyond the target gets beta = 0."""
    target = declared * margin * sd
    gap = delta0 - target if declared == -1 else target - delta0
    return max(0.0, gap / k), target


def apply_correction(panel: dict, beta: float, length, declared: int = -1, side=None) -> np.ndarray:
    """Shift K/R-gain rows (native residue uncharged) by -/+ (beta/2) decay(d) on the cytosolic/outer side."""
    side = panel["side"] if side is None else side
    shift = np.where(side == 0, -1.0, 1.0) * (beta / 2) * decay(panel["distance"], length)
    if declared == 1:
        shift = -shift
    shift = np.where(np.isin(panel["ref_aa"], NEUTRAL_INDEX), shift, 0.0)
    out = panel["scores"].copy()
    out[:, KR_INDEX] += shift[:, None]
    return out


def shuffle_side(panel: dict, rng: np.random.Generator) -> np.ndarray:
    """Permute side labels within protein x native-residue groups, keeping every group's side counts."""
    group = panel["protein"] * 20 + panel["ref_aa"]
    in_order = np.argsort(group, kind="stable")
    at_random = np.lexsort((rng.random(len(group)), group))
    side = panel["side"].copy()
    side[in_order] = panel["side"][at_random]
    return side


def fidelity(original: np.ndarray, corrected: np.ndarray) -> dict:
    from scipy.stats import spearmanr
    finite = np.isfinite(original) & np.isfinite(corrected)
    a, b = original[finite], corrected[finite]
    return {"spearman": float(spearmanr(a, b).statistic), "mean_abs_change": float(np.abs(a - b).mean()),
            "fraction_changed": float((a != b).mean()), "n_rows": int(finite.sum())}


def protein_contrasts(window: dict, min_strata: int = 3) -> dict:
    """Each protein's own weighted K/R-gain contrast, for proteins with enough two-sided strata."""
    contrast = P.gain_contrast(window)
    keep = contrast["eligible"] & np.isfinite(contrast["penalty"])
    ids, weight, cyto, outer = _strata(window, contrast["penalty"], keep)
    proteins, inverse = np.unique(ids // 20, return_inverse=True)
    numerator = np.bincount(inverse, weights=weight * (cyto - outer), minlength=len(proteins))
    denominator = np.bincount(inverse, weights=weight, minlength=len(proteins))
    count = np.bincount(inverse, minlength=len(proteins))
    enough = count >= min_strata
    return {"protein": proteins[enough], "contrast": numerator[enough] / denominator[enough],
            "weight": denominator[enough], "n_strata": count[enough]}


def protein_parts(window: dict, length) -> dict:
    """Per-protein sums behind the contrast: D (weighted side difference), K (weighted mean decay), W (weight)."""
    contrast = P.gain_contrast(window)
    keep = contrast["eligible"] & np.isfinite(contrast["penalty"])
    ids, weight, cyto, outer = _strata(window, contrast["penalty"], keep)
    _, _, decay_cyto, decay_outer = _strata(window, decay(window["distance"], length), keep)
    proteins, inverse = np.unique(ids // 20, return_inverse=True)
    return {"protein": proteins,
            "D": np.bincount(inverse, weights=weight * (cyto - outer), minlength=len(proteins)),
            "K": np.bincount(inverse, weights=weight * (decay_cyto + decay_outer) / 2, minlength=len(proteins)),
            "W": np.bincount(inverse, weights=weight, minlength=len(proteins))}


def cluster_interval(numerator, denominator, cluster, draws: int = PANEL_DRAWS, seed: int = PANEL_SEED) -> list:
    """Percentile interval for sum(numerator)/sum(denominator), resampling whole clusters."""
    codes = np.unique(cluster, return_inverse=True)[1]
    size = codes.max() + 1
    num = np.bincount(codes, weights=numerator, minlength=size)
    den = np.bincount(codes, weights=denominator, minlength=size)
    rng = np.random.default_rng(seed)
    samples = []
    for start in range(0, draws, 64):
        counts = rng.multinomial(size, np.full(size, 1 / size), size=min(64, draws - start))
        samples.extend((counts @ num / (counts @ den)).tolist())
    return [float(x) for x in np.percentile(samples, [2.5, 97.5])]


def run_report_card(tool, panel: dict, scores: np.ndarray, out_dir: Path) -> dict:
    """Write the (possibly corrected) scores in the mapped schema and run the working-copy report card."""
    npz = out_dir / "scores.npz"
    lib_sites = {name: panel[name] for name in GEOMETRY}
    write_scores_npz(npz, lib_sites, scores, panel["accessions"])
    summary = tool.run_audit(TOPOLOGY, npz, False, out_dir, max_distance=REPORT_CARD["max_distance"],
                             folds=REPORT_CARD["folds"], draws=REPORT_CARD["draws"], seed=REPORT_CARD["seed"])
    npz.unlink()
    return summary


def _leaves(obj, prefix=()):
    if isinstance(obj, dict):
        for key, value in obj.items():
            yield from _leaves(value, prefix + (key,))
    else:
        yield prefix, (tuple(obj) if isinstance(obj, list) else obj)


def card_diff(tool, before: dict, after: dict, target: str = "positive_charge_gain") -> dict:
    """Which report-card entries changed; entries that involve the corrected rule are expected to."""
    skip = {"as_of_utc", "scores", "topology"}
    involved = {name for name, pair in tool.CONTRASTS.items() if target in pair} | {target}
    a = {path: value for path, value in _leaves({k: v for k, v in before.items() if k not in skip})}
    b = {path: value for path, value in _leaves({k: v for k, v in after.items() if k not in skip})}
    changed = [path for path in a.keys() | b.keys() if a.get(path) != b.get(path)]
    expected = [p for p in changed if involved & set(p)]
    unexpected = [p for p in changed if not (involved & set(p))]
    return {"entries": len(a.keys() | b.keys()), "changed": len(changed), "expected_changes": len(expected),
            "unexpected_changes": ["/".join(p) for p in sorted(unexpected)[:50]],
            "n_unexpected": len(unexpected), "involved_names": sorted(involved)}

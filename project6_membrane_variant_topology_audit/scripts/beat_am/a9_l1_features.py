"""A9: L1 feature matrix for the development pool.

Builds one row per usable measured single substitution across the fitting-pool
proteins, with:

  predictors            AM pathogenicity, ESM1b damage (both damage-oriented)
  substitution          BLOSUM62, Grantham distance and its three components,
                        hydrophobicity / volume / charge / helix / strand deltas
  topology              side, distance to nearest TM boundary, in-TM flag,
                        TM index, loop length, TM count  (membrane proteins only)
  assay                 endpoint-type indicator, is_membrane
  target                within-assay rank-normalised damage

Coordinates come through route_a_common.align_construct exactly as in a8 -- no
global offset -- and the >=0.99 native-residue gate is applied.

Scope discipline: reads ONLY features/fitting_pool/. Confirmation proteins are
absent from this table entirely.

Creates new files only.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

INV = ROOT / "results" / "beat_am_20261016" / "inventory"
FEAT = ROOT / "results" / "beat_am_20261016" / "features"
POOL = FEAT / "fitting_pool"

PGYM_REF = ROOT / "data" / "raw" / "proteingym" / "reference_github_main" / "DMS_substitutions.csv"
ASSAY_DIR = ROOT / "data" / "raw" / "proteingym" / "v1.3_extracted" / "DMS_ProteinGym_substitutions"
UNIPROT_MAP = ROOT / "data" / "raw" / "proteingym" / "uniprot_entry_map_v1.3.json"
CENSUS = ROOT / "data" / "raw" / "topology_viewer.html"

MIN_AGREEMENT = 0.99
MIN_VARIANTS = 50
CONFIRM_ROLES = {"confirmation_membrane", "confirmation_allhuman"}

AA = "ACDEFGHIKLMNPQRSTVWY"

# Kyte-Doolittle hydropathy
HYDRO = {"A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8, "G": -0.4,
         "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8, "M": 1.9, "N": -3.5,
         "P": -1.6, "Q": -3.5, "R": -4.5, "S": -0.8, "T": -0.7, "V": 4.2,
         "W": -0.9, "Y": -1.3}
# Grantham's composition, polarity and molecular volume
GR_C = {"A": 0.00, "C": 2.75, "D": 1.38, "E": 0.92, "F": 0.00, "G": 0.74,
        "H": 0.58, "I": 0.00, "K": 0.33, "L": 0.00, "M": 0.00, "N": 1.33,
        "P": 0.39, "Q": 0.89, "R": 0.65, "S": 1.42, "T": 0.71, "V": 0.00,
        "W": 0.13, "Y": 0.20}
GR_P = {"A": 8.1, "C": 5.5, "D": 13.0, "E": 12.3, "F": 5.2, "G": 9.0,
        "H": 10.4, "I": 5.2, "K": 11.3, "L": 4.9, "M": 5.7, "N": 11.6,
        "P": 8.0, "Q": 10.5, "R": 10.5, "S": 9.2, "T": 8.6, "V": 5.9,
        "W": 5.4, "Y": 6.2}
GR_V = {"A": 31.0, "C": 55.0, "D": 54.0, "E": 83.0, "F": 132.0, "G": 3.0,
        "H": 96.0, "I": 111.0, "K": 119.0, "L": 111.0, "M": 105.0, "N": 56.0,
        "P": 32.5, "Q": 85.0, "R": 124.0, "S": 32.0, "T": 61.0, "V": 84.0,
        "W": 170.0, "Y": 136.0}
CHARGE = {a: 0 for a in AA} | {"D": -1, "E": -1, "K": 1, "R": 1, "H": 1}
# Chou-Fasman helix / strand propensity
HELIX = {"A": 1.42, "C": 0.70, "D": 1.01, "E": 1.51, "F": 1.13, "G": 0.57,
         "H": 1.00, "I": 1.08, "K": 1.16, "L": 1.21, "M": 1.45, "N": 0.67,
         "P": 0.57, "Q": 1.11, "R": 0.98, "S": 0.77, "T": 0.83, "V": 1.06,
         "W": 1.08, "Y": 0.69}
STRAND = {"A": 0.83, "C": 1.19, "D": 0.54, "E": 0.37, "F": 1.38, "G": 0.75,
          "H": 0.87, "I": 1.60, "K": 0.74, "L": 1.30, "M": 1.05, "N": 0.89,
          "P": 0.55, "Q": 1.10, "R": 0.93, "S": 0.75, "T": 1.19, "V": 1.70,
          "W": 1.37, "Y": 1.47}

# Grantham constants, mean-distance normalised as in the original paper
GR_ALPHA, GR_BETA, GR_GAMMA, GR_SCALE = 1.833, 0.1018, 0.000399, 50.723


def grantham(a: str, b: str) -> float:
    dc, dp, dv = GR_C[a] - GR_C[b], GR_P[a] - GR_P[b], GR_V[a] - GR_V[b]
    return GR_SCALE * (GR_ALPHA * dc * dc + GR_BETA * dp * dp
                       + GR_GAMMA * dv * dv) ** 0.5


def blosum62():
    from Bio.Align import substitution_matrices
    m = substitution_matrices.load("BLOSUM62")
    return {(a, b): float(m[a, b]) for a in AA for b in AA}


def endpoint_type(dms_id: str, ref: dict) -> str:
    low = (dms_id + " " + (ref.get("selection_assay") or "")).lower()
    if "abundance" in low or "expression" in low or "proteolysis" in low:
        return "abundance"
    if "trafficking" in low:
        return "trafficking"
    if "activity" in low or "function" in low or "growth" in low:
        return "function"
    if "binding" in low:
        return "binding"
    return "unspecified"


def topology_features(record: dict | None) -> dict[int, dict]:
    """position -> {side, distance, in_tm, tmd_index, loop_len}."""
    if not record:
        return {}
    out: dict[int, dict] = {}
    tmds = record.get("tmds") or []
    for i, t in enumerate(tmds, start=1):
        for p in range(int(t["start"]), int(t["end"]) + 1):
            out[p] = {"side": 0, "distance": 0, "in_tm": 1,
                      "tmd_index": i, "loop_len": 0}
    loops = record.get("loops") or {}
    for loop_id, loop in loops.items():
        idx = int(loop_id)
        side_raw = (loop.get("location") or "").strip().lower()
        side = {"inside": 1, "outside": 2}.get(side_raw, 0)
        start, end = int(loop["start"]), int(loop["end"])
        loop_len = end - start + 1
        bounds = []
        if idx > 1 and idx - 2 < len(tmds):
            bounds.append(("left", int(tmds[idx - 2]["end"])))
        if idx <= len(tmds):
            bounds.append(("right", int(tmds[idx - 1]["start"])))
        for p in range(start, end + 1):
            ds = [p - b if d == "left" else b - p for d, b in bounds]
            ds = [x for x in ds if x > 0]
            if not ds:
                continue
            out[p] = {"side": side, "distance": min(ds), "in_tm": 0,
                      "tmd_index": idx, "loop_len": loop_len}
    return out


def read_am(acc: str):
    path = POOL / f"{acc}.am.tsv"
    if not path.exists():
        return {}, ""
    scores, native_at = {}, {}
    with path.open(encoding="utf-8") as fh:
        next(fh, None)
        for line in fh:
            v, s = line.rstrip("\n").split("\t")
            try:
                pos = int(v[1:-1])
            except ValueError:
                continue
            scores[(pos, v[-1])] = float(s)
            native_at[pos] = v[0]
    if not native_at:
        return scores, ""
    n = max(native_at)
    return scores, "".join(native_at.get(i, "X") for i in range(1, n + 1))


def read_esm(acc: str):
    path = POOL / f"{acc}.esm1b.tsv"
    if not path.exists():
        return {}
    out = {}
    with path.open(encoding="utf-8") as fh:
        next(fh, None)
        for line in fh:
            v, s = line.rstrip("\n").split("\t")
            try:
                out[(int(v[1:-1]), v[-1])] = float(s)
            except ValueError:
                continue
    return out


def load_assay(dms_id: str):
    path = ASSAY_DIR / f"{dms_id}.csv"
    if not path.exists():
        return []
    out = []
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames:
            return []
        if "mutant" not in reader.fieldnames or "DMS_score" not in reader.fieldnames:
            return []
        for row in reader:
            tok = str(row["mutant"]).split(":")
            if len(tok) != 1 or len(tok[0].strip()) < 3:
                continue
            m = tok[0].strip()
            try:
                out.append((m[0], int(m[1:-1]), m[-1], float(row["DMS_score"])))
            except ValueError:
                continue
    return out


def main() -> None:
    from route_a_common import align_construct
    from run_gate1 import read_topology

    manifest = {}
    with (INV / "group_manifest.csv").open(encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh):
            if r["role"] not in CONFIRM_ROLES and r["fitting_eligible"] == "True":
                manifest[r["accession"]] = r

    uni = json.loads(UNIPROT_MAP.read_text(encoding="utf-8"))["entries"]
    refs: dict[str, list[dict]] = {}
    with PGYM_REF.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("taxon") != "Human":
                continue
            acc = (uni.get(row["UniProt_ID"]) or {}).get("accession")
            if acc in manifest:
                refs.setdefault(acc, []).append(row)

    census = read_topology(CENSUS)
    bl = blosum62()

    cols: dict[str, list] = {}
    def push(**kw):
        for k, v in kw.items():
            cols.setdefault(k, []).append(v)

    n_assays = 0
    for acc in sorted(refs):
        am_scores, am_seq = read_am(acc)
        if not am_scores or not am_seq:
            continue
        esm = read_esm(acc)
        topo = topology_features(census.get(acc))
        n_tmds = len((census.get(acc) or {}).get("tmds") or [])
        is_memb = 1 if census.get(acc) else 0

        for ref in refs[acc]:
            dms_id = ref["DMS_id"]
            data = load_assay(dms_id)
            if len(data) < MIN_VARIANTS:
                continue
            construct = (ref.get("target_seq") or "").strip().upper()
            keep = []
            if construct:
                m = align_construct(am_seq, construct)
                placed = 0
                for n, p, a, s in data:
                    if 1 <= p < len(m) and int(m[p]):
                        placed += 1
                        rp = int(m[p])
                        if (rp, a) in am_scores and am_seq[rp - 1] == n:
                            keep.append((n, rp, a, s))
                agreement = len(keep) / placed if placed else 0.0
            else:
                agreement = 0.0
            if agreement < MIN_AGREEMENT or len(keep) < MIN_VARIANTS:
                continue

            dms = np.array([s for *_, s in keep])
            target = rankdata(-dms) / len(dms)
            etype = endpoint_type(dms_id, ref)
            n_assays += 1

            for (native, pos, alt, _), t in zip(keep, target):
                tf = topo.get(pos, {"side": -1, "distance": -1, "in_tm": -1,
                                    "tmd_index": -1, "loop_len": -1})
                push(
                    accession=acc, gene=manifest[acc]["gene"],
                    group_id=manifest[acc]["group_id"], dms_id=dms_id,
                    endpoint_type=etype, is_membrane=is_memb,
                    position=pos, native=native, alt=alt,
                    am=am_scores[(pos, alt)],
                    esm1b=esm.get((pos, alt), np.nan),
                    blosum62=bl.get((native, alt), np.nan),
                    grantham=grantham(native, alt),
                    d_composition=GR_C[alt] - GR_C[native],
                    d_polarity=GR_P[alt] - GR_P[native],
                    d_volume=GR_V[alt] - GR_V[native],
                    d_hydro=HYDRO[alt] - HYDRO[native],
                    d_charge=CHARGE[alt] - CHARGE[native],
                    d_helix=HELIX[alt] - HELIX[native],
                    d_strand=STRAND[alt] - STRAND[native],
                    native_is_kr=1 if native in "KR" else 0,
                    alt_is_kr=1 if alt in "KR" else 0,
                    side=tf["side"], distance=tf["distance"],
                    in_tm=tf["in_tm"], tmd_index=tf["tmd_index"],
                    loop_len=tf["loop_len"], n_tmds=n_tmds,
                    target=t,
                )

    out = FEAT / "l1_table.npz"
    arrays = {}
    for k, v in cols.items():
        a = np.array(v)
        arrays[k] = a
    np.savez_compressed(out, **arrays)

    n = len(cols["target"])
    n_prot = len(set(cols["accession"]))
    esm_arr = np.array(cols["esm1b"], dtype=float)
    (FEAT / "l1_table_manifest.json").write_text(
        json.dumps(
            {
                "rows": n,
                "assays": n_assays,
                "proteins": n_prot,
                "membrane_proteins": len({a for a, m in
                                          zip(cols["accession"], cols["is_membrane"])
                                          if m == 1}),
                "columns": sorted(cols),
                "esm1b_missing_rows": int(np.isnan(esm_arr).sum()),
                "target": "within-assay rank-normalised damage "
                          "(rankdata(-DMS_score)/n)",
                "mapping": "route_a_common.align_construct onto the AM-derived "
                           "reference; no global offset; >=0.99 native gate",
                "scope": "fitting pool only; no confirmation protein present",
                "l1_objective_note": "L1 regresses the within-assay "
                                     "rank-normalised target. A true pairwise "
                                     "ranking loss within assay is an L2 feature, "
                                     "not approximated here.",
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"rows {n:,}  assays {n_assays}  proteins {n_prot}")
    print(f"ESM1b missing on {int(np.isnan(esm_arr).sum()):,} rows "
          f"({np.isnan(esm_arr).mean():.2%})")
    print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

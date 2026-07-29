"""Route B, step B1a: build and verify the ProtT5 chunk plan before any GPU time is spent.

Why this runs first
-------------------
ProtT5 attention is quadratic and the census holds proteins up to 14,507 residues, so sequences are
embedded in overlapping chunks. Only chunks containing panel sites need to be run: 726,444 sites across
4,837 proteins, which is roughly 5,100 chunks rather than 4,863 whole proteins.

The failure mode is silent. An off-by-one in the chunk offset misaligns every embedding to the wrong
residue, the model still trains, and nothing downstream complains. SERT is the precedent: a best-fit
global offset mismapped 215 residues at 66.8% agreement and only an explicit reference-residue check
caught it.

The check here is free and exact. The panel stores `ref_aa` for every site, so for each site the residue
the census sequence carries at that position must equal the residue the panel says is native. Any
mismatch is a mapping bug, and the plan is refused rather than written.

Output
------
`chunk_plan.npz` holds, per panel site, the chunk it belongs to and its offset inside that chunk, plus the
chunk table itself (accession, start, end). The embedding job reads the chunk table; the trainer reads the
per-site index.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import form_ab_lib as lib
from run_gate1 import read_topology

CHUNK = 1022          # ProtT5 context budget, leaving room for the terminator
OVERLAP = 128         # so every panel site sits at least OVERLAP/2 from a chunk edge where possible


def plan_chunks(length: int) -> list[tuple[int, int]]:
    """Half-open [start, end) windows covering 0..length, striding by CHUNK - OVERLAP."""
    if length <= CHUNK:
        return [(0, length)]
    step = CHUNK - OVERLAP
    starts = list(range(0, max(length - CHUNK, 0) + 1, step))
    if starts[-1] + CHUNK < length:
        starts.append(length - CHUNK)
    return [(s, min(s + CHUNK, length)) for s in starts]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, default=lib.ROOT / "results/route_b_20261012")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    census = read_topology(lib.TOPOLOGY)
    panel = lib.P.load_panel(lib.AM_ARCHIVE, 60)
    accessions = [str(a) for a in panel["accessions"]]
    protein = panel["protein"].astype(np.int64)
    position = panel["position"].astype(np.int64)
    ref_aa = panel["ref_aa"].astype(np.int64)

    # --- the correctness check, before anything is built -------------------------------------------
    mismatch = []
    missing = []
    for p in np.unique(protein):
        acc = accessions[p]
        if acc not in census:
            missing.append(acc)
            continue
        seq = census[acc]["sequence"]
        rows = np.flatnonzero(protein == p)
        pos, ref = position[rows], ref_aa[rows]
        inside = (pos >= 1) & (pos <= len(seq))
        if not inside.all():
            mismatch.append({"accession": acc, "reason": "position outside sequence",
                             "n": int((~inside).sum()), "length": len(seq),
                             "max_position": int(pos.max())})
            continue
        native = np.array([lib.ALPHABET.index(seq[q - 1]) if seq[q - 1] in lib.ALPHABET else -1
                           for q in pos])
        bad = native != ref
        if bad.any():
            mismatch.append({"accession": acc, "reason": "native residue differs",
                             "n": int(bad.sum()), "of": int(len(pos)),
                             "first_position": int(pos[bad][0]),
                             "census_residue": seq[pos[bad][0] - 1],
                             "panel_residue": lib.ALPHABET[ref[bad][0]]})

    report = {
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "note": "exploratory chunk planning",
        "chunk": CHUNK, "overlap": OVERLAP,
        "n_panel_sites": int(len(protein)), "n_panel_proteins": int(len(np.unique(protein))),
        "n_census_proteins": len(census),
        "accessions_missing_from_census": missing,
        "reference_residue_mismatches": mismatch,
        "topology_sha256": lib.P.sha256(lib.TOPOLOGY),
        "archive_sha256": lib.P.sha256(lib.AM_ARCHIVE),
        "implementation_sha256": lib.P.sha256(Path(__file__)),
    }

    if mismatch or missing:
        (args.out / "chunk_plan_refused.json").write_text(json.dumps(report, indent=2) + "\n",
                                                          encoding="utf-8")
        raise SystemExit(f"Refusing to write a chunk plan: {len(mismatch)} proteins with reference "
                         f"mismatches, {len(missing)} missing from the census. "
                         f"See {args.out / 'chunk_plan_refused.json'}")
    print(f"reference residues verified for all {len(protein):,} panel sites", flush=True)

    # --- build the plan ---------------------------------------------------------------------------
    chunk_acc, chunk_start, chunk_end, chunk_seq = [], [], [], []
    site_chunk = np.full(len(protein), -1, dtype=np.int64)
    site_offset = np.full(len(protein), -1, dtype=np.int64)

    for p in np.unique(protein):
        acc = accessions[p]
        seq = census[acc]["sequence"]
        windows = plan_chunks(len(seq))
        rows = np.flatnonzero(protein == p)
        pos0 = position[rows] - 1            # 0-based residue index
        # For each site pick the chunk whose centre it is closest to, so sites avoid chunk edges.
        centres = np.array([(s + e) / 2 for s, e in windows])
        starts = np.array([s for s, _ in windows])
        ends = np.array([e for _, e in windows])
        contains = (pos0[:, None] >= starts[None, :]) & (pos0[:, None] < ends[None, :])
        distance = np.abs(pos0[:, None] - centres[None, :])
        distance = np.where(contains, distance, np.inf)
        pick = distance.argmin(axis=1)
        if not np.isfinite(distance[np.arange(len(pos0)), pick]).all():
            raise RuntimeError(f"{acc}: a panel site falls in no chunk")
        used = {}
        for w in np.unique(pick):
            used[int(w)] = len(chunk_acc)
            chunk_acc.append(acc)
            chunk_start.append(int(windows[w][0]))
            chunk_end.append(int(windows[w][1]))
            # Carry the subsequence itself so the embedding job needs only numpy, not the analysis stack.
            chunk_seq.append(seq[windows[w][0]:windows[w][1]])
        site_chunk[rows] = [used[int(w)] for w in pick]
        site_offset[rows] = pos0 - np.array([windows[w][0] for w in pick])

    assert (site_chunk >= 0).all() and (site_offset >= 0).all()
    lengths = np.array(chunk_end) - np.array(chunk_start)
    assert (site_offset < lengths[site_chunk]).all(), "offset past the end of its chunk"

    # The embedding job reads this file with numpy alone; `alphabet` lets it check its own gather.
    np.savez_compressed(
        args.out / "chunk_plan.npz",
        chunk_accession=np.array(chunk_acc), chunk_start=np.array(chunk_start),
        chunk_end=np.array(chunk_end), chunk_sequence=np.array(chunk_seq),
        site_chunk=site_chunk, site_offset=site_offset,
        site_protein=protein, site_position=position, site_ref_aa=ref_aa,
        alphabet=np.array(list(lib.ALPHABET)))

    report["n_chunks"] = len(chunk_acc)
    report["tokens_total"] = int(lengths.sum())
    report["chunk_length"] = {"min": int(lengths.min()), "median": int(np.median(lengths)),
                              "max": int(lengths.max())}
    report["embedding_bytes_fp16_1024"] = int(len(protein) * 1024 * 2)
    (args.out / "chunk_plan.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"chunks          : {len(chunk_acc):,}")
    print(f"tokens total    : {lengths.sum():,}")
    print(f"chunk length    : min {lengths.min()} median {int(np.median(lengths))} max {lengths.max()}")
    print(f"embedding output: {len(protein) * 1024 * 2 / 1e9:.2f} GB at fp16 x 1024 dims")
    print(f"wrote           : {args.out / 'chunk_plan.npz'}")


if __name__ == "__main__":
    main()

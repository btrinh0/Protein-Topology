"""Per-residue embeddings for every measured site in the fitting pool.

The cached panel embeddings cover flank sites only (726,444 sites at census
distance <=60), which reach just 9.1% of the L1 table's rows and none of the
83,820 non-membrane rows. This embeds the sites that actually carry measured
outcomes.

Recipe matched EXACTLY to the cached panel so new and old vectors are
interchangeable (`scripts/route_b_embed_modal.py`):
  Rostlab/prot_t5_xl_uniref50, T5Tokenizer(do_lower_case=False, legacy=True),
  T5EncoderModel in float16, residues space-separated, U/Z/O/B -> X per the
  model card, add_special_tokens=True, and row r of the output is residue r
  because the tokenizer appends </s> AFTER the last residue.
Chunking also matches: window 1022, overlap 128.

Frontend choice. ProtT5 is the default here and ESM2 is opt-in behind
--model esm2, because of the project's own documented reasoning in
route_b_embed_modal.py: ESM1b already satisfies the K/R-gain rule
(comparator-robust at -0.073), so building on ESM representations risks rule
compliance arriving from the frontend rather than from the rule loss term, which
would make the central claim unfalsifiable. ESM2 is the same lineage and the
same risk. ProtT5 is a different lineage.

Coordinates. Reference sequences are rebuilt from the extracted AlphaMissense
wild-type residues, exactly as a8/a9 did, so site positions line up with the L1
table with no offset.

Sharding and resume. --shard k --n-shards n takes every n-th accession, so one
process per GPU pinned with CUDA_VISIBLE_DEVICES. One output file per accession,
skipped if present, so a rerun only does what is missing.

Label-free: reads positions, never outcomes.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

ROOT = Path("/home/ubuntu/project6_membrane_variant_topology_audit")
FEAT = ROOT / "results" / "beat_am_20261016" / "features"
POOL = FEAT / "fitting_pool"

CHUNK, OVERLAP = 1022, 128
PROT_T5 = "Rostlab/prot_t5_xl_uniref50"
ESM2 = "facebook/esm2_t33_650M_UR50D"


def reference_sequence(acc: str) -> str:
    """Rebuild AM's own reference sequence from the extracted wild-type residues."""
    path = POOL / f"{acc}.am.tsv"
    if not path.exists():
        return ""
    native_at: dict[int, str] = {}
    with path.open() as fh:
        next(fh, None)
        for line in fh:
            v = line.split("\t", 1)[0]
            try:
                native_at[int(v[1:-1])] = v[0]
            except ValueError:
                continue
    if not native_at:
        return ""
    return "".join(native_at.get(i, "X") for i in range(1, max(native_at) + 1))


def sites_needed() -> dict[str, list[int]]:
    t = np.load(FEAT / "l1_table.npz", allow_pickle=True)
    acc = t["accession"].astype(str)
    pos = t["position"].astype(int)
    out: dict[str, set[int]] = {}
    for a, p in zip(acc, pos):
        out.setdefault(a, set()).add(int(p))
    return {a: sorted(v) for a, v in out.items()}


def chunks_for(length: int) -> list[tuple[int, int]]:
    """[start, end) windows of CHUNK with OVERLAP, covering 1..length."""
    if length <= CHUNK:
        return [(0, length)]
    step = CHUNK - OVERLAP
    spans = []
    start = 0
    while start < length:
        end = min(start + CHUNK, length)
        spans.append((start, end))
        if end == length:
            break
        start += step
    return spans


def assign(positions: list[int], spans: list[tuple[int, int]]) -> dict[int, int]:
    """Each 1-based position -> index of the span giving it the most context."""
    out = {}
    for p in positions:
        i0 = p - 1
        best, best_margin = None, -1
        for k, (s, e) in enumerate(spans):
            if s <= i0 < e:
                margin = min(i0 - s, e - 1 - i0)
                if margin > best_margin:
                    best, best_margin = k, margin
        if best is not None:
            out[p] = best
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", choices=["prot_t5", "esm2"], default="prot_t5")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--out", default=None)
    ap.add_argument("--all-positions", action="store_true",
                    help="embed every residue, not only measured sites, so L2 can "
                         "build a local window around each position; the GPU cost "
                         "is identical because the same chunks are forwarded either "
                         "way, only the gather widens")
    args = ap.parse_args()

    import torch

    out_dir = Path(args.out) if args.out else FEAT / f"embeddings_{args.model}"
    out_dir.mkdir(parents=True, exist_ok=True)

    needed = sites_needed()
    accs = sorted(needed)
    mine = [a for i, a in enumerate(accs) if i % args.n_shards == args.shard]
    print(f"model={args.model} shard {args.shard}/{args.n_shards}: "
          f"{len(mine)} of {len(accs)} accessions", flush=True)

    if args.model == "prot_t5":
        from transformers import T5EncoderModel, T5Tokenizer
        tok = T5Tokenizer.from_pretrained(PROT_T5, do_lower_case=False, legacy=True)
        model = T5EncoderModel.from_pretrained(PROT_T5, torch_dtype=torch.float16)
        dims = 1024
    else:
        from transformers import AutoModel, AutoTokenizer
        tok = AutoTokenizer.from_pretrained(ESM2)
        model = AutoModel.from_pretrained(ESM2, torch_dtype=torch.float16)
        dims = 1280
    model = model.cuda().eval()
    print(f"  model loaded, dims={dims}", flush=True)

    t0 = time.time()
    done = skipped = 0
    for acc in mine:
        out_path = out_dir / f"{acc}.npz"
        if out_path.exists():
            skipped += 1
            continue
        seq = reference_sequence(acc)
        if not seq:
            print(f"  {acc}: no reference sequence, skipped", flush=True)
            continue
        positions = (list(range(1, len(seq) + 1)) if args.all_positions
                     else [p for p in needed[acc] if 1 <= p <= len(seq)])
        spans = chunks_for(len(seq))
        amap = assign(positions, spans)
        if not amap:
            continue

        vectors = np.zeros((len(positions), dims), dtype=np.float16)
        index = {p: i for i, p in enumerate(positions)}

        by_span: dict[int, list[int]] = {}
        for p, k in amap.items():
            by_span.setdefault(k, []).append(p)

        with torch.no_grad():
            for k, ps in by_span.items():
                s, e = spans[k]
                sub = seq[s:e]
                if args.model == "prot_t5":
                    spaced = " ".join(
                        "X" if c in "UZOB" else c for c in sub
                    )
                    enc = tok(spaced, return_tensors="pt", padding=True,
                              add_special_tokens=True)
                    enc = {kk: vv.cuda() for kk, vv in enc.items()}
                    rep = model(**enc).last_hidden_state[0]
                    # row r is residue r; </s> follows the last residue
                    offset = 0
                else:
                    enc = tok(sub, return_tensors="pt", add_special_tokens=True)
                    enc = {kk: vv.cuda() for kk, vv in enc.items()}
                    rep = model(**enc).last_hidden_state[0]
                    # ESM2 prepends <cls>, so residue r is row r+1
                    offset = 1
                rep = rep.float().cpu().numpy()
                for p in ps:
                    row = (p - 1 - s) + offset
                    if 0 <= row < rep.shape[0]:
                        vectors[index[p]] = rep[row].astype(np.float16)

        np.savez_compressed(out_path, position=np.array(positions, dtype=np.int32),
                            embedding=vectors, model=args.model,
                            seq_len=len(seq), n_chunks=len(spans))
        done += 1
        if done % 5 == 0:
            print(f"  {done} done, {time.time()-t0:.0f}s elapsed", flush=True)

    meta = {
        "model": PROT_T5 if args.model == "prot_t5" else ESM2,
        "dims": dims, "chunk": CHUNK, "overlap": OVERLAP,
        "shard": args.shard, "n_shards": args.n_shards,
        "accessions_written": done, "accessions_skipped_existing": skipped,
        "seconds": round(time.time() - t0, 1),
        "recipe_matches_cached_panel": args.model == "prot_t5",
        "coordinates": "AM-derived reference sequence, same as a8/a9",
        "label_free": True,
    }
    (out_dir / f"_meta_shard{args.shard}.json").write_text(json.dumps(meta, indent=2))
    print(f"done: wrote {done}, skipped {skipped}, {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

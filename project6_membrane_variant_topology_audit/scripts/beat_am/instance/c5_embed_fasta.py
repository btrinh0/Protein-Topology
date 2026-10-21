"""Embed every residue of every sequence in a FASTA with ProtT5-XL.

Standalone: needs only the FASTA and the model, not the project bundle. Used to
embed the confirmation proteins without re-uploading 4.13 GB.

Recipe matched exactly to route_b_embed_modal.py and c1_embed_pool.py, so these
vectors are interchangeable with the cached panel (validated at cosine 1.0000):
Rostlab/prot_t5_xl_uniref50, T5Tokenizer(do_lower_case=False, legacy=True),
T5EncoderModel in float16, residues space-separated, U/Z/O/B -> X, row r is
residue r because the tokenizer appends </s> after the last residue.
Chunking: window 1022, overlap 128; each position takes the chunk giving it the
most context.

--shard k --n-shards n splits sequences across GPUs. One npz per accession,
skipped if present, so reruns only do what is missing.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

CHUNK, OVERLAP = 1022, 128
MODEL = "Rostlab/prot_t5_xl_uniref50"


def read_fasta(path: Path) -> dict[str, str]:
    out, name, buf = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].strip(), []
        elif line.strip():
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


def chunks_for(length: int) -> list[tuple[int, int]]:
    if length <= CHUNK:
        return [(0, length)]
    step, spans, start = CHUNK - OVERLAP, [], 0
    while start < length:
        end = min(start + CHUNK, length)
        spans.append((start, end))
        if end == length:
            break
        start += step
    return spans


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    args = ap.parse_args()

    import torch
    from transformers import T5EncoderModel, T5Tokenizer

    seqs = read_fasta(Path(args.fasta))
    names = sorted(seqs)
    mine = [n for i, n in enumerate(names) if i % args.n_shards == args.shard]
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"shard {args.shard}/{args.n_shards}: {len(mine)} of {len(names)} "
          f"sequences, {sum(len(seqs[n]) for n in mine):,} residues", flush=True)

    tok = T5Tokenizer.from_pretrained(MODEL, do_lower_case=False, legacy=True)
    model = T5EncoderModel.from_pretrained(MODEL, torch_dtype=torch.float16)
    model = model.cuda().eval()
    print("  model loaded", flush=True)

    t0, done = time.time(), 0
    for name in mine:
        out_path = out_dir / f"{name}.npz"
        if out_path.exists():
            continue
        seq = seqs[name]
        spans = chunks_for(len(seq))
        positions = list(range(1, len(seq) + 1))
        vectors = np.zeros((len(positions), 1024), dtype=np.float16)

        assign: dict[int, int] = {}
        for p in positions:
            i0, best, margin = p - 1, None, -1
            for k, (s, e) in enumerate(spans):
                if s <= i0 < e:
                    m = min(i0 - s, e - 1 - i0)
                    if m > margin:
                        best, margin = k, m
            if best is not None:
                assign[p] = best
        by_span: dict[int, list[int]] = {}
        for p, k in assign.items():
            by_span.setdefault(k, []).append(p)

        with torch.no_grad():
            for k, ps in by_span.items():
                s, e = spans[k]
                spaced = " ".join("X" if c in "UZOB" else c for c in seq[s:e])
                enc = tok(spaced, return_tensors="pt", padding=True,
                          add_special_tokens=True)
                enc = {kk: vv.cuda() for kk, vv in enc.items()}
                rep = model(**enc).last_hidden_state[0].float().cpu().numpy()
                for p in ps:
                    row = p - 1 - s
                    if 0 <= row < rep.shape[0]:
                        vectors[p - 1] = rep[row].astype(np.float16)

        np.savez_compressed(out_path, position=np.array(positions, dtype=np.int32),
                            embedding=vectors, seq_len=len(seq),
                            n_chunks=len(spans), model="prot_t5")
        done += 1
        print(f"  {done}/{len(mine)} {name} len={len(seq)} "
              f"({time.time()-t0:.0f}s)", flush=True)

    (out_dir / f"_meta_fasta_shard{args.shard}.json").write_text(json.dumps({
        "model": MODEL, "dims": 1024, "chunk": CHUNK, "overlap": OVERLAP,
        "shard": args.shard, "n_shards": args.n_shards, "written": done,
        "seconds": round(time.time() - t0, 1),
        "recipe_matches_cached_panel": True, "label_free": True}, indent=2))
    print(f"done: {done} written in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

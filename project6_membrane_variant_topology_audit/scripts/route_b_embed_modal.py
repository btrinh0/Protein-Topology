"""Route B, step B1b: ProtT5-XL per-residue embeddings for the panel sites, on Modal.

Reads `results/route_b_20261012/chunk_plan.npz`, which `route_b_chunk_plan.py` wrote after verifying every
one of the 726,444 panel sites against the census sequence. Embeds the 5,043 chunks, gathers the
per-residue vector for each panel site, and writes one fp16 array per batch to a Modal Volume.

The gather happens inside the GPU function on purpose. Returning whole chunk embeddings would move
2,605,296 x 1024 x 2 bytes, about 5.3 GB; gathering first reduces it to 726,444 x 1024 x 2, about 1.49 GB.

Frontend choice. ProtT5-XL-UniRef50, not ESM. ESM1b already satisfies the K/R-gain rule
(comparator-robust at -0.073, 2026-09-29 entry), so building on ESM representations risks rule compliance
arriving from the frontend rather than from the loss term, which would make the central claim
unfalsifiable. ESM2 is the same family and the same risk. ProtT5 is a different lineage, and the plain
student control measures its baseline contrast rather than assuming anything.

Run
---
    pip install modal && modal setup          # one time, authenticates to your account
    modal run scripts/route_b_embed_modal.py

Then fetch the output and assemble it:

    modal volume get route-b-embeddings / ./results/route_b_20261012/embeddings
    python scripts/route_b_embed_modal.py --assemble

Cost and GPU choice
-------------------
The weights are baked into the image by `_bake_model`, so the ~11 GB pull happens once at image build and
no container ever downloads them. Chunks are length-bucketed into padded forwards under TOKEN_BUDGET
tokens, which turns 5,043 single-sequence forwards into a few hundred batched ones; without that the job
is launch-overhead-bound and a faster GPU buys nothing.

With both fixes, 2,605,296 tokens is roughly 1-2 minutes of compute plus about a minute of model load:

    A100 40GB  $2.10/hr   ~4 min container   ~$0.14
    H100 SXM5  $3.95/hr   ~3 min container   ~$0.20
    L4 24GB    $0.80/hr   ~9 min container   ~$0.12

The fixed overhead (cold start, model load to GPU) does not shrink with a faster card, so H100 saves about
a minute over A100 and costs a few cents more. Any of the three is fine; pick whichever has capacity.
Modal bills per second and the batch files are idempotent, so a failed or partial run is cheap and a rerun
only redoes missing batches.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

# One task for the whole job. The model is loaded inside the function, so each task costs another load
# (about a minute); with only ~2 minutes of compute in total, splitting the work would cost more in
# repeated loads than it saves in parallelism. 5,043 chunks pickle to ~8 MB, which Modal moves as a blob.
BATCH = 6000
TOKEN_BUDGET = 32768             # padded tokens per GPU forward, length-bucketed to limit padding waste
MAX_CONTAINERS = 1
GPU = "A100"                     # "H100" and "L4" both work; see the note in the module docstring
MODEL = "Rostlab/prot_t5_xl_uniref50"
VOLUME = "route-b-embeddings"
PLAN = "results/route_b_20261012/chunk_plan.npz"
OUT = "results/route_b_20261012"

try:
    import modal
except ImportError:                                  # --assemble runs without modal installed
    modal = None


ROOT = Path(__file__).resolve().parents[1]


def _bake_model() -> None:
    """Pull the ~11 GB of ProtT5 weights at image build time so no container ever downloads them."""
    from transformers import T5EncoderModel, T5Tokenizer
    T5Tokenizer.from_pretrained(MODEL, do_lower_case=False, legacy=True)
    T5EncoderModel.from_pretrained(MODEL)


def load_chunks() -> list[dict]:
    """Read the verified plan with numpy alone, and re-check the gather before sending anything.

    This module deliberately does not import `form_ab_lib`, which keeps the Modal CLI environment
    independent from the local analysis environment.
    """
    import numpy as np

    plan = np.load(ROOT / PLAN, allow_pickle=False)
    if "chunk_sequence" not in plan:
        raise SystemExit("chunk_plan.npz predates sequence export; rerun scripts/route_b_chunk_plan.py")
    sequences = plan["chunk_sequence"]
    site_chunk, site_offset, site_ref = plan["site_chunk"], plan["site_offset"], plan["site_ref_aa"]
    alphabet = "".join(plan["alphabet"].tolist())

    by_chunk: dict[int, list[int]] = {}
    for row, c in enumerate(site_chunk):
        by_chunk.setdefault(int(c), []).append(row)

    chunks, bad = [], 0
    for c in range(len(sequences)):
        seq = str(sequences[c])
        rows = by_chunk.get(c, [])
        offsets = [int(site_offset[r]) for r in rows]
        for r, o in zip(rows, offsets):
            if o >= len(seq) or seq[o] != alphabet[int(site_ref[r])]:
                bad += 1
        chunks.append({"chunk_id": c, "sequence": seq, "site_rows": rows, "site_offsets": offsets})
    if bad:
        raise SystemExit(f"{bad} sites would read the wrong residue; the plan is not usable")
    gathered = sum(len(c["site_rows"]) for c in chunks)
    if gathered != len(site_chunk):
        raise SystemExit(f"plan covers {gathered} sites but has {len(site_chunk)} rows")
    print(f"plan verified: {len(chunks):,} chunks, {gathered:,} sites, 0 alignment mismatches")
    return chunks


if modal is not None:
    app = modal.App("route-b-prott5")
    volume = modal.Volume.from_name(VOLUME, create_if_missing=True)
    image = (modal.Image.debian_slim(python_version="3.12")
             .pip_install("torch==2.4.1", "transformers==4.44.2", "sentencepiece==0.2.0",
                          "numpy==2.1.1", "protobuf==5.28.0", "hf_transfer==0.1.8")
             .env({"HF_HUB_ENABLE_HF_TRANSFER": "1"})
             .run_function(_bake_model))

    @app.function(image=image, gpu=GPU, volumes={"/out": volume}, timeout=3600,
                  max_containers=MAX_CONTAINERS, retries=modal.Retries(max_retries=2))
    def embed_batch(batch_id: int, chunks: list[dict]) -> dict:
        """Embed a batch of chunks and gather one vector per panel site inside them."""
        import numpy as np
        import torch
        from transformers import T5EncoderModel, T5Tokenizer

        out_path = Path(f"/out/batch_{batch_id:05d}.npz")
        if out_path.exists():                        # idempotent across reruns
            return {"batch": batch_id, "skipped": True}

        tokenizer = T5Tokenizer.from_pretrained(MODEL, do_lower_case=False, legacy=True)
        model = T5EncoderModel.from_pretrained(MODEL, torch_dtype=torch.float16).cuda().eval()

        # Length-bucket so each forward pads to roughly its own longest sequence, not to 1022.
        order = sorted(range(len(chunks)), key=lambda i: len(chunks[i]["sequence"]))
        groups, current = [], []
        for i in order:
            longest = max(len(chunks[j]["sequence"]) for j in current + [i])
            if current and (len(current) + 1) * (longest + 1) > TOKEN_BUDGET:
                groups.append(current)
                current = [i]
            else:
                current.append(i)
        if current:
            groups.append(current)

        rows, vectors = [], []
        with torch.no_grad():
            for group in groups:
                seqs = [chunks[i]["sequence"] for i in group]
                # ProtT5 wants space-separated residues; U/Z/O/B map to X per the model card.
                spaced = [" ".join("X" if r in "UZOB" else r for r in s) for s in seqs]
                encoded = tokenizer(spaced, return_tensors="pt", padding=True, add_special_tokens=True)
                hidden = model(input_ids=encoded["input_ids"].cuda(),
                               attention_mask=encoded["attention_mask"].cuda()).last_hidden_state
                for slot, i in enumerate(group):
                    seq = chunks[i]["sequence"]
                    # Row r is residue r; the tokenizer appends </s> after the last residue.
                    residues = hidden[slot, : len(seq)].float().cpu().numpy()
                    if residues.shape[0] != len(seq):
                        raise RuntimeError(f"chunk {chunks[i]['chunk_id']}: {residues.shape[0]} rows "
                                           f"for {len(seq)} residues")
                    for site_row, offset in zip(chunks[i]["site_rows"], chunks[i]["site_offsets"]):
                        rows.append(site_row)
                        vectors.append(residues[offset])
        np.savez(out_path, site_row=np.array(rows, dtype=np.int64),
                 embedding=np.asarray(vectors, dtype=np.float16))
        volume.commit()
        return {"batch": batch_id, "chunks": len(chunks), "groups": len(groups), "sites": len(rows)}

    @app.local_entrypoint()
    def main() -> None:
        chunks = load_chunks()
        batches = [(i // BATCH, chunks[i:i + BATCH]) for i in range(0, len(chunks), BATCH)]
        print(f"{len(chunks):,} chunks in {len(batches)} task(s) on {GPU}")
        done = list(embed_batch.starmap(batches))
        sites = sum(d.get("sites", 0) for d in done)
        print(f"finished {len(done)} task(s), {sites:,} site vectors written to volume {VOLUME!r}")
        print(f"next:  modal volume get {VOLUME} / ./{OUT}/embeddings")


def assemble() -> None:
    """Stitch the downloaded batch files into one (n_sites, 1024) fp16 array, in panel-row order."""
    import numpy as np

    root = ROOT / OUT
    plan = np.load(root / "chunk_plan.npz", allow_pickle=False)
    n_sites = len(plan["site_chunk"])
    # `modal volume get <vol> / <dest>` writes a lone volume file *as* <dest> rather than into a
    # directory of that name, so accept both shapes.
    target = root / "embeddings"
    if target.is_dir():
        files = sorted(target.glob("batch_*.npz"))
    elif target.is_file():
        files = [target]
    else:
        files = []
    if not files:
        raise SystemExit(f"No batch files at {target}; run the modal job and fetch first")

    table = np.zeros((n_sites, 1024), dtype=np.float16)
    filled = np.zeros(n_sites, dtype=bool)
    for path in files:
        with np.load(path, allow_pickle=False) as blob:
            rows, vectors = blob["site_row"], blob["embedding"]
        table[rows] = vectors
        filled[rows] = True
    if not filled.all():
        raise SystemExit(f"{(~filled).sum():,} of {n_sites:,} sites have no embedding; "
                         f"{len(files)} batch files present. Fetch the rest and rerun --assemble.")

    np.save(root / "site_embeddings_fp16.npy", table)
    meta = {"n_sites": int(n_sites), "dims": 1024, "dtype": "float16",
            "batch_files": len(files), "model": MODEL,
            "bytes": int(table.nbytes),
            "note": "row i corresponds to row i of the panel used by route_b_chunk_plan.py"}
    (root / "site_embeddings.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print(f"assembled {n_sites:,} x 1024 fp16 ({table.nbytes / 1e9:.2f} GB) "
          f"-> {root / 'site_embeddings_fp16.npy'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--assemble", action="store_true",
                        help="stitch downloaded batch files into one array (no GPU, no modal needed)")
    args = parser.parse_args()
    if args.assemble:
        assemble()
    else:
        raise SystemExit("Run the GPU job with:  modal run scripts/route_b_embed_modal.py\n"
                         "Then:                  modal volume get route-b-embeddings / "
                         "./results/route_b_20261012/embeddings\n"
                         "Then:                  python scripts/route_b_embed_modal.py --assemble")

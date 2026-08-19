"""Generate synthetic inputs for the audit command-line interfaces."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
NEUTRAL = "ACFGILMNPQSTVWY"
POSITIONS = {"A": (7, 25), "N": (9, 27), "Q": (11, 29), "V": (13, 31),
             "K": (8, 32), "R": (10, 33), "D": (12, 34), "E": (14, 35)}


def generate(out: Path) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    topology, rows = [], []
    arrays = {k: [] for k in ("protein", "position", "side", "distance", "tmd_index", "ref_aa", "alt_aa", "score")}
    for protein in range(2):
        accession = f"TOY{protein + 1:03d}"
        sequence = list("A" * 50)
        sequence[14:24] = list("L" * 10)
        for ref, positions in POSITIONS.items():
            for position in positions:
                sequence[position - 1] = ref
        topology.append({"uniprot": accession, "gene": accession, "sequence": "".join(sequence),
                         "tmds": [{"start": 15, "end": 24}], "numTMDs": 1,
                         "loops": {"1": {"start": 1, "end": 14, "location": "inside"},
                                   "2": {"start": 25, "end": 50, "location": "outside"}}})
        for ref, positions in POSITIONS.items():
            for side, position in enumerate(positions):
                distance = 15 - position if side == 0 else position - 24
                for alt_index, alt in enumerate(ALPHABET):
                    if alt == ref:
                        continue
                    score = 0.3 + 0.15 * protein + (0.5 if side == 0 else 0.0)
                    score += 0.01 * (alt_index - 9.5)
                    if ref in NEUTRAL and alt in "KR":
                        score += -0.10 if side == 0 else 0.20
                    if ref in "KR" and alt in NEUTRAL:
                        score += 0.15 if side == 0 else 0.05
                    rows.append({"uniprot_id": accession, "protein_variant": f"{ref}{position}{alt}", "score": score})
                    for key, value in {"protein": protein, "position": position, "side": side,
                        "distance": distance, "tmd_index": 0, "ref_aa": ALPHABET.index(ref),
                        "alt_aa": alt_index, "score": score}.items():
                        arrays[key].append(value)
    (out / "toy_topology.json").write_text(json.dumps(topology, indent=2) + "\n", encoding="utf-8")
    with (out / "toy_scores.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    np.savez_compressed(out / "toy_mapped_scores.npz", **{k: np.asarray(v) for k, v in arrays.items()},
                        accessions=np.array([x["uniprot"] for x in topology]), aa_alphabet=np.array(list(ALPHABET)))
    return {"synthetic": True, "variants": len(rows), "proteins": len(topology),
            "expected_raw_KR_gain_side_contrast": -0.30, "out": str(out.resolve())}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path(__file__).parent / "toy")
    args = parser.parse_args()
    print(json.dumps(generate(args.out), indent=2))

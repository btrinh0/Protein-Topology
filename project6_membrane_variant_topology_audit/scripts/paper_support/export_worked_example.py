"""Export one topology-selected illustration from saved prediction scores."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from membrane_rule_audit.audit import flank_map

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--out", type=Path, default=ROOT / "results/worked_example")
OUT = parser.parse_args().out
OUT.mkdir(parents=True, exist_ok=True)
SCORES = ROOT / "results/standalone_multirule_20261025/label_free_audit_20261025/standalone_position_scores.npz"
BANK = ROOT / "results/standalone_multirule_20261025/annotation_bank.json"
bank = json.loads(BANK.read_text(encoding="utf-8"))
with np.load(SCORES, allow_pickle=False) as blob:
    acc, pos, native = (blob[k].copy() for k in ("accession", "position", "native"))
    alphabet = blob["alphabet"].astype(str).tolist()
    score = blob["guarded_oof"].copy()

candidates = []
for accession in sorted(bank):
    coordinates = flank_map(bank[accession]["topology"], 10)
    for i in np.flatnonzero((acc == accession) & (native == "G")):
        coordinate = coordinates.get(int(pos[i]))
        if coordinate is not None and coordinate[0] == 0 and coordinate[1] == 6:
            assert np.isfinite(score[i]).sum() == 19
            candidates.append((accession, int(pos[i]), int(i), coordinate))
assert candidates, "No topology-matched illustrative position was found."
accession, position, row, coordinate = sorted(candidates)[0]
assert bank[accession]["sequence"][position - 1] == "G"
uncharged = "ACFGILMNPQSTVWY"
alternatives = {aa: float(s) for aa, s in zip(alphabet, score[row]) if aa != "G"}
kr_mean = float(np.mean([alternatives[a] for a in "KR"]))
comparators = [a for a in uncharged if a != "G"]
comparator_mean = float(np.mean([alternatives[a] for a in comparators]))
assert len(alternatives) == 19 and len(comparators) == 14
rows = [dict(accession=accession, position=position, native="G", alt=a,
             variant=f"G{position}{a}", score=s,
             arm="K/R" if a in "KR" else "comparator" if a in comparators else "excluded D/E/H")
        for a, s in sorted(alternatives.items())]
with (OUT / "worked_example_scores.csv").open("w", encoding="utf-8-sig", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
summary = dict(accession=accession, position=position, native="G", side="cytosolic",
               distance_to_nearest_TM_boundary=coordinate[1], tmd_index=coordinate[2],
               score_source="unchanged standalone guarded group-out predictions",
               selection="First accession/position with native G, cytosolic side and distance 6; scores did not enter selection.",
               eligible_positions=len(candidates),
               K_score=alternatives["K"], R_score=alternatives["R"], KR_mean=kr_mean,
               comparator_residues="".join(comparators), comparator_count=len(comparators),
               comparator_mean=comparator_mean, z=kr_mean-comparator_mean,
               score_sha256=hashlib.sha256(SCORES.read_bytes()).hexdigest(),
               annotation_sha256=hashlib.sha256(BANK.read_bytes()).hexdigest())
(OUT / "worked_example.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
print(json.dumps(summary, indent=2))

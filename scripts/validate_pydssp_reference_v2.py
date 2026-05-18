from __future__ import annotations

"""Validate the NumPy PyDSSP sensitivity implementation on its bundled reference set."""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(1, str(Path(__file__).parent.parent / ".deps" / "pydssp-source"))

import pydssp


def read_reference(path: Path) -> np.ndarray:
    values = []
    reading = False
    conversion = {" ": 0, "S": 0, "T": 0, "H": 1, "G": 1, "I": 1, "E": 2, "B": 2}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "!" in line:
            continue
        if line.startswith("  # "):
            reading = True
            continue
        if reading:
            values.append(conversion[line[16]])
    return np.asarray(values, dtype=int)


def main() -> int:
    root = Path(__file__).parent.parent / ".deps" / "pydssp-source" / "tests" / "testset" / "TS50"
    targets = [line.strip() for line in (root / "list").read_text(encoding="utf-8").splitlines() if line.strip()]
    correlations = []
    for target in targets:
        reference = read_reference(root / "dssp" / f"{target}.dssp")
        coordinates, sequence = pydssp.read_pdbtext((root / "pdb" / f"{target}.pdb").read_text(encoding="utf-8"), return_sequence=True)
        predicted = pydssp.assign(coordinates, donor_mask=sequence != "PRO", out_type="index")
        correlations.append(float((reference == predicted).mean()))
    summary = {
        "reference_set": "PyDSSP TS50",
        "structures": len(correlations),
        "mean_residue_state_agreement": float(np.mean(correlations)),
        "minimum_residue_state_agreement": float(np.min(correlations)),
        "threshold": 0.97,
        "passes_mean_threshold": bool(np.mean(correlations) > 0.97),
        "assignment_role": "sensitivity_only",
    }
    output = Path("data/processed/v2/pydssp_reference_validation_v2.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Validated {len(correlations)} PyDSSP reference structures at mean agreement {summary['mean_residue_state_agreement']:.5f}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

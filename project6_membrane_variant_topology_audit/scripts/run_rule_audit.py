"""Demonstrate the reusable tool on two models using the earlier raw settings.

Only the original raw-score, all-uncharged K/R-gain comparison is run.
No comparator, scale, subgroup or simulation robustness analyses are invoked.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"tools"))
from membrane_rule_audit.panel import run_panel_audit, sha256


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    base = ROOT/"data"/"processed"/"gate0"
    parser.add_argument("--am", type=Path, default=base/"am_full_flank_scores.npz")
    parser.add_argument("--esm", type=Path, default=base/"esm1b_flank_scores.npz")
    parser.add_argument("--out", type=Path, default=ROOT/"results"/"tool_demo_20261003")
    parser.add_argument("--draws", type=int, default=1999)
    parser.add_argument("--seed", type=int, default=20261003)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    results = {}
    for model, path in (("alphamissense", args.am), ("esm1b", args.esm)):
        print(f"Running packaged raw-score comparison: {model}", flush=True)
        # Both stored archives already orient scores as higher-is-damaging.
        results[model] = run_panel_audit(path, args.out/model, model,
                                        "higher-is-damaging", 10, args.draws, args.seed)
    summary = {
        "scope": "Reusable tool and two-model demonstration using earlier settings",
        "robustness_analyses_run": False,
        "implementation_sha256": sha256(Path(__file__)), "results": results,
    }
    (args.out/"two_model_demo.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    for model, output in results.items():
        value = output["result"]["estimate_raw"]
        text = f"{value:+.6f}" if value is not None else "not estimable"
        print(f"{model}: {text}; {output['result']['n_sites']:,} eligible sites", flush=True)


if __name__ == "__main__":
    main()

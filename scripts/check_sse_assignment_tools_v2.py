from __future__ import annotations

"""Record the available secondary-structure assignment tools for Gate 2."""

import argparse
import json
import shutil
from importlib.util import find_spec
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Check DSSP and STRIDE availability.")
    parser.add_argument("--output", type=Path, default=Path("data/processed/v2/sse_tool_gate_v2.json"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    executables = {name: shutil.which(name) or "" for name in ("mkdssp", "dssp", "stride")}
    result = {
        "executables": executables,
        "biopython_dssp_module": bool(find_spec("Bio.PDB.DSSP")),
        "primary_assignment_available": bool(executables["mkdssp"] or executables["dssp"]),
        "secondary_assignment_available": bool(executables["stride"]),
        "current_pipeline_assignment": "author_annotation",
        "gate_status": "blocked_pending_external_sse_tool" if not any(executables.values()) else "tool_available",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

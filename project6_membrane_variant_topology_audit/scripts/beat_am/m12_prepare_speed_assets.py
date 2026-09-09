"""Discover and pin public inference assets for the speed comparison."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import urllib.request

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/external/speed_assets_20261026"
sys.path.insert(0,str(Path(__file__).resolve().parent))
from m10_fetch_external_baselines import request, dump


def retrieve(url, path):
    if path.exists():
        payload = path.read_bytes()
    else:
        with request(url) as response:
            payload = response.read()
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_bytes(payload)
    return dict(path=str(path.relative_to(ROOT)),url=url,bytes=len(payload),
                sha256=hashlib.sha256(payload).hexdigest())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--weights",action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    manifest = dict(files=[],repositories={},errors=[])
    for label,repo in (("vespa","Rostlab/VESPA"),("vespag","JSchlensok/VespaG"),("alphamissense","google-deepmind/alphamissense")):
        directory = OUT/label
        directory.mkdir(exist_ok=True)
        commit_path = directory/"commit.json"
        if not commit_path.exists():
            with request(f"https://api.github.com/repos/{repo}/commits/main") as response:
                commit_path.write_bytes(response.read())
        commit = json.loads(commit_path.read_text(encoding="utf-8"))["sha"]
        tree_path = directory/"tree.json"
        if not tree_path.exists():
            with request(f"https://api.github.com/repos/{repo}/git/trees/{commit}?recursive=1") as response:
                tree_path.write_bytes(response.read())
        tree = json.loads(tree_path.read_text(encoding="utf-8"))["tree"]
        manifest["repositories"][label] = dict(repository=repo,commit=commit)
        weights = []
        for item in tree:
            if item["type"] != "blob":
                continue
            name = item["path"]
            if name.startswith(("model_weights/","models/")) and name.endswith((".pt",".pth",".pkl",".onnx")):
                weights.append(dict(path=name,bytes=item.get("size")))
            selected = name in ("README.md","params.yaml","hubconf.py","requirements.txt","pyproject.toml")
            selected |= label in ("vespa","vespag") and name.endswith(".py") and not name.startswith("tests/")
            selected |= args.weights and label == "vespa" and name.startswith("models/") and name.endswith((".pt",".pkl"))
            selected |= args.weights and label == "vespag" and name.startswith("model_weights/v1/") and name.endswith((".pt",".pth")) and item.get("size",0) < 20_000_000
            if not selected:
                continue
            url = f"https://raw.githubusercontent.com/{repo}/{commit}/{name}"
            try:
                manifest["files"].append(retrieve(url,directory/name))
            except Exception as exc:
                manifest["errors"].append(dict(url=url,error=str(exc)))
        manifest["repositories"][label]["available_weight_files"] = weights
        print(label,json.dumps(dict(commit=commit,weights=weights)),flush=True)
        dump(OUT/"asset_manifest.json",manifest)
    try:
        path = OUT/"vespag_zenodo_record.json"
        manifest["files"].append(retrieve("https://zenodo.org/api/records/11085958",path))
        record = json.loads(path.read_text(encoding="utf-8"))
        print("ZENODO_FILES",json.dumps([dict(name=x["key"],bytes=x["size"],url=x["links"]["self"])
                                               for x in record.get("files",[])]),flush=True)
    except Exception as exc:
        manifest["errors"].append(dict(source="vespag_zenodo_record",error=str(exc)))
        print("ZENODO_METADATA_UNAVAILABLE",str(exc),flush=True)
    dump(OUT/"asset_manifest.json",manifest)
    print("Downloaded/pinned",len(manifest["files"]),"files; errors",len(manifest["errors"]),flush=True)


if __name__ == "__main__":
    main()

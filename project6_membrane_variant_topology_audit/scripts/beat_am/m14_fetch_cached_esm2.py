"""Fetch the required public ESM2 cache members using HTTP range requests."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(Path(__file__).resolve().parent))
from m10_fetch_external_baselines import HTTPRangeFile, request, dump

OUT = ROOT/"data/external/speed_assets_20261026/cached_esm2"
URL = "https://zenodo.org/records/11085958/files/proteingym_217_esm2.h5?download=1"


class H5RangeFile(HTTPRangeFile):
    def readinto(self, buffer):
        data = self.read(len(buffer))
        buffer[:len(data)] = data
        return len(data)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inspect-only",action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    with np.load(ROOT/"results/beat_am_20261016/features/l1_table.npz",allow_pickle=False) as table:
        membrane = table["is_membrane"] == 1
        requested = sorted(set(table["dms_id"][membrane].astype(str)))
    remote = H5RangeFile(URL)
    manifest = dict(url=URL,source_bytes=remote.length,etag=remote.etag,requested_assays=requested,datasets=[])
    with h5py.File(remote,"r") as source:
        keys = list(source)
        manifest["all_dataset_names"] = keys
        print("H5_DATASETS",len(keys),"requested",len(requested),flush=True)
        if not args.inspect_only:
            for assay in requested:
                choices = [assay,assay+".csv"]
                key = next((k for k in choices if k in source),None)
                if key is None:
                    continue
                path = OUT/f"{assay}.npz"
                if path.exists():
                    with np.load(path,allow_pickle=False) as cached:
                        shape = cached["embedding"].shape
                else:
                    values = np.array(source[key],dtype=np.float32)
                    np.savez_compressed(path,embedding=values,source_key=np.array(key))
                    shape = values.shape
                manifest["datasets"].append(dict(assay=assay,key=key,shape=list(shape),path=str(path.relative_to(ROOT))))
                print("CACHED",assay,shape,flush=True)
    manifest["http_requests"] = remote.requests
    manifest["downloaded_bytes"] = remote.bytes_read
    dump(OUT/"download_manifest.json",manifest)
    print("TRANSFER",remote.bytes_read,"bytes",remote.requests,"requests",flush=True)


if __name__ == "__main__":
    main()

"""CPU timings for cached-score lookup and prediction with cached inputs.

No training, experimental labels, model selection, foundation-model download,
or sequence-encoder timing. Published weights are pinned by m12; VespaG cached
inputs are read by m14. The distinctions between timing stages are explicit.
"""
from __future__ import annotations

import argparse
import ast
from copy import deepcopy
import ctypes
from datetime import datetime, timezone
import gc
import hashlib
import itertools
import json
from pathlib import Path
import pickle
import platform
import sys
import time
import typing
import warnings

import joblib
import numpy as np
import pandas as pd
import sklearn
import torch
from threadpoolctl import threadpool_limits, threadpool_info

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(Path(__file__).resolve().parent))
from m2_standalone_refinement_scorer import FEATURES, build_tabular, load_bundle
from m10_fetch_external_baselines import dump
from m11_compare_external_baselines import CANDIDATE

OUT = ROOT/"results/speed_comparison_20261026"
ASSETS = ROOT/"data/external/speed_assets_20261026"
TABLE = ROOT/"results/beat_am_20261016/features/l1_table.npz"
PREP = ROOT/"results/standalone_multirule_20261025"
BENCH = PREP/"benchmark_20261025"
MODEL = ROOT/"results/beat_am_20261016/speed_benchmark/predictor_free_model.joblib"
EMB = ROOT/"results/beat_am_20261016/_from_instance_gsq37aoc/results/beat_am_20261016/features/embeddings_prot_t5_all"
BEST = "positive_inside+disulfide_loss"
AA = "ACDEFGHIKLMNPQRSTVWY"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def hardware():
    result = dict(platform=platform.platform(),processor=platform.processor(),
                  python=sys.version,logical_cpus=__import__("os").cpu_count(),
                  torch=torch.__version__,numpy=np.__version__,sklearn=sklearn.__version__,
                  torch_cuda_available=torch.cuda.is_available())
    try:
        import winreg
        with winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE,r"HARDWARE\DESCRIPTION\System\CentralProcessor\0") as key:
            result["cpu_name"] = winreg.QueryValueEx(key,"ProcessorNameString")[0].strip()
        class MemoryStatus(ctypes.Structure):
            _fields_ = [("length",ctypes.c_ulong),("load",ctypes.c_ulong)] + [
                (name,ctypes.c_ulonglong) for name in ("total","free","page_total","page_free","virtual_total","virtual_free","extended_free")]
        status = MemoryStatus()
        status.length = ctypes.sizeof(status)
        ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status))
        result["physical_ram_gib"] = status.total/2**30
    except Exception as exc:
        result["hardware_detail_error"] = str(exc)
    return result


def load_trusted_definitions(path, names, namespace):
    """Use only the reviewed published class/function definitions, not its CLI."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    definitions = [node for node in tree.body if isinstance(node,(ast.ClassDef,ast.FunctionDef)) and node.name in names]
    if set(node.name for node in definitions) != set(names):
        raise ValueError("Missing published model definition")
    source = "from __future__ import annotations\n" + ast.unparse(ast.Module(body=definitions,type_ignores=[]))
    exec(compile(source,str(path),"exec"),namespace)
    return namespace


def published_heads():
    nn = torch.nn
    namespace = dict(torch=torch,nn=nn,deepcopy=deepcopy,itertools=itertools,Type=typing.Type)
    load_trusted_definitions(ASSETS/"vespag/vespag/models/utils.py",["construct_fnn"],namespace)
    load_trusted_definitions(ASSETS/"vespag/vespag/models/fnn.py",["FNN"],namespace)
    # The published v1 state uses net.0 and net.2: no dropout module at inference.
    vespag = namespace["FNN"]([256],input_dim=2560,dropout_rate=None)
    weights = torch.load(ASSETS/"vespag/model_weights/v1/esm2.pt",map_location="cpu",weights_only=True)
    vespag.load_state_dict(weights,strict=True)
    vespag.eval()
    load_trusted_definitions(ASSETS/"vespa/vespa/predict/conspred.py",["ConsCNN"],namespace)
    conservation = namespace["ConsCNN"]()
    checkpoint = torch.load(ASSETS/"vespa/models/ProtT5cons_checkpoint.pt",map_location="cpu",weights_only=True)
    conservation.load_state_dict(checkpoint["state_dict"],strict=True)
    conservation.eval()
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        with (ASSETS/"vespa/models/VESPAl-10LR_Cons_Blsm.pkl").open("rb") as handle:
            vespal = pickle.load(handle)
    compatibility = list(dict.fromkeys(str(item.message) for item in captured))
    # Check the old published binary classifier against its coefficient formula.
    real_check = np.zeros((2,10),dtype=float)
    for classifier in vespal.values():
        actual = classifier.predict_proba(real_check)[:,1]
        expected = 1/(1+np.exp(-(real_check@classifier.coef_.T+classifier.intercept_)))[:,0]
        np.testing.assert_allclose(actual,expected,atol=1e-12,rtol=0)
    return vespag,conservation,vespal,compatibility


def prepare():
    with np.load(TABLE,allow_pickle=False) as table:
        # Do not read target / experimental-outcome arrays for timing.
        names = ["accession","position","native","alt","dms_id","is_membrane"]
        data = {name:table[name].copy() for name in names}
        data["tabular"] = np.column_stack([table[name] for name in FEATURES]).astype(np.float32)
    selected = np.flatnonzero(data["is_membrane"] == 1)
    frame = pd.DataFrame({name:data[name][selected] for name in names})
    tabular = data["tabular"][selected]
    with np.load(BENCH/"guarded_loocv/out_of_fold_predictions.npz",allow_pickle=False) as archive:
        if not np.array_equal(archive["accession"],data["accession"]):
            raise ValueError("Cached standalone score order mismatch")
        model_names = archive["models"].astype(str).tolist()
        scores = archive["predictions"][selected].copy()
    with np.load(ROOT/"results/external_predictor_comparison_20261026/aligned_external_scores.npz",allow_pickle=False) as archive:
        if not np.array_equal(archive["accession"],data["accession"]):
            raise ValueError("Cached external score order mismatch")
        external_names = archive["models"].astype(str).tolist()
        external = archive["scores"][selected].copy()
        frame["source_mutant"] = archive["source_mutant"][selected]
    caches = {"Standalone_raw":scores[:,model_names.index("raw")],
              "Standalone_PI_disulfide":scores[:,model_names.index(BEST)],
              "Standalone_PI_glyco_disulfide":scores[:,model_names.index(CANDIDATE)]}
    caches.update({name:external[:,j] for j,name in enumerate(external_names)})
    if not all(np.isfinite(values).all() for values in caches.values()):
        raise ValueError("Nonfinite cached score")
    bank = json.loads((PREP/"annotation_bank.json").read_text(encoding="utf-8"))
    embeddings = {}
    for accession in sorted(set(frame.accession)):
        with np.load(EMB/f"{accession}.npz",allow_pickle=False) as archive:
            pos = archive["position"].astype(int)
            values = archive["embedding"].astype(np.float32)
        record = bank[accession]
        if not np.array_equal(pos,np.arange(1,len(record["sequence"])+1)):
            raise ValueError("Full canonical embeddings are required")
        embeddings[accession] = values
        rows = np.flatnonzero(frame.accession.to_numpy() == accession)
        verified = build_tabular(record,frame.position.to_numpy()[rows],frame.native.to_numpy()[rows],frame.alt.to_numpy()[rows])
        np.testing.assert_array_equal(verified,tabular[rows])
    with np.load(PREP/"prepared_variant_features.npz",allow_pickle=False) as archive:
        feature_names = archive["feature_names"].astype(str).tolist()
        lookup = {(str(a),int(p),str(v)):f for a,p,v,f in zip(archive["accession"],archive["position"],archive["alt"],archive["features"])}
    F = np.array([lookup[str(row.accession),int(row.position),str(row.alt)] for row in frame.itertuples()],dtype=float)
    B = {name:np.zeros_like(F) for name in (BEST,CANDIDATE)}
    for path in sorted((BENCH/"guarded_loocv").glob("fold_*.json")):
        fold = json.loads(path.read_text(encoding="utf-8"))
        mask = np.isin(frame.accession,fold["evaluation_accessions"])
        for name in B:
            coefficients = fold["correction_fits"][name]["coefficients"]
            B[name][mask,:] = [coefficients[rule] for rule in feature_names]
    for name in B:
        delta = np.sum(F*B[name],axis=1)
        np.testing.assert_allclose(scores[:,model_names.index(name)]-scores[:,model_names.index("raw")],delta,atol=1e-10,rtol=0)
    return frame,tabular,caches,bank,embeddings,F,B


def esm2_queries(frame):
    """Align published cached embeddings through the already verified source tokens."""
    directory = ASSETS/"cached_esm2"
    matrix, source_rows, offset, provenance = [],np.empty(len(frame),dtype=int),0,[]
    for assay,part in frame.groupby("dms_id",sort=True):
        path = directory/f"{assay}.npz"
        if not path.exists():
            raise FileNotFoundError(f"Published cached embedding missing: {assay}")
        with np.load(path,allow_pickle=False) as archive:
            values = archive["embedding"].astype(np.float32)
        if values.ndim != 2 or values.shape[1] != 2560 or not np.isfinite(values).all():
            raise ValueError("Wrong cached ESM2 representation")
        positions = np.array([int(token[1:-1])-1 for token in part.source_mutant])
        if np.any(positions<0) or np.any(positions>=len(values)):
            raise ValueError("Embedding/source-mutation coordinate mismatch")
        unique,inverse = np.unique(positions,return_inverse=True)
        matrix.append(values[unique])
        source_rows[part.index] = offset+inverse
        offset += len(unique)
        provenance.append(dict(assay=assay,source_shape=list(values.shape),queried_sites=len(unique),sha256=sha(path)))
    return np.concatenate(matrix),source_rows,provenance


def statistics(samples, rows, inner_loops):
    values = np.asarray(samples)/inner_loops
    median = float(np.median(values))
    return dict(median_seconds=median,p10_seconds=float(np.quantile(values,.1)),
        p90_seconds=float(np.quantile(values,.9)),min_seconds=float(values.min()),
        rows=rows,variants_per_second=rows/median,inner_loops=inner_loops,
        repeats=len(values),samples_seconds=values.tolist())


def time_methods(methods, *, repeats, rows, seed=20261026):
    # Each method is warmed before interleaved, fixed-seed execution order.
    for function,loops in methods.values():
        value = function()
        if np.asarray(value).shape != (rows,) or not np.isfinite(value).all():
            raise ValueError("Timed method returned invalid scores")
    rng = np.random.default_rng(seed)
    samples = {name:[] for name in methods}
    was_enabled = gc.isenabled()
    gc.disable()
    try:
        for _ in range(repeats):
            for name in rng.permutation(list(methods)):
                function,loops = methods[name]
                start = time.perf_counter()
                for _ in range(loops):
                    value = function()
                elapsed = time.perf_counter()-start
                samples[name].append(elapsed)
    finally:
        if was_enabled:
            gc.enable()
    return {name:statistics(values,rows,methods[name][1]) for name,values in samples.items()}


def batch_functions(indices, frame,tabular,caches,bank,embeddings,F,B,bundle,heads,esm2):
    part = frame.iloc[indices]
    grouped = [(str(a),np.array(rows,dtype=int)) for a,rows in part.groupby("accession",sort=True).indices.items()]
    # part indices are positions in this requested batch, not global row indices.
    positions = part.position.to_numpy().astype(int)
    natives,alts = part.native.to_numpy(),part.alt.to_numpy()
    projection = {}
    for accession,local in grouped:
        unique,inverse = np.unique(positions[local]-1,return_inverse=True)
        projection[accession] = (local,embeddings[accession][unique],inverse)
    def make_x(rebuild=False):
        tab = tabular[indices].copy()
        pcs = np.empty((len(indices),64),dtype=np.float32)
        for accession,(local,E,inverse) in projection.items():
            pcs[local] = bundle["pca"].transform(E)[inverse]
            if rebuild:
                tab[local] = build_tabular(bank[accession],positions[local],natives[local],alts[local])
        return np.column_stack([tab,pcs]).astype(np.float32)
    X = make_x()
    baseline = bundle["model"].predict(X)
    coefficients = {name:B[name][indices] for name in B}
    features = F[indices]
    methods = {}
    for model,values in caches.items():
        # Same indexed in-memory lookup implementation for every predictor.
        methods["lookup/"+model] = (lambda values=values:np.take(values,indices),1000 if len(indices)==1 else 10)
    methods["head_only/Standalone_raw"] = (lambda:bundle["model"].predict(X),1)
    methods["cached_embedding/Standalone_raw"] = (lambda:bundle["model"].predict(make_x()),1)
    methods["cached_embedding_rebuild_features/Standalone_raw"] = (lambda:bundle["model"].predict(make_x(True)),1)
    for name,label in ((BEST,"Standalone_PI_disulfide"),(CANDIDATE,"Standalone_PI_glyco_disulfide")):
        beta = coefficients[name]
        methods["cached_embedding/"+label] = (lambda beta=beta:bundle["model"].predict(make_x())+np.sum(features*beta,axis=1),1)
        methods["correction_only/"+label] = (lambda beta=beta:baseline+np.sum(features*beta,axis=1),100 if len(indices)==1 else 1)
    vespag,conservation,vespal = heads
    if conservation is not None:
        tensors = {a:torch.from_numpy(embeddings[a]).unsqueeze(0) for a,_ in grouped}
        def light_predict():
            probabilities = np.empty((len(indices),9),dtype=np.float32)
            for accession,local in grouped:
                probability = conservation.extract_probabilities(conservation(tensors[accession]))[0].T.numpy()
                probabilities[local] = probability[positions[local]-1]
            # Nine conservation probabilities and the real variant's BLOSUM62.
            inputs = np.column_stack([probabilities,tabular[indices,FEATURES.index("blosum62")]])
            return np.mean([model.predict_proba(inputs)[:,1] for model in vespal.values()],axis=0)
        methods["cached_embedding/VESPAl"] = (light_predict,1)
    if vespag is not None and esm2 is not None:
        E,source_rows = esm2
        unique,inverse = np.unique(source_rows[indices],return_inverse=True)
        tensor = torch.from_numpy(E[unique])
        alternate = np.array([AA.index(str(a)) for a in alts])
        def gemme_student_predict():
            landscape = vespag(tensor).numpy()
            return landscape[inverse,alternate]
        methods["cached_embedding/VespaG"] = (gemme_student_predict,1)
    return methods


def write_report(result):
    lines = ["# CPU speed benchmark with cached inputs", "", "Date: 26 October 2026 (America/Los_Angeles).", "",
        "Scope selected by the user: prediction with cached inputs, plus lookup for every model. No foundation-model encoder or MSA-search runtime is included. No new training or coefficient fitting was performed.", "",
        "The main workload is the same 43,991 membrane-panel measurement rows used in the accuracy comparison (13 proteins, 15 assays). Duplicated variant requests across assays remain duplicated for every method. Smaller batches use the same fixed-seed requests.", "",
        "Model initialization, disk loading, imports and tensor construction are excluded from warm timings. Single-query lookups are repeated within each timed sample to reduce clock noise. Methods are warmed and interleaved in random fixed-seed order; results are medians of seven measurements with 10th/90th percentiles retained in JSON.", "",
        "## Hardware", "", "```json",json.dumps(result["hardware"],indent=2),"```", "",
        "## Prediction after cached embeddings", "",
        "This includes each measured prediction head and its output selection. For the standalone model, PCA and supplied biochemical/topology features are included; rebuilding those features is also timed separately. VESPA-light includes its published conservation CNN plus its ten logistic-regression models. VespaG uses its published v1 ESM2 head. Cached embeddings have different dimensions and were produced by each method's own encoder.", "",
        "| Threads | Stage / model | Batch size | Median ms | Variants / second |", "|---:|---|---:|---:|---:|"]
    for row in result["measurements"]:
        if row["rows"] == result["primary_rows"] and not row["method"].startswith("lookup/"):
            lines.append(f"| {row['threads']} | {row['method']} | {row['rows']:,} | {row['median_seconds']*1000:.3f} | {row['variants_per_second']:,.0f} |")
    lines += ["", "## Lookup for every model", "",
        "All models use the same in-memory indexed-array retrieval implementation on the same requests. These timings measure cached lookup, not original model inference; differences here are memory/access noise and do not establish a model speed advantage.", "",
        "| Threads | Model | Requests | Median ms |", "|---:|---|---:|---:|"]
    for row in result["measurements"]:
        if row["rows"] == result["primary_rows"] and row["method"].startswith("lookup/"):
            lines.append(f"| {row['threads']} | {row['method'].split('/',1)[1]} | {row['rows']:,} | {row['median_seconds']*1000:.3f} |")
    lines += ["", "## Stages that could not be measured", ""]
    for model,reason in result["unavailable_stages"].items():
        lines.append(f"- **{model}:** {reason}")
    lines += ["", "## Interpretation", "",
        "A fast prediction head does not establish a fast complete predictor. ProtT5 generation, ESM2 generation, ESM1b masked-sequence inference, MSA search, annotation acquisition, and loading an entire public proteome archive are outside this selected benchmark. Lookup timings must not be compared to another method's sequence-to-score timing.", "",
        "The saved standalone head is the existing full-development 18+64-feature bundle, with the same fixed architecture used for the accuracy comparison. Runtime is not measured by refitting the head or mixing in AlphaMissense features. The additional layer's measured coefficients and features come from the existing guarded folds; their effects reproduce the saved refined-minus-raw predictions.", "",
        "Published model classes are loaded from pinned upstream source definitions with strict state-dictionary matching. No synthetic embeddings or random model weights are timed as model predictions. The LR compatibility check is recorded in verification. Numeric outputs from these runtime adapters are not used to replace the published scores in the accuracy comparison.", "",
        "## Reproduce", "", "```powershell", ".\\.venv_ml\\Scripts\\python.exe scripts\\beat_am\\m13_local_speed_benchmark.py", "```", "",
        "Sources: [AlphaMissense availability](https://github.com/google-deepmind/alphamissense), [VESPA](https://github.com/Rostlab/VESPA), [VespaG](https://github.com/JSchlensok/VespaG), [published cached ESM2 inputs](https://zenodo.org/records/11085958).", ""]
    text = "\n".join(lines)
    (OUT/"REPORT.md").write_text(text,encoding="utf-8")
    (ROOT/"docs/STANDALONE_SPEED_COMPARISON_20261026.md").write_text(text,encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--threads",type=int,nargs="+",default=[1,8])
    parser.add_argument("--repeats",type=int,default=7)
    args = parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    torch.set_num_threads(1)
    frame,tabular,caches,bank,embeddings,F,B = prepare()
    bundle = load_bundle(MODEL)
    vespag,conservation,vespal,compatibility = published_heads()
    unavailable = {
        "AlphaMissense fresh inference":"Official trained weights are not released; lookup was timed.",
        "ESM1b fresh inference":"Masked language-model forward passes are outside the selected cached-input scope; lookup was timed.",
        "GEMME cached-alignment prediction":"No local R/Python2/Docker GEMME runtime; the official download host was unreachable. Lookup was timed; no inferred or literature runtime is substituted.",
        "VESPA full cached-input prediction":"No cached ProtT5 mutation log-odds inputs are available. Its complete pipeline is not replaced by VESPA-light or by synthetic feature values; lookup was timed."}
    try:
        E,source_rows,esm2_provenance = esm2_queries(frame)
        esm2 = E,source_rows
    except Exception as exc:
        esm2,esm2_provenance = None,[]
        unavailable["VespaG cached-input prediction"] = str(exc)
    result = dict(scope="user-selected cached inputs and all-model lookup",hardware=hardware(),
        primary_rows=len(frame),primary_proteins=int(frame.accession.nunique()),
        repeats=args.repeats,threads=args.threads,warmup=1,
        no_fit=True,no_encoder_inference=True,no_experimental_outcome_reads=True,
        source_hashes={"table":sha(TABLE),"saved_standalone_head":sha(MODEL),"source_code":sha(__file__)},
        lr_compatibility_warnings=compatibility,esm2_inputs=esm2_provenance,
        unavailable_stages=unavailable,measurements=[])
    dump(OUT/"PLAN.json",{key:value for key,value in result.items() if key != "measurements"})
    rng = np.random.default_rng(20261026)
    requests = {1:np.array([0]),1000:np.sort(rng.choice(len(frame),1000,replace=False)),len(frame):np.arange(len(frame))}
    for threads in args.threads:
        torch.set_num_threads(threads)
        with threadpool_limits(limits=threads),torch.inference_mode():
            for count,indices in requests.items():
                methods = batch_functions(indices,frame,tabular,caches,bank,embeddings,F,B,bundle,
                                          (vespag,conservation,vespal),esm2)
                measured = time_methods(methods,repeats=args.repeats,rows=count)
                for name,row in measured.items():
                    result["measurements"].append(dict(method=name,threads=threads,**row))
                dump(OUT/"results.json",result)
                print(f"Completed threads={threads}, rows={count}, methods={len(measured)}",flush=True)
    pd.DataFrame(result["measurements"]).drop(columns=["samples_seconds"]).to_csv(OUT/"timings.csv",index=False)
    dump(OUT/"verification.json",dict(feature_parity=True,saved_refinement_effects_reproduced=True,
        native_published_head_state_dicts_strict=True,lr_probability_formula_checked=True,
        synthetic_embedding_benchmarks=False,experimental_targets_read=False,
        model_hash_unchanged=sha(MODEL)==result["source_hashes"]["saved_standalone_head"],
        thread_libraries=threadpool_info(),asset_manifest_sha256=sha(ASSETS/"asset_manifest.json")))
    write_report(result)
    print("PRIMARY WARM TIMINGS",flush=True)
    for row in result["measurements"]:
        if row["rows"] == len(frame) and row["threads"] == 1:
            print(row["method"],round(row["median_seconds"]*1000,3),"ms",flush=True)
    print("Report:",ROOT/"docs/STANDALONE_SPEED_COMPARISON_20261026.md",flush=True)


if __name__ == "__main__":
    main()

"""Benchmark cached-feature prediction and correction for the standalone scorer."""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "scripts" / "beat_am"))

FEAT = ROOT / "results" / "beat_am_20261016" / "features"
EMB = (ROOT / "results" / "beat_am_20261016" / "_from_instance_gsq37aoc"
       / "results" / "beat_am_20261016" / "features" / "embeddings_prot_t5_all")
OUT = ROOT / "results" / "beat_am_20261016" / "speed_benchmark"

from b1_predictor_free import PREDICTOR_FREE, KEEP_IDX, load_dev  # noqa: E402
from a17_general_scorer import _backend, GENERAL  # noqa: E402

N_PC = 64
SEED = 20261016

# The 4.7-second ProtT5 encoding time is the recorded 39-protein run on two L40 GPUs.
PRIOR_EMBEDDING = {"hardware": "2x NVIDIA L40",
                   "proteins": 39, "seconds": 4.7,
                   "source": "beat_am run log, route B embedding stage"}


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def timeit(fn, repeats=5):
    fn()
    best = float("inf")
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    from sklearn.decomposition import PCA
    import joblib

    make_model, backend = _backend()
    log(f"backend: {backend}")

    log("loading development features")
    data = load_dev()
    X, E, y = data["X"], data["E"], data["y"]
    n = len(y)
    log(f"  {n:,} variants, {X.shape[1]} tabular + {N_PC} PCs "
        f"= {X.shape[1] + N_PC} features")

    log("fitting the predictor-free scorer (and persisting it this time)")
    t0 = time.perf_counter()
    pca = PCA(n_components=N_PC, random_state=SEED).fit(E)
    Xf = np.hstack([X, pca.transform(E)]).astype(np.float32)
    model = make_model({"max_depth": 6, "learning_rate": 0.05,
                        "n_estimators": 500}, SEED)
    model.fit(Xf, y)
    fit_s = time.perf_counter() - t0
    log(f"  fitted in {fit_s:.1f}s")
    joblib.dump({"model": model, "pca": pca, "features": PREDICTOR_FREE,
                 "n_pc": N_PC, "backend": backend},
                OUT / "predictor_free_model.joblib")
    size_mb = (OUT / "predictor_free_model.joblib").stat().st_size / 1048576
    log(f"  artifact {size_mb:.1f} MB")

    results = {}


    for batch in (1, 1_000, 100_000):
        Xb = Xf[:batch] if batch <= n else Xf
        secs = timeit(lambda: model.predict(Xb))
        rate = len(Xb) / secs
        results[f"tree_inference_batch_{batch}"] = {
            "rows": int(len(Xb)), "seconds": secs,
            "variants_per_second": rate}
        log(f"  tree inference, batch {len(Xb):>7,}: {secs * 1e3:9.3f} ms "
            f"-> {rate:>12,.0f} variants/s")


    Eb = E[:100_000]
    secs = timeit(lambda: pca.transform(Eb))
    results["pca_1024_to_64"] = {"rows": int(len(Eb)), "seconds": secs,
                                 "rows_per_second": len(Eb) / secs}
    log(f"  PCA 1024->64, {len(Eb):,} rows: {secs * 1e3:.1f} ms "
        f"-> {len(Eb) / secs:,.0f} rows/s")


    t = np.load(FEAT / "l1_table.npz", allow_pickle=True)
    d = {k: t[k] for k in t.files}

    def build_tabular():
        return np.column_stack(
            [d[f].astype(np.float32) for f in GENERAL])[:, KEEP_IDX]

    secs = timeit(build_tabular, repeats=3)
    results["tabular_features"] = {"rows": n, "seconds": secs,
                                   "rows_per_second": n / secs}
    log(f"  tabular features, {n:,} rows: {secs * 1e3:.1f} ms "
        f"-> {n / secs:,.0f} rows/s")


    sizes = {}
    for label, path in (
            ("predictor_free_model", OUT / "predictor_free_model.joblib"),
            ("am_containing_model", ROOT / "results" / "beat_am_20261016"
             / "final_model" / "model.joblib"),
            ("alphamissense_released_scores", ROOT / "data" / "raw"
             / "AlphaMissense_aa_substitutions.tsv.gz"),
            ("esm1b_released_scores", ROOT / "data" / "raw"
             / "ESM1b_ALL_hum_isoforms_LLR.zip"),
            ("prot_t5_site_embeddings", ROOT / "results" / "route_b_20261012"
             / "site_embeddings_fp16.npy")):
        if path.exists():
            sizes[label] = round(path.stat().st_size / 1048576, 1)
    results["artifact_sizes_mb"] = sizes
    log("  artifact sizes (MB): " + json.dumps(sizes))

    summary = {
        "question": "is the DMS-trained, predictor-free scorer faster than "
                    "AlphaMissense and ESM1b?",
        "honest_answer": "not in any way a user would feel, and the framing is "
                         "misleading. See accounting below.",
        "accounting": {
            "the_scorer_is_not_small": "it is a gradient-boosted tree head on "
                                       "ProtT5-XL, a 3B-parameter transformer. "
                                       "The trees are ~MB; the frontend is not.",
            "vs_esm1b_end_to_end": "ProtT5-XL at 3B parameters is LARGER than "
                                   "ESM1b at 650M, so on a novel sequence this "
                                   "scorer needs more compute, not less",
            "vs_alphamissense_in_practice": "AlphaMissense and ESM1b ship as "
                                            "precomputed proteome-wide score "
                                            "files. For any human protein a "
                                            "user's cost is a lookup, not a "
                                            "forward pass, so there is no "
                                            "inference time to beat",
            "where_speed_could_matter": "sequences outside those released files "
                                        "-- engineered proteins, other organisms "
                                        "-- and there this scorer is the heavier "
                                        "option, not the lighter one",
        },
        "measured_here": results,
        "not_measured_here": {
            "prot_t5_forward_pass": "no GPU in this session; the project's own "
                                    "prior measurement is cited instead",
            "prior_measurement": PRIOR_EMBEDDING,
            "alphamissense_inference": "weights are not runnable without "
                                       "AlphaFold-scale infrastructure; only "
                                       "released predictions were used anywhere "
                                       "in this project",
        },
        "fit_seconds": fit_s,
        "backend": backend,
    }
    (OUT / "speed_benchmark.json").write_text(
        json.dumps(summary, indent=2, default=str), encoding="utf-8")
    log("")
    log(f"wrote {OUT / 'speed_benchmark.json'}")


if __name__ == "__main__":
    main()

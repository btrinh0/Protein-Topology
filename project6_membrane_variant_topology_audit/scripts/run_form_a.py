"""Train Form A student models on planted rules and audit held-out predictions.

The pilot stage uses separate seeds to verify the pipeline and model capacity.
Synthetic targets use observed membrane geometry without
reading AlphaMissense scores. `summarize_form_a.py` reads the outputs.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import form_ab_lib as lib
from run_gate1 import read_topology

FIELDS = ["scenario", "student", "replicate", "oracle_raw", "independent_raw", "estimate_raw", "ci_lo", "ci_hi",
          "penalty_sd", "n_sites", "n_strata", "n_proteins", "train_seed", "noise_seed", "best_iter",
          "val_mse_at_best", "input_sha256", "seconds"]
NAIVE_FIELDS = ["scenario", "student", "replicate", "oracle_raw", "matched", "charged_arm_only", "pooled_paired",
                "sd_z0"]
DIAG_FIELDS = ["scenario", "student", "replicate", "oracle_validation", "estimate_validation", "ratio"]
TABLES = (("estimates", FIELDS, "estimates.csv"), ("naive", NAIVE_FIELDS, "naive_comparators.csv"),
          ("diag", DIAG_FIELDS, "validation_capacity.csv"))
STUDENT_ORDER = ("S0", "S0n", "S1", "S2", "S3")
STAGES = {
    "pilot": dict(seed=lib.PILOT_SEED, replicates=3, train=40_000, validation=20_000,
                  out=lib.ROOT / "results/form_ab_pilot_20261008/form_a"),
    # The full-stage cap exceeds the available training positions, so all are used.
    "full": dict(seed=lib.MASTER_SEED, replicates=20, train=1_000_000, validation=30_000,
                 out=lib.ROOT / "results/form_a_20261008"),
}


def prepare(args) -> dict:
    census = read_topology(lib.TOPOLOGY)
    sites, accessions = lib.load_sites()
    split = lib.cluster_split(census)
    pools = lib.site_pools(sites, accessions, split["label"])
    window = sites["distance"] <= 10
    ctx = {
        "census": census, "sites": sites, "accessions": accessions, "split": split, "pools": pools,
        "base": lib.blosum_base(), "context": lib.context_counts(sites, accessions, census),
        "test": pools["test"][window[pools["test"]]],
        "validation_window": pools["validation"][window[pools["validation"]]],
    }
    ctx["columns"] = lib.site_columns(sites, ctx["context"])
    ctx["sd_z0"] = lib.sd_z0(lib.take(sites, ctx["test"]), ctx["base"])
    return ctx


def job_streams(args, scenario: str, replicate: int) -> dict:
    return {name: lib.rng_for(args.seed, scenario, replicate, name) for name in lib.STREAMS}


def training_data(ctx, args, spec, amp, scenario, replicate):
    streams = job_streams(args, scenario, replicate)
    out = {}
    for name, pool, size in (("train", ctx["pools"]["train"], args.train),
                             ("validation", ctx["pools"]["validation"], args.validation)):
        index = streams["subsample"].choice(pool, size=min(size, len(pool)), replace=False)
        sub = lib.take(ctx["sites"], index)
        parts = lib.draw_target(sub, ctx["base"], amp, spec["lam"], spec["side_main"], streams[name])
        full = parts["clean"] + parts["site"] + parts["noise"]
        out[name] = (index, lib.row_targets(full, sub["ref_aa"]))
    return out, streams


def audit_official(ctx, args, scenario, student, replicate, test_sites, matrix, oracle, extra) -> dict:
    """Write the mapped archive and run the packaged estimator on it, as the command line does."""
    out = args.out / scenario / student / f"{replicate:02d}"
    npz = args.out / "predictions" / f"{scenario}_{student}_{replicate:02d}.npz"
    lib.write_scores_npz(npz, test_sites, matrix, ctx["accessions"])
    report = lib.P.run_panel_audit(npz, out, f"{scenario}_{student}_{replicate:02d}", "higher-is-damaging",
                                   10, lib.PANEL_DRAWS, lib.PANEL_SEED)
    result = report["result"]
    if not args.keep_predictions:
        npz.unlink()
    return {"scenario": scenario, "student": student, "replicate": replicate, "oracle_raw": oracle,
            "independent_raw": lib.independent_delta(test_sites, matrix),
            "estimate_raw": result["estimate_raw"], "ci_lo": result["ci95_raw"][0], "ci_hi": result["ci95_raw"][1],
            "penalty_sd": result["penalty_sd"], "n_sites": result["n_sites"], "n_strata": result["n_strata"],
            "n_proteins": result["n_proteins"], "input_sha256": report["input_sha256"], **extra}


def run_job(ctx, args, scenario: str, replicate: int) -> dict:
    spec = lib.SCENARIOS[scenario]
    amp = spec["rho"] * ctx["sd_z0"]
    started = time.time()
    test_sites = lib.take(ctx["sites"], ctx["test"])
    streams = job_streams(args, scenario, replicate)
    target = lib.draw_target(test_sites, ctx["base"], amp, spec["lam"], spec["side_main"], streams["test"])
    oracle, _ = lib.oracle_delta(test_sites, target["h"])
    if spec["lam"] is None and abs(oracle - amp) > 1e-12:
        raise RuntimeError("Flat planted rule must give an oracle equal to its amplitude")
    noise_seed = f"{args.seed}.{list(lib.SCENARIOS).index(scenario)}.{replicate}.{lib.STREAMS['test']}"
    payload = {"estimates": [], "naive": [], "diag": []}

    for student, matrix in (("S0", lib.mask_native(target["clean"] + target["site"], test_sites["ref_aa"])),
                            ("S0n", lib.mask_native(target["clean"] + target["site"] + target["noise"],
                                                    test_sites["ref_aa"]))):
        payload["estimates"].append(audit_official(
            ctx, args, scenario, student, replicate, test_sites, matrix, oracle,
            {"train_seed": "", "noise_seed": noise_seed, "best_iter": "", "val_mse_at_best": "", "seconds": ""}))

    data, _ = training_data(ctx, args, spec, amp, scenario, replicate)
    (train_index, y_train), (val_index, y_val) = data["train"], data["validation"]
    val_window = lib.take(ctx["sites"], ctx["validation_window"])
    oracle_val, _ = lib.oracle_delta(val_window, lib.planted_h(val_window, amp, spec["lam"]))
    model_seed = int(streams["model"].integers(0, 2 ** 31 - 1))
    for student, names in lib.STUDENT_FEATURES.items():
        fit_started = time.time()
        X_train = lib.design_matrix(ctx["columns"], train_index, names)
        X_val = lib.design_matrix(ctx["columns"], val_index, names)
        model, best, curve = lib.fit_student(X_train, y_train, X_val, y_val, model_seed)
        del X_train, X_val
        X_test = lib.design_matrix(ctx["columns"], ctx["test"], names)
        matrix = lib.predictions_to_matrix(lib.predict_at(model, X_test, best), test_sites["ref_aa"])
        row = audit_official(ctx, args, scenario, student, replicate, test_sites, matrix, oracle,
                             {"train_seed": model_seed, "noise_seed": noise_seed, "best_iter": best,
                              "val_mse_at_best": float(curve[best - 1]),
                              "seconds": round(time.time() - fit_started, 1)})
        payload["estimates"].append(row)
        if student in ("S1", "S2") and scenario in ("N0", "S_main"):
            payload["naive"].append(
                {"scenario": scenario, "student": student, "replicate": replicate, "oracle_raw": oracle,
                 "matched": row["estimate_raw"], **lib.naive_comparators(test_sites, matrix),
                 "sd_z0": ctx["sd_z0"]})
        if args.diagnose_validation:
            X_diag = lib.design_matrix(ctx["columns"], ctx["validation_window"], names)
            diag = lib.predictions_to_matrix(lib.predict_at(model, X_diag, best), val_window["ref_aa"])
            estimate = lib.audit_matrix(val_window, diag, draws=2)[0]["estimate_raw"]
            payload["diag"].append(
                {"scenario": scenario, "student": student, "replicate": replicate, "oracle_validation": oracle_val,
                 "estimate_validation": estimate, "ratio": estimate / oracle_val if oracle_val else ""})
    print(f"{scenario} rep {replicate}: {time.time() - started:.0f}s", flush=True)
    return payload


def job_files(args) -> list:
    folder = args.out / "jobs"
    return sorted(folder.glob("*.json")) if folder.exists() else []


def save_job(args, scenario: str, replicate: int, payload: dict) -> None:
    """Atomic: a job exists only when all of its rows were written."""
    folder = args.out / "jobs"
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f"{scenario}_{replicate:02d}.json"
    temp = path.with_suffix(".tmp")
    temp.write_text(json.dumps(payload), encoding="utf-8")
    os.replace(temp, path)


def completed(args) -> set:
    return {(p.stem.rsplit("_", 1)[0], int(p.stem.rsplit("_", 1)[1])) for p in job_files(args)}


def write_tables(args) -> None:
    """Rebuild the CSV tables from the job files, so a restart can never duplicate or drop rows."""
    merged = {name: [] for name, _, _ in TABLES}
    for path in job_files(args):
        payload = json.loads(path.read_text(encoding="utf-8"))
        for name in merged:
            merged[name].extend(payload[name])
    order = list(lib.SCENARIOS)
    for name in merged:
        merged[name].sort(key=lambda r: (r["replicate"], order.index(r["scenario"]),
                                         STUDENT_ORDER.index(r["student"])))
    for name, fields, filename in TABLES:
        if merged[name]:
            with (args.out / filename).open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=fields)
                writer.writeheader()
                writer.writerows(merged[name])


def run(args) -> None:
    args.out.mkdir(parents=True, exist_ok=True)
    if job_files(args) and not args.resume:
        raise SystemExit("Jobs already exist in the output folder; pass --resume to continue them")
    ctx = prepare(args)
    split = ctx["split"]
    test_sites = lib.take(ctx["sites"], ctx["test"])
    _, weight = lib.oracle_delta(test_sites, np.zeros(len(ctx["test"])))
    manifest = {
        "stage": args.stage, "as_of_utc": datetime.now(timezone.utc).isoformat(), "seed": args.seed,
        "replicates": args.replicates, "train_positions_cap": args.train, "validation_positions": args.validation,
        "scenarios": args.scenarios, "split": {k: split[k] for k in ("kind", "seed", "counts", "assignment_sha256")},
        "train_sites_available": int(len(ctx["pools"]["train"])),
        "test_sites_in_window": int(len(ctx["test"])), "sum_of_weights": weight, "sd_z0": ctx["sd_z0"],
        "equivalence_margin": lib.EQUIVALENCE * ctx["sd_z0"], "panel_draws": lib.PANEL_DRAWS,
        "panel_seed": lib.PANEL_SEED, "tool": lib.validated_tool(), "packages": lib.package_versions(),
        "archive_sha256": lib.P.sha256(lib.AM_ARCHIVE), "topology_sha256": lib.P.sha256(lib.TOPOLOGY),
        "implementation_sha256": {"form_ab_lib": lib.P.sha256(Path(lib.__file__)),
                                  "run_form_a": lib.P.sha256(Path(__file__))},
        "score_array_read": False, "clinical_labels_used": False,
    }
    manifest_path = args.out / "run_manifest.json"
    if manifest_path.exists():
        manifest_path = args.out / f"run_manifest_resume_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (args.out / "split.json").write_text(json.dumps(
        {k: split[k] for k in ("kind", "seed", "counts", "assignment_sha256")}, indent=2) + "\n", encoding="utf-8")
    done = completed(args)
    for replicate in range(args.replicates):
        for scenario in args.scenarios:
            if (scenario, replicate) in done:
                continue
            save_job(args, scenario, replicate, run_job(ctx, args, scenario, replicate))
            write_tables(args)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("pilot", "full"), required=True)
    parser.add_argument("--resume", action="store_true", help="skip jobs that already finished")
    parser.add_argument("--scenarios", nargs="+", choices=list(lib.SCENARIOS), default=list(lib.SCENARIOS))
    parser.add_argument("--replicates", type=int)
    parser.add_argument("--train", type=int, help="cap on training positions per replicate (19 alternatives each)")
    parser.add_argument("--validation", type=int)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--diagnose-validation", action="store_true")
    parser.add_argument("--keep-predictions", action="store_true",
                        help="keep each prediction archive; by default it is deleted after its audit, "
                             "leaving the audit's input_sha256")
    args = parser.parse_args()
    defaults = STAGES[args.stage]
    for name in ("replicates", "train", "validation", "out"):
        if getattr(args, name) is None:
            setattr(args, name, defaults[name])
    args.seed = defaults["seed"]
    run(args)


if __name__ == "__main__":
    main()

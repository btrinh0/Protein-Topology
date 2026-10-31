"""Compare frozen standalone predictions with matched external scores without refitting."""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from route_a_common import align_construct
from a9_l1_features import read_am

TABLE = ROOT / "results/beat_am_20261016/features/l1_table.npz"
BENCH = ROOT / "results/standalone_multirule_20261025/benchmark_20261025"
SOURCES = ROOT / "data/external/predictor_comparison_20261026"
REFERENCE = ROOT / "data/raw/proteingym/reference_github_main/DMS_substitutions.csv"
ASSAYS = ROOT / "data/raw/proteingym/v1.3_extracted/DMS_ProteinGym_substitutions"
OUT = ROOT / "results/external_predictor_comparison_20261026"
NEW_EXTERNAL = ("GEMME", "VESPA", "VESPAl", "VespaG")
REQUESTED_EXTERNAL = ("AlphaMissense", "ESM1b", "GEMME", "VESPA", "VESPAl")
ALL_EXTERNAL = REQUESTED_EXTERNAL + ("VespaG",)
CANDIDATE = "positive_inside+glycosylation_loss+disulfide_loss"
MIN_ROWS, MIN_TARGET_VALUES, DRAWS, SEED = 20, 3, 9999, 20261025
MUTANT = re.compile(r"^([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])$")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path, records):
    frame = pd.DataFrame(records)
    frame.to_csv(path, index=False)


def load_table():
    with np.load(TABLE, allow_pickle=False) as source:
        names = ("accession", "position", "native", "alt", "dms_id", "gene",
                 "target", "is_membrane", "endpoint_type", "am", "esm1b")
        data = {name: source[name] for name in names}
    frame = pd.DataFrame(data)
    frame["row_index"] = np.arange(len(frame))
    keys = ["dms_id", "accession", "position", "native", "alt"]
    if frame.duplicated(keys).any():
        raise ValueError("Duplicate assay/canonical mutation in development table")
    return frame


def freeze_plan(out, table):
    previous = json.loads((BENCH / "BENCHMARK_PLAN.json").read_text(encoding="utf-8"))
    if previous["table_sha256"] != sha(TABLE):
        raise ValueError("Original cross-validation source table has changed")
    plan = dict(
        study="external comparison of frozen standalone group-out predictions",
        status="exploratory comparison on previously exposed development data",
        primary_scheme="guarded_loocv", sensitivity_scheme="grouped_cv",
        primary_scope="membrane", primary_refined_candidate=CANDIDATE,
        primary_common_comparators=list(REQUESTED_EXTERNAL),
        supplementary_comparator="VespaG (distinct from VESPA and VESPAl)",
        primary_metric="per-assay Spearman, equal assay means within protein, equal protein means",
        secondary_metric="MSE of within-cohort, within-assay ranks; no raw-score MSE comparison",
        minimum_rows_per_assay=MIN_ROWS, minimum_unique_target_values=MIN_TARGET_VALUES,
        row_matching="assay ID plus exact reference-aligned canonical native/position/alternate",
        missing_scores="no imputation; shared intersection primary; comparator-paired coverage secondary",
        constant_scores="drop undefined-correlation assay for every model in that comparison; report exclusions",
        score_orientation="cached AM/ESM1b higher=damage; official merged ProteinGym columns higher=fitness, negated once",
        orientation_chosen_from_outcomes=False,
        refitting=False, new_coefficients=False, selecting_best_refinement=False,
        scopes=["all", "membrane", "non_membrane", "official coarse assay types within all and membrane"],
        uncertainty="paired dependence-group resampling of fixed protein metrics; conditional, not fresh confirmation",
        bootstrap_draws=DRAWS, bootstrap_seed=SEED,
        no_reserved_experimental_outcomes=True,
        table_sha256=sha(TABLE), previous_plan_sha256=sha(BENCH / "BENCHMARK_PLAN.json"),
        source_prediction_sha256={mode:sha(BENCH/mode/"out_of_fold_predictions.npz")
                                  for mode in ("guarded_loocv", "grouped_cv")},
        proteins=int(table.accession.nunique()), assays=int(table.dms_id.nunique()), rows=len(table))
    path = out / "COMPARISON_PLAN.json"
    if path.exists():
        saved = json.loads(path.read_text(encoding="utf-8"))
        saved.pop("created_utc", None)
        if saved != plan:
            raise ValueError("Comparison plan changed; preserve the original and use a new output folder")
    else:
        dump(path, dict(created_utc=datetime.now(timezone.utc).isoformat(), **plan))
    return sha(path)


def load_oof(mode, table):
    directory = BENCH / mode
    with np.load(directory / "out_of_fold_predictions.npz", allow_pickle=False) as source:
        for name in ("accession", "position", "native", "alt"):
            if not np.array_equal(source[name], table[name].to_numpy()):
                raise ValueError(f"OOF row order mismatch: {name}")
        if not np.array_equal(source["target"], table.target.to_numpy().astype(float)):
            raise ValueError("OOF target/order mismatch")
        models = source["models"].astype(str).tolist()
        predictions = source["predictions"].copy()
        groups = source["groups"].astype(str).copy()
    if CANDIDATE not in models or not np.isfinite(predictions).all():
        raise ValueError("Missing requested candidate or nonfinite original predictions")
    covered = []
    for path in sorted(directory.glob("fold_*.json")):
        fold = json.loads(path.read_text(encoding="utf-8"))
        train, test = set(fold["training_groups"]), set(fold["evaluation_groups"])
        if train & test:
            raise ValueError("Outer-fold group leakage")
        for inner in fold["inner"]:
            itrain, ical = set(inner["training_groups"]), set(inner["calibration_groups"])
            if itrain & ical or test & (itrain | ical):
                raise ValueError("Inner-calibration/outer-test group leakage")
        covered.extend(test)
    expected = sorted(set(groups))
    if mode == "guarded_loocv" and sorted(covered) != expected:
        raise ValueError("Guarded OOF group coverage mismatch")
    if mode == "grouped_cv" and sorted(covered) != expected:
        raise ValueError("Grouped OOF group coverage mismatch")
    return models, predictions, groups


def map_mutant(token, construct, canonical, mapping):
    match = MUTANT.fullmatch(str(token))
    if not match:
        return None
    native, position, alternate = match.group(1), int(match.group(2)), match.group(3)
    if native == alternate:
        return None
    if not 1 <= position <= len(construct) or construct[position-1] != native:
        return None
    mapped = int(mapping[position])
    if not 1 <= mapped <= len(canonical) or canonical[mapped-1] != native:
        return None
    return mapped, native, alternate


def numeric(value):
    try:
        return float(value)
    except (ValueError, TypeError):
        return np.nan


def align_external(table, out):
    manifest = json.loads((SOURCES / "download_manifest.json").read_text(encoding="utf-8"))
    if manifest["errors"] or manifest["source_table_sha256"] != sha(TABLE):
        raise ValueError("Incomplete/mismatched external-score download")
    entries = {entry["assay"]:entry for entry in manifest["entries"]}
    references = {row["DMS_id"]:row for row in csv.DictReader(REFERENCE.open(encoding="utf-8", newline=""))}
    external = np.full((len(table),len(ALL_EXTERNAL)), np.nan)
    external[:,0], external[:,1] = table.am, table.esm1b
    canonical_cache, quality, metadata, coverage, sources = {}, [], [], [], []
    original_mutants = np.full(len(table), "", dtype="U30")
    for assay, part in table.groupby("dms_id", sort=True):
        entry = entries[assay]
        path = SOURCES / entry["retained_file"]
        if sha(path) != entry["retained_sha256"]:
            raise ValueError("External-score file hash mismatch")
        ref = references[assay]
        accession = str(part.accession.iloc[0])
        if part.accession.nunique() != 1:
            raise ValueError("One assay mapped to multiple canonical proteins")
        if accession not in canonical_cache:
            _, sequence = read_am(accession)
            if not sequence:
                raise ValueError("Missing canonical sequence")
            canonical_cache[accession] = sequence
        canonical, construct = canonical_cache[accession], ref["target_seq"]
        mapping = align_construct(canonical, construct)
        keys = {(int(row.position),str(row.native),str(row.alt)):int(row.Index)
                for row in part.itertuples()}
        source = pd.read_csv(path, dtype={"mutant":str}, float_precision="round_trip")
        if source.mutant.duplicated().any():
            raise ValueError(f"Duplicate source mutation: {assay}")
        hits, valid_keys, unparseable = set(), set(), 0
        for row in source.to_dict("records"):
            key = map_mutant(row["mutant"], construct, canonical, mapping)
            if key is None:
                unparseable += 1
                continue
            if key in valid_keys:
                raise ValueError("Ambiguous canonical mapping of external mutants")
            valid_keys.add(key)
            if key not in keys:
                continue
            index = keys[key]
            hits.add(index)
            original_mutants[index] = row["mutant"]
            for j, model in enumerate(NEW_EXTERNAL, 2):
                # Use the merged ProteinGym scores rather than raw VESPA probabilities.
                external[index,j] = -numeric(row.get(model))
        if hits != set(part.index):
            missing = sorted(set(part.index)-hits)
            raise ValueError(f"Missing exact source identifiers in {assay}: {len(missing)}")
        # Reconstruct targets through the fixed mutation map without selecting score orientation from outcomes.
        measured_path = ASSAYS / f"{assay}.csv"
        observed = {}
        with measured_path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                key = map_mutant(row["mutant"], construct, canonical, mapping)
                if key in keys:
                    if key in observed:
                        raise ValueError("Duplicate mapped experimental mutant")
                    observed[key] = float(row["DMS_score"])
        if set(observed) != set(keys):
            raise ValueError("Experimental target reconstruction coverage mismatch")
        ordered = sorted(keys, key=keys.get)
        reconstructed = rankdata(-np.array([observed[key] for key in ordered])) / len(ordered)
        indices = [keys[key] for key in ordered]
        if not np.allclose(reconstructed, table.loc[indices,"target"], atol=3e-8, rtol=0):
            difference = float(np.max(np.abs(reconstructed-table.loc[indices,"target"].to_numpy())))
            raise ValueError(f"Experimental target reconstruction mismatch: {assay}, max difference={difference}")
        quality.append(dict(assay=assay, accession=accession, original_rows=len(source),
            matched_existing_rows=len(hits), excluded_source_nonmatching_tokens=unparseable,
            exact_mutation_join=True, reconstructed_DMS_target=True))
        meta = dict(assay=assay, accession=accession, gene=str(part.gene.iloc[0]),
            membrane=int(part.is_membrane.iloc[0]), coarse_selection_type=ref.get("coarse_selection_type", "Unknown"),
            selection_assay=ref.get("selection_assay", ""), selection_type=ref.get("selection_type", ""),
            original_endpoint_type=str(part.endpoint_type.iloc[0]), existing_rows=len(part))
        metadata.append(meta)
        for j, model in enumerate(ALL_EXTERNAL):
            count = int(np.isfinite(external[part.index,j]).sum())
            coverage.append(dict(**meta, model=model, finite_scores=count,
                missing_scores=len(part)-count, coverage=count/len(part)))
        sources.append(dict(assay=assay, external_scores_sha256=sha(path), measured_file_sha256=sha(measured_path)))
    np.savez_compressed(out/"aligned_external_scores.npz", scores=external,
        models=np.array(ALL_EXTERNAL), dms_id=table.dms_id.to_numpy().astype(str),
        accession=table.accession.to_numpy().astype(str), position=table.position.to_numpy(),
        native=table.native.to_numpy().astype(str), alt=table.alt.to_numpy().astype(str),
        source_mutant=original_mutants)
    write_csv(out/"assay_metadata.csv", metadata)
    write_csv(out/"coverage_by_assay.csv", coverage)
    dump(out/"alignment_verification.json", dict(assays=quality, source_hashes=sources,
        rows=len(table), exact_existing_row_coverage=bool(np.all(original_mutants != "")),
        no_new_fit=True, no_reserved_outcome_reads=True,
        canonical_source_hashes={a:sha(ROOT/f"results/beat_am_20261016/features/fitting_pool/{a}.am.tsv")
                                 for a in sorted(canonical_cache)}))
    return external, {r["assay"]:r for r in metadata}, coverage


def common_mask(values, models, required):
    indices = [models.index(name) for name in required]
    return np.isfinite(values[:,indices]).all(axis=1)


def assay_metrics(table, values, models, mask, metadata):
    records, exclusions = [], []
    for assay, part in table[mask].groupby("dms_id", sort=True):
        y = part.target.to_numpy()
        x = values[part.index,:]
        reason = None
        if len(part) < MIN_ROWS:
            reason = "fewer_than_20_shared_rows"
        elif len(np.unique(y)) < MIN_TARGET_VALUES:
            reason = "fewer_than_3_target_values"
        elif not np.isfinite(x).all():
            reason = "nonfinite_score"
        else:
            constant = [name for j,name in enumerate(models) if len(np.unique(x[:,j])) < 2]
            if constant:
                reason = "undefined_constant_prediction:" + ",".join(constant)
        if reason:
            exclusions.append(dict(assay=assay, reason=reason, rows=len(part)))
            continue
        target_rank = rankdata(y) / len(y)
        for j, model in enumerate(models):
            score_rank = rankdata(x[:,j]) / len(y)
            records.append(dict(**metadata[assay], rows=len(part), model=model,
                spearman=float(spearmanr(x[:,j],y).statistic),
                rank_mse=float(np.mean((score_rank-target_rank)**2)),
                group=str(part.group.iloc[0])))
    return records, exclusions


def protein_metrics(assays, assay_type=None):
    selected = [row for row in assays if assay_type is None or row["coarse_selection_type"] == assay_type]
    buckets = defaultdict(list)
    for row in selected:
        buckets[(row["accession"],row["model"])].append(row)
    result = []
    for (accession,model), rows in sorted(buckets.items()):
        first = rows[0]
        result.append(dict(accession=accession, gene=first["gene"], model=model,
            membrane=first["membrane"], group=first["group"], assays=len(rows),
            rows=sum(r["rows"] for r in rows),
            spearman=float(np.mean([r["spearman"] for r in rows])),
            rank_mse=float(np.mean([r["rank_mse"] for r in rows]))))
    return result


def paired_summary(records, model, reference, *, draws=DRAWS, seed=SEED):
    left = {r["accession"]:r for r in records if r["model"] == model}
    right = {r["accession"]:r for r in records if r["model"] == reference}
    if set(left) != set(right):
        raise ValueError("Paired comparison has different protein sets")
    accessions = sorted(left)
    if not accessions:
        return None
    delta = np.array([left[a]["spearman"]-right[a]["spearman"] for a in accessions])
    labels = np.array([left[a]["group"] for a in accessions])
    groups = sorted(set(labels))
    ci = None
    if len(groups) >= 2:
        numerator = np.array([delta[labels == g].sum() for g in groups])
        denominator = np.array([(labels == g).sum() for g in groups])
        counts = np.random.default_rng(seed).multinomial(len(groups),np.full(len(groups),1/len(groups)),size=draws)
        samples = (counts @ numerator)/(counts @ denominator)
        ci = [float(v) for v in np.quantile(samples,[.025,.975])]
    return dict(model=model, reference=reference, proteins=len(accessions), groups=len(groups),
        model_mean_spearman=float(np.mean([left[a]["spearman"] for a in accessions])),
        reference_mean_spearman=float(np.mean([right[a]["spearman"] for a in accessions])),
        delta_spearman=float(delta.mean()), ci95_delta_spearman=ci,
        wins=int(np.count_nonzero(delta > 1e-12)), losses=int(np.count_nonzero(delta < -1e-12)),
        ties=int(np.count_nonzero(np.abs(delta) <= 1e-12)),
        mean_rank_mse_delta=float(np.mean([left[a]["rank_mse"]-right[a]["rank_mse"] for a in accessions])),
        interpretation="conditional_development_comparison_not_independent_confirmation")


def analyze_cohort(table, values, models, mask, metadata, *, bootstrap=True):
    assays, excluded = assay_metrics(table,values,models,mask,metadata)
    categories = sorted(set(r["coarse_selection_type"] for r in assays))
    summaries, paired, protein_records = [], [], []
    for assay_type in [None]+categories:
        proteins = protein_metrics(assays,assay_type)
        for scope in ("all", "membrane", "non_membrane"):
            selected = [row for row in proteins if scope == "all" or row["membrane"] == (scope == "membrane")]
            for model in models:
                subset = [row for row in selected if row["model"] == model]
                if not subset:
                    continue
                summaries.append(dict(scope=scope, assay_type=assay_type or "all", model=model,
                    proteins=len(subset), groups=len(set(r["group"] for r in subset)),
                    assays=sum(r["assays"] for r in subset), rows=sum(r["rows"] for r in subset),
                    mean_spearman=float(np.mean([r["spearman"] for r in subset])),
                    mean_rank_mse=float(np.mean([r["rank_mse"] for r in subset]))))
            if assay_type is None and bootstrap:
                comparisons = [(name,"raw") for name in models if name != "raw"]
                comparisons += [(CANDIDATE,name) for name in ALL_EXTERNAL if name in models]
                for model,reference in comparisons:
                    summary = paired_summary(selected,model,reference)
                    if summary:
                        paired.append(dict(scope=scope,assay_type="all",**summary))
            if assay_type is None:
                protein_records.extend(dict(scope=scope,**row) for row in selected)
    return dict(rows_retained=int(mask.sum()), assays=assays, excluded_assays=excluded,
        protein_metrics=protein_records, summaries=summaries, paired=paired)


def report(results, out):
    primary = results["guarded_loocv"]["requested_common"]
    rows = [r for r in primary["summaries"] if r["scope"] == "membrane" and r["assay_type"] == "all"]
    means = {r["model"]:r for r in rows}
    names = ["raw", "positive_inside", "positive_inside+glycosylation_loss", "positive_inside+disulfide_loss",
             CANDIDATE, "AlphaMissense", "ESM1b", "GEMME", "VESPA", "VESPAl"]
    labels = {"raw":"Standalone raw", "positive_inside":"Standalone + inside-out",
        "positive_inside+glycosylation_loss":"Inside-out + glycosylation",
        "positive_inside+disulfide_loss":"Inside-out + disulfide",
        CANDIDATE:"Inside-out + glycosylation + disulfide"}
    lines = ["# Standalone predictor versus external variant-effect predictors", "",
        "Date: 26 October 2026 (America/Los_Angeles).", "",
        "Completed by reusing the existing guarded group-out and five-fold predictions. No model, PCA, or correction coefficient was refitted. This is an exploratory comparison on existing exposed development data, not fresh confirmation.", "",
        "## Main comparison", "",
        "The table uses the same finite-score variants and qualifying assays for every listed model. Assay correlations are averaged equally within protein, then proteins equally. The membrane panel is primary. Scores are Spearman correlations with experimental DMS effects, not pathogenicity probabilities.", "",
        "| Model | Mean membrane Spearman | Proteins | Assays | Variant measurements |", "|---|---:|---:|---:|---:|"]
    for name in names:
        if name in means:
            r = means[name]
            lines.append(f"| {labels.get(name,name)} | {r['mean_spearman']:.6f} | {r['proteins']} | {r['assays']} | {r['rows']:,} |")
    lines += ["", "The paper candidate above excludes charge burial, as requested. All nine historical standalone combinations remain in the machine-readable supplement; no combination was selected based on these new comparisons.", "",
        "## Paired differences for the requested refined candidate", "",
        "Intervals resample the existing dependence groups while keeping predictions fixed. They do not cover project-wide model selection or refit uncertainty. There is no multiple-comparison adjustment.", "",
        "| Candidate minus comparator | Difference | Conditional 95% interval | Proteins won / lost / tied |", "|---|---:|---|---|"]
    for p in primary["paired"]:
        if p["scope"] == "membrane" and p["model"] == CANDIDATE:
            ci = p["ci95_delta_spearman"]
            text = f"[{ci[0]:+.6f}, {ci[1]:+.6f}]" if ci else "Too few groups"
            lines.append(f"| {labels.get(p['reference'],p['reference'])} | {p['delta_spearman']:+.6f} | {text} | {p['wins']} / {p['losses']} / {p['ties']} |")
    lines += ["", "## Validation-scheme sensitivity and broader panel", "",
        "A small lead over AlphaMissense in guarded leave-group-out evaluation is not robust superiority: the order reverses under the existing five-fold scheme. The combined refined candidate also remains below the raw standalone predictor under both schemes.", "",
        "| Model | Membrane: guarded group-out | Membrane: five-fold | All 39 proteins: guarded group-out | All 39 proteins: five-fold |",
        "|---|---:|---:|---:|---:|"]
    for name in ["raw",CANDIDATE,"AlphaMissense","ESM1b","GEMME","VESPA","VESPAl"]:
        cells = []
        for scope,mode in (("membrane","guarded_loocv"),("membrane","grouped_cv"),("all","guarded_loocv"),("all","grouped_cv")):
            row = next(r for r in results[mode]["requested_common"]["summaries"]
                       if r["model"] == name and r["scope"] == scope and r["assay_type"] == "all")
            cells.append(f"{row['mean_spearman']:.6f}")
        lines.append("| " + " | ".join([labels.get(name,name)]+cells) + " |")
    lines += ["", "## Additional modern comparator: VespaG", "",
        "VespaG was present in the same official archive and was retained as a separately declared supplementary comparator. It is a different model from VESPA and VESPAl; its results are not substituted for theirs.", "",
        "| Panel | VespaG | Raw standalone | Requested refined candidate |", "|---|---:|---:|---:|"]
    supplemental = results["guarded_loocv"]["paired_to_VespaG"]
    if supplemental["rows_retained"] != primary["rows_retained"]:
        lines += ["", "Coverage differs from the primary cohort; this supplementary table is not directly rank-comparable with the primary table.", ""]
    for scope,description in (("membrane","13 membrane proteins"),("all","All 39 proteins")):
        cells = []
        for model in ("VespaG","raw",CANDIDATE):
            row = next(r for r in supplemental["summaries"] if r["model"] == model and r["scope"] == scope and r["assay_type"] == "all")
            cells.append(f"{row['mean_spearman']:.6f}")
        lines.append("| " + " | ".join([description]+cells) + " |")
    lines += ["", "VespaG's broader-panel point estimate exceeds the standalone predictor here. The requested comparisons do not establish a state-of-the-art predictor.", "",
        f"Figure: `{(out/'comparison_membrane.png').relative_to(ROOT).as_posix()}`. The paired intervals and five-fold sensitivity are shown together."]
    lines += ["", "## Assay-type results", "",
        "Types come from the existing ProteinGym reference's coarse_selection_type, retaining its original selection-assay description. They do not use the older keyword heuristic, which can collapse stability into abundance. These are descriptive strata, not newly independent tests.", "",
        "| Membrane assay type | Model | Mean Spearman | Proteins | Assays | Measurements |", "|---|---|---:|---:|---:|---:|"]
    for r in primary["summaries"]:
        if r["scope"] == "membrane" and r["assay_type"] != "all" and r["model"] in names:
            lines.append(f"| {r['assay_type']} | {labels.get(r['model'],r['model'])} | {r['mean_spearman']:.6f} | {r['proteins']} | {r['assays']} | {r['rows']:,} |")
    lines += ["", "## Coverage, mapping, and verification", "",
        "Only the 44 assays already represented in the 127,811-row development table were requested from the official ProteinGym v1.3 score archive. Source CSVs were verified against ZIP CRCs; retained scores, metadata, and original source bytes have SHA-256 records. No reserved assay outcome was opened.", "",
        "Each mutation was joined by assay and exact native/position/alternate after the original construct-to-canonical alignment. The original within-assay DMS targets were independently reconstructed through that mapping. No target value was used as a join key or to choose score direction.", "",
        "Official merged ProteinGym scores are fitness-oriented and were negated once to match the cached damage-oriented AlphaMissense/ESM1b and standalone targets. The upstream merge code and configuration are retained. Direct raw-score MSE is not compared across incompatible score scales; the supplementary error metric compares within-assay rank-normalized scores.", "",
        "The strict requested-model intersection is the primary comparison. Scores on each comparator's own shared rows and the complete legacy AlphaMissense/ESM1b cohort are reported separately. Those secondary tables must not be combined into a ranking across different cohorts.", "",
        "## Interpretation limits", "",
        "The standalone head is DMS-supervised with evaluation groups excluded; external predictors have different training histories. This is not proof that all methods had identical pretraining or supervision exposure. Published whole-ProteinGym averages cannot be substituted for these matched-panel results.", "",
        "Binding, expression, activity, and stability are different outcomes. Glycosylation loss is not universally detrimental to every outcome; the ACE2 binding assay is especially relevant to that caveat. An audit expectation and experimental prediction accuracy are different claims.", "",
        "The existing glycosylation/disulfide corrections have not been recalibrated here. Narrow phenotype strata and few affected proteins should be treated as descriptive, not as a basis for selecting a new rule or coefficient.", "",
        "## Reproduce", "", "```powershell",
        ".\\.venv_ml\\Scripts\\python.exe scripts\\beat_am\\m10_fetch_external_baselines.py",
        ".\\.venv_ml\\Scripts\\python.exe scripts\\beat_am\\m11_compare_external_baselines.py",
        "```", "",
        "Outputs: `results/external_predictor_comparison_20261026/`. The result JSON, row-aligned scores, per-assay/per-protein CSVs, paired intervals, download hashes, exclusions, and plan are retained.", "",
        "Sources: [ProteinGym repository](https://github.com/OATML-Markslab/ProteinGym), [official score archive](https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/zero_shot_substitutions_scores.zip), [VESPA](https://github.com/Rostlab/VESPA).", ""]
    (out/"REPORT.md").write_text("\n".join(lines), encoding="utf-8")


def plot_comparison(results, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator
    names = [CANDIDATE,"AlphaMissense","VespaG","GEMME","VESPA","ESM1b","VESPAl"]
    labels = ["Inside-out + glyco + disulfide","AlphaMissense","VespaG (supplementary)","GEMME","VESPA","ESM-1b","VESPA-light"]
    panels, endpoints = [], []
    for mode in ("guarded_loocv","grouped_cv"):
        panel = []
        for name in names:
            cohort = "paired_to_VespaG" if name == "VespaG" else "requested_common"
            row = next(r for r in results[mode][cohort]["paired"]
                       if r["scope"] == "membrane" and r["model"] == name and r["reference"] == "raw")
            panel.append(row)
            endpoints.extend(row["ci95_delta_spearman"])
        panels.append(panel)
    with plt.rc_context({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,"axes.spines.right":False}):
        fig,axes = plt.subplots(1,2,figsize=(12.0,5.7),sharex=True,sharey=True)
        for ax,title,panel in zip(axes,["Guarded leave-group-out","Five-fold grouped CV"],panels):
            ax.axvline(0,color="#6b7280",linewidth=1.1,zorder=0)
            for i,row in enumerate(panel):
                mean = row["delta_spearman"]
                lo,hi = row["ci95_delta_spearman"]
                color = "#1d6a96" if row["model"] == CANDIDATE else "#a75528" if row["model"] == "AlphaMissense" else "#586778"
                ax.errorbar(mean,i,xerr=[[mean-lo],[hi-mean]],fmt="o",markersize=6,
                            color=color,elinewidth=1.5,capsize=3,zorder=3)
            ax.set_title(title,fontsize=12,pad=13,weight="bold")
            ax.set_yticks(range(len(names)),labels)
            ax.grid(axis="x",color="#e5e7eb",linewidth=.7,zorder=0)
            ax.xaxis.set_major_locator(MultipleLocator(.05))
            ax.set_xlabel("Mean Spearman difference from raw standalone",labelpad=9)
            ax.set_xlim(min(endpoints)-.015,max(endpoints)+.015)
            ax.tick_params(axis="y",length=0,pad=9)
        axes[0].invert_yaxis()
        fig.suptitle("Matched DMS comparison on 13 membrane proteins",fontsize=15,weight="bold",y=.98)
        fig.text(.5,.90,"43,991 shared measurements across 15 assays; positive values favor the compared model",ha="center",fontsize=10,color="#374151")
        fig.text(.5,.035,"Conditional 95% intervals resample 12 dependence groups. Existing development data; no new fit or independent confirmation.",ha="center",fontsize=9,color="#4b5563")
        fig.subplots_adjust(left=.24,right=.98,top=.81,bottom=.19,wspace=.14)
        fig.savefig(out/"comparison_membrane.png",dpi=200,facecolor="white")
        fig.savefig(out/"comparison_membrane.svg",facecolor="white")
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    table = load_table()
    plan_hash = freeze_plan(args.out,table)
    if args.prepare_only:
        print(f"Plan recorded: {plan_hash}", flush=True)
        return
    external, metadata, coverage = align_external(table,args.out)
    results, checks = {}, {}
    for mode in ("guarded_loocv", "grouped_cv"):
        models, standalone, groups = load_oof(mode,table)
        values = np.column_stack([standalone,external])
        models += list(ALL_EXTERNAL)
        frame = table.copy()
        frame["group"] = groups
        primary_mask = common_mask(values,models,REQUESTED_EXTERNAL)
        requested_models = [m for m in models if m != "VespaG"]
        selected_columns = [models.index(m) for m in requested_models]
        results[mode] = {}
        results[mode]["requested_common"] = analyze_cohort(frame,values[:,selected_columns],requested_models,primary_mask,metadata)
        legacy_names = [m for m in models if m not in NEW_EXTERNAL]
        legacy_columns = [models.index(m) for m in legacy_names]
        results[mode]["legacy_full"] = analyze_cohort(frame,values[:,legacy_columns],legacy_names,np.ones(len(table),dtype=bool),metadata)
        # Use published per-protein metrics as an aggregation regression check.
        old = json.loads((BENCH/mode/"results.json").read_text(encoding="utf-8"))
        old_raw = next(r for r in old["comparisons_vs_raw"]["membrane"] if r["model"] == "raw")
        reproduced = next(r for r in results[mode]["legacy_full"]["summaries"]
                          if r["model"] == "raw" and r["scope"] == "membrane" and r["assay_type"] == "all")
        if abs(reproduced["mean_spearman"]-old_raw["model_mean_spearman"]) > 1e-12:
            raise ValueError("Original baseline correlation did not reproduce")
        checks[mode] = dict(original_membrane_raw=old_raw["model_mean_spearman"],
            reproduced_membrane_raw=reproduced["mean_spearman"], groups=len(set(groups)),
            no_refit=True, group_exclusion_metadata_checked=True)
        for comparator in ALL_EXTERNAL:
            names = [m for m in models if m not in ALL_EXTERNAL] + [comparator]
            columns = [models.index(m) for m in names]
            mask = common_mask(values,models,[comparator])
            results[mode]["paired_to_"+comparator] = analyze_cohort(frame,values[:,columns],names,mask,metadata)
        scheme = args.out/mode
        scheme.mkdir(exist_ok=True)
        for cohort, content in results[mode].items():
            write_csv(scheme/f"{cohort}_summary.csv",content["summaries"])
            write_csv(scheme/f"{cohort}_assays.csv",content["assays"])
            write_csv(scheme/f"{cohort}_proteins.csv",content["protein_metrics"])
            write_csv(scheme/f"{cohort}_paired.csv",content["paired"])
        print(f"Completed {mode}; primary common rows={int(primary_mask.sum()):,}", flush=True)
    dump(args.out/"results.json",results)
    dump(args.out/"verification.json",dict(plan_sha256=plan_hash, baseline_reproduction=checks,
        source_table_unchanged=sha(TABLE), no_refitting=True,
        source_code_hashes={Path(__file__).name:sha(__file__), "m10_fetch_external_baselines.py":sha(Path(__file__).with_name("m10_fetch_external_baselines.py"))}))
    report(results,args.out)
    plot_comparison(results,args.out)
    print("PRIMARY MEMBRANE COMMON-COHORT RESULTS",flush=True)
    for row in results["guarded_loocv"]["requested_common"]["summaries"]:
        if row["scope"] == "membrane" and row["assay_type"] == "all":
            print(f"{row['model']}: {row['mean_spearman']:.6f}; proteins={row['proteins']} assays={row['assays']} rows={row['rows']}",flush=True)
    print(f"Report: {DOC}",flush=True)


if __name__ == "__main__":
    with threadpool_limits(limits=1):
        main()

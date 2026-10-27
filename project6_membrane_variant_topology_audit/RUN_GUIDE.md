Keep the enclosed project directory intact so that the scripts' relative paths resolve. Use the existing frozen files for the paper's exact cohorts. Several runners protect frozen plans or write to dated output paths. For a new run, choose a fresh output directory where the CLI supports one. Other paths are constants and require an explicitly recorded path adjustment before recomputation.

The delivered evidence is enough for checking saved tables and redrawing the Matplotlib figures. Retraining and regeneration from public sources require the inputs listed in REQUIRED_INPUTS.csv. Original model-cache files and ProtT5 embeddings remain separate because they are much larger than the scripts.

The recorded analysis environment is Python 3.12.13 with versions in requirements-analysis.txt. These versions were inspected during bundle assembly. The original GPU environment and the full environment used to typeset the final PDF were not recovered.

```text
cd project6_membrane_variant_topology_audit
python -m pip install -r requirements-analysis.txt
```

For the original research runners that import form_ab_lib, install the exact validated wheel. That module verifies the installed package against RELEASE_MANIFEST.json. The later working-tree modules are included under tools and are loaded by the refinement scripts through their own paths.

```text
python -m pip install --no-deps --force-reinstall dist/membrane-rule-audit-v0.1.0-validated-20261004/wheel/membrane_rule_audit-0.1.0-py3-none-any.whl
```

Installing the current package with python -m pip install ./tools provides the newer public CLI. The original Form A manifest check expects its earlier validated package, so use the wheel command above when reproducing that experiment.

The initial sequence prepares the raw data, canonical mappings, and historical grouping. a1 and a4 require the saved inventory metadata enclosed in the bundle. a12 consumes the saved MMseqs comparison. Its optional instance runner is included for provenance and requires a Linux root-path adjustment plus an MMseqs2 installation if rerun.

```text
python scripts/route_a_fetch_proteingym.py
python scripts/route_a_resolve_uniprot.py
python scripts/route_a_a0_feasibility.py
python scripts/beat_am/a4_group_manifest.py
python scripts/beat_am/a12_curate_pairs.py
python scripts/beat_am/a15_shuffle_test.py
python scripts/beat_am/a22_freeze_manifest.py
python scripts/beat_am/a6_extract_predictors.py
python scripts/beat_am/a9_l1_features.py
```

These historical scripts are not a turnkey replacement for the enclosed frozen cohort. Some depend on earlier inventory decisions and snapshots. Retain the recorded group_manifest_frozen.csv and table hash when reproducing the paper's 39-protein pool.

For the embeddings, use the original c1_embed_pool.py recipe with ProtT5 and all positions. It runs on CUDA and has ROOT fixed to /home/ubuntu/project6_membrane_variant_topology_audit. Change that location explicitly or run in that directory layout. Its --out option can point at the embedding directory expected by m4. The alternate c5_embed_fasta.py accepts a FASTA and output directory. requirements-embedding.txt names the additional packages. Their exact original versions were unavailable.

```text
python scripts/beat_am/instance/c1_embed_pool.py --model prot_t5 --all-positions --out EMBEDDING_DIRECTORY
```

Download or restore the cached UniProt snapshots before preparing the rule features. c1_disulfide_sources.py and p1_proline_sources.py contain the original acquisition code. The p1 runner also fetches structural metadata used by a separate proline experiment, so its full output exceeds what the final paper needs. The final paper uses its cached UniProt records through m1.

```text
python scripts/beat_am/m1_prepare_multirule_refinement.py
python scripts/beat_am/m4_benchmark_multirule.py --stage all --threads 8
python scripts/beat_am/m5_report_multirule_benchmark.py
```

m4 is the current benchmark entry point. It trains the 34 guarded outer folds and the five-fold appendix check, with three inner GroupKFold partitions. It creates the outer model cache used by the label-free full-position audits. The saved full-data head used for additional characterisation and the speed benchmark is created by b9_speed_benchmark.py. Its backend helper can select LightGBM if available, while the paper's m4 benchmark explicitly uses HistGradientBoostingRegressor.

```text
python scripts/beat_am/b9_speed_benchmark.py
python scripts/beat_am/m7_audit_standalone_rules.py
python scripts/beat_am/m8_audit_existing_refinements.py
python scripts/beat_am/m10_disulfide_internal_consistency.py
python scripts/beat_am/m17_positive_inside_standalone.py
python scripts/beat_am/m33_glycosylation_constraint_repair.py --stage audit --out results/glycosylation_constraint_repair_20261031/run2
python scripts/beat_am/m34_verify_glycosylation_repair.py --out results/glycosylation_constraint_repair_20261031/run2
python scripts/beat_am/m33_glycosylation_constraint_repair.py --stage accuracy --out results/glycosylation_constraint_repair_20261031/run2
python scripts/beat_am/m35_verify_current_disulfide.py
```

The listed run2 directory identifies the final paper's existing constrained calibration. A rerun with changed code or inputs requires a new output directory and corresponding source-path changes in consumers such as m35. Source hashes in the saved plans are deliberate checks against mixing versions.

The external score comparison and paired statistical calculations follow the main benchmark.

```text
python scripts/beat_am/m10_fetch_external_baselines.py
python scripts/beat_am/m11_compare_external_baselines.py
python scripts/beat_am/m20_dms_comparator_pvalues.py
```

For the full AlphaMissense/ESM1b census, restore the source score distributions and topology HTML, extract the mapped flanks, then run the unchanged panel estimator. run_rule_audit supplies the main Table 6 values and intervals.

```text
python scripts/extract_am_full_flank_scores.py
python scripts/extract_esm1b_flank_scores.py
python scripts/run_rule_audit.py
```

For the public GEMME census, the original candidate_members.json index and download manifest are included as small source metadata. Restore vespag_fasta_training_datasets.zip and the score matrices from the original workspace or public source. m28 uses the index and sequence archive before it can download or audit scores. The initial source-indexing runner was unavailable, so preserve the enclosed frozen candidate index.

```text
python scripts/beat_am/m28_gemme_membrane_census.py --fetch
python scripts/beat_am/m29_verify_positive_inside_expansion.py results/gemme_membrane_census_20261031
```

m29 accepts the result directory as its positional argument. m28a and m28b provide recovery paths for interrupted score transfers. The saved matched result and independent verification are included in the bundle.

For Tier 1 validation, regenerate the six full-stage scenarios with the original validated package installed, summarise them, and apply the documented width-criterion amendment if reproducing the final amended criteria file.

```text
python scripts/run_form_a.py --stage full
python scripts/summarize_form_a.py --results results/form_a_20261008
python scripts/apply_form_a_amendment.py
python scripts/validate_method_benchmark.py
python scripts/run_method_label_controls.py
```

For the runtime measurements, restore the published-model head assets and caches first. m12 --weights downloads the needed weights. m14 retrieves published ESM-2 inputs. These timings depend on the processor and numerical-library settings.

```text
python scripts/beat_am/m12_prepare_speed_assets.py --weights
python scripts/beat_am/m14_fetch_cached_esm2.py
python scripts/beat_am/m13_local_speed_benchmark.py --threads 1 8 --repeats 7
python scripts/beat_am/b10_refinement_cost.py
```

m15_finalize_speed_report.py verifies the retained provenance/runner_at_measurement_20261026.py against the original measured-source hash before amending its presentation. It preserves the timing samples and historical test record, and records publication-source hashes separately. Use copied measurement outputs for report regeneration.

For the Matplotlib figures, figures_f4p_trial.py is the script whose styling matches the supplied PDF. Figure 1 uses its --arch option. Figures 3, 4, and 5 use trial_accuracy, trial_refinement, and trial_forest. They read the enclosed saved tables. The final Figure 2 requires the missing TikZ block.

```text
python scripts/beat_am/figures_f4p_trial.py --arch
python scripts/beat_am/figures_f4p_trial.py
python scripts/paper_support/export_worked_example.py --out results/worked_example
```

For the bundle's existing software tests and saved-result checks, use the following commands. The three timing tests additionally require the separate published-model assets. The Windows sandbox used during assembly required standard-permission temporary directories and single-thread settings, recorded in TEST_RESULTS.json.

```text
python -m unittest discover -s tests -p "test_*.py" -v
python scripts/paper_support/recheck_final_paper.py --out ../MY_SAVED_RESULT_CHECKS.json
```

To repeat the added score-array checks against a full original workspace, supply its project root. This reconstructs the ten external comparator point estimates and the native K/R frequencies and checks them against the archived values. It carries forward the archived comparator confidence intervals without asserting that their missing original bootstrap runner has been recovered.

```text
python scripts/paper_support/recheck_final_paper.py --project-root FULL_PROJECT_DIRECTORY --score-audits --out ../MY_SCORE_RECONSTRUCTION_CHECKS.json
```

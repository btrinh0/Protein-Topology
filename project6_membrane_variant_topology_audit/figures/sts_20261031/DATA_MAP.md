# Figure data map: STS paper figures (31 Oct 2026)

## v2 redesign: names and changes

| v2 file | v1 equivalent | data change |
|---|---|---|
| `fig1_overview` | `fig1_pipeline` | none (illustration) |
| `fig2_validation` | same | (a) side-blind dropped; (b) production points dropped, label `>10^15x`; (c) dots instead of bars |
| `fig3_forest` | — | built from the one saved four-comparator run (see Fig 3 below) |
| `fig4_casestudy` | — | new; built by `scripts/beat_am/figures_casestudy.py` with `.venv_render` (see below) |
| `fig5_accuracy` | `fig4_accuracy` | raw column dropped (table only); mean CI = 95% t-interval across 13 proteins (computed) |
| `fig6_refinement` | `fig5_refinement` (a) | one standardized, oriented axis (see below) |
| `figS1_robustness` | `fig5_refinement` (b) | same 20 values; 2 significant figures |

**Fig 6 standardization.** Each audit saved its estimate divided by the SD of the
predictor's scores at that audit's eligible sites. The script reads that SD and
divides the estimate and interval by it:

| row | SD source (per model) |
|---|---|
| positive-inside | `all_audits.csv` `penalty_sd` |
| disulfide | `internal_consistency.json` `results.<model>.S.measurement.sd_of_values` |
| glycosylation A, B | `glycosylation_constraint_repair_20261031/run2/audit_results.json` `...original_all_motifs.{A,B}_noncytosolic.legacy_whole_protein.sd_of_z` |
| charge burial | `existing_refinement_audits.json` `charge_burial.legacy_results.z.sd_of_z_all_sites` (asserted equal to estimate / saved `delta_standardized_within_this_scorer`) |

Positive-inside and glycosylation S↔T are multiplied by −1 so that positive always
agrees with the rule. The SD is rule-specific, not one SD per predictor.

**Orders of magnitude (Fig 2b).** The smallest side-aware |Δ̂| divided by the largest
side-blind |Δ̂| is 5.785 × 10¹⁵ (log10 = 15.76). "Fifteen orders of magnitude" is
right in the sense of the ratio being more than 10¹⁵; v1's "≈16" rounded the
logarithm. The v2 label reads ">10¹⁵×", which agrees with the paper.

Written before any plotting. Every plotted value comes from a file listed here.
Paths are relative to `project6_membrane_variant_topology_audit/`.

**Location note.** The figure guide cites `results/standalone_refinement_20261016/`,
`docs/REFINEMENT_RUN_LOG.md`, `manifest.json` and
`scripts/standalone_refinement/figures_paper.py`. None of these exist in this folder
or in `btrinh0/Protein-Topology` `main`. The LaTeX paper (`latex/` on
`btrinh0/project-search` branch `claude/awesome-turing-y45vqs`) does not exist either:
that repository has only branches `claude/hopeful-archimedes-5kxtl0` and
`claude/zen-mayer-4e24cr`, each holding `PAPER_DRAFT_20261027.md`. Cross-checks below
therefore use (a) the values quoted in the figure guide and (b)
`PAPER_DRAFT_20261027.md`. The existing figure script is
`scripts/beat_am/figures_paper.py`; the new script imports its paths from it.

Tier 1 model types (`scripts/form_ab_lib.py`, `scripts/run_form_a.py`):
S0 = noise-free generator, S0n = noisy generator, S1 = full model (native, replacement,
side, distance), S2 = side-blind model (native, replacement only), S3 = side-free
chemistry model (native, replacement, K/R, D/E, hydrophobicity).

---

## Fig 1 `fig1_pipeline` (schematic)

No data. Box labels follow `PAPER_DRAFT_20261027.md` §3.1 (ProtT5-XL frozen; 64 PCA
components + 18 tabular features; histogram gradient-boosted head; PCA, head and
correction coefficients fitted on training partitions only).

## Fig 2 `fig2_validation`

| panel | file | columns / keys | row filter | missing? |
|---|---|---|---|---|
| A points | `results/form_a_20261008/estimates.csv` | `scenario, student, replicate, oracle_raw, estimate_raw` | `scenario` in P_neg05, P_neg10, P_pos10, P_decay; all 5 students; 20 replicates each | none |
| A slope | `results/form_a_20261008/criteria.json` | `S1.slope_through_origin` (0.8192481) | — | none |
| A band | `results/form_a_20261008/criteria.json` | `gates.slope` = [0.7, 1.1] | — | none |
| B side-blind strip | `estimates.csv` | `estimate_raw` | `student == "S2"` (all 6 scenarios, 120 runs) | none |
| B side-blind max | `criteria.json` | `S2_identity.max_abs_estimate` (3.6211e-18) | — | none |
| B real-model strip | `estimates.csv` | `estimate_raw` | `student == "S1"`, planted scenarios (80 runs) | none |
| B production | `results/method_completion_20261004/figure_manifest.json` | `main_figure_models[].estimate_raw` | alphamissense, esm1b | none |
| C bars | `results/method_completion_20261004/benchmark_results.json` | `cells[].zero_rejection`, `zero_rejection_mc95`, `bias` | `scenario` in shared_position_confound ("position sensitivity varies"), reference_composition_confound ("composition differs by side"); `method` packaged_within_site_stratified = full method, naive_charge_arm_only = without Stage 2, naive_unstratified_paired = without Stage 3 | none |

Computed, not read: per-(scenario, model) replicate mean and its 95% t-interval
(Fig 2A); an independent slope-through-origin from S1 rows as a check on
`criteria.json`; the orders-of-magnitude gap in Fig 2B.

Stage mapping evidence: `docs/PAPER_DRAFT_IMRAD_20261020.md` lines 432–433 name the
Stage 2 shortcut's failing scenario as position sensitivity (bias ≈ +0.600) and the
Stage 3 shortcut's as composition (≈ +0.543); `benchmark_results.json` gives bias
+0.5999 for naive_charge_arm_only and +0.5429 for naive_unstratified_paired in exactly
those scenarios.

## Fig 3 `fig3_forest`: built 31 Oct 2026 from the saved run

**Decision (31 Oct 2026, project owner): keep Fig 3 with both models, built from
`results/comparator_sensitivity_20261001/comparator_sensitivity.json` alone.** All
five comparator rows (`results.<model>.<comparator>.estimate_raw`, `ci95_raw`;
window 10) come from that one run, so they share one bootstrap. Its declared point
estimates equal `figure_manifest.json` to 1e-12; its declared intervals differ from
the paper's (separate bootstrap run): AlphaMissense [+0.0106, +0.0136] vs paper
[+0.0105, +0.0138], ESM1b [−0.170, −0.124] vs paper [−0.172, −0.121]. The paper text
must be updated to this run's intervals and alternative ranges. Contributing-site
counts for the declared comparator come from `figure_manifest.json`.

The earlier record of why the figure was first blocked follows.

| element | file | value in file | value in paper/guide | status |
|---|---|---|---|---|
| AM declared point | `results/method_completion_20261004/figure_manifest.json` | +0.012089 | +0.0121 | match |
| AM declared CI | same, `ci95_sd × normalizing_penalty_sd` | [+0.01047, +0.01377] | [+0.0105, +0.0138] | match |
| ESM1b declared point / CI | same | −0.14670 [−0.17183, −0.12087] | −0.147 [−0.172, −0.121] | match |
| site / protein counts | same, `contributing_sites`, `contributing_proteins` | 173,413 · 4,517 / 173,622 · 4,534 | same | match |
| AM four alternatives | `results/comparator_sensitivity_20261001/comparator_sensitivity.json` | −0.00121 to +0.01118 | **−0.051 to +0.059** | **MISMATCH** |
| ESM1b four alternatives | same | −0.36350 to −0.08377 | **−0.175 to −0.023** | **MISMATCH** |

The paper's alternative-comparator ranges first appear in
`docs/FULL_PAPER_DRAFT_20261013.md` Table 7 with no source file. A scan of every
JSON/CSV/MD under `results/` and `docs/` found no output that contains those four
endpoints together (the only co-occurrences are coincidental matches in raw DMS score
tables). `results/step2_20260929/step2_arms.json` has a different comparator
sensitivity (slopes, not contrasts) that does not match either. The only saved
four-alternative output (`comparator_sensitivity.json`, w = 10) disagrees with the
paper: under it, AlphaMissense's sign is negative only for the polar set (−0.0012,
interval spans zero), and ESM1b's range is −0.364 to −0.084. **Figure not made.**

**Decision (31 Oct 2026, project owner): Fig 3 is left out of the paper.** The §5.6
text still quotes the −0.051 to +0.059 and −0.175 to −0.023 ranges, which no saved
output supports.

## Fig 4 `fig4_accuracy`

| panel | file | columns | row filter | missing? |
|---|---|---|---|---|
| A dots and means | `results/dms_comparator_pvalues_20261028/per_protein.csv` | `accession, gene, model, spearman` | `scheme == "guarded_loocv"`; `model` in GEMME, VespaG, raw, positive_inside+disulfide_loss | none |
| A open circles (constrained three-rule) | `results/glycosylation_constraint_repair_20261031/run2/accuracy_per_protein.csv` | `accession, model, spearman`, `membrane == 1` | `model == PI_DS_constraint_glyco`; joined on accession; its `raw` and `existing_PI_DS` columns are checked equal to the comparator file | none |
| A aggregation check | `results/dms_comparator_pvalues_20261028/per_assay.csv` | `assay, accession, model, spearman` | same scheme/models; recompute mean across assays per protein | none |
| B differences | `per_protein.csv` (as above) | refined minus VespaG per accession | — | none |
| B test | `results/dms_comparator_pvalues_20261028/highlighted_configuration_main_comparisons.csv` | `delta_spearman, ci95_delta_spearman, p_two_sided` | `reference == "VespaG"` | none |
| raw vs refined Δ | `results/dms_comparator_pvalues_20261028/refinements_vs_raw.csv` | `delta_spearman, ci95_delta_spearman` | guarded_loocv, positive_inside+disulfide_loss | none |

Computed, not read: column means (Fig 4A), per-protein differences (Fig 4B).
Optional AlphaMissense/ESM1b columns not drawn: their per-protein files
(`external_predictor_comparison_20261026/guarded_loocv/paired_to_*_proteins.csv`) cover
a different protein set (scope "all", 780 rows), not the 13-protein panel.

## Fig 5 `fig5_refinement`

| mini-panel | file | key | before | after |
|---|---|---|---|---|
| positive-inside | `results/positive_inside_standalone_20261028/all_audits.csv` | `window == 10`, `comparator_name == "all_uncharged"`; `estimate_raw, group_ci_low, group_ci_high` | `guarded_oof_raw` | `guarded_oof__positive_inside+disulfide_loss` |
| disulfide gap | `results/standalone_multirule_20261025/disulfide_internal_consistency_20261025/internal_consistency.json` | `results.<model>` | `guarded_oof_PI` | `guarded_oof_PI_plus_disulfide` |
| glycosylation A, B | `results/glycosylation_constraint_repair_20261031/run2/audit_results.json` (updated 1 Nov 2026) | `guarded_group_out.<model>.original_all_motifs.{A,B}_noncytosolic.group_resampling` | `existing_PI_DS` | `PI_DS_constraint_glyco` |
| charge burial | `results/standalone_multirule_20261025/label_free_audit_20261025/existing_refinement_audits.json` | `results.<model>.charge_burial.grouped_endpoint` | — (no active coefficient; one point) | `positive_inside` (equals raw and the current three-rule model per the paper; run2 leaves every TM score unchanged) |
| B heatmap | `all_audits.csv` | `estimate_raw` by `comparator_name` × `window` (5, 10, 20, 60) | — | `guarded_oof__positive_inside+disulfide_loss` (20 cells) |

All "before" values are the leave-group-out (guarded) heads, matching the paper's
before/after tables. No values missing.

## Fig 6 `fig6_distance`: SKIPPED (optional)

Distance-binned K/R-gain contrasts exist only from older report-card runs:
AlphaMissense in `results/beat_am_20261016/rule/_card_work/alphamissense/membrane_rule_audit.json`
(held-out split, its own comparator) and ESM1b in a separate Form B report card. They
use different site populations, there is no binned output for the standalone
predictor, and `PAPER_DRAFT_20261027.md` cites none of them. Plotting them would add
values the paper text cannot be checked against, so the figure is not made.

## Fig 4 `fig4_casestudy` and supplement `supplement/figS_casestudy_<UNIPROT>.pdf`

Built by `scripts/beat_am/figures_casestudy.py`, run with
`.venv_render/Scripts/python.exe` (PyMOL open-source 3.2 from the
`pymol-open-source-whl` package; conda was not installed).

| element | file | key |
|---|---|---|
| DMS (experiment) | `results/standalone_multirule_20261025/benchmark_20261025/guarded_loocv/out_of_fold_predictions.npz` | `target` (within-assay rank damage, higher = damaging) |
| Ours | same file | `predictions[:, models == "positive_inside+disulfide_loss"]` (out of fold) |
| AlphaMissense, VespaG, GEMME | `results/external_predictor_comparison_20261026/aligned_external_scores.npz` | `scores` (stored higher = damaging; rows asserted identical to the file above); ESM1b is in the same file but not drawn (left out by the project owner, 31 Oct 2026) |
| topology, sequence | `results/standalone_multirule_20261025/annotation_bank.json` via `tools/membrane_rule_audit.refinement.topology_map` (the audit's own loader) | side 1 = cytosol, 2 = outside; `distance`; `tmds` |
| structures | `data/external/afdb_sasa_20261025/cif/` (4) and `data/external/afdb_casestudy_20261031/` (9, downloaded 31 Oct 2026 from alphafold.ebi.ac.uk, model v6) | AF-<acc>-F1-model_v6 |

Selection (`case_study_selection.csv`, written before any colouring): count of
positions with uncharged native residue, side known, 1–10 residues from a boundary,
and at least one K or R DMS row; most wins, ties to shorter sequence. Uses presence
of measurements only, no values. Chosen: P31645 (119); next P51681 (70).

Computed, not read: per-residue mean of K and R scores (averaging repeated assays),
within-source ranks, label choice (largest pairwise rank difference > 15, at most 6),
membrane frame (mean TM-segment direction, oriented by the annotated flanking side;
boundary heights = median projections of TM ends).

Panels were drawn only for proteins with at least 5 ranked residues: P31645,
P51681, Q9BQB6, P11712. The other 9 panel proteins have 0–2 eligible residues with
K/R measurements, so no panel.

**Main panel changed to rhodopsin (31 Oct 2026, project owner's choice).** P08100 has
only 2 residues under the audit eligibility, so the main panel uses "measured"
eligibility: any uncharged native residue with at least one K/R DMS row, at any depth
(29 residues: 21 membrane, 5 outside, 3 cytosol). Assay: OPSD_HUMAN_Wan_2019,
expression by flow cytometry (`assay_metadata.csv`). The rule-based pick (P31645)
and its table stay in `case_study_selection.csv`; SERT, CCR5, Q9BQB6 and P11712 are
in `supplement/` with the audit eligibility. ERBB2 (P04626) was considered: its DMS
covers residues 651–674 only (19 of 21 K/R positions inside its single TM helix).

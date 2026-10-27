# Standalone predictor versus external variant-effect predictors

Date: 26 October 2026 (America/Los_Angeles).

Completed by reusing the existing guarded group-out and five-fold predictions. No model, PCA, or correction coefficient was refitted. This is an exploratory comparison on existing exposed development data, not fresh confirmation.

## Main comparison

The table uses the same finite-score variants and qualifying assays for every listed model. Assay correlations are averaged equally within protein, then proteins equally. The membrane panel is primary. Scores are Spearman correlations with experimental DMS effects, not pathogenicity probabilities.

| Model | Mean membrane Spearman | Proteins | Assays | Variant measurements |
|---|---:|---:|---:|---:|
| Standalone raw | 0.462330 | 13 | 15 | 43,991 |
| Standalone + inside-out | 0.462194 | 13 | 15 | 43,991 |
| Inside-out + glycosylation | 0.459865 | 13 | 15 | 43,991 |
| Inside-out + disulfide | 0.462596 | 13 | 15 | 43,991 |
| Inside-out + glycosylation + disulfide | 0.460260 | 13 | 15 | 43,991 |
| AlphaMissense | 0.455125 | 13 | 15 | 43,991 |
| ESM1b | 0.402594 | 13 | 15 | 43,991 |
| GEMME | 0.427989 | 13 | 15 | 43,991 |
| VESPA | 0.404522 | 13 | 15 | 43,991 |
| VESPAl | 0.308022 | 13 | 15 | 43,991 |

The paper candidate above excludes charge burial, as requested. All nine historical standalone combinations remain in the machine-readable supplement; no combination was selected based on these new comparisons.

## Paired differences for the requested refined candidate

Intervals resample the existing dependence groups while keeping predictions fixed. They do not cover project-wide model selection or refit uncertainty. There is no multiple-comparison adjustment.

| Candidate minus comparator | Difference | Conditional 95% interval | Proteins won / lost / tied |
|---|---:|---|---|
| Standalone raw | -0.002071 | [-0.007711, +0.000836] | 4 / 4 / 5 |
| AlphaMissense | +0.005134 | [-0.038072, +0.048488] | 6 / 7 / 0 |
| ESM1b | +0.057666 | [+0.002789, +0.108847] | 7 / 6 / 0 |
| GEMME | +0.032271 | [-0.012814, +0.082034] | 7 / 6 / 0 |
| VESPA | +0.055737 | [+0.008565, +0.097713] | 10 / 3 / 0 |
| VESPAl | +0.152238 | [+0.073940, +0.239346] | 12 / 1 / 0 |

## Validation-scheme sensitivity and broader panel

A small lead over AlphaMissense in guarded leave-group-out evaluation is not robust superiority: the order reverses under the existing five-fold scheme. The combined refined candidate also remains below the raw standalone predictor under both schemes.

| Model | Membrane: guarded group-out | Membrane: five-fold | All 39 proteins: guarded group-out | All 39 proteins: five-fold |
|---|---:|---:|---:|---:|
| Standalone raw | 0.462330 | 0.445566 | 0.491563 | 0.484768 |
| Inside-out + glycosylation + disulfide | 0.460260 | 0.443578 | 0.490873 | 0.484105 |
| AlphaMissense | 0.455125 | 0.455125 | 0.492011 | 0.492011 |
| ESM1b | 0.402594 | 0.402594 | 0.421522 | 0.421522 |
| GEMME | 0.427989 | 0.427989 | 0.483680 | 0.483680 |
| VESPA | 0.404522 | 0.404522 | 0.466083 | 0.466083 |
| VESPAl | 0.308022 | 0.308022 | 0.410077 | 0.410077 |

## Additional modern comparator: VespaG

VespaG was present in the same official archive and was retained as a separately declared supplementary comparator. It is a different model from VESPA and VESPAl; its results are not substituted for theirs.

| Panel | VespaG | Raw standalone | Requested refined candidate |
|---|---:|---:|---:|
| 13 membrane proteins | 0.450895 | 0.462330 | 0.460260 |
| All 39 proteins | 0.504338 | 0.491563 | 0.490873 |

VespaG's broader-panel point estimate exceeds the standalone predictor here. The requested comparisons do not establish a state-of-the-art predictor.

Figure: `results/external_predictor_comparison_20261026/comparison_membrane.png`. The paired intervals and five-fold sensitivity are shown together.

## Assay-type results

Types come from the existing ProteinGym reference's coarse_selection_type, retaining its original selection-assay description. They do not use the older keyword heuristic, which can collapse stability into abundance. These are descriptive strata, not newly independent tests.

| Membrane assay type | Model | Mean Spearman | Proteins | Assays | Measurements |
|---|---|---:|---:|---:|---:|
| Activity | Standalone raw | 0.476745 | 4 | 4 | 17,866 |
| Activity | Standalone + inside-out | 0.476636 | 4 | 4 | 17,866 |
| Activity | Inside-out + glycosylation | 0.476636 | 4 | 4 | 17,866 |
| Activity | Inside-out + disulfide | 0.476836 | 4 | 4 | 17,866 |
| Activity | Inside-out + glycosylation + disulfide | 0.476836 | 4 | 4 | 17,866 |
| Activity | AlphaMissense | 0.509101 | 4 | 4 | 17,866 |
| Activity | ESM1b | 0.455525 | 4 | 4 | 17,866 |
| Activity | GEMME | 0.492774 | 4 | 4 | 17,866 |
| Activity | VESPA | 0.473936 | 4 | 4 | 17,866 |
| Activity | VESPAl | 0.398655 | 4 | 4 | 17,866 |
| Binding | Standalone raw | 0.375345 | 3 | 3 | 14,502 |
| Binding | Standalone + inside-out | 0.375398 | 3 | 3 | 14,502 |
| Binding | Inside-out + glycosylation | 0.365306 | 3 | 3 | 14,502 |
| Binding | Inside-out + disulfide | 0.375323 | 3 | 3 | 14,502 |
| Binding | Inside-out + glycosylation + disulfide | 0.365197 | 3 | 3 | 14,502 |
| Binding | AlphaMissense | 0.419180 | 3 | 3 | 14,502 |
| Binding | ESM1b | 0.373599 | 3 | 3 | 14,502 |
| Binding | GEMME | 0.406006 | 3 | 3 | 14,502 |
| Binding | VESPA | 0.355625 | 3 | 3 | 14,502 |
| Binding | VESPAl | 0.311158 | 3 | 3 | 14,502 |
| Expression | Standalone raw | 0.530247 | 6 | 6 | 10,160 |
| Expression | Standalone + inside-out | 0.530047 | 6 | 6 | 10,160 |
| Expression | Inside-out + glycosylation | 0.530047 | 6 | 6 | 10,160 |
| Expression | Inside-out + disulfide | 0.530822 | 6 | 6 | 10,160 |
| Expression | Inside-out + glycosylation + disulfide | 0.530822 | 6 | 6 | 10,160 |
| Expression | AlphaMissense | 0.486065 | 6 | 6 | 10,160 |
| Expression | ESM1b | 0.436907 | 6 | 6 | 10,160 |
| Expression | GEMME | 0.445886 | 6 | 6 | 10,160 |
| Expression | VESPA | 0.425054 | 6 | 6 | 10,160 |
| Expression | VESPAl | 0.274247 | 6 | 6 | 10,160 |
| OrganismalFitness | Standalone raw | 0.073956 | 1 | 1 | 224 |
| OrganismalFitness | Standalone + inside-out | 0.073956 | 1 | 1 | 224 |
| OrganismalFitness | Inside-out + glycosylation | 0.073956 | 1 | 1 | 224 |
| OrganismalFitness | Inside-out + disulfide | 0.073956 | 1 | 1 | 224 |
| OrganismalFitness | Inside-out + glycosylation + disulfide | 0.073956 | 1 | 1 | 224 |
| OrganismalFitness | AlphaMissense | 0.210238 | 1 | 1 | 224 |
| OrganismalFitness | ESM1b | 0.141651 | 1 | 1 | 224 |
| OrganismalFitness | GEMME | 0.143973 | 1 | 1 | 224 |
| OrganismalFitness | VESPA | 0.186027 | 1 | 1 | 224 |
| OrganismalFitness | VESPAl | 0.180728 | 1 | 1 | 224 |
| Stability | Standalone raw | 0.680860 | 1 | 1 | 1,239 |
| Stability | Standalone + inside-out | 0.680860 | 1 | 1 | 1,239 |
| Stability | Inside-out + glycosylation | 0.680860 | 1 | 1 | 1,239 |
| Stability | Inside-out + disulfide | 0.680860 | 1 | 1 | 1,239 |
| Stability | Inside-out + glycosylation + disulfide | 0.680860 | 1 | 1 | 1,239 |
| Stability | AlphaMissense | 0.567372 | 1 | 1 | 1,239 |
| Stability | ESM1b | 0.472680 | 1 | 1 | 1,239 |
| Stability | GEMME | 0.627412 | 1 | 1 | 1,239 |
| Stability | VESPA | 0.554557 | 1 | 1 | 1,239 |
| Stability | VESPAl | 0.515195 | 1 | 1 | 1,239 |

## Coverage, mapping, and verification

Only the 44 assays already represented in the 127,811-row development table were requested from the official ProteinGym v1.3 score archive. Source CSVs were verified against ZIP CRCs; retained scores, metadata, and original source bytes have SHA-256 records. No reserved assay outcome was opened.

Each mutation was joined by assay and exact native/position/alternate after the original construct-to-canonical alignment. The original within-assay DMS targets were independently reconstructed through that mapping. No target value was used as a join key or to choose score direction.

Official merged ProteinGym scores are fitness-oriented and were negated once to match the cached damage-oriented AlphaMissense/ESM1b and standalone targets. The upstream merge code and configuration are retained. Direct raw-score MSE is not compared across incompatible score scales; the supplementary error metric compares within-assay rank-normalized scores.

The strict requested-model intersection is the primary comparison. Scores on each comparator's own shared rows and the complete legacy AlphaMissense/ESM1b cohort are reported separately. Those secondary tables must not be combined into a ranking across different cohorts.

## Interpretation limits

The standalone head is DMS-supervised with evaluation groups excluded; external predictors have different training histories. This is not proof that all methods had identical pretraining or supervision exposure. Published whole-ProteinGym averages cannot be substituted for these matched-panel results.

Binding, expression, activity, and stability are different outcomes. Glycosylation loss is not universally detrimental to every outcome; the ACE2 binding assay is especially relevant to that caveat. An audit expectation and experimental prediction accuracy are different claims.

The existing glycosylation/disulfide corrections have not been recalibrated here. Narrow phenotype strata and few affected proteins should be treated as descriptive, not as a basis for selecting a new rule or coefficient.

## Reproduce

```powershell
.\.venv_ml\Scripts\python.exe scripts\beat_am\m10_fetch_external_baselines.py
.\.venv_ml\Scripts\python.exe scripts\beat_am\m11_compare_external_baselines.py
```

Outputs: `results/external_predictor_comparison_20261026/`. The result JSON, row-aligned scores, per-assay/per-protein CSVs, paired intervals, download hashes, exclusions, and plan are retained.

Sources: [ProteinGym repository](https://github.com/OATML-Markslab/ProteinGym), [official score archive](https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/zero_shot_substitutions_scores.zip), [VESPA](https://github.com/Rostlab/VESPA).

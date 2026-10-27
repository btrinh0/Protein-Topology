# Original DMS cross-validation: VESPA, VespaG and GEMME p-values

Computed 28 October 2026 (America/Los_Angeles). This report excludes AlphaMissense and ESM1b as requested. It reuses the original cross-validated scores without fitting a model, changing a coefficient or constructing a new split.

**Evaluation naming:** the primary scores are guarded leave-group-out predictions. The original five-fold grouped predictions are a sensitivity analysis. Both procedures used unshuffled group splits; calling these predictions "shuffled scoring" would be inaccurate. Sign flips below are statistical null calculations, not new train/test splits.

The primary panel contains 13 membrane proteins, 15 assays and 43,991 measurements across 12 dependence groups. Correlations are averaged equally across assays within each protein, then equally across proteins. All displayed models use the same rows.

## Score table

| Model | Mean Spearman |
|---|---:|
| raw | 0.462330 |
| positive_inside+disulfide_loss | 0.462596 |
| positive_inside+glycosylation_loss+disulfide_loss | 0.460260 |
| positive_inside+charge_burial+glycosylation_loss+disulfide_loss | 0.459003 |
| GEMME | 0.427989 |
| VESPA | 0.404522 |
| VespaG | 0.450895 |
| VESPAl | 0.308022 |

## Previously highlighted PI plus disulfide configuration

| Comparator | Difference in mean Spearman | Conditional 95% interval | Two-sided exact p | Holm p: three main comparisons | Holm p: all 27 main/configuration comparisons |
|---|---:|---|---:|---:|---:|
| GEMME | +0.034608 | [-0.009904, +0.083573] | 0.211426 | 0.422852 | 1.000000 |
| VESPA | +0.058074 | [+0.011003, +0.099848] | 0.048340 | 0.145020 | 1.000000 |
| VespaG | +0.011702 | [-0.027806, +0.050386] | 0.585938 | 0.585938 | 1.000000 |

Supplementary VESPA-light contrast: difference +0.154575, exact p = 0.004395, Holm p across all 36 configuration/comparator contrasts = 0.158203.

## Unchanged standalone predictor

| Comparator | Difference in mean Spearman | Conditional 95% interval | Two-sided exact p | Holm p: three main comparisons | Holm p: all 27 main/configuration comparisons |
|---|---:|---|---:|---:|---:|
| GEMME | +0.034342 | [-0.010077, +0.083772] | 0.215332 | 0.430664 | 1.000000 |
| VESPA | +0.057808 | [+0.010902, +0.099649] | 0.049805 | 0.149414 | 1.000000 |
| VespaG | +0.011436 | [-0.027991, +0.050276] | 0.592285 | 0.592285 | 1.000000 |

Supplementary VESPA-light contrast: difference +0.154308, exact p = 0.004395, Holm p across all 36 configuration/comparator contrasts = 0.158203.

## PI plus glycosylation plus disulfide configuration

| Comparator | Difference in mean Spearman | Conditional 95% interval | Two-sided exact p | Holm p: three main comparisons | Holm p: all 27 main/configuration comparisons |
|---|---:|---|---:|---:|---:|
| GEMME | +0.032271 | [-0.012814, +0.082034] | 0.248535 | 0.497070 | 1.000000 |
| VESPA | +0.055737 | [+0.008565, +0.097713] | 0.055176 | 0.165527 | 1.000000 |
| VespaG | +0.009365 | [-0.029765, +0.047885] | 0.660645 | 0.660645 | 1.000000 |

Supplementary VESPA-light contrast: difference +0.152238, exact p = 0.005371, Holm p across all 36 configuration/comparator contrasts = 0.177246.

## All four corrections

| Comparator | Difference in mean Spearman | Conditional 95% interval | Two-sided exact p | Holm p: three main comparisons | Holm p: all 27 main/configuration comparisons |
|---|---:|---|---:|---:|---:|
| GEMME | +0.031015 | [-0.016075, +0.082528] | 0.282227 | 0.564453 | 1.000000 |
| VESPA | +0.054481 | [+0.003207, +0.098530] | 0.077148 | 0.231445 | 1.000000 |
| VespaG | +0.008109 | [-0.032357, +0.048184] | 0.716309 | 0.716309 | 1.000000 |

Supplementary VESPA-light contrast: difference +0.150981, exact p = 0.006348, Holm p across all 36 configuration/comparator contrasts = 0.184082.

## Does refinement improve the unchanged standalone?

| Configuration | Difference | Exact p | Holm p: eight refinements |
|---|---:|---:|---:|
| positive_inside | -0.000136 | 0.562500 | 1.000000 |
| positive_inside+charge_burial | -0.001343 | 0.859375 | 1.000000 |
| positive_inside+glycosylation_loss | -0.002465 | 0.281250 | 1.000000 |
| positive_inside+disulfide_loss | +0.000266 | 0.656250 | 1.000000 |
| positive_inside+charge_burial+glycosylation_loss | -0.003706 | 0.429688 | 1.000000 |
| positive_inside+charge_burial+disulfide_loss | -0.000956 | 0.906250 | 1.000000 |
| positive_inside+glycosylation_loss+disulfide_loss | -0.002071 | 0.687500 | 1.000000 |
| positive_inside+charge_burial+glycosylation_loss+disulfide_loss | -0.003327 | 0.453125 | 1.000000 |

## Existing five-fold sensitivity

| Comparator | PI+D mean | Comparator mean | Difference | Exact p | Holm p: three comparisons |
|---|---:|---:|---:|---:|---:|
| GEMME | 0.445561 | 0.427989 | +0.017573 | 0.522461 | 1.000000 |
| VESPA | 0.445561 | 0.404522 | +0.041039 | 0.127441 | 0.382324 |
| VespaG | 0.445561 | 0.450895 | -0.005333 | 0.794434 | 1.000000 |

## Statistical interpretation

The two-sided exact block sign-flip test enumerates all 4,096 sign assignments across 12 dependence groups. OPSD and CCR5 share a group and always receive the same sign. The statistic preserves the equal-protein mean; groups are not incorrectly reweighted as one equal observation when they contain different numbers of proteins. The analysis does not treat 43,991 mutations as independent replicates.

The null calculation assumes sign-symmetric/exchangeable group differences. Because cross-validated heads share training data, this is a conditional exploratory comparison of existing predictions, not an independent repeated-training significance result. The percentile-bootstrap confidence interval and exact test are different statistical constructions and need not agree at a 0.05 threshold. Report the actual test result rather than choosing the method that is favorable.

Holm adjustment across three comparisons applies to a fixed displayed configuration. It does not correct selecting that configuration after viewing earlier results. More conservative columns adjust all 27 main-comparator/configuration contrasts, and the CSV/JSON also include adjustment over 36 contrasts when VESPA-light is included. Both schemes are saved in full. No claim of model-selection-adjusted or fresh confirmatory significance follows from a nominal p-value.

The external methods have different training objectives and supervision histories. The chosen DMS comparators are relevant to this outcome; omitting other predictors does not establish matched pretraining or supervision across the remaining methods.

Fixed hold-out experiments were archived at the user's request with byte-for-byte verification. Their historical outcome exposure remains in archive/fixed_holdout_20261028/. No archived hold-out scores were used to fit or select anything in this calculation.

## Reproduce

```powershell
.\.venv_ml\Scripts\python.exe scripts\beat_am\m20_dms_comparator_pvalues.py
```

Full per-assay/per-protein metrics, every configuration comparison, analysis specification, source hashes and independent exact-p verification are saved in results/dms_comparator_pvalues_20261028/.

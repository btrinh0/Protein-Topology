# Do the standalone refinements beat their baseline?

Date: 25 October 2026 (user's local date).

**Both requested evaluations are complete: five-fold grouped cross-validation, followed by leave-one-protein-out evaluation with relatives excluded together.** These are new refits on existing exposed development data, not fresh reserved-protein confirmation.

In guarded LOOCV, the full combination changes mean membrane-protein Spearman correlation from **0.462330 to 0.459003**, a change of **-0.003327** with conditional group-resampling interval **[-0.011392, +0.002683]**. A negative mean change; the conditional stability interval overlaps zero.

Five-fold grouped CV gives a full-combination change of **-0.002874**, interval **[-0.010038, +0.002766]**. A negative mean change; the conditional stability interval overlaps zero.

## Fair baseline and split design

The baseline is the same standalone architecture: 18 substitution/topology features plus 64 ProtT5 principal components, followed by a fixed gradient-boosted tree head. AlphaMissense and ESM scores are not inputs. The head has depth 6, learning rate 0.05 and 500 iterations; no model setting was selected using these new results. The older reported 0.4873 mean over 39 proteins came from another evaluation design and is not the comparator here.

All 127,811 existing measured rows from 44 assays and 39 proteins were used. The main comparison is the 13 membrane proteins to which these annotation features apply. The 26 other proteins receive zero correction and supply an exact scope check. Relatedness uses the conservative union of the two existing dependence-group registries: 34 groups overall, 12 among the membrane proteins. These automated historical groups are not claimed to be a newly validated biological family classification.

For each outer test group, three-fold inner group-out-of-fold predictions on the remaining proteins produce calibration residuals. Coefficients are fitted on those residuals and the outer-training membrane assay data only. PCA and the supervised head are refitted within every applicable training partition. The final outer head is refitted on all outer-training data; no outer test group participates in its PCA, supervised training or coefficient calibration. Each rule subset is independently calibrated with inactive terms fixed at zero, using the previously implemented nonnegative ridge calibration (ridge 0.001). This measures that specific calibration procedure; it does not test every possible method of enforcing the same chemistry.

Guarded LOOCV produces one held-out prediction for every protein, removing all known related-group members together. There are 34 distinct outer fits rather than 39 independent fits because some proteins share an excluded group. No old fitted correction coefficient is reused.

## All declared combinations

Values are the mean of protein-level Spearman correlations; assays are first averaged equally within protein. Deltas are relative to the raw baseline under the same evaluation scheme.

| Model | Five-fold mean rho | Change | Guarded LOOCV mean rho | Change | LOOCV delta interval |
|---|---:|---:|---:|---:|---|
| Raw standalone baseline | 0.445566 | +0.000000 | 0.462330 | +0.000000 | [+0.000000, +0.000000] |
| Positive-inside only | 0.445556 | -0.000010 | 0.462194 | -0.000136 | [-0.000493, +0.000079] |
| PI + Charge burial | 0.444667 | -0.000899 | 0.460988 | -0.001343 | [-0.007468, +0.002653] |
| PI + Glycosylation loss | 0.443585 | -0.001981 | 0.459865 | -0.002465 | [-0.007733, +0.000058] |
| PI + Disulfide loss | 0.445561 | -0.000005 | 0.462596 | +0.000266 | [-0.000393, +0.000967] |
| PI + Charge burial + Glycosylation loss | 0.442705 | -0.002861 | 0.458624 | -0.003706 | [-0.011555, +0.002161] |
| PI + Charge burial + Disulfide loss | 0.444665 | -0.000901 | 0.461374 | -0.000956 | [-0.007304, +0.003134] |
| PI + Glycosylation loss + Disulfide loss | 0.443578 | -0.001988 | 0.460260 | -0.002071 | [-0.007711, +0.000836] |
| Positive-inside + all three additions | 0.442692 | -0.002874 | 0.459003 | -0.003327 | [-0.011392, +0.002683] |

## Which addition helps most?

The first comparison adds one rule to positive-inside. The second measures the loss on removing that rule from the full combination, with the remaining weights refitted. A positive number favors including the rule. These two views distinguish an individual gain from a contribution that depends on another component.

| Addition | LOOCV: add to PI | Interval | LOOCV: full minus without this rule | Interval |
|---|---:|---|---:|---|
| Disulfide loss | +0.000402 | [-4.24e-07, +0.000993] | +0.000379 | [-0.000029, +0.000948] |
| Charge burial | -0.001207 | [-0.007389, +0.002854] | -0.001256 | [-0.007422, +0.002805] |
| Glycosylation loss | -0.002329 | [-0.007569, +0.000000] | -0.002371 | [-0.007705, +0.000000] |

The largest add-one point estimate is **Disulfide loss (+0.000402)**. This ranking is descriptive; comparison of several additions on this exposed data does not establish a uniquely best generalizable rule.

The disulfide addition has the best observed average, but its positive-inside-plus-disulfide model improves over the raw baseline by only +0.000266, with interval [-0.000393, +0.000967]. The gain does not establish a reliable improvement. Relative to positive-inside alone, disulfide improves 2 proteins, worsens 3 and leaves 8 unchanged; charge burial improves 8 but worsens 2 enough to give a negative mean; glycosylation worsens 1 and leaves 12 unchanged. These counts explain why the ranking by average size is not the same as the ranking by number of proteins helped.

## Experimental support and limits

| Feature | Measured variant rows | Distinct sites | Proteins | Dependence groups |
|---|---:|---:|---:|---:|
| Positive-inside only | 935 | 447 | 7 | 6 |
| Charge burial | 2318 | 599 | 10 | 9 |
| Glycosylation loss | 115 | 11 | 3 | 3 |
| Disulfide loss | 120 | 8 | 3 | 2 |

The number actually receiving a nonzero direct term can be lower because a calibration fit may return zero. For the full combination:

| Feature | Five-fold: nonzero direct term on measured test rows | LOOCV: nonzero direct term on measured test rows |
|---|---:|---:|
| Positive-inside only | 931 | 935 |
| Charge burial | 645 | 2318 |
| Glycosylation loss | 38 | 38 |
| Disulfide loss | 0 | 120 |

In five-fold CV, all measured disulfide examples occur in the same outer test fold. Its remaining calibration data have no disulfide-feature variation, so that fit returns zero. Consequently, the five-fold run supplies no direct held-out test of a nonzero disulfide term. Tiny differences between the disulfide-containing and disulfide-free combinations there arise through re-estimation of other weights, not a demonstrated bond-loss benefit. Guarded LOOCV separates the transporter and GPCR groups and addresses this particular coverage problem.

Disulfide and glycosylation coverage is particularly narrow. The disulfide result uses the requested control-free annotation feature, not a bonded-versus-verified-free cysteine experiment. It can evaluate predictive utility, not isolate an intrinsic bond-loss penalty. Bond endpoints are not forced to have equal scores. A direction-constrained feature alone does not guarantee total-score compliance with the corresponding rule; this benchmark is about functional prediction accuracy, not a new whole-census rule-compliance verdict.

All candidate features, assay outcomes and the model architecture had prior project exposure. The frozen new splits prevent direct evaluation-group leakage inside this run, but do not turn the project history into untouched confirmation. The reported 95% intervals resample dependence groups while keeping fitted cross-validated predictions fixed; they are conditional stability summaries, do not refit the entire learning procedure in every resample, and are not corrected for selecting a winner among several additions. Very small effects and few affected groups require particular caution.

## Individual-protein effects and secondary scopes

For the full combination in guarded LOOCV: **8 membrane proteins improved, 3 worsened, and 2 were unchanged**, using a 1e-12 numerical tolerance. Mean squared-error change was +0.000564 (negative is better).

- all: 39 proteins; raw mean rho 0.491563, full 0.490454, change -0.001109.
- non membrane: 26 proteins; raw mean rho 0.506180, full 0.506180, change +0.000000.

Full per-protein changes and every fold's correction weights are retained. A favorable mean must not be described as no individual protein worsening.

### Membrane-protein changes in guarded LOOCV

Every change below is relative to that protein's raw held-out baseline. These are paired protein results, not additional independent variant replications.

| Protein | Raw rho | PI + depth: change | PI + glyco: change | PI + disulfide: change | Full: change |
|---|---:|---:|---:|---:|---:|
| GLPA (P02724) | 0.514573 | +0.003567 | +0.000000 | +0.000000 | +0.003513 |
| ERBB2 (P04626) | 0.618568 | +0.002328 | -0.001910 | -0.001910 | +0.002328 |
| OPSD (P08100) | 0.641274 | +0.000743 | +0.000091 | +0.004747 | +0.005201 |
| MET (P08581) | 0.570600 | +0.000000 | +0.000000 | +0.000000 | +0.000000 |
| CP2C9 (P11712) | 0.568942 | +0.001552 | +0.000031 | +0.000030 | +0.001558 |
| LYAM1 (P14151) | 0.428857 | -0.009043 | +0.000000 | +0.000000 | -0.008977 |
| EPHB2 (P29323) | 0.680860 | +0.000000 | +0.000000 | +0.000000 | +0.000000 |
| SC6A4 (P31645) | 0.480575 | +0.003011 | +0.000186 | +0.000986 | +0.003863 |
| CCR5 (P51681) | 0.320307 | +0.002180 | +0.000108 | -0.000116 | +0.001956 |
| KCNH2 (Q12809) | 0.489747 | +0.004305 | -0.000536 | -0.000536 | +0.004157 |
| SCN5A (Q14524) | 0.073956 | -0.032938 | +0.000000 | +0.000000 | -0.032933 |
| VKOR1 (Q9BQB6) | 0.390075 | +0.006841 | +0.000260 | +0.000258 | +0.006906 |
| ACE2 (Q9BYF1) | 0.231958 | +0.000000 | -0.030277 | +0.000000 | -0.030821 |

These assays mix experimental endpoints. For example, the glycosylation-supported set includes ACE2_HUMAN_Chan_2020 (binding), OPSD_HUMAN_Wan_2019 (abundance) and SC6A4_HUMAN_Young_2021. An annotated chemical event is not by itself evidence that every assay readout should decrease. The benchmark tests one transferable directional correction against those measured readouts, not clinical pathogenicity.

## Reproducibility and integrity

The benchmark plan was frozen before opening experimental targets for this run. It records code/input/embedding hashes, all outer and inner group assignments, fixed model settings, every ablation, endpoints and interpretation limits. Post-run verification checked **39 outer folds**, including all inner training/calibration separations, inactive zero coefficients, finite predictions and the exact non-membrane scope check. Six known-answer split/metric tests passed in addition to the earlier 107 refinement/regression checks.

- Runner: `scripts/beat_am/m4_benchmark_multirule.py`.
- Frozen plan and fitted fold models: `results/standalone_multirule_20261025/benchmark_20261025/`.
- Five-fold results: `grouped_cv/results.json` and `grouped_cv/per_protein.csv` under that folder.
- Guarded LOOCV: `guarded_loocv/results.json` and `guarded_loocv/per_protein.csv`.
- Post-run verification: `BENCHMARK_VERIFICATION.json`.

No reserved experimental outcome was opened and no old model artifact or manuscript was overwritten. No new final full-data predictor was selected or released from these comparisons.

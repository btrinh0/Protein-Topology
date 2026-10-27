# Glycosylation refinement: a completed score-rule repair

31 October 2026. The new constrained refinement supports the expected directions
of both established glycosylation audits on the existing held-out membrane-panel
outputs. The original additive feature is unchanged; its coefficient now has a
label-free constraint floor. Mean DMS Spearman is slightly lower than PI+disulfide
alone and slightly higher than the previous three-rule configuration.

## What changed

The earlier coefficient-fitting objective minimized experimental residual error.
Its glycosylation weight moved the audits toward the expected directions but did
not complete the repair. The new procedure chooses the smallest nonnegative
coefficient satisfying the supported calibration proteins' two matched
glycosylation contrasts. No DMS outcome is used in this coefficient calculation.

For each protein, the original B8 matching produces affine contrasts
`A(beta) = A(0) + beta * a` and `B(beta) = B(0) + beta * b`.
The constraints are `A(beta) >= 0` and `B(beta) <= 0`, with zero margin.
The shared minimum coefficient is the maximum of their individual lower bounds.
Intervals and accuracy scores are not used to select that coefficient.

The original feature remains limited to directly supported glycosylation sites
in resolved mature outer loops. It adds a penalty to mutations that destroy an
annotated N-X-S/T sequon. It does not alter the sequence encoder or supervised
head, and does not add AlphaMissense or other predictor scores as inputs.

## Held-out evaluation

The existing guarded leave-group-out heads were reused. For each of the 12 outer
groups containing membrane proteins, calibration proteins were scored with that
outer head. The held-out group and known related group members were absent from
both its supervised training and the new coefficient calculation. Each existing
head reproduced the cached held-out scores to numerical precision.

The primary candidate holds the existing PI+disulfide coefficients fixed and
adds the constrained glycosylation term. A secondary candidate starts from the
raw standalone head. No new head or PCA was fitted. This is an exploratory
extension on the existing development panel, rather than a new confirmation set.

## Original motif-wide audits

Contrasts and group-resampling intervals use the unchanged B8 instrument and
the same sequence/topology panel, matched controls, weighting and resampling seed.

| Configuration | A: extra Pro penalty at X, expected positive | B: conserving S/T exchange, expected negative |
|---|---|---|
| Raw standalone | -0.01150 [-0.01869, -0.00387] | +0.01805 [+0.00847, +0.03080] |
| Existing PI+disulfide | -0.01140 [-0.01862, -0.00379] | +0.01811 [+0.00858, +0.03096] |
| Previous PI+glycosylation+disulfide | -0.00191 [-0.01351, +0.01654] | +0.00878 [-0.01844, +0.03005] |
| **PI+disulfide+constrained glycosylation** | **+0.04696 [+0.01473, +0.08041]** | **-0.03977 [-0.06310, -0.01592]** |

A contains 54 motif sites and 488 controls across 10 contributing proteins and
9 dependence groups. B contains 57 motif sites and 410 controls across 11
proteins and 10 groups. Both repaired intervals exclude zero in the expected
direction. Cytosolic contrasts are exactly unchanged by the glycosylation addition.

## Directly annotated correction scope

| Contrast | Repaired estimate | 95% group interval | Target sites | Proteins |
|---|---:|---|---:|---:|
| A | +0.12351 | [+0.09231, +0.14635] | 23 | 7 |
| B | -0.11157 | [-0.12208, -0.10141] | 26 | 9 |

All **16 eligible protein-by-contrast estimates** have the expected sign after
the repair. Previously 5 had the expected sign and 11 were opposite. These are
individual point-estimate directions, not 16 separate significance claims.
Four opposite-direction estimates remain among the broader potential-motif
comparisons outside these eligible matched comparisons.

The new glycosylation coefficient ranges from **0.09575 to 0.15978** across the
primary candidate's outer folds. These weights were fixed before the subsequent
accuracy calculation.

## Prediction accuracy

After the complete profiles and audit outputs were written and hashed, the
existing development benchmark predictions and targets were opened. Spearman
is calculated within each assay, averaged over assays within a protein, and
then averaged equally over the 13 membrane proteins.

| Configuration | Mean protein Spearman |
|---|---:|
| Raw standalone | 0.462330 |
| Existing PI+disulfide | 0.462596 |
| Previous PI+glycosylation+disulfide | 0.460260 |
| Raw + constrained glycosylation | 0.460555 |
| **PI+disulfide+constrained glycosylation** | **0.460834** |

The primary repair changes mean Spearman by **-0.001762** relative to
PI+disulfide, with paired group-resampling interval [-0.008324, +0.001968].
It changes mean Spearman by +0.000575 relative to the old three-rule version,
with interval [-0.000730, +0.002023]. These results support a biological
score-consistency repair; they do not establish an accuracy improvement.

Only three proteins' Spearman values change relative to PI+disulfide:

| Protein | Before | After | Change |
|---|---:|---:|---:|
| RHO | 0.646021 | 0.655862 | +0.009842 |
| SLC6A4 | 0.481561 | 0.482385 | +0.000824 |
| ACE2 | 0.231958 | 0.198387 | -0.033572 |

ACE2's assay measures binding. This is the principal accuracy cost of applying
the common glycosylation-loss penalty across this mixed-endpoint panel.
No coefficient was changed after inspecting these accuracy results.

## Other rules and software checks

- The original positive-inside result is exactly unchanged:
  -0.0202865 [-0.0283027, -0.0039752].
- Every C-to-S score is exactly unchanged, preserving the existing disulfide
  magnitude and pairing measurements that use this substitution.
- Every TM score is exactly unchanged, preserving the primary charge-burial
  centre/interface comparison.
- Native/missing masks and all scores outside the glycosylation feature support
  are unchanged.
- Eight known-answer tests passed, including directionality, unsupported
  constraints, held-out group rejection, and a synthetic matched-site example.
- Independent pandas reconstruction reproduced 30 audit rows, including their
  intervals and contributing counts, with maximum error 2.78e-17.
- The exported full-head coefficient artifact reproduces the repair through
  the existing `refine_position_scores` API exactly.

## Saved full-head profile

For the saved standalone raw head, the minimum coefficient is **0.1408278556**.
Its repaired motif-wide A is +0.04472 [+0.01304, +0.07491] and B is -0.04276
[-0.06974, -0.01728]. The coefficient artifact has all other terms set to zero,
so it applies specifically to that raw head. The held-out PI+disulfide candidate
uses the separate outer-fold profiles listed above.

Use the exported `coefficients` with the existing standalone refinement API,
starting from raw scores and applying refinement once. Recalibrate the floor if
the head, annotation scope, score scale, or other rule coefficients change.

## Files and reproduction

Results: `results/glycosylation_constraint_repair_20261031/run2/`.

- `audit_results.json`: both audit scopes and all configurations.
- `group_out_profiles.json`: coefficients, constraints and excluded groups.
- `per_protein_contrasts.csv`: individual protein directions.
- `accuracy_results.json`, `accuracy_per_protein.csv`: accuracy and paired intervals.
- `INDEPENDENT_VERIFICATION.json`: independent reconstruction and collateral checks.
- `saved_full_glycosylation_coefficients.json`: existing-API-compatible coefficients.

From the scientific project directory, reproduce into a new output folder:

```powershell
.\.venv_ml\Scripts\python.exe scripts/beat_am/m33_glycosylation_constraint_repair.py --stage audit --out results/glycosylation_constraint_repair_reproduction
.\.venv_ml\Scripts\python.exe scripts/beat_am/m33_glycosylation_constraint_repair.py --stage accuracy --out results/glycosylation_constraint_repair_reproduction
.\.venv_ml\Scripts\python.exe scripts/beat_am/m34_verify_glycosylation_repair.py --out results/glycosylation_constraint_repair_reproduction
.\.venv_ml\Scripts\python.exe -m unittest discover -s tests -p test_glycosylation_constraints.py -v
```

Earlier model and result artifacts were preserved. The completed audit's runner
snapshot matches its frozen source hash; the later accuracy-stage change only
avoided repeatedly decompressing a score matrix. No manuscript was edited.

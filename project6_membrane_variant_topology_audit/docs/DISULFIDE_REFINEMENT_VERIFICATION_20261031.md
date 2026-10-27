# Disulfide loss after refinement: pass confirmed

31 October 2026. The current PI+disulfide+constrained-glycosylation configuration
already passes the established annotation-defined disulfide-loss criterion.
All seven contributing proteins have positive C-to-S endpoint-minus-background
point estimates. No further coefficient or model adjustment was necessary.

## Established loss criterion

The comparison uses the identical C-to-S substitution at eligible annotated
disulfide endpoints and background cysteines, matched within protein and
compartment. A positive gap is expected. The original pass criterion requires
the 95% whole-protein resampling interval to lie above zero. The same
annotations select the correction and its evaluation targets, so this is an
internal consistency check; background sites are unannotated cysteines.

| Output | Matched gap | Original 95% interval | Result |
|---|---:|---|---|
| Unrefined guarded standalone | +0.059857 | [-0.000929, +0.104381] | Inconclusive |
| Existing PI+disulfide | +0.125719 | [+0.071391, +0.170225] | Pass |
| Current stack with constrained glycosylation | +0.125719 | [+0.071391, +0.170225] | Pass |

There are 44 eligible endpoint sites and 55 background sites in eight matched
strata across seven proteins. Resampling the six dependence groups instead
also gives a positive interval: [+0.070874, +0.172083]. The original seven-protein
interval remains the primary criterion; the resampling unit was not changed
to obtain a pass.

## Individual proteins

| Protein | Before disulfide refinement | Current refined stack |
|---|---:|---:|
| RHO | +0.086733 | +0.112959 |
| MET | +0.086006 | +0.142743 |
| SELL | +0.083808 | +0.172629 |
| EPHB2 | -0.052370 | +0.035941 |
| CCR5 | +0.139826 | +0.166053 |
| SCN5A | -0.004465 | +0.034910 |
| ACE2 | +0.161135 | +0.237068 |

EPHB2 and SCN5A, the two originally opposite-direction proteins, now have
positive estimates. These are per-protein directions, not seven separate
significance claims.

## Replacement sensitivities

| Replacement summary, current stack | Gap | Original 95% interval |
|---|---:|---|
| C-to-S | +0.125719 | [+0.071391, +0.170225] |
| C-to-A | +0.148753 | [+0.106114, +0.201336] |
| Mean of 19 non-native replacements | +0.123684 | [+0.057855, +0.174514] |

All three summaries support the expected positive association.

## What the previous inconclusive result referred to

The domain-matched bond-pairing comparison is a separate measurement: whether
scores at the two endpoints of a true bond are more similar than scores at
matched non-bond pairs. For PI+disulfide its difference is -0.037839 with interval
[-0.101008, +0.010870], which remains inconclusive. A common additive loss penalty
does not enforce equality of scores at two different positions. This pairing
result is not the pass criterion for the bond-loss penalty.

The new glycosylation correction leaves every C-to-S score exactly unchanged
relative to PI+disulfide. Therefore the established C-to-S loss and pairing
results remain unchanged in the current stack.

## Verification and reproduction

The current complete score matrices were rechecked without fitting a model,
coefficient, margin or threshold and without reading an experimental outcome.
The old PI+disulfide point estimate and original interval were reproduced to
1e-12 tolerance; per-protein weighted reconstruction matched the instrument;
source hashes remained unchanged.

Run from the scientific project directory:

```powershell
.\.venv_ml\Scripts\python.exe scripts/beat_am/m35_verify_current_disulfide.py
```

Saved outputs:

- `results/disulfide_refinement_verification_20261031/verification.json`
- `results/disulfide_refinement_verification_20261031/per_protein.csv`

The separate pairing figures above are retained in
`results/standalone_multirule_20261025/label_free_audit_20261025/existing_refinement_audits.json`.
The predictor, correction coefficients, accuracy results and manuscript were
not changed during this verification.

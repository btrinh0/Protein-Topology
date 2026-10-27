# Label-free rule audits of the standalone predictor

Date: 25 October 2026 (user's local date).

**Completed without a new fit, new data, or experimental outcome reads.** Existing B8/B9/B10/B10D audit functions were given scores from the standalone predictor itself. The three rules do not have the same status:

| Rule | Saved standalone head | Existing leave-group-out heads | Interpretation |
|---|---|---|---|
| Charge burial | Primary raw curve, endpoint, charge splits and scale check meet the expected directions | Same primary directions; scaled strict monotonicity has weaker resampling support | The aggregate depth trend is already present. |
| Glycosylation | Both non-cytosolic motif contrasts support the opposite direction | Both wrong directions repeat; cytosolic placebo intervals include zero | A reproducible score-behavior deficit under the established audit. |
| Disulfide | Pairing is lower than the random-pair null; proximity/domain controls remain inconclusive | Pairing signal repeats; primary controlled comparison remains inconclusive | Pairing is not loss magnitude. No clean certification of the full loss rule. |

The saved full-fitting-pool head is the primary object. The existing guarded-LOOCV heads are a sensitivity analysis using already fitted models; they exclude each scored dependence group. A score audit on fitted proteins measures behavior, not fresh generalization. Neither audit proves every individual variant score is correct.

## Inputs and execution

The available cached full-length embeddings cover 13 membrane proteins and 10,185 canonical positions. Each head supplied all 193,515 possible missense scores (19 per position); native identity cells were NaN. All three audits use this same available protein panel, not an asserted complete-proteome standalone panel. The scores were regenerated from frozen weights with the corrected 18-feature standalone construction and 64 ProtT5 components. No AlphaMissense or ESM prediction value entered the scorer. The earlier cached standalone census scores with the side-code/native-mask problems were not used.

The audit read cached sequence/topology, annotation and model/embedding artifacts. It did not read a DMS target table, ClinVar, or an experimental confirmation outcome. Model, PCA and correction fit methods were patched to raise an error during execution; all stages completed with those gates active. Existing input hashes were unchanged. No old model, coefficient or manuscript was edited.

## Glycosylation: opposite contrasts on the standalone output

The original B8 estimator compares the same substitution at motif and non-motif sites, matched within protein and native residue, separately by side. Its point estimator and controls were retained. Both the original whole-protein intervals and a dependence-group resampling sensitivity are saved. The table uses dependence-group intervals.

A: at N-X-S/T's X position, Pro breaks the motif. The excess Pro penalty should be larger at the motif than at matched non-motif positions, so the expected difference is positive.

B: at its S/T position, S↔T preserves the motif. The relative benefit of that swap should be stronger at the motif than at matched non-motif positions, so the expected difference is negative. A wrong-signed difference is not a claim that the preserving swap necessarily has a larger absolute damage score than every breaking mutation.

| Contrast | Expected sign | Saved full-fit head | Leave-group-out heads |
|---|---|---|---|
| A_noncytosolic | positive | -0.015512 [-0.027756, -0.007257] | -0.011496 [-0.018686, -0.003873] |
| B_noncytosolic | negative | +0.018143 [+0.007511, +0.028243] | +0.018048 [+0.008474, +0.030798] |
| A_cytosolic | placebo, near zero | +0.001068 [-0.007989, +0.014167] | -0.000614 [-0.003693, +0.005869] |
| B_cytosolic | placebo, near zero | +0.000851 [-0.019700, +0.008026] | +0.004187 [-0.014890, +0.011075] |

The non-cytosolic A arm contains 54 motif sites, 488 controls, 10 proteins and 9 dependence groups. B contains 57 motif sites, 410 controls, 11 proteins and 10 groups. The placebo intervals include zero; this is absence of a detected placebo effect, not a formal equivalence result.

Restricting the target arm to directly supported occupied glycosites—the refinement's actual scope—retains both opposite directions in both head sets. The eligible A arm uses 23 motif sites in 7 proteins; B uses 26 in 9 proteins. Thus the main finding is not solely a mismatch between potential motifs and annotated glycosites. This remains a conditional score comparison, with the original audit's structural and multiplicative-scale caveats; it is not an assay or clinical causal experiment.

### Do the glycosylation features actually move this audit?

A fixed +0.05 diagnostic coefficient was declared before scoring. It is not a fitted or validated repair. On the feature-eligible sites it changes A by exactly +0.05 and B by exactly −0.05, as intended. In the leave-group-out heads, both feature-eligible probe intervals support the expected direction:

| Eligible contrast | Raw leave-group-out | Fixed diagnostic probe |
|---|---|---|
| A_noncytosolic | -0.013303 [-0.031403, -0.004737] | +0.036697 [+0.018597, +0.045263] |
| B_noncytosolic | +0.022603 [+0.004629, +0.040165] | -0.027397 [-0.045371, -0.009835] |

The previously fitted experimental-residual coefficients also move both original contrasts toward the expected direction, but do not produce a supported pass:

| Original contrast | Existing PI-only scores | Existing PI + glycosylation scores |
|---|---|---|
| A_noncytosolic | -0.011404 [-0.018621, -0.003792] | -0.001944 [-0.013508, +0.016414] |
| B_noncytosolic | +0.018113 [+0.008579, +0.030961] | +0.008817 [-0.018337, +0.030046] |

The actionable finding is a calibration/objective gap, not a demonstrated arithmetic wiring bug: the feature can move the intended statistic, but the old coefficient-fitting objective minimized experimental residual error and did not enforce these audit inequalities. No coefficient was adjusted in this audit. An accuracy delta near zero therefore does not establish that the feature fails to represent the rule.

## Charge burial: already present

The B9 primary contrast is charge-introduction penalty at helix centres minus interfaces, matched within protein/native residue. Original minimum TM length is 10; the refinement-eligible length-17 sensitivity is also retained. The primary matched endpoint uses 235 centre sites and 260 interface sites in 52 strata and 8 contributing dependence groups.

| Head | Raw matched centre-minus-interface | Scaled matched endpoint (whole-protein interval) |
|---|---|---|
| saved_full_fit | +0.053299 [+0.043500, +0.075679] | +0.221214 [+0.156908, +0.320220] |
| guarded_oof | +0.039452 [+0.028827, +0.051985] | +0.146936 [+0.123535, +0.187746] |

Both raw point curves rise through all five bins, both charge signs support the gradient, and the endpoint remains positive after position-spread normalization. The leave-group-out scaled curve's strict five-bin ordering is reproduced in only about 43% of resamples, so its exact shape is less firm than its positive endpoint.

The previously fitted charge correction increases the leave-group-out depth endpoint from +0.039452 [+0.028827, +0.051985] to +0.101388 [+0.075243, +0.126967]. It amplifies a trend already present. This supports keeping charge burial as a compliance control; it does not show that every individual charge variant is already optimally predicted.

## Disulfide: separate pairing from magnitude

The B10 P2 audit compares true bond endpoints' C→S score differences with random re-pairing of the same endpoints within protein. It is label-free and runnable without free-Cys controls. That null preserves marginal endpoint scores; it does not remove every co-location confound, so the existing proximity and domain controls were also run.

| Head | True-minus-random mean absolute difference | Permutation p | Domain-matched residual and interval |
|---|---|---:|---|
| saved_full_fit | -0.014446 | 0.0045 | -0.005525 [-0.016593, +0.004743] |
| guarded_oof | -0.015802 | 0.0015 | -0.030102 [-0.082096, +0.010870] |

The original pairing set has 55 bonds in 7 proteins. The primary domain-matched comparison has 22 bonds in only 4 proteins, and its interval includes zero for both head sets. The separation-matched intervals also include zero. A position-order domain sensitivity supports the expected direction for the leave-group-out heads, but the primary domain control is inconclusive; this is not a robust clean pass after all controls.

On the directly eligible refinement bond subset, the saved head's pairing permutation p is 0.056 and the leave-group-out result is 0.0215 (47 bonds, 6 proteins). This dependence on head/annotation scope is another reason to avoid an unconditional disulfide pass/fail claim.

### The correction cannot repair a pairing statistic

Adding the same bond-loss offset to both endpoints leaves |score(endpoint 1) − score(endpoint 2)| unchanged. The within-protein re-pairing null is unchanged as well when all eligible endpoints receive that offset. This is a mathematical mismatch between a magnitude correction and a pairing endpoint, not evidence of a numerical defect.

Both the +0.05 diagnostic probe and the previously fitted disulfide coefficients leave the eligible pairing gap unchanged to numerical precision. For the leave-group-out heads, the gap stays −0.009899. The feature does change its intended magnitude association: a +0.05 probe increases the matched bonded-versus-unannotated C→S gap by exactly +0.05. The old fitted disulfide coefficients increase that association from +0.059857 to +0.125719.

That background set is unannotated, not experimentally verified free thiols. The legacy fields named `free` retain their old schema but must not be read as affirmative free-Cys evidence. P1's chemically controlled loss-magnitude question remains unidentified with the current bank. Pairing compliance, magnitude association and experimentally accurate loss penalties are distinct measurements.

## What the near-zero accuracy deltas mean now

- Charge burial: already satisfies the aggregate directional audit; adding more of that penalty is not a demonstrated rule repair.
- Glycosylation: a real opposite-direction audit result on this panel, reproduced with previously held-out heads. The feature has the correct direction, but its learned calibration does not complete the repair. A rule-enforcing objective would be a separate future fit, requiring its own accuracy evaluation.
- Disulfide: pairing signal is present but the stronger controls are inconclusive. A uniform loss feature cannot change P2. Its magnitude effect is implemented, but verified-free controls are absent for a clean P1 interpretation.

Thus the proposed two-way explanation needs qualification: audit behavior and accuracy benefit are different endpoints. An audit failure plus a near-zero accuracy delta alone does not prove an implementation bug. The no-fit perturbation checks identify whether the selected feature actually acts on the chosen audit endpoint.

The earlier five-fold/LOOCV accuracy findings and exploratory targeted-MSE result are unchanged. No new accuracy claim is made here. The positive-inside loop contrast from B9 is retained as a secondary reference but is inconclusive in this small panel; it neither validates nor overturns a whole-census claim from an older model.

## Reproducibility

- Raw score/audit runner: `scripts/beat_am/m7_audit_standalone_rules.py`.
- Existing-coefficient movement check: `scripts/beat_am/m8_audit_existing_refinements.py`.
- Output folder: `results/standalone_multirule_20261025/label_free_audit_20261025/`.
- Frozen plan, cached standalone scores, full legacy results, grouping sensitivities and probe checks are retained there.
- Five known-answer adapter tests passed. Point estimates reproduce the original B8/B9 estimators exactly; added grouping changes resampling units only. Input hashes and all probe identities were verified.

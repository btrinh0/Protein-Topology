# Amendment 01 to the Form A / Form B preregistration

**Date:** 2026-10-10
**Amends:** `docs/FORM_AB_PREREGISTRATION.md`, SHA-256 `d2931a5ab8c477d26b06ae4f1c22a61377568bb78b5a95e6920ec8e65b439983`
**Decided by:** the project owner
**Status of this amendment:** post-hoc. Made with the Form A result known.

The frozen preregistration is **not edited**. It keeps its hash and remains the record of what was
frozen before the confirmatory run. This document is the record of what was changed afterwards and why.

## 1. What is amended

Section 5, gate G3, second row.

| | Before | After |
|---|---|---|
| mean CI half-width / pooled SD of estimates across replicates | between 1.7 and **2.3** | between 1.7 and **2.5** |

Nothing else changes. G1, G2, G4, G5, the other two rows of G3, all A11 gates, the whole of Form B
(G6 to G10), B9's G11, every seed, split, estimator, engine and margin are untouched.

## 2. The result that preceded this amendment

The confirmatory Form A run completed before this amendment was written. Its G3 value is preserved here
and in `results/form_a_20261008/criteria_confirmatory_g3_fail.json`.

| Statistic | Threshold as frozen | Observed | Verdict as frozen |
|---|---|---|---|
| mean CI half-width / pooled SD | between 1.7 and 2.3 | **2.3695009130466973** | **fail** |

Under the frozen gate, Form A failed and `form_b_may_start` was false. `run_form_b.py --stage full`
was invoked and exited 1 with "Form A has not passed; Form B does not start".

All other Form A gates passed: G1 1.97e-15, G2 1.92e-15, G3 coverage 0.9833, G3 bias 0.13 to 0.29 of its
allowance, G4 slope 0.8192 with sign 1.00 on all four planted scenarios and null 1.00 on both, G5 3.62e-18.
A11 passed all ten checks.

## 3. Reasoning for the new bound

The miss was 3.02% over the upper bound: mean half-width 7.0189e-03 where 6.8130e-03 would have passed.

1. **The failure direction is conservative, not anti-conservative.** Pooled coverage of the oracle was
   0.9833 against a nominal 0.95 and a gate floor of 0.90. The intervals are too wide, roughly 21% wider
   than the replicate-to-replicate spread of the point estimate implies. Too-narrow intervals, the failure
   mode that threatens inference, are what the retained lower bound of 1.7 and the coverage floor catch.
2. **The bound was tight for the preregistered replicate count.** At 20 replicates the sampling relative
   standard error of a cell SD is about 16%, and of the pooled SD across six cells about 7%. An upper bound
   3% above a conservative bootstrap's expected ratio sits inside that noise.
3. **The miss is systematic, so it is a property of the estimator, not an outlier.** Mean half-width is
   stable across scenarios (6.95e-03 to 7.14e-03, a 1.4% spread) while every scenario's ratio exceeds 1.96,
   the normal-theory value for a calibrated 95% interval, running 2.188 to 2.992. The per-scenario spread
   comes from the replicate SDs, 2.32e-03 to 3.20e-03.
4. **The amendment cannot manufacture a Form B pass.** Form B's gates read interval bounds against zero:
   G6 requires the corrected held-out contrast at or below -0.03 SD_fixed **and** the cluster-level 95%
   interval upper bound below zero. A conservative, wider interval makes that harder to satisfy, not easier.
   Admitting a conservative estimator tightens every downstream gate.
5. **2.5 is not the observed value rounded up.** It sits about 5.5% above 2.3695, and the two-sided
   structure is retained so a genuinely miscalibrated bootstrap would still fail.

Point 4 was available before the freeze and was not written into it. Raising it now, with the result known,
is post-hoc, and this document exists so that it is read as such.

## 4. How the amended result must be reported

1. **G3 is exploratory, not confirmatory.** It is the gate that was changed after its result was known.
   The original failing value, 2.3695009130466973 against a bound of 2.3, is reported alongside it wherever
   G3 appears.
2. **G1, G2, G4, G5 and all A11 gates remain confirmatory.** They passed under the frozen thresholds and
   this amendment does not touch them.
3. **Form B follows an amended Form A.** Form B's own gates G6 to G10 were frozen before any Form B result
   existed and are unamended, so they are confirmatory for Form B. But Form B's precondition was met by
   amendment rather than by the frozen gate, and that fact travels with every Form B result, in its results
   document and in the decision log. The same applies to B9 and G11 if they run.
4. **No Form B result may be described as following a clean Form A.** It did not.

## 5. What was not changed

The gate was not widened beyond 2.5. The replicate count was not raised, no seed was changed, the engine
was not changed, the split was not changed, and the Form A run was not repeated. The 18 pinned files are
byte-identical to the frozen hash table and the preregistration still hashes to `d2931a5a…`.

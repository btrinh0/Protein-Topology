# Positive-inside audit of the standalone predictor

28 October 2026. Completed on existing standalone predictions, with no new fit, new coefficient calibration, network data, or experimental-target reads.

## Primary result

The primary statistic is the cytosolic-minus-outer difference in within-position K/R-introduction penalty, matched within protein and native residue, at distances 1–10 residues from a membrane boundary. Negative supports the specified positive-inside expectation. The comparator is the original 15 uncharged residues, excluding the native identity cell.

| Output | Raw contrast | 95% group-resampling interval | Result |
|---|---:|---|---|
| Saved full-fit raw | +0.011962 | [-0.000936, +0.017798] | inconclusive |
| Group-out raw | +0.008541 | [-0.007188, +0.015172] | inconclusive |
| Group-out: PI | -0.020391 | [-0.029660, -0.005042] | expected direction supported |
| Group-out: PI+charge | -0.020577 | [-0.029775, -0.005321] | expected direction supported |
| Group-out: PI+glyco | -0.020391 | [-0.029659, -0.005042] | expected direction supported |
| Group-out: PI+disulfide | -0.020287 | [-0.029555, -0.004843] | expected direction supported |
| Group-out: PI+charge+glyco | -0.020577 | [-0.029775, -0.005321] | expected direction supported |
| Group-out: PI+charge+disulfide | -0.020466 | [-0.029663, -0.005244] | expected direction supported |
| Group-out: PI+glyco+disulfide | -0.020286 | [-0.029554, -0.004842] | expected direction supported |
| Group-out: PI+charge+glyco+disulfide | -0.020466 | [-0.029663, -0.005244] | expected direction supported |

The primary comparison contains 596 contributing positions (306 cytosolic and 290 outer), 77 matched strata, 11 proteins and 10 dependence groups. The broader score archive contains 10,185 positions and 193,515 missense scores per raw head across 13 membrane proteins; that is not the independent sample size.

The saved full-fit head describes behavior on the fitted panel. Existing group-out heads exclude each scored dependence group from supervised fitting and correction calibration. Both use historically exposed project data; neither is new prospective confirmation.

## Correction movement

| Existing correction | Paired change versus group-out raw | 95% group-resampling interval |
|---|---:|---|
| Group-out: PI | -0.028932 | [-0.039358, -0.016249] |
| Group-out: PI+charge | -0.029118 | [-0.039633, -0.016706] |
| Group-out: PI+glyco | -0.028932 | [-0.039358, -0.016249] |
| Group-out: PI+disulfide | -0.028828 | [-0.039255, -0.016087] |
| Group-out: PI+charge+glyco | -0.029118 | [-0.039633, -0.016706] |
| Group-out: PI+charge+disulfide | -0.029008 | [-0.039510, -0.016545] |
| Group-out: PI+glyco+disulfide | -0.028828 | [-0.039255, -0.016087] |
| Group-out: PI+charge+glyco+disulfide | -0.029007 | [-0.039509, -0.016545] |

The coefficients were read unchanged from the existing guarded-fold metadata. A negative movement is movement in the declared audit direction; it is not an experimental accuracy gain.

## Comparator sensitivity at the primary window

| Output | Comparator | Contrast | 95% group-resampling interval | Result |
|---|---|---:|---|---|
| Saved full-fit raw | all_uncharged | +0.011962 | [-0.000936, +0.017798] | inconclusive |
| Saved full-fit raw | polar_neutral | +0.012860 | [+0.002327, +0.017627] | opposite direction supported |
| Saved full-fit raw | hydrophobic | +0.011368 | [-0.009808, +0.021364] | inconclusive |
| Saved full-fit raw | small | +0.014641 | [-0.000185, +0.022517] | inconclusive |
| Saved full-fit raw | aromatic | +0.006583 | [-0.003696, +0.010950] | inconclusive |
| Group-out raw | all_uncharged | +0.008541 | [-0.007188, +0.015172] | inconclusive |
| Group-out raw | polar_neutral | +0.008516 | [-0.005071, +0.014383] | inconclusive |
| Group-out raw | hydrophobic | +0.007985 | [-0.012024, +0.017070] | inconclusive |
| Group-out raw | small | +0.011215 | [-0.007322, +0.019054] | inconclusive |
| Group-out raw | aromatic | +0.004858 | [-0.002933, +0.008925] | inconclusive |
| Group-out: PI | all_uncharged | -0.020391 | [-0.029660, -0.005042] | expected direction supported |
| Group-out: PI | polar_neutral | -0.020416 | [-0.029406, -0.004386] | expected direction supported |
| Group-out: PI | hydrophobic | -0.020947 | [-0.032761, -0.003798] | expected direction supported |
| Group-out: PI | small | -0.017717 | [-0.028702, -0.000643] | expected direction supported |
| Group-out: PI | aromatic | -0.024074 | [-0.032490, -0.009719] | expected direction supported |

## Window sensitivity with the original comparator

| Output | Window | Contrast | 95% group-resampling interval | Result |
|---|---:|---:|---|---|
| Saved full-fit raw | 5 | +0.015467 | [-0.003625, +0.023598] | inconclusive |
| Saved full-fit raw | 10 | +0.011962 | [-0.000936, +0.017798] | inconclusive |
| Saved full-fit raw | 20 | +0.007880 | [-0.007571, +0.014974] | inconclusive |
| Saved full-fit raw | 60 | +0.001800 | [-0.009425, +0.008368] | inconclusive |
| Group-out raw | 5 | +0.010564 | [-0.007490, +0.018709] | inconclusive |
| Group-out raw | 10 | +0.008541 | [-0.007188, +0.015172] | inconclusive |
| Group-out raw | 20 | +0.006047 | [-0.013031, +0.014417] | inconclusive |
| Group-out raw | 60 | +0.000025 | [-0.013742, +0.007119] | inconclusive |
| Group-out: PI | 5 | -0.024256 | [-0.033605, -0.005747] | expected direction supported |
| Group-out: PI | 10 | -0.020391 | [-0.029660, -0.005042] | expected direction supported |
| Group-out: PI | 20 | -0.015620 | [-0.027882, -0.005529] | expected direction supported |
| Group-out: PI | 60 | -0.012578 | [-0.023134, -0.004411] | expected direction supported |

## Scope and reproducibility

The primary interval resamples dependence groups 1,999 times with predictions fixed. Whole-protein intervals and within-scorer standardization are secondary columns in `all_audits.csv`; standardized magnitudes must not be used as cross-model effect-size rankings. Comparator and window sensitivities were recorded before this run and are all retained. A confidence interval spanning zero is inconclusive, not proof of no biological preference.

Canonical sequences, all native-cell masks, explicit side-code conversion and cached-generation source hashes were checked. A plain-Python estimator independently reproduces every primary-comparator point from the packaged audit. Known-answer direction, null, offset-cancellation and side-reversal checks passed. The original inputs and model/coefficient artifacts are unchanged.

Runner: `scripts/beat_am/m17_positive_inside_standalone.py`. Results: `results/positive_inside_standalone_20261028/`. The plan, full sensitivity table, per-stratum contributions and verification record are retained there.

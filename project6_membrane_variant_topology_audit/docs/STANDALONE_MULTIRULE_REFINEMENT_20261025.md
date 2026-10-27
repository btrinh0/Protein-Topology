# Standalone predictor: three additional refinement components

Date: 25 October 2026.

**This document records the implementation/preparation stage, when real fitting and benchmarking were deferred as requested.** The subsequently requested [five-fold and guarded LOOCV comparison](STANDALONE_MULTIRULE_BENCHMARK_20261025.md) is now complete. Disulfide loss had the largest small positive add-one estimate; a reliable accuracy gain was not established, and the full combination had lower mean correlation.

This is a continuation of the standalone DMS-trained ProtT5 predictor and positive-inside work. It takes no AlphaMissense or ESM1b prediction as an input. The historical predictor-free results and manuscript are preserved. Statements below about deferred fitting describe this original preparation stage.

The new components are charge burial, N-linked glycosylation loss and disulfide loss. Their implementation does not depend on finding an AlphaMissense failure first. A component can be useful on the standalone predictor even if AlphaMissense already behaves in the expected direction.

## What the refinements add

The prepared architecture is an additive correction layer on the raw standalone ranking score:

`refined = raw + b_inside*f_inside + b_depth*f_depth + b_glyco*f_glyco + b_bond*f_bond`

Each coefficient is nonnegative, separately named and inspectable. Zero disables a component. Coefficients have **not** been fitted to real data in this stage; synthetic coefficients in the software-validation directory are test fixtures, never a usable fitted model.

| Component | Implemented feature | Scope and limit |
|---|---|---|
| Positive-inside | Signed K/R-gain adjustment, negative cytosol and positive outer side, decaying with distance from a TM boundary | Neutral native residues, known side, distances 1–60; unknown side is not treated as inside. The old fitted coefficient is not reused. |
| Charge burial | Charge introduction multiplied by fractional helix depth | Neutral native to D/E/K/R in annotated TM segments of at least 17 residues. Histidine is excluded. Depth is an annotation proxy, and channels or buried functional polar sites can be exceptions. |
| Glycosylation loss | A feature when a substitution destroys a supported N-X-S/T glycosylation site | N replacement, X-to-Pro, or S/T-to-a-residue-outside-S/T. S↔T is preserved. The default requires direct evidence for N-linked glycosylation and an outer mature-loop context; an unannotated sequon is not assumed occupied. |
| Disulfide loss | A feature for any ordinary non-Cys replacement at either endpoint of a supported bond | Direct annotation evidence by default; exact canonical endpoints. Known competing-function/signal sites, ambiguous partners and explicitly redox/transient/interchain bonds are excluded. No verified-free-cysteine bank is required. |

These are directional refinement features. They **do not guarantee that the complete predictor satisfies each aggregate rule**: the underlying model and other correction components can counteract a contribution. The final stage must measure both complete-score rule behavior and functional accuracy before a component is described as a successful repair. Scores remain ranking scores; they are not clipped into probabilities.

## The control-free disulfide alternative

The earlier controlled audit remains unchanged and still requires positively evidenced free cysteines. The new route deliberately asks a different, answerable question: **does adding an annotated bond-loss feature improve predictions of measured variant effects?**

Its nonnegative weight can later be learned against experimental residuals on calibration families. The calibration interface uses within-assay centering and equal family/protein/assay weights. It compares predictive behavior across ordinary measured missense variants; it does not require a chemically matched free-Cys subset. If there is no within-assay identifying variation, the interface returns an unsupported zero coefficient instead of claiming a rule pass.

This route cannot isolate an intrinsic disulfide-specific effect from all conservation, structural importance and assay-context differences. A held-out improvement would support a useful annotation-based refinement, not causal proof of bond-loss magnitude. Both partners receive the same additive term; their scores are not forced equal. The software explicitly checks that the correction does not manufacture the earlier endpoint-pairing result.

## Standalone integration and foundation fixes

`scripts/beat_am/m2_standalone_refinement_scorer.py` implements the successor feature builder and exposes `score_positions`, `refine_position_scores` and `same_model_probe`. Importing it does not load or run the real model.

It uses the supervised table's side codes: 1 cytosolic, 2 outer, 0 TM. Native identity cells are excluded from missense scoring and remain NaN. A correction probe starts from the supplied standalone scores. It never imports an AlphaMissense or ESM score archive. Missing embeddings fail explicitly rather than becoming zero vectors.

The new builder exactly reproduced **all 18 original tabular features on the existing supervised rows for each of the 13 membrane proteins**. This is a feature-consistency check, not a model-performance result. It corrects the census-construction problems documented in `archive/predictor_free_review_20261022/REVIEW.md`; the invalidated old enforcement verdict is not reinstated.

## Prepared annotation data

The preparation used cached UniProt JSON and canonical topology/sequence records, with their source hashes retained. It produced 193,515 possible single-missense feature rows over 10,185 positions in the existing 13 fitted membrane proteins:

| Annotation/feature coverage | Count |
|---|---:|
| Supported N-linked glycosylation sites | 26 |
| Supported disulfide bonds | 50 |
| Variant rows with positive-inside feature | 3,796 |
| Variant rows with charge-burial feature | 5,452 |
| Variant rows with glycosylation-loss feature | 988 |
| Variant rows with disulfide-loss feature | 1,900 |

These are engineering-coverage counts. Variants at one site, endpoints of one bond, and related proteins are not independent replications. This bank consists of existing fitted proteins and is not a new holdout. It demonstrates that all components can be constructed without the free-cysteine control bottleneck.

Preparation read identifier/topology/tabular columns only. It did not read the table's experimental target, AlphaMissense or ESM score arrays. No real model was loaded or scored, no real coefficient fitted, and no performance benchmark run.

## Final stage, after the additions

The implementation prepares nine comparisons: the raw baseline and positive-inside combined with every subset of the three additions. The full combination is one declared comparison; an apparent best combination must not be selected and then presented as independent confirmation on the same data.

Before fitting or evaluating:

1. Freeze shared family roles and a historical-exposure ledger. Separate supervised fitting, correction calibration and final evaluation; use verified group-out-of-fold predictions where needed. The two old overlapping train/test assignments cannot establish independence for one shared coefficient.
2. Fit preprocessing inside each applicable training partition. Use the corrected standalone construction throughout. Preserve the established experimental target orientation and scale.
3. Calibrate components only on their declared calibration rows, without tuning on final evaluation. A passing rule can be locked at zero. The implementation rejects evaluation-role rows and nonfinite outcomes during calibration.
4. Run the full ablation set and every rule audit at the end. Report functional accuracy, affected-region results, per-protein harm, uncertainty by family, and collateral changes across rules. A low aggregate cost cannot be described as no individual protein worsening.

The generic final-stage interfaces are available from the source tree:

```text
python -m tools.membrane_rule_audit.refinement fit --calibration CALIBRATION.npz --target-scale within_assay_rank_damage --label-provenance EXPERIMENTAL_SOURCE --out coefficients.json
python -m tools.membrane_rule_audit.refinement apply --input RAW_VARIANTS.npz --coefficients coefficients.json --score-direction higher-is-damaging --out refined_variants.npz
```

These commands are not a performance benchmark. The calibration NPZ must supply `features`, `feature_names` in the four-rule order, `baseline`, `target`, `accession`, `assay`, `family` and `role=calibration`. The application NPZ must supply `features`, `feature_names`, `baseline`, `score_state=raw` and aligned unique `accession`, `position`, `native`, `alt` identifiers. Feature/reference truth and cohort separation remain the caller's upstream responsibilities. The standalone integration builds the features directly from the prepared annotation bank.

## Files and verification

- Core refinement features, combination layer and calibration/application CLI: `tools/membrane_rule_audit/refinement.py`.
- Prepared bank, named feature archive, source hashes and exact training-feature parity: `results/standalone_multirule_20261025/`.
- Preparation: `scripts/beat_am/m1_prepare_multirule_refinement.py`.
- Standalone integration: `scripts/beat_am/m2_standalone_refinement_scorer.py`.
- Software verification: `scripts/beat_am/m3_verify_multirule_software.py`.

**107 software tests passed**, comprising 37 new feature/calibration/integration checks plus 70 existing disulfide, proline and charge/placebo regressions. The command interface recovers known synthetic coefficients and preserves synthetic variant identifiers through application. No synthetic result is biological validation. Logs and `verification.json` are in `results/standalone_multirule_20261025/software_validation/`.

The first preparation manifest is retained as written. Subsequent interface-validation changes are pinned separately in the verification source hashes; the prepared bank and feature archive retain their initial hashes. No old model artifact or manuscript was overwritten.

# membrane_rule_audit

Reusable within-site comparisons of mutation scores near membrane boundaries. Supply released predictor scores and a topology census, or an already mapped flank-score archive. The tool reports the specified contrast, contributing counts and machine-readable outputs.

An additional **disulfide-loss audit** uses a separate, matched-site estimator.
Its software is implemented; real-data validation is currently blocked by the
absence of a qualified native-context free-cysteine control set. See
[the implementation and blocker report](../../docs/DISULFIDE_RULE_STATUS_20261025.md).

The **Pro-gain structural audit** has been run on real AlphaMissense/ESM1b
scores. It found no defensible AlphaMissense opposite-direction failure;
ESM1b supports the conditional direction. Functional validation remains narrow.
See [the proline result and limitations](../../docs/PROLINE_TM_STATUS_20261025.md).

Each rule carries a declared directional prediction and is reported against neutral placebo rules, so the output is a verdict per rule rather than only a number. See [Placebo rules and verdicts](#placebo-rules-and-verdicts).

## Install and run

From the project root:

```text
python -m pip install ./tools
membrane-rule-audit --topology census.json --scores scores.tsv --score-direction higher-is-damaging --out results/my_model
```

Without installation, the general command is:

```text
python -m tools.membrane_rule_audit --topology census.json --scores scores.tsv --score-direction higher-is-damaging --out results/my_model
```

Only NumPy and Python 3.10 or newer are required. Score direction is mandatory: `higher-is-damaging` or `higher-is-tolerated`. A direction change negates the input before estimating contrasts.

## Disulfide-loss extension

The new command is `disulfide-rule-audit`. It compares the same Cys->Ser
replacement at matched, experimentally supported structural-bond endpoints and
verified free cysteines. Cys->Ala is a separately reported sensitivity analysis.
Every non-Cys replacement removes the original bond, so the membrane audit's
within-position centering is deliberately not reused.

Inputs are a canonical accession-to-sequence JSON map, an evidence/context site
CSV and a canonical variant-score CSV. The site schema and frozen study design
are in `docs/DISULFIDE_RULE_PROTOCOL_20261025.json` and clarification 01. Each
site carries sequence SHA-256, native Cys identity, evidence/source, native-context
review, family/data role, structural domain, compartment, secondary structure
and relative solvent accessibility. Bonded sites require the canonical Cys
partner and a reviewed structural role. Free controls require affirmative
experimental evidence and review of competing functions. Unknown is not free.

Matching is within accession/domain/compartment/secondary-structure strata,
smallest RSA distance first with a 0.15 caliper and no control reuse. Differences
are averaged within homology families and then equally across families; intervals
resample entire families. Six families is the frozen confirmation feasibility
minimum, not a power guarantee. Annotation/matching/exposure gates run before
the score file is opened. A blocked run writes its reasons and exits with code 2.

Run from the project root without installing:

```text
python -m tools.membrane_rule_audit.disulfide --sites curated_sites.csv --sequences canonical_sequences.json --scores canonical_scores.csv --score-direction higher-is-damaging --out results/my_disulfide_audit
```

Scores require `accession,position,reference,alternate,score`; coordinates are
one-based canonical coordinates, with explicit mapping before import. Assay
readout direction is mandatory. Missing pairs are reported, never imputed.
An experimental-validation run additionally declares `--kind experimental_validation`,
`--assay-endpoint`, and whether the assay context was reviewed.

For a declared synthetic known answer:

```text
python examples/generate_disulfide_example.py
python -m tools.membrane_rule_audit.disulfide --sites examples/disulfide_toy/sites.csv --sequences examples/disulfide_toy/sequences.json --scores examples/disulfide_toy/scores.csv --score-direction higher-is-damaging --out results/disulfide_rule_20261025/synthetic_demo
```

Both contrasts should equal +0.30 across six synthetic families. This tests
software behavior and does not demonstrate an AlphaMissense error.

The Python API is in `membrane_rule_audit.disulfide`: `eligible_sites`,
`match_sites`, `audit_pairs`, `fit_correction`, `apply_correction`, and
`run_disulfide_audit`. Fitting requires calibration-only pairs, an explicit
positive target on the predictor's own scale, and an endpoint-aligned experimental
direction established on development/calibration data. Confirmation outcomes
cannot select or gate the correction. A passing rule receives a zero correction.
Offsets apply only to verified eligible endpoints and the calibrated replacement;
corrected scores are not calibrated clinical probabilities.

Joint evaluation must include raw, positive-inside only, disulfide only and both.
A Cys->Ser shift can change the neutral comparator of the K/R-gain audit even
though K/R scores themselves are unchanged. The test suite explicitly covers
this interaction. No real two-rule model has been fitted while eligibility is
blocked.

## Proline introduction in regular TM helices

`proline-rule-audit` compares a within-position Pro-versus-comparator penalty
between independently defined regular TM windows and coil controls, within
protein/native-residue/packing strata. Native Pro loss is a different mutation
class and is excluded. A caller supplies a JSON structural site bank and a named
NPZ score matrix. Score direction is mandatory. NumPy is the only runtime
dependency of this audit; reproducing the structural site builder additionally
requires Biopython and the pinned experimental mmCIF sources.

```text
python -m tools.membrane_rule_audit.proline --sites results/proline_rule_20261025/audit_sites.json --scores results/proline_rule_20261025/alphamissense_site_scores.npz --score-direction higher-is-damaging --model AlphaMissense --out results/my_proline_audit
```

The NPZ schema is `scores` (positions x 20), `alphabet` (the declared
ACDEFGHIKLMNPQRSTVWY order), and per-row `accession`, one-based canonical
`position` and `native`. Archive rows are aligned by key, duplicate rows and
native mismatches are rejected, and missing alternatives remain unmeasured.
The JSON bank declares `accession`, `position`, `native`, `context` (`tm_helix`,
`coil`, `non_tm_helix`), `packing_bin` and uncertainty `family` group. Source
truth/geometry must be established upstream; a JSON context label is not
experimental verification by this command.

The primary comparator is A/I/L/M/V excluding native; Ser/Thr is a sensitivity.
Pro-minus-Gly is diagnostic rather than a universal ordering. All comparator
analyses share complete positions. Protein- and supplied-group resampling,
position-spread normalization and a nondirectional TM-versus-non-TM-helix
diagnostic are reported. Automated sequence linkage should not be assumed to
produce certified biological families; the delivered study records its grouping
problem and a separately labeled post-score curated-annotation sensitivity.

`scripts/beat_am/p1_proline_sources.py` through `p5_proline_report.py` reproduce
source retrieval, structural mapping, the frozen predictor/assay analysis,
grouping review and report. The two eligibility amendments preceded score
access. Added structural-panel proteins had no experimental holdout outcomes
opened. The measured check uses previously exposed Phase C assays only.

This candidate currently supplies an audit/positive-control result, not a
validated nonzero correction. The standalone ProtT5 model was not scored or
retrained for this candidate. The full report preserves all comparator results,
grouping sensitivities and sample-size limits.

## Small known-answer example

Generate the included synthetic input from the project root:

```text
python examples/generate_toy_example.py
membrane-rule-audit --topology examples/toy/toy_topology.json --scores examples/toy/toy_scores.tsv --score-direction higher-is-damaging --out results/toy_general
membrane-rule-panel --scores examples/toy/toy_mapped_scores.npz --model SyntheticExample --score-direction higher-is-damaging --out results/toy_panel
```

Both gain outputs should be −0.30 in raw score units, and both neutral placebo rules should be exactly 0.00, since the generator plants a charge effect only. The input has two artificial proteins and 608 variants; it contains no patient or real predictor data. Changing declared score direction flips the contrast.

The general command's legacy console label `sites kept` counts retained scored variant rows. This toy has 32 native positions and 16 eligible gain sites; the mapped panel reports eligible and contributing position counts separately. Use those position counts when describing the gain estimator's sample, rather than the 608 variant rows.

The [method validation](../../docs/METHOD_VALIDATION_20261006.md) records an actual wheel installation into an isolated environment, both installed commands, the API and source-byte checks, as well as the larger predeclared known-answer benchmark and shuffled-label controls. Plotting those benchmark results requires Matplotlib; the delivered tool requires only NumPy.

## Placebo rules and verdicts

Four charge rules are audited against two neutral placebo rules:

| Rule | Arm | Declared prediction for cytosolic minus outer |
| --- | --- | --- |
| `positive_charge_loss` | K/R removed | `+1` — losing a cytosolic K/R should score worse inside |
| `positive_charge_gain` | K/R introduced | `-1` — gaining a K/R should score worse outside |
| `negative_charge_loss` | D/E removed | none declared |
| `negative_charge_gain` | D/E introduced | none declared |
| `neutral_placebo_gain_vi` | V/I introduced | `0` |
| `neutral_placebo_gain_lm` | L/M introduced | `0` |
| `neutral_placebo_gain_nq` | N/Q introduced | `0` |
| `neutral_placebo_gain_q` | Q introduced | `0` |

The first two placebo groups are uncharged, conservative and hydrophobic, with no
side asymmetry of their own, so whatever they register is a side baseline the charge
rules also inherit. There are two because a placebo reading should not rest on one
choice of residues; report both.

The last two are better described as **polar probes** than as placebos. They test
whether a side effect tracks polarity rather than charge: if it does, the polar
groups carry it and the hydrophobic pair does not. Note that N is the
N-glycosylation residue, so gaining it on the outer side can create a sequon, which
makes `nq` ambiguous between polarity and glycosylation. `q` is the
non-glycosylatable isostere of N and isolates polarity on its own, so read the two
together and prefer `q` when they disagree. A placebo that excludes zero does not invalidate
a charge result, but it does mean the raw side difference is not wholly
charge-specific, which is what the contrasts below are for.

`verdict` reads each estimate against its declared prediction: `consistent`,
`inconsistent`, `indeterminate: interval includes zero`, `placebo consistent with
zero`, `placebo excludes zero`, or `no declared prediction`. A rule with no
declared prediction is described, never graded.

### Differences in differences

`rule_audit` reports contrasts between rules, each cancelling a side baseline the
two rules share:

| Contrast | Reads as |
| --- | --- |
| `positive_loss_minus_gain` | the rule's signed directional prediction, since it wants opposite signs on the two arms |
| `positive_gain_minus_placebo_vi` / `_lm` | the part of the K/R-gain side difference that is specific to charge |
| `negative_loss_minus_gain` | the same directional contrast for D/E |
| `negative_gain_minus_placebo_vi` | the charge-specific part of D/E gain |
| `positive_gain_minus_placebo_nq` / `_q` | the K/R-gain side difference net of a polar baseline |
| `placebo_nq_minus_placebo_vi` / `placebo_q_minus_placebo_vi` | polar versus hydrophobic: large means the side effect tracks polarity, not charge |
| `placebo_lm_minus_placebo_vi` | whether two hydrophobic pairs agree |

Every rule in a scope is resampled under **one shared protein draw** over the
union of contributing proteins, so the two rules in a contrast keep their
correlation. Resampling each rule on its own draw would treat them as independent
and widen the contrast interval. A rule whose protein set is a strict subset of
the union contributes a varying number of clusters per draw, which is ordinary
behaviour for a shared cluster bootstrap over a common population.

Each scope is reported twice: `all` over the retained distance range, and
`within_10_residues` over the near window.

### Matched-comparator audit

Holding each rule's group out of its own comparator leaves the rules on *different*
comparators: a charge rule keeps all 14 uncharged alternatives, because K/R/D/E were
never in that pool, while a placebo rule drops to 12. A contrast between two gain
rules then differs in its comparator as well as its numerator, so it is not a clean
numerator contrast.

`matched_gain_audit` fixes that. Every gain rule is run against one shared pool,
`ACFGPSTWY` — the uncharged alternatives with every gain group removed at once — and
then restricted to the sites eligible for all of them, so each contrast differs only
in which residues the numerator introduces. Two consequences to report honestly:

- the pool is narrower (9 residues, 8 where the native residue is one of them), so
  these estimates are **not** the headline definition and should be read alongside it,
  not in place of it;
- the site population is the intersection across rules, which is smaller than any
  single rule's. A singleton group such as `q` has no numerator where Q is native, so
  without the intersection its contrasts would compare two different populations.

Use `rule_audit` for the headline contrast and `matched_gain_audit` when the question
is whether one residue class behaves differently from another.

### Arm disjointness

The rule's residue group is held out of the opposing arm. For a charged group that
was already the case, since the group sits outside the uncharged comparator. For an
uncharged placebo group it is what keeps the two arms from sharing terms and
diluting the contrast toward zero, so a placebo added without it reads as a
spuriously small baseline. Charge-rule outputs are unchanged by this; see
`tests/test_rule_audit_placebo.py`, which asserts the charge arms are bit-identical
to the previous construction.

## General input format

CSV/TSV, optionally gzipped, with these columns:

```text
uniprot_id    protein_variant    score
P51681       R123W              0.84
```

Alternatively supply `ref`, `position`, `alt` columns. Column names can be overridden by the command options. NPZ inputs can contain `accessions`, `protein`, `position`, `ref_aa`, `alt_aa`, `score` (or `am_pathogenicity`), and `aa_alphabet`. Protein and amino-acid indices are zero-based; native positions are one-based.

The topology input is the project's viewer HTML with embedded `PROTEINS_INLINE`, or a JSON array such as:

```json
[{"uniprot":"P51681","gene":"CCR5","sequence":"MDYQ...","tmds":[{"start":55,"end":75}],"loops":{"1":{"start":1,"end":54,"location":"inside"}}}]
```

Loop location is `inside` for cytosolic or `outside` for outer. References must agree with the census sequence. The general interface computes positive-charge gain/loss and negative-charge gain/loss, excluding charge reversals and histidine. Distance bins and contributing/dropped counts are included.

Its composition output is a descriptive pooled calibration using gene groups, rather than the study's cross-fitted calibration using sequence-similarity groups. It is not an exact reproduction of the study's fitted slope or a validated model-correctness test. A matching shared-site key includes the native position so different residues in one protein remain distinct.

## Precomputed panel interface

```text
membrane-rule-panel --scores flank_scores.npz --model MyPredictor --score-direction higher-is-damaging --out results/my_model
```

The mapped NPZ format additionally requires per-variant `side` (0 cytosolic, 1 outer), `distance`, and `tmd_index`. The alphabet is `ACDEFGHIKLMNPQRSTVWY`. Identical duplicates are combined; conflicting scores or site metadata are rejected.

This interface packages the earlier raw K/R-gain comparison with the all-uncharged comparator `ACFGILMNPQSTVWY`. Synonymous alternatives are excluded. It contrasts charged and uncharged substitutions at each eligible uncharged reference, then computes a harmonic-size-weighted cytosolic-minus-outer difference within protein/reference strata. One-sided strata do not contribute. The default window is ten residues.

The saved JSON includes source and implementation hashes, direction, settings, counts and the point estimate. With at least two contributing proteins, the default 1,999 draws supply a protein-resampling stability interval. Such an interval does not by itself establish independent-family or clinical generalization.

## Python API and two-model demonstration

```python
from membrane_rule_audit import run_panel_audit

summary = run_panel_audit(
    "flank_scores.npz", "results/my_model", "MyPredictor",
    "higher-is-damaging",
)
```

Run `python scripts/run_rule_audit.py` from the project root to apply this same function to AlphaMissense and ESM1b. The two-model demonstration uses raw scores, the original comparator and the original ten-residue window. No alternative-comparator, transformed-scale, context or simulation checks are invoked.

## Scope

Within-site subtraction removes additive position effects on the selected score scale. It does not isolate charge from every property of an amino-acid substitution. Report the comparator, score direction and units with the result.

The package demonstrates reusable score analysis. It does not determine an internal learned mechanism, insertion energy, clinical harm, or which predictor is biologically correct. The scientific history is retained separately in the project's decision log.

## Standalone multi-rule refinement components

The additive correction implementation in `membrane_rule_audit.refinement`
supplies named positive-inside, charge-burial, glycosylation-loss and disulfide-loss
features. `prepare_record` requires matching canonical annotation/sequence;
`variant_features` constructs missense-only features; `apply_refinement` combines
explicit nonnegative coefficients. `fit_coefficients` is a later calibration API
with family/protein/assay weighting and evaluation-role rejection. The
`multi-rule-refine` entry point exposes calibration and application to supplied
tables. The later [development benchmark](../../docs/STANDALONE_MULTIRULE_BENCHMARK_20261025.md)
evaluated fitted coefficients with grouped CV and guarded LOOCV. It found no
established reliable accuracy gain; the full combination reduced mean correlation.

The control-free disulfide option evaluates an annotation-based predictive
refinement without a verified-free-cysteine bank. It does not replace the
controlled disulfide estimator or identify a causal bond-loss magnitude. These
features do not guarantee total-score rule compliance. See the
[implementation and final-stage specification](../../docs/STANDALONE_MULTIRULE_REFINEMENT_20261025.md)
for schemas, integration, limitations and verification.


# Retained source-method notes

These module descriptions and extended rationale were retained from the original source before the maintained Python introductions were shortened. They describe historical implementations and settings; use REPRODUCING.md for current entry points and known limitations. Original byte hashes are recorded in the source manifest.

## examples/generate_toy_example.py

Create small, fully synthetic inputs for both installed command interfaces.

These examples are software fixtures, not measurements or real predictor data.
The mapped K/R-gain contrast has a known cytosolic-minus-outer value of -0.30.

## scripts/apply_form_a_amendment.py

Recompute Form A's criteria under Amendment 01, preserving the confirmatory verdict.

Amendment 01 (`docs/FORM_AB_AMENDMENT_01.md`, owner decision, 2026-10-10) raises G3's
width-ratio upper bound from 2.3 to 2.5. It is post-hoc: the Form A result was known when
it was made.

This script does not edit any pinned file and does not edit the preregistration. It reuses
`summarize_form_a.criteria` unchanged, with only `GATES["width_ratio"]` replaced in memory,
so the amended verdict is computed by exactly the frozen code path.

Before writing anything it copies the confirmatory verdict to
`criteria_confirmatory_g3_fail.json`. That file is the record of the frozen gate's outcome and
is never overwritten once it exists. The amended `criteria.json` carries explicit provenance
fields naming the amendment, the original G3 value, the original verdict and G3's exploratory
status, so no consumer of that file can read it as a clean confirmatory pass.

## scripts/beat_am/a12_curate_pairs.py

A12: curate the 22 homologous pairs, then rebuild the split from the survivors.

MMseqs2 (instance) and the local Biopython audit disagreed badly: 4 pairs vs 22,
with 0 found by MMseqs2 alone. Neither is simply right. MMseqs2's k-mer
prefilter and e-value cutoff lose remote homologs near 30% identity; exhaustive
local Smith-Waterman with no low-complexity masking invents pairs from
compositionally biased segments. BLAST and MMseqs2 both mask low complexity by
default; a4 did not, which is the defect being corrected here.

Three independent lines of evidence per pair, none of them my own alignment
score:

  1. Curated family membership from the Oxford census annotation
     (isGPCR / isIonChannel / isSLC / isABC / isATPase / isE3Ligase).
  2. Whether MMseqs2 independently found the pair.
  3. Composition of the aligned segments: Shannon entropy and the largest
     single-residue fraction. Real homologous segments look like normal protein
     (entropy near 4.0 bits over 20 letters); low-complexity artifacts do not.

A pair is dropped only when it fails ALL supporting tests and shows positive
evidence of compositional bias, because over-grouping merely costs training
data while under-grouping risks leaking a confirmation protein's relative into
a fitted stage.

Creates new files only.

## scripts/beat_am/a15_shuffle_test.py

A15: composition-preserving shuffle test on the UNRESOLVED homology pairs.

a12 left 9 pairs UNRESOLVED: 30-35% identity over long, compositionally normal
alignments. They are not low-complexity artifacts, but nor is there positive
evidence they are homologous. One of them decides a real question:

  G003 chains KCNQ4 --[SUPPORTED]-- KCNQ1 --[UNRESOLVED]-- AMFR
                                       |--[UNRESOLVED]-- PPM1D

KCNQ4 is the only branch-(b) candidate that adds genuine two-sided K/R-gain
flank coverage (69 cytosolic / 20 outer), so admitting it is what would give the
flank claim its missing margin (n=6 -> 7, and the smallest attainable exact p
0.031 -> 0.016). But admitting it bars every G003 member from fitting, and
**PPM1D was already fitted in L1**. If the two UNRESOLVED links are chance
alignments, G003 splits and KCNQ4 becomes nearly free.

The test. For each pair, align, then re-align the same pair 300 times with one
sequence SHUFFLED, preserving its exact amino-acid composition. Composition is
the confound -- two proteins rich in the same residues score well by chance --
so holding composition fixed while destroying order isolates whether the
observed score reflects shared ancestry or shared composition. The empirical
p-value is the fraction of shuffles scoring at least the real alignment.

Real homology sits far above its shuffled distribution (p ~ 0). A chance
alignment sits inside it. The 8 SUPPORTED pairs are included as positive
controls: if the test does not separate them from the UNRESOLVED ones, the test
itself is uninformative and nothing should be concluded from it.

Also recorded: whether the pair's curated census families CONFLICT. a12 used
shared family only as support, never conflicting family as counter-evidence,
which was a gap -- KCNQ1 is isIonChannel while AMFR is isE3Ligase, and PPM1D is
not a membrane protein at all.

Sequences only. No outcome label is read. Creates new files only.

## scripts/beat_am/a17_general_scorer.py

A17: a GENERAL scorer, applicable to any missense variant with no assay metadata.

Why this exists. The best development model consumes two features a census panel
site cannot supply:

  * `endpoint_type` one-hot -- a panel site belongs to no assay, so it is undefined
  * the local-window embedding mean -- needs neighbours, and the panel stops at
    distance 60, so the window is truncated at the edges

The audit needs scores for all 726,444 panel sites, and the deliverable requires
`score_variants.py` to score *any* missense variant list. Both demand a scorer
with no assay-dependent and no neighbour-dependent inputs. This trains that
scorer and measures what the generality costs, so the price is visible rather
than buried.

Feature set: AM, ESM1b (damage-oriented), the substitution descriptors, the
per-site topology, and a PCA projection of the ProtT5 vector AT the position.
No endpoint, no window.

PCA is fit on outer-training rows only, per locked decision 5.

Scope: fitting pool only. Confirmation proteins are absent from l1_table.npz.

## scripts/beat_am/a18_panel_score_and_correct.py

A18: score the census panel, fit the one-parameter correction, measure its cost.

Section 5 step 6(a), the prescribed cheapest route.

Pipeline:
  1. fit the general scorer (a17 feature set) on the whole fitting pool
  2. score all 726,444 panel sites x 20 alternatives, in site chunks
  3. build the census split; remove confirmation proteins AND their homology-group
     relatives from the TRAIN side, rather than reusing cluster_split unchanged
  4. on census train: audit -> delta0 and penalty SD; probe the sensitivity k
     empirically; fit_beta
  5. apply_correction to the full panel
  6. on HELD-OUT census clusters: audit the corrected scores -- the headline
     rule number -- and run the full 8-rule report card, with the diff vs AM
  7. fidelity(original, corrected), and the accuracy cost on development
     proteins by applying the identical shift to their rows

Unit note. `fit_beta` computes `target = declared * margin * sd`, so the -0.03
audit target is in **SD units** and the raw target is -0.03 * penalty_sd. That
resolves the ambiguity flagged in session 14; it is read from the code, not
assumed.

Label-free throughout on the panel side. Confirmation outcomes are never opened.

## scripts/beat_am/a19_report_card_and_cost.py

A19: the rest of the rule deliverable -- report card, splits, and the accuracy cost.

a18 produced the headline rule number. Four required pieces were still missing:

  1. the full 8-rule report card for the corrected model, on held-out census
     clusters, plus the diff against the SAME model without the correction --
     which is what "every other report-card rule keeps its verdict" means, and
     what shows the placebo rules still read null
  2. the report card for AlphaMissense on the same held-out clusters, for the
     diff the deliverable asks for
  3. the **symbol-family split**, the one Form B missed
  4. the accuracy cost: per-protein Spearman on development proteins before and
     after the identical correction

Side-convention hazard, handled explicitly. `apply_correction` treats panel
`side == 0` as cytosolic. The L1 feature table uses a DIFFERENT encoding from
a9: 1 = inside/cytosolic, 2 = outside/outer, 0 = in-TM or unknown. Applying the
panel rule directly to L1 rows would invert the shift for every cytosolic row.
The mapping is made explicit below and the fraction of changed rows is checked
against the panel's 7.35%.

Label-free on the panel side. Development outcomes only; confirmation outcomes
are never opened.

## scripts/beat_am/a1_inventory_feasibility.py

A1: outcome-blind pool table, exposure register and confirmation feasibility.

Section 5 step 1 of the beat-AlphaMissense prompt. Reads only outcome-blind
inventory metadata (results/route_a_20261012/a0_feasibility.json, which records
"DMS_score not read; mutant column only") plus development-only paired
differences on the four already-opened proteins (ALL_COMPARISONS.csv).

No confirmation-set outcome label is read anywhere in this file.

Writes into results/beat_am_<date>/inventory/. Creates new files only.

### Retained implementation rationale for scripts/beat_am/a1_inventory_feasibility.py

Proteins whose outcomes or scores a previous project analysis has already touched. Accessions, because gene symbols differ between sources (OPSD/RHO, SC6A4/SERT, VKOR1/VKORC1, CP2C9/CYP2C9, TPOR/MPL).

## scripts/beat_am/a21_beta_both_splits.py

A21: fit a single beta that satisfies BOTH census splits, not just the sequence one.

a19 found the correction passes the sequence split (standardized -0.0543, interval
upper -0.0026) and FAILS the symbol-family split (-0.0212, interval upper +0.0003,
so the interval includes zero). The prompt asks for the symbol-family split to be
reported and for a genuine attempt to pass it.

Cause: beta was fitted on the sequence split's TRAIN partition only. The symbol
split starts from a larger violation (+0.1277 standardized vs +0.0938), so the
beta that just clears the sequence target undershoots the symbol one.

Fix, declared before running: fit beta on each split's own train partition, then
take the LARGER of the two. A single beta must satisfy the stricter requirement,
and over-correcting is the safe direction because the accuracy cost measured in
a20 was zero to slightly positive (+0.0001 overall, +0.0017 on the flank).
Both splits are then re-audited on their own held-out clusters with that one beta,
and the accuracy cost is re-measured so the trade is visible.

Label-free. Confirmation outcomes never opened.

### Retained implementation rationale for scripts/beat_am/a21_beta_both_splits.py

Deliberate over-correction on TRAIN only. Declared before seeing whether it passes: fit to 1.5x the real target so held-out clears 0.03 with headroom. Justified because a20 measured the accuracy cost of the correction at zero to slightly POSITIVE (+0.0001 overall, +0.0017 on the flank), so over-correcting is the cheap direction. The real target stays 0.03 for every verdict.

## scripts/beat_am/a22_freeze_manifest.py

A22: the frozen split. Admit KCNQ4, apply the shuffle-test filter, verify no leak.

Three changes, each with its justification stated:

1. **Homology pairs are filtered by the shuffle test, uniformly.** Any pair whose
   local alignment score is indistinguishable from composition-matched chance
   (empirical p > 0.05 over 300 shuffles) is dropped; every pair with p <= 0.05
   is kept. The 8 positive controls all have p = 0.0000 and z from 12.1 to 103.4,
   so the instrument separates signal from chance before it is used to decide
   anything.

   This is applied as a RULE rather than pair by pair on purpose. Dropping the
   PPM1D-KCNQ1 link happens to be convenient -- PPM1D is the one protein in G003
   that was actually fitted in L1 -- and a convenient edit to a split is exactly
   the kind of thing that should not be made ad hoc. The shuffle test was run in
   session 11, before this decision existed, and the filter applies to all 17
   non-artifact pairs equally.

2. **KCNQ4 is admitted to the confirmation set.** Branch (b) recovery qualified
   it on every hard gate (reference agreement 1.000000 on the measured mutant
   column, zero non-finite, 197 flank positions at <=60) and it is the only
   candidate adding genuine two-sided K/R-gain coverage (69 cytosolic / 20 outer).
   It lifts the flank-substantial confirmation set from 6 to 7, taking the
   smallest attainable exact sign-flip p from 0.031 to 0.016.

3. **beta must be refit.** KCNQ4's group retains AMFR (shuffle z = 4.70, kept
   conservatively), and AMFR sits in the census sequence-TRAIN partition that
   beta was fitted on. Admitting KCNQ4 therefore contaminates the existing beta,
   which has to be refitted with AMFR excluded. KCNQ4 itself is in census
   *validation* and KCNQ1 in *test*, so neither was ever in a fitting partition.

Then the gate: no member of any confirmation-containing group may appear in
l1_table.npz, i.e. may have been fitted. Checked exhaustively, not spot-checked.

Writes new files only. No outcome is read.

## scripts/beat_am/a4_group_manifest.py

A4: group manifest and homology overlap audit (Section 5 step 2).

One manifest, one split, for every fitted stage. Groups cover protein, gene and
homologs; all assays of one protein stay together.

Sequences only -- no outcome label is read anywhere in this file. Sequences come
from two metadata sources: the Oxford census embedded in topology_viewer.html
(field `sequence`, 4863 membrane proteins) and ProteinGym's reference table
(field `target_seq`). No assay score file is opened.

Clustering: MMseqs2 is not installed locally, so this falls back to the
Biopython aligner already used for construct mapping in route_a_common.py --
same BLOSUM62 matrix and same gap scores -- but in LOCAL mode rather than
global. Local is the right tool for "do these two share a homologous region",
and over-grouping is the safe direction for leakage. Recorded in the manifest.

Grouping criterion, from the prompt's overlap audit: any pair at >=30% identity
over >=50% coverage joins the same group, by single linkage.

Coverage uses min(len_a, len_b) as the denominator. That inflates coverage
relative to a max() denominator, so more pairs qualify and more merging happens
-- again the safe direction.

Creates new files only.

### Retained implementation rationale for scripts/beat_am/a4_group_manifest.py

---- resolution rule "exclude them AND THEIR RELATIVES from every fitted stage": any group that contains a confirmation protein (membrane or all-human) has its entire membership barred from supervised training, rule fitting, correction coefficients, learned normalisation and hyperparameter selection.

## scripts/beat_am/a6_extract_predictors.py

A6: extract AlphaMissense and ESM1b scores for the manifest proteins.

One streaming pass over AlphaMissense_aa_substitutions.tsv.gz (1.2 GB gzipped)
pulls every needed accession at once, so the file is never scanned twice.
ESM1b comes from per-accession CSVs inside ESM1b_ALL_hum_isoforms_LLR.zip by
random access.

Scope discipline. Confirmation-protein scores are written to a SEPARATE
directory, `features/confirmation_held/`, which no development or fitting
script reads. Extracting a predictor score is not looking at an outcome and is
not a model-versus-AM comparison, but keeping the two sets physically apart
makes an accidental peek take a deliberate act.

Orientation (locked decision 1): AM pathogenicity is already damage-oriented
(higher = more pathogenic). ESM1b LLR is fitness-oriented (alternative minus
wild type), so the damage score written here is NEGATED: damage = -LLR.

AlphaMissense data is (c) 2023 DeepMind Technologies Limited. The file header
declares CC BY-NC-SA 4.0 (note: the project prompt says CC BY 4.0; the file
itself is the more restrictive NC-SA). Not validated for, and not approved for,
any clinical use.

Creates new files only.

## scripts/beat_am/a9_l1_features.py

A9: L1 feature matrix for the development pool.

Builds one row per usable measured single substitution across the fitting-pool
proteins, with:

  predictors            AM pathogenicity, ESM1b damage (both damage-oriented)
  substitution          BLOSUM62, Grantham distance and its three components,
                        hydrophobicity / volume / charge / helix / strand deltas
  topology              side, distance to nearest TM boundary, in-TM flag,
                        TM index, loop length, TM count  (membrane proteins only)
  assay                 endpoint-type indicator, is_membrane
  target                within-assay rank-normalised damage

Coordinates come through route_a_common.align_construct exactly as in a8 -- no
global offset -- and the >=0.99 native-residue gate is applied.

Scope discipline: reads ONLY features/fitting_pool/. Confirmation proteins are
absent from this table entirely.

Creates new files only.

## scripts/beat_am/b10_disulfide_rule.py

B10: a fourth named rule for the audit -- DISULFIDE BOND LOSS.

A structural disulfide covalently links two cysteine side chains. Substituting
either endpoint destroys that linkage. A scorer that has learned disulfide
chemistry should penalise substitutions at bonded cysteines more than at
comparable cysteines that are free under native conditions.

============================================================================
WHY THE WITHIN-POSITION ESTIMATOR OF B8 AND B9 IS NOT USED HERE
============================================================================

B8 (glycosylation sequon) and B9 (charge burial) both measure a WITHIN-POSITION
contrast: at one site, one class of substitution against another class at that
same site. That design is why those audits are credible. It cancels exactly
whatever the model applies to the position as a whole -- a conservation offset,
a domain prior, a per-protein calibration -- because both arms are measured at
the same residue.

That estimator CANNOT be used for this rule, and reusing it would silently
destroy the signal rather than fail loudly.

The reason is that the bond-loss effect has no within-position comparator. At a
bonded cysteine, EVERY ordinary non-Cys replacement removes the linkage. There
is no "motif-preserving" substitution to subtract against, the way ->S/->T
preserves a sequon in B8 or an uncharged alternative preserves the environment
in B9. The within-position mean at a bonded cysteine already contains the
bond-loss penalty in every one of its terms. Subtracting it would subtract the
very quantity being measured, and the audit would report approximately zero for
a scorer that had learned the chemistry perfectly.

No comparator was invented to paper over this. Inventing one -- treating some
substitution as bond-preserving when it is not -- would be worse than using a
weaker design honestly.

So this rule uses a MATCHED-SITE design: the same substitution compared BETWEEN
different cysteine sites, bonded against free, matched within protein and
compartment. That is a genuinely weaker instrument than B8's and B9's, and the
cost is stated rather than hidden:

  - It gives up the exact cancellation of per-position effects.
  - It therefore reintroduces the position-sensitivity confound that the
    within-position contrast was introduced to remove.
  - Bonded cysteines are, on average, more buried and more conserved than free
    ones. A positive result here is consistent with bond loss AND with generic
    burial or conservation, and matching on annotation alone does not separate
    them. See LIMITATIONS and prediction S3.

The second test below (P2, the pairing test) is the part of this script that is
largely immune to that confound, and it may be the more interesting result.

============================================================================
DECLARED PREDICTIONS -- fixed before any scorer output was read
============================================================================

Recorded verbatim in the output JSON under "declared_before_running". Any of
these that fails is reported AS A FAILURE, not reinterpreted.

The eligibility census -- how many bonded and free cysteines exist, their
evidence codes and their compartments -- WAS computed before these predictions
were written, because the predictions have to be answerable by the sample that
exists. That census reads UniProt feature annotation only. No AlphaMissense or
ESM1b value was read until after the predictions below were fixed in this file.
One prediction (S5) is declared as unavailable on the basis of that census, and
the counts that make it unavailable are reported.

  P1  PRIMARY, magnitude. score(C->S) is higher -- more damaging -- at bonded
      cysteines than at matched free cysteines, with a 95% interval excluding
      zero. Matched = same protein, same compartment.

  P2  PRIMARY, pairing. If a scorer represents SPECIFIC bonds rather than
      "cysteines are important", the two endpoints of the same bond should
      receive more similar C->S penalties than two cysteines of the same
      protein paired at random. So the mean within-bond |difference| is LOWER
      than the mean under within-protein permuted re-pairing.

      P2 is the sharper test. It is computed from exactly the same set of
      scores as its own null -- the permutation only changes WHICH cysteine is
      paired with which -- so a protein-level, position-level or burial-level
      difference in score magnitude cannot produce it.

  S1  SECONDARY, substitution sensitivity. C->A shows the same direction as
      P1. Run separately and reported separately; C->S and C->A are NOT
      averaged, because they are not interchangeable probes (alanine removes
      the side-chain hydroxyl as well as the thiol).

  S2  SECONDARY, evidence quality. P1 survives restricting the bonded set to
      endpoints with direct experimental evidence (ECO:0000269, or
      ECO:0007744 with a PDB source). Experimental and inferred annotation are
      never pooled silently; the breakdown is reported.

  S3  SECONDARY, burial-proxy matching. P1 survives adding a local
      hydrophobicity stratum to the match, so that bonded and free cysteines
      are compared only within the same protein, the same compartment AND the
      same quartile of local sequence hydrophobicity. Local hydrophobicity is a
      sequence-only proxy for burial. If P1 disappears under this match, the
      primary result is substantially a burial artefact and must be read that
      way.

  S4  SECONDARY, free-set quality. P1 survives restricting free controls to
      proteins carrying at least two annotated disulfide bonds -- proteins
      whose disulfide pattern a curator has clearly characterised, so that an
      unannotated cysteine there is better supported as genuinely free.

  S5  CONTROL, compartment placebo. Structural disulfides form in the
      oxidising ER lumen and extracellular space, not in the reducing cytosol.
      A cytosolic arm would be a natural placebo, as the cytosolic sequons are
      in B8. DECLARED UNAVAILABLE: the eligibility census finds 16,517
      non-cytosolic bonded cysteines against 16 cytosolic ones in the census
      proteins, so no cytosolic placebo can be run. The scarcity is itself
      consistent with the chemistry, but it is an annotation count and not a
      test of any scorer. Reported as NOT AVAILABLE, not as a pass.

  X1  CROSS-SCORER. AlphaMissense and ESM1b are reported separately. Only
      DIRECTION and the pass/fail pattern are compared between them. Their
      scales and spreads differ, so no standardised (SD-unit) quantity is
      compared across scorers anywhere in this file.

============================================================================
DESIGN
============================================================================

1. PRIMARY STATISTIC. The RAW damage score for C->S at each eligible cysteine.
   Not a contrast. C->A is run separately as S1.

2. ELIGIBILITY. Cysteines in the Oxford topology census, restricted to
   accessions whose census sequence is byte-identical to the UniProt sequence,
   so that annotation positions and scored positions cannot disagree. Every
   annotated disulfide position is checked to be a cysteine in that sequence
   and mismatches are counted, not dropped silently.

   BONDED: carries a DISULFID feature with a defined intrachain partner.
   Interchain bonds are excluded -- their partner is in another molecule, so a
   whole-protein resample would not keep the pair together. Endpoints with
   conflicting alternative partners, and singleton features with no partner,
   are excluded as uncertain.

   FREE: carries NO disulfide annotation. This is the hard half and the
   assumption is stated plainly: ABSENCE OF A DISULFIDE ANNOTATION IS AN
   UNKNOWN, NOT EVIDENCE THAT A CYSTEINE IS FREE. The annotation source
   supplies no usable set of positively evidenced free thiols -- the count of
   sites that could be supported that way is reported, and it is approximately
   zero. An independent earlier pass over the same source reached the same
   conclusion (results/disulfide_rule_20261025/inventory_summary.json,
   "verified_free_controls_in_import": 0).

   One thing does partly mitigate this, and it is a property of the matched
   design rather than a claim about the annotation: because every stratum
   requires a bonded arm and a free arm in the SAME protein, free controls are
   drawn almost entirely from proteins whose disulfides a curator did annotate.
   An unannotated cysteine in such a protein is better supported as genuinely
   free than an unannotated cysteine in a protein nobody characterised.
   Prediction S4 pushes that further.

   EXCLUDED FROM BOTH ARMS: cysteines carrying a binding-site, active-site or
   modified-residue annotation, whose role is something other than structural,
   and cysteines inside a signal peptide, which is cleaved and cannot supply a
   mature-protein control. Following the control-pool rule used elsewhere in
   this project, these positions are removed from EVERY arm, not only from the
   arm being scored.

3. COMPARTMENT. Taken from the census topology: a loop with location "inside"
   is cytosolic, "outside" is non-cytosolic, and positions inside an annotated
   TM segment form their own compartment. Compartment is a stratum variable, so
   a bonded cysteine is only ever compared against free cysteines on the same
   side of the membrane.

4. AGGREGATION. Algorithm 1 as B8 and B9 implement it: the bonded-minus-free
   difference is formed within protein x compartment strata, strata are weighted
   by effective sample size w = n1*n0/(n1+n0), and the bootstrap resamples
   WHOLE PROTEINS -- 1999 draws, seed 20261016. All quantities in one draw come
   from the same protein resample.

5. BOND PARTNERS STAY DEPENDENT. The two cysteines of one bond are not
   independent observations. Every retained bond is intrachain, so both of its
   endpoints belong to the same protein, and a whole-protein resample therefore
   carries both endpoints together automatically -- either both are in a draw
   or neither is. That is why interchain bonds are excluded rather than kept:
   their endpoints sit in different molecules and no protein-level resample
   could preserve their dependence. Bonds with exactly one eligible endpoint
   (the partner excluded as a binding site, signal peptide and so on) are
   counted and reported as orphaned endpoints; the eligible endpoint is kept in
   the bonded arm, since the rule is about that residue's own linkage.

6. THE PAIRING TEST (P2). For each protein with at least two bonds whose BOTH
   endpoints are eligible and scored, the observed statistic is the mean over
   bonds of |score(C->S, p) - score(C->S, q)|. The null re-pairs the same
   endpoints within the same protein: the 2k endpoint scores are shuffled and
   paired consecutively, which draws a uniformly random perfect matching. 1999
   permutation draws, same seed.

   The null uses the identical multiset of scores, so it holds every level
   effect fixed -- protein, compartment, burial, conservation -- and varies only
   the pairing. A protein with exactly two bonds has three possible matchings,
   one of them the true one, so the null contains the observed pairing with
   positive probability; that makes the test conservative rather than liberal.

============================================================================
LABEL-FREE
============================================================================

Inputs are the Oxford topology census (sequence, TM boundaries, loop sides),
UniProt sequence and feature annotation, and the scorers' own outputs. No
ClinVar, no humsavar, no gnomAD, no HGMD, no deep mutational scan, no outcome
label of any kind is read. UniProt feature annotation is sequence and structure
annotation, not clinical data.

The UniProt features are read from the response already cached on disk at
data/external/disulfide_rule_20261025/human_reviewed_features.tsv, fetched by
c1_disulfide_sources.py from the REST API with its URL and SHA-256 recorded in
that directory's source_manifest.json. Nothing is re-fetched on a rerun.

============================================================================
LIMITATIONS -- also written into the output JSON
============================================================================

1. THE BURIAL AND CONSERVATION CONFOUND. This is the main one. Bonded
   cysteines are on average more buried and more conserved than free ones, so a
   positive P1 is consistent with the scorer having learned bond loss AND with
   it having learned only that these positions matter. Matching within protein
   and compartment helps and is NOT sufficient. S3 bounds the burial part with
   a sequence-hydrophobicity proxy. Using annotation alone, the matched-site
   design CANNOT fully separate bond loss from burial, and a positive P1 on its
   own should not be read as evidence that a scorer represents disulfide
   chemistry. P2 is the test that does not have this problem.

2. THE FREE SET IS AN UNKNOWN SET, not a verified-free set. See eligibility
   above. Any genuinely bonded cysteine sitting in the free arm because nobody
   annotated it biases P1 TOWARD zero, so this cost works against the declared
   direction rather than manufacturing it.

3. NO STRUCTURE WAS READ. No solvent accessibility, no PDB coordinates, no
   predicted structure. Burial enters only through the hydrophobicity proxy in
   S3, which is weak.

4. EVIDENCE MIXTURE. Most annotated endpoints are inferred from family rules
   rather than observed directly. The breakdown is reported and S2 is the
   experimental-only arm.

5. NO CYTOSOLIC PLACEBO. See S5. This rule lacks the built-in placebo that
   makes B8's sequon audit self-checking.

This script MEASURES. It fits nothing, trains nothing and tunes no correction
parameter. Writes only into results/beat_am_20261016/disulfide_rule/ and
creates no file anywhere else.

AlphaMissense data (c) 2023 DeepMind Technologies Limited, CC BY-NC-SA 4.0 per
the distributed file header. Not validated for, and not approved for, clinical use.

### Retained implementation rationale for scripts/beat_am/b10_disulfide_rule.py

--------------------------------------------------------------------------- Score extraction -- one streaming AlphaMissense pass, random access for ESM1b. ---------------------------------------------------------------------------

## scripts/beat_am/b10_refinement_cost.py

B10: what the refinement layer costs computationally, across the full panel.

B9 showed the "my scorer is faster than AlphaMissense" framing does not survive
its own accounting: the predictor-free scorer sits on ProtT5-XL at 3B
parameters, and both baselines ship as precomputed proteome-wide files, so
there is no inference time to beat.

Under the refinement framing the useful question is different and the tool is
the subject: what does it cost to AUDIT a scorer for a named rule, and to
ENFORCE that rule on its output? Those are the two operations a user of this
tool actually performs, and neither has ever been timed.

Four stages, measured separately so no single number hides the expensive one:

  1. load the panel                      -- I/O, paid once
  2. audit                               -- does this scorer obey the rule?
  3. enforce                             -- apply the correction to every cell
  4. fidelity                            -- how much did enforcement disturb?

The audit is timed at several bootstrap depths, because the resampling
dominates it and the honest number depends on how many draws a user wants.

Nothing here reads an outcome label or fits anything; it re-runs operations
already performed elsewhere in the project purely to time them.

Writes only into results/beat_am_20261016/speed_benchmark/.

## scripts/beat_am/b10d_pairing_domain_control.py

B10D: the DOMAIN-matched control for B10's disulfide pairing test (P2).

B10C established that B10's claim -- that its permutation null holds "burial and
conservation ... all held fixed" so that "a level confound cannot produce" P2 -- is
wrong, and that matching control pairs on a constraint proxy leaves only ~21% of the
AlphaMissense P2 gap and none of the ESM1b gap.

That control had one structural weakness, stated in B10C's own output: the constraint
proxy is derived from the audited scorer's OWN outputs, and at a bonded cysteine every
substitution removes the linkage, so the proxy still contains bond-loss signal.
Matching on it therefore over-corrects, which makes 21% a LOWER bound rather than an
estimate. B10's separation-matched control leaves ~82%, an upper bound, because
sequence separation is a weak proxy for co-location.

THIS SCRIPT REPLACES THE PROXY WITH ANNOTATION. The confound is that the two cysteines
of one disulfide usually sit in the same domain and therefore the same structural
environment, so their scores are correlated for reasons having nothing to do with the
bond. UniProt DOMAIN annotation states domain membership directly, is independent of
both scorers, and contains no bond-loss signal. So it can bound the effect from neither
side -- it simply measures it.

----------------------------------------------------------------------------
THE DESIGN
----------------------------------------------------------------------------

A pair of cysteines is INTRA-DOMAIN when both positions fall inside the same annotated
DOMAIN interval of that protein. Control pairs are drawn, as in B10 and B10C, from the
protein's own BONDED endpoints and never from a true bond, so the marginal distribution
of scores is held fixed exactly as P2's permutation does.

  PRIMARY. Restrict to intra-domain true bonds. For each one, draw a control pair from
  the same protein's intra-domain non-bond pairs, PREFERRING a pair inside the same
  domain instance as the true bond, falling back to any intra-domain pair in that
  protein. Among eligible candidates the one with the closest sequence separation is
  taken, so proximity is balanced as a side effect rather than left free; ties break by
  position order and each control is used once where possible. The realised separations
  and the share of same-instance matches are both reported, so a poor match is visible.

  SENSITIVITY. The same comparison with no separation preference at all (first by
  position order), to show the tie-break rule does not drive the result.

  CROSS-DOMAIN ARM. Declared NOT AVAILABLE before running. Only 26 of 7688 true bonds
  have endpoints in different annotated domains, against 5500 intra-domain. That is far
  too thin to carry a comparison and no number is computed from it. This is reported as
  a structural limit of the rule, not worked around.

Because the primary runs on a subset, the P2 gap it is compared against is recomputed
on that same subset -- observed and within-protein permutation null over the
intra-domain bonds only -- so the "share surviving" figure has a matching denominator.

----------------------------------------------------------------------------
DECLARED PREDICTIONS -- fixed before any score was read
----------------------------------------------------------------------------

  E1  PREMISE. True bond partners are far more likely to share an annotated domain than
      randomly re-paired bonded cysteines of the same protein. Confound-friendly and
      expected to hold. Established from annotation alone, with no scorer involved,
      which is the point of this control.

  E2  THE TEST. Among intra-domain pairs only, true bonds still show a LOWER mean
      |difference| in C->S score than intra-domain non-partner controls, with a 95%
      interval excluding zero.

  E3  The share of the subset P2 gap surviving domain matching falls BETWEEN B10C's
      constraint-matched 21% and B10's separation-matched 82% for AlphaMissense.
      Reasoning: domain membership is a better co-location proxy than separation, so it
      should remove more than 18%; and unlike the constraint proxy it cannot
      over-correct, so it should remove less than 79%.

  E4  No prediction is declared for ESM1b. B10C found its P2 fully explained by shared
      constraint; whether annotation reproduces that is reported either way.

  Misses are reported as misses and are not reinterpreted.

----------------------------------------------------------------------------

Label-free: topology census, the UniProt feature response B10 had already cached, and
the scorers' own cached outputs. No ClinVar, humsavar, gnomAD, HGMD or DMS data read.
No network request, no new AlphaMissense pass. Nothing fitted, trained or tuned.
Writes only into results/beat_am_20261016/conservation_confound/.

AlphaMissense data (c) 2023 DeepMind Technologies Limited, CC BY-NC-SA 4.0 per the
distributed file header. Not validated for, and not approved for, clinical use.

## scripts/beat_am/b1_predictor_free.py

B1: a scorer that consumes NO audited predictor, and enforcement on it.

Why this exists. Every enforcement result in the paper so far is measured on a
model whose input #1 is AlphaMissense and whose input #2 is ESM1b. A reader can
reasonably object that the correction is patching the audited predictor's own
output rather than enforcing a rule on an independent model. This builds the
scorer that removes the objection.

Declared before running:

  * **Both** audited predictors are dropped, not just AlphaMissense. ESM1b
    already satisfies the rule on its own (-0.073 SD), so leaving it in would
    let compliance arrive through the frontend. Dropping both supports the
    stronger claim: a scorer containing no audited predictor.
  * Everything else is held identical to a17/a18/a21 -- same rows, same groups,
    same grid, same seed, same CV protocol, same beta-fitting rule (fit on each
    split's own train partition, take the larger, 1.5x margin), same decay
    length, same 0.03 verdict target. The only change is the feature set.
  * ProtT5 embeddings stay. They are a language-model representation and
    section 3.6(b) showed such representations already carry membrane topology,
    so this scorer may well show a contrast of its own. That is the point: if it
    does, enforcement on it is a genuine test; if it does not, the violation was
    attributable to the audited predictors, which is itself a result.

Outcome is not predicted here. Either direction is reported as found.

Scope: fitting pool only. The 7 reserved membrane proteins, the 24 reserved
non-membrane proteins and the 6 pending candidates are absent from
`l1_table.npz` and are never read. No confirmation outcome is opened. Accuracy
figures below are development-set figures and are labelled as such.

Writes only into results/beat_am_20261016/predictor_free/.

### Retained implementation rationale for scripts/beat_am/b1_predictor_free.py

---- 5. accuracy cost on development rows --------------------------- L1 encodes 1 = cytosolic, 2 = outer; the panel uses 0 for cytosolic. Applying the panel rule here would invert every cytosolic shift (a20).

## scripts/beat_am/b8_sequon_rule.py

B8: a second named rule for the audit -- the N-glycosylation sequon.

The paper currently enforces one rule. A refinement tool is judged by how many
named rules it can audit and enforce, so this adds a second, built to the same
design as Algorithm 1 and declared before any contrast was computed.

THE RULE. N-linked glycosylation requires the motif N-X-S/T with X != P. The
asparagine is glycosylated only when that motif is intact, and only on the
non-cytosolic face, because the oligosaccharyltransferase active site faces the
ER lumen. A scorer that has learned this should treat motif-breaking
substitutions differently from motif-preserving ones AT THE SAME POSITION -- and
should do so on the non-cytosolic side only.

TWO CONTRASTS, both within-position, both declared with a direction.

  A. The X position. Substituting proline destroys the motif; other residues do
     not. Proline is disruptive everywhere, so the raw penalty is uninformative;
     what matters is whether the ->P penalty is LARGER at sequon X positions
     than at matched non-sequon positions of the same native residue.
         z_i = score(i->P) - mean(score(i->other non-P))
         prediction: Delta = z(sequon) - z(matched) > 0 on the non-cytosolic side

  B. The S/T position. S and T are interchangeable here -- either sustains the
     motif -- while anything else destroys it. S<->T is chemically conservative
     everywhere, so again only the DIFFERENCE against matched non-sequon S/T
     positions is informative.
         z_i = score(i->partner) - mean(score(i->other))
         prediction: Delta = z(sequon) - z(matched) < 0 on the non-cytosolic side

THE PLACEBO IS BUILT IN. Cytosolic sequons are never glycosylated: the motif is
there, the chemistry is not. A scorer that has learned glycosylation should read
near zero on the cytosolic side. A scorer that has merely learned "N-X-S/T looks
conserved" will show the same effect on both sides. That control costs nothing
and no label-based version of this test can have it.

WHY LABEL-FREE. A version of this rule built on ClinVar or humsavar would face
two problems this design does not have. Those labels increasingly cite
computational predictions under ACMG PP3/BP4, so auditing AlphaMissense against
them is partly circular. And restricting to labelled variants shrinks the sample
from ~18,000 sequons to ~100. Nothing here reads an outcome label.

Structure follows Algorithm 1 exactly: within-position contrast, then
comparison within protein x native-residue strata (removes composition
imbalance), then a whole-protein bootstrap (preserves within-protein
dependence).

On what the within-position contrast removes: it cancels anything the model
ADDS uniformly to every substitution at a position. It does not cancel a purely
multiplicative position effect, so a difference of scores at a highly conserved
position can still carry some of that conservation. This limit is inherited
from Algorithm 1 rather than introduced here, and the matched controls -- same
protein, same side, same native residue -- are what bound it.

This measures. It fits nothing. The enforcement step -- fitting a correction
parameter against a declared target -- is deliberately not here.

Writes only into results/beat_am_20261016/sequon_rule/.

### Retained implementation rationale for scripts/beat_am/b8_sequon_rule.py

Every position belonging to ANY N-X-S/T motif, whether or not that sequon has a known side. A sequon residue is never a neutral control. The first version excluded only X and S/T, which let sequon N positions serve as controls for the 3.7% of sequons whose X is itself an N. That affected contrast A only; contrast B draws its controls from S/T positions and was never exposed to it.

## scripts/beat_am/b9_charge_burial.py

B9: a third named rule for the audit -- CHARGE BURIAL DEPTH.

This one is deliberately a POSITIVE CONTROL. The first two rules in this audit
(positive-inside, B8's glycosylation sequon) were chosen because a scorer might
plausibly fail them. An instrument that only ever reports violations cannot be
trusted: it has never been shown to say "pass" when the physics is obeyed. This
rule is one that any scorer with real membrane physics in it should obey, and
obey strongly. If both AlphaMissense and ESM1b show it, the instrument detects
compliance as well as violation. If they diverge, that is itself a finding.

THE RULE. Moving a charged side chain (D, E, K, R) out of water and into the
hydrocarbon core of a lipid bilayer costs free energy, because the core has no
water to solvate the charge. The cost is not flat across the bilayer: it rises
toward the centre, where the dielectric constant is lowest and water penetration
is least. The biological-hydrophobicity scales (Hessa/von Heijne) measure this
directly -- the per-residue insertion penalty for K, R, D and E is near zero at
the lipid-water interface and largest at the helix midpoint. So a scorer that
has learned membrane physics should penalise charge-INTRODUCING substitutions
more at deep positions in a transmembrane helix than at shallow ones.

Note what this does and does not assume. It does not assume the scorer knows
about membranes explicitly. It assumes only that if the scorer's outputs track
the physics, the depth gradient must be visible in them. That is the whole
design of this audit: the rule is named, the direction is declared, and the
scorer's own outputs are the only measurement.

----------------------------------------------------------------------------
DECLARED PREDICTIONS -- written before any score was read
----------------------------------------------------------------------------

Recorded verbatim in the output JSON under "declared_before_running". Any of
these that fails is reported AS A FAILURE, not reinterpreted.

  P1  PRIMARY, dose-response. Mean z rises monotonically across the five depth
      bins, from the lipid-water interface (bin 0) to the helix centre (bin 4).

  P2  PRIMARY, endpoint gap. z at the centre bin exceeds z at the interface
      bin, and their 95% intervals do not overlap. The paired, stratified form
      of the same statement -- Delta(bin 4 - bin 0) > 0 with a 95% interval
      excluding zero -- is also required.

  P3  POSITIVE CONTROL. BOTH AlphaMissense and ESM1b satisfy P1 and P2. A
      scorer that fails them has not learned the depth dependence of charge
      burial. Divergence between the two scorers is a reportable finding.

  S1  SECONDARY, loop anchor. z at the helix centre greatly exceeds z at loop
      positions outside the membrane. Secondary because loops differ from TM
      interiors in many ways at once; the within-helix gradient is the cleaner
      comparison, since it holds the structural environment roughly constant
      and varies only depth.

  S2  SECONDARY, charge sign. Both the K/R split and the D/E split show the
      depth gradient. No prediction is made about which is larger.

  S3  SECONDARY, positive-inside in loops. Among loop positions, introducing
      K or R is penalised LESS on the cytosolic side than on the non-cytosolic
      side. This is the project's headline rule, measured here for free on a
      pool this script already has to score.

  S4  SECONDARY, scale control. The depth gradient survives dividing each
      position's z by that same position's own spread across all 19 substitutions
      (see the limitation on multiplicative effects below). If the gradient
      disappears under that normalisation, the raw gradient is substantially a
      conservation-scale artefact and P1/P2 should be read with that in mind.

----------------------------------------------------------------------------
DESIGN
----------------------------------------------------------------------------

1. ELIGIBLE SITES. Positions inside an annotated transmembrane segment of
   length >= 10 whose native residue is one of the fifteen uncharged residues
   ACFGILMNPQSTVWY. Charged natives are excluded because the contrast is about
   INTRODUCING charge. Histidine is excluded from both arms throughout: its
   protonation state in a low-dielectric environment is ambiguous, so it is
   neither a clean "charged" nor a clean "uncharged" alternative.

2. WITHIN-POSITION CONTRAST. At each eligible position i,

       z_i = mean(score(i->D), score(i->E), score(i->K), score(i->R))
             - mean(score(i->u) for u in the uncharged residues, u != native)

   Both arms are measured at the SAME position, so anything the model adds
   uniformly to every substitution at that position -- a conservation offset, a
   domain-level prior, a per-protein calibration -- cancels exactly.

   Two further splits, same reference arm:
       z_KR = mean(score(i->K), score(i->R)) - (same uncharged mean)
       z_DE = mean(score(i->D), score(i->E)) - (same uncharged mean)

   WHAT THIS DOES NOT CANCEL. The contrast removes an ADDITIVE per-position
   effect. It does not remove a MULTIPLICATIVE one. If the scorer's whole
   output range is simply stretched at conserved positions, and helix centres
   are more conserved than helix ends, then z would rise with depth even if the
   scorer knew nothing about charge burial. This limitation is inherited from
   Algorithm 1 and is real here. Prediction S4 is the bound on it: the same
   curve computed on z divided by the position's own spread across all 19
   substitutions. That normalisation removes a purely multiplicative scale.

3. FRACTIONAL DEPTH. For a TM segment spanning residues s..e of length L,

       d = 1 - |(p - s) - (L-1)/2| / ((L-1)/2)

   so d = 0 at either lipid-water interface and d = 1 at the helix centre.
   Binned into five equal-width bins on d: bin 0 = interface, bin 4 = centre.
   Depth is a position's distance from the membrane midplane expressed in
   fractions of a half-helix; it is a proxy for burial depth, not a measured
   z-coordinate, and a tilted or kinked helix will be mis-assigned somewhat.
   That mis-assignment is noise, and noise flattens a gradient rather than
   manufacturing one.

4. STRATIFICATION AND AGGREGATION -- Algorithm 1, exactly as B8 implements it.
   Every comparison between two groups is made within protein x native-residue
   strata, so a bin is never compared against a different mix of proteins or a
   different mix of amino acids. Strata are weighted by the effective sample
   size w = n1*n0/(n1+n0), and the bootstrap resamples WHOLE PROTEINS (1999
   draws, seed 20261016), which preserves within-protein dependence.

   This matters more here than it might look. Helix centres are richer in
   L/I/V/A and helix ends in W/Y/S/T/G/P, and z depends strongly on the native
   residue. An unadjusted bin-to-bin comparison would partly be measuring that
   composition shift. So the curve is reported twice:

     - RAW per-bin mean z, the plain and directly interpretable number, with a
       whole-protein bootstrap interval. Per-bin native-residue composition is
       written alongside so the reader can see the imbalance.
     - ADJUSTED per-bin Delta against the interface bin, which IS the Algorithm
       1 contrast: same protein, same native residue, differing only in depth.
       Delta(bin 0) is zero by construction. This is the composition-controlled
       curve and it is the one P1 and P2 are judged on.

   All quantities in a given bootstrap draw are computed from the SAME protein
   resample, so the five bins, the monotonicity statistic and the endpoint gap
   are mutually consistent within a draw.

5. REFERENCE ARM. Loop positions -- outside EVERY annotated TM segment, with a
   known side -- act as a shallow anchor. Following the control-pool rule used
   elsewhere in this project, positions in any TM segment are excluded from the
   loop pool regardless of that segment's length, including the six segments
   shorter than 10 residues that are not themselves eligible sites. Loop
   positions are capped at 10 per protein x native residue, drawn with the
   seeded RNG, which is B8's MAX_CONTROLS convention; the cap is unbiased for a
   stratum mean and keeps the scored pool near a million positions instead of
   2.1 million. Loop residues immediately adjacent to a TM boundary are kept,
   which makes this anchor conservative rather than generous.

----------------------------------------------------------------------------
LABEL-FREE
----------------------------------------------------------------------------

Inputs are the Oxford topology census (sequence and segment boundaries) and the
scorers' own outputs. No ClinVar, no humsavar, no gnomAD, no deep mutational
scan, no outcome label of any kind is read. The rule being tested is a physical
one, and compliance with it is measurable without ever knowing whether any
particular variant is pathogenic.

----------------------------------------------------------------------------
CROSS-SCORER COMPARISON
----------------------------------------------------------------------------

AlphaMissense pathogenicity and ESM1b damage (= -LLR) live on different scales
with different spreads. Standardised values -- anything expressed in SD units --
are therefore NOT comparable between the two scorers, and no such comparison is
made or implied. Standardised values appear in the output only to make effect
sizes within one scorer readable. The cross-scorer claim this script does make
is purely about DIRECTION and MONOTONICITY, which are scale-free.

ESM1b orientation: the distributed LLR is fitness-oriented (lower = more
damaging), so it is NEGATED to damage = -LLR before any contrast, matching
locked decision 1 and a6_extract_predictors.py. The flank cache at
data/processed/gate0/esm1b_flank_scores.npz is deliberately NOT used: it covers
only positions within 60 residues of a membrane boundary and would miss TM
interiors.

----------------------------------------------------------------------------

This script MEASURES. It fits nothing, trains nothing and tunes no correction
parameter. Writes only into results/beat_am_20261016/charge_burial/ and creates
no file anywhere else.

AlphaMissense data (c) 2023 DeepMind Technologies Limited, CC BY-NC-SA 4.0 per
the distributed file header. Not validated for, and not approved for, clinical use.

## scripts/beat_am/b9_speed_benchmark.py

B9: honest speed accounting for the predictor-free scorer vs AM and ESM1b.

The claim being tested is the grad student's suggestion: a smaller model that
is slightly less accurate but much faster is a publishable result. This measures
what can be measured here and states plainly what cannot.

THE ACCOUNTING PROBLEM, STATED FIRST. The predictor-free scorer is not a small
model. It is a gradient-boosted tree head on top of ProtT5-XL, a 3-billion
parameter transformer. Any "it is faster" claim depends entirely on whether the
embedding forward pass is counted:

  * amortised (embeddings precomputed): the tree head alone is trivially fast
  * end to end (novel sequence): a 3B forward pass dominates everything, and is
    LARGER than ESM1b's 650M

And the comparison is further confounded because AlphaMissense and ESM1b are
distributed as precomputed proteome-wide score files. For any human protein, a
user's "inference cost" for either is a file lookup, not a forward pass. Speed
only becomes a real axis for sequences outside those files.

So this benchmarks the three stages separately rather than reporting one
misleading end-to-end number:

  1. tree inference, given features          (the amortised case)
  2. tabular feature construction            (chemistry + geometry, no embedding)
  3. PCA projection, 1024 -> 64
  4. artifact sizes

The ProtT5 forward pass is NOT benchmarked here -- this session has no GPU. The
project's own prior measurement is cited instead.

Writes only into results/beat_am_20261016/speed_benchmark/.

## scripts/beat_am/c1_disulfide_sources.py

Fetch public annotation/identifier data for the frozen disulfide protocol.

No variant-effect score or clinical-label table is requested. Cache responses
atomically with hashes and source URLs. Run inventory before any score audit.

## scripts/beat_am/figures_f4p_trial.py

Trial: the STS data figures redrawn in the figures4papers house style.

A look-and-feel comparison only, not for the paper. Same saved inputs and loaders as
figures_sts.py; only the styling follows ChenLiu-1996/figures4papers
(scientific-figure-making/references/design-theory.md and api.md):

  - Arial/Helvetica sans-serif, font 24 for large bar panels, 16 for compact plots
  - ultra-wide canvases (width 3-4x height), shrunk to fit the page afterwards
  - bars with black edges, values printed above them, no x ticks (legend names them)
  - blue for the proposed method, red/pink for comparators, neutral gray for baselines
  - y-limits tightened toward min - std, as the guide recommends
  - dedicated legend panel, tight_layout(pad=2), 300 dpi

Run from the project root:

    python scripts/beat_am/figures_f4p_trial.py

Writes PNG and PDF into figures/sts_20261031/trial_figures4papers/.

## scripts/beat_am/figures_paper.py

Paper figures for PAPER_DRAFT_IMRAD_20261020.md.

Four figures, one per numbered result. Every value is read from a saved output
file; nothing is retyped from a draft. Run from the project root:

    python scripts/beat_am/figures_paper.py

Outputs PNG (viewing) and SVG (typesetting) into figures/paper_20261020/.
Creates new files only; overwrites nothing outside its own output folder.

## scripts/beat_am/figures_sts.py

STS paper figures, redesign (figure guide v2, 31 Oct 2026).

Every plotted value is read from a saved output file; see
figures/sts_20261031/DATA_MAP.md for the file behind each panel. Each number the
paper also prints is checked against the paper value, and a figure with any
mismatch is not saved. Run from the project root:

    python scripts/beat_am/figures_sts.py

Writes PDF (paper) and PNG (slides) plus one VALUES_<name>.txt per figure into
figures/sts_20261031/. Creates new files only.

The comparator forest (Fig 3) takes all five comparator rows from one saved run,
results/comparator_sensitivity_20261001; the paper's earlier alternative ranges are
not in any saved output. Not built here: the structure case study (built separately
by figures_casestudy.py with PyMOL) and the optional distance figure (no matching
binned outputs).

### Retained implementation rationale for scripts/beat_am/figures_sts.py

(3) rule corrections, grouped by how each coefficient is set: PI and disulfide fitted to DMS residuals; glycosylation set by a label-free score constraint; charge burial available but off in the final model (gray)

## scripts/beat_am/instance/b1_mmseqs_manifest.py

MMseqs2 grouping and overlap audit, compared against the Biopython fallback.

Section 5 step 2 as the prompt specifies it, which is TWO stages with DIFFERENT
thresholds, and they must not be conflated:

  1. Cluster with MMseqs2: `--min-seq-id 0.3 -c 0.8 --cov-mode 0`
  2. Overlap audit: an all-against-all search, where "any pair at or above 30%
     identity with at least 50% coverage joins the same group"

Stage 2 is the more permissive of the two (0.5 coverage vs 0.8), so it can only
add merges on top of stage 1. The local Biopython run (a4) applied stage 2's
criterion alone, so its groups should already be a superset of stage 1 and
close to stage 1+2 combined. This script tests that rather than assuming it.

Sequences come from a4's own `build_universe`, imported directly, so the input
is byte-identical to what the local run clustered. No outcome label is read.

Writes new files only, under results/beat_am_20261016/inventory_mmseqs/.

## scripts/beat_am/instance/c1_embed_pool.py

Per-residue embeddings for every measured site in the fitting pool.

The cached panel embeddings cover flank sites only (726,444 sites at census
distance <=60), which reach just 9.1% of the L1 table's rows and none of the
83,820 non-membrane rows. This embeds the sites that actually carry measured
outcomes.

Recipe matched EXACTLY to the cached panel so new and old vectors are
interchangeable (`scripts/route_b_embed_modal.py`):
  Rostlab/prot_t5_xl_uniref50, T5Tokenizer(do_lower_case=False, legacy=True),
  T5EncoderModel in float16, residues space-separated, U/Z/O/B -> X per the
  model card, add_special_tokens=True, and row r of the output is residue r
  because the tokenizer appends </s> AFTER the last residue.
Chunking also matches: window 1022, overlap 128.

Frontend choice. ProtT5 is the default here and ESM2 is opt-in behind
--model esm2, because of the project's own documented reasoning in
route_b_embed_modal.py: ESM1b already satisfies the K/R-gain rule
(comparator-robust at -0.073), so building on ESM representations risks rule
compliance arriving from the frontend rather than from the rule loss term, which
would make the central claim unfalsifiable. ESM2 is the same lineage and the
same risk. ProtT5 is a different lineage.

Coordinates. Reference sequences are rebuilt from the extracted AlphaMissense
wild-type residues, exactly as a8/a9 did, so site positions line up with the L1
table with no offset.

Sharding and resume. --shard k --n-shards n takes every n-th accession, so one
process per GPU pinned with CUDA_VISIBLE_DEVICES. One output file per accession,
skipped if present, so a rerun only does what is missing.

Label-free: reads positions, never outcomes.

## scripts/beat_am/instance/c5_embed_fasta.py

Embed every residue of every sequence in a FASTA with ProtT5-XL.

Standalone: needs only the FASTA and the model, not the project bundle. Used to
embed the confirmation proteins without re-uploading 4.13 GB.

Recipe matched exactly to route_b_embed_modal.py and c1_embed_pool.py, so these
vectors are interchangeable with the cached panel (validated at cosine 1.0000):
Rostlab/prot_t5_xl_uniref50, T5Tokenizer(do_lower_case=False, legacy=True),
T5EncoderModel in float16, residues space-separated, U/Z/O/B -> X, row r is
residue r because the tokenizer appends </s> after the last residue.
Chunking: window 1022, overlap 128; each position takes the chunk giving it the
most context.

--shard k --n-shards n splits sequences across GPUs. One npz per accession,
skipped if present, so reruns only do what is missing.

## scripts/beat_am/m10_disulfide_internal_consistency.py

Circular annotation-consistency evaluation, explicitly not biological validation.

Use the existing standalone scores and fitted coefficients. The annotation
that selects a positive score offset also defines the evaluation target.
No coefficient, model, target margin or ranking threshold is fitted here.

## scripts/beat_am/m10_fetch_external_baselines.py

Fetch only the existing development assays' public ProteinGym predictions.

Uses HTTP byte ranges instead of downloading the entire 1.9 GB archive. No
models are trained and no reserved assay is requested. Original CSV hashes and
archive CRCs are recorded before retaining the requested score columns.

## scripts/beat_am/m11_compare_external_baselines.py

Matched, no-refit external comparison of the existing standalone predictions.

Primary: guarded group-out predictions, 13 membrane proteins, intersection of
available requested comparators. All-protein, assay-type, individual coverage,
and five-fold results are secondary. Existing exposed development data only.

## scripts/beat_am/m12_prepare_speed_assets.py

Discover and pin public inference assets for the speed comparison.

Downloads source and selected small published head weights, never foundation
model weights or new experimental outcomes. All files stay inside this project.

## scripts/beat_am/m13_local_speed_benchmark.py

CPU timings for cached-score lookup and prediction with cached inputs.

No training, experimental labels, model selection, foundation-model download,
or sequence-encoder timing. Published weights are pinned by m12; VespaG cached
inputs are read by m14. The distinctions between timing stages are explicit.

## scripts/beat_am/m1_prepare_multirule_refinement.py

Prepare three standalone refinements without scoring or benchmarking a model.

Only cached annotations and metadata/tabular features are opened. No target,
AlphaMissense/ESM score array, model.predict call or correction fit is used.
The existing experiment and its invalidated enforcement verdict are preserved.

## scripts/beat_am/m20_dms_comparator_pvalues.py

Paired block sign-flip tests on the original membrane CV scores.

No predictor is fitted. No test split is shuffled or replaced. AlphaMissense
and ESM1b are excluded from the analysis; GEMME, VESPA and VespaG are the
requested comparators, with VESPA-light retained as a supplementary comparator.

## scripts/beat_am/m28_gemme_membrane_census.py

Positive-inside audit of the authors' public full-protein GEMME scores.

Uses exact sequence matches between the public human prediction corpus and the
existing topology census. No training, DMS outcomes, or choice based on effects.
The published ProteinGym comparison remains a separate cohort.

## scripts/beat_am/m28a_resume_gemme_download.py

Resume the frozen GEMME audit with fewer requests and explicit backoff.

Only the HTTP transfer grouping changes. Source bytes still pass the original
ZIP member CRC checks, and the scientific protocol and estimator are unchanged.

## scripts/beat_am/m28b_finish_gemme_download.py

Finish the same public GEMME transfer in four bounded concurrent batches.

Fewer, larger requests avoid the small-request rate limit; explicit Retry-After
backoff is retained. The original scientific protocol and source checks remain.

## scripts/beat_am/m2_standalone_refinement_scorer.py

Correct feature construction for the saved predictor-free head.

Importing this file does not load a model or score a variant. The preparation
stage checks only feature parity. Real scoring is deliberately deferred.

## scripts/beat_am/m33_glycosylation_constraint_repair.py

Label-free glycosylation repair with group-disjoint coefficient calibration.

Audit stage: cached standalone outputs, annotations, frozen heads and coefficient
metadata only. No experimental outcomes and no model/PCA fit. The optional
accuracy stage reads the existing exposed development benchmark only after
the complete repair profiles and audit results have been written and hashed.

## scripts/beat_am/m35_verify_current_disulfide.py

Recheck the established disulfide-loss criterion on the current score stack.

No model, coefficient, threshold or assay outcome is fitted or opened. The
annotation-selected endpoint/background design is the existing internal
consistency test; the separate bond-pairing criterion is not redefined.

## scripts/beat_am/m4_benchmark_multirule.py

Compare the standalone refinements by grouped CV, then guarded LOOCV.

The existing development table is reused. No reserved experimental outcome is
opened. Each outer evaluation group is absent from every supervised fit, PCA,
and correction calibration used to predict it. Correction calibration uses
inner group-out-of-fold residuals, not in-sample residuals of the saved head.

## scripts/beat_am/m7_audit_standalone_rules.py

Label-free B8/B9/B10 audits of already fitted standalone predictor heads.

No fitting, label table, predictor-score lookup, download or outcome evaluation.
Primary: saved full-fitting-pool standalone head. Sensitivity: existing guarded
LOOCV heads. Native identity cells remain missing; old buggy census scores are
not reused. Fixed coefficient probes diagnose feature/audit alignment only.

## scripts/export_blinded_flank_events.py

Export the Gate 1 flank cohort as a per-event table with no clinical label.

The ClinVar classification is used only as an inclusion filter, exactly as in
run_gate1.py: an event must carry a high-confidence P/LP or B/LB aggregate
classification at two stars or better. The label itself is never written out, so
this table can be used to plan the analysis (cluster sizes, power, covariate
structure) while the outcome stays blind.

Labels are attached later, after the preregistration is frozen, by joining on
variation_id.

## scripts/extract_am_flank_scores.py

Extract AlphaMissense scores at TMD-flank positions of the human ER membrane census.

Label-free: this script never reads ClinVar. It produces the substrate for Gate 0,
which asks whether AlphaMissense itself encodes the positive-inside rule.

Output is a compact .npz with one row per (canonical protein position, alternate
amino acid) inside `--max-distance` residues of an annotated TMD edge.

## scripts/extract_am_full_flank_scores.py

Extract AlphaMissense scores at TMD flanks from the ALL-substitutions table.

The genomic hg38 table only contains substitutions reachable by a single nucleotide
change, about 6 of the 19 per residue, and which 6 depends on the codon. That is a
real limitation of the earlier audit: the charge-preserving comparator K<->R needs
an AGA/AGG arginine codon, so CGN arginines were dropped entirely.

AlphaMissense_aa_substitutions.tsv.gz carries all 19 substitutions per residue, so
this removes the codon selection and puts AlphaMissense on the same alphabet as
ESM1b and a saturating deep mutational scan.

Output matches the schema of extract_am_flank_scores.py so every downstream
estimator runs unchanged.

## scripts/extract_esm1b_flank_scores.py

Extract ESM1b scores at TMD flanks, in the same shape as the AlphaMissense extract.

Source: ALL_hum_isoforms_ESM1b_LLR.zip from the esm_variants portal (Brandes et al.,
Nat Genet 2023), one CSV per UniProt accession. Each file is a matrix: columns are
"<wild-type residue> <position>", rows are the 20 possible residues, values are
log-likelihood ratios.

Two conventions differ from AlphaMissense and are normalised here:

* **Sign.** ESM1b's LLR is negative for damaging substitutions; AlphaMissense's
  pathogenicity is positive for them. This writes `score = -LLR`, so in the output
  higher always means more damaging and every downstream estimator is unchanged.
* **Alphabet coverage.** ESM1b scores all 19 substitutions per residue, while the
  AlphaMissense genomic table only covers those reachable by one nucleotide change.
  Everything is kept here; restricting to a common set is done at comparison time.

## scripts/form_ab_lib.py

Shared implementation for the Form A and Form B studies.

Form A reads site geometry without loading AlphaMissense scores. Estimation uses
the installed validated panel module. The working-copy audit module supplies the
report card and is pinned by SHA-256.

## scripts/form_b9.py

Train the B9 neural student with a K/R-gain consistency term.

The model uses native residue, alternate residue, side, and distance. Rule
contrasts are exact linear functions of its 20 x 20 x 2 x 60 prediction table and
the report card's gene-by-native-residue weights. The pilot stage evaluates on
validation clusters, and the full stage evaluates on test clusters. Derived
AlphaMissense weights are not published.

## scripts/paper_support/recheck_final_paper.py

Check the final paper against saved outputs, optionally revisiting score arrays.

Added on 10 November 2026 for the support bundle. This was not an original
experiment runner. Saved-result checks require only the evidence in this bundle.
The optional score checks require the original mapped NPZ files and recompute
comparator point estimates and the descriptive GEMME-cohort K/R frequencies.

## scripts/phase_c_engine.py

Shared Phase C estimator for the analysis and design simulation.

Implements corrections F2 through F4 from the frozen correction plan. The same
stratified estimator and weights are used for point estimates, percentile
bootstrap draws, wild bootstrap draws, and design simulations. Confidence
intervals are obtained by inverting the preregistered test.

## scripts/phase_c_power.py

Phase C feasibility, computed from Phase B counts only.

Reads the DMS inventory, which records how many variants fall in each flank cell.
It does not read any measured effect, so running this does not break the Phase C
blind on DMS outcomes.

The binding constraint in Phase C is not the variant count, it is the number of
independent proteins. Three GPCRs assayed by surface expression contribute most of
the data, so this script reports the leave-one-family-out picture alongside the
full one.

## scripts/plan_power.py

Design quantities for the preregistration, computed without unblinding outcomes.

Reads the blinded flank-event table. The only outcome information it uses is the
pooled P/LP fraction across all flank events, which Gate 1 already reported and
which is not split by side or charge class.

Reports, per preregistered cell: cluster structure (events per gene), the number of
conditional strata that carry information, and the difference a two-sided test can
detect at 80% power once gene clustering is taken into account.

## scripts/plot_gate0.py

Figure for Gate 0: side asymmetry in AlphaMissense scores versus distance to the TMD.

Reads results/gate0_*/gate0_summary.json and writes gate0_side_asymmetry.png.
The matching table view is gate0_estimates.csv in the same directory.

## scripts/plot_hero_grid.py

Hero figure: evolution's footprint, two predictors, and the lab data.

Rows are sources of evidence, columns are the two positive-charge rules. Every
predictor curve is divided by that model's own spread of within-site penalties,
because AlphaMissense is bounded 0-1 and ESM1b is an unbounded log-likelihood
ratio; raw numbers across the two mean nothing.

The composition row is a residue count, not a model output, so it keeps its own
units. In the gain column it is drawn sign-flipped, because the rule's prediction
flips: positive residues are depleted from outer flanks, so *adding* one there is
the costly direction.

The fourth row uses the finalized Phase C summary and reports its sample sizes.

## scripts/plot_within_site.py

Figure: within-site charge penalties by membrane side and distance.

Reads results/gate0_*/gate0_within_site.json and writes within_site_charge_rules.png.

## scripts/recover_unlocated_score_sets.py

Recover DMS score sets that Phase B could not place on a census sequence.

Phase B located a deposit by searching for its target sequence as an exact
substring of the census protein. That fails whenever the deposit uses a different
isoform, and eleven targets were lost that way.

Exact substring matching is not needed. Every scored variant already carries its
own reference residue and position, so the offset can be recovered from the
variants themselves: scan every offset that keeps the variants inside the census
sequence and keep the one that reconciles the most reference residues. This is the
rule already frozen in PHASE_B_INCLUSION_RULES_FROZEN.md section 2, applied
without the substring precondition.

Reads counts only. No measured effect is examined.

## scripts/route_a_a0_feasibility.py

Route A, step A0: the feasibility count, which is the decision gate.

Three numbers, then a power simulation, then the decision floor. The DMS_score column is never
read; every count uses variant identity and position only.

Counts
  1. ProteinGym assays whose protein maps onto the topology census (accepted mapping).
  2. Of those, membrane proteins with at least one measured flank site inside the 10-residue window.
  3. K/R-gain substitutions (single mutants, alternative K or R, native residue in ACFGILMNPQSTVWY)
     at in-window flank sites.

Floor, fixed in the build plan before looking: at least 8 assays with at least 30 K/R-gain flank
substitutions each. Below the floor, Route A's capability test is too thin to answer.

The power simulation uses the real K/R-gain flank variants and their real uncorrected AlphaMissense
scores, with a hypothetical DMS outcome whose side effect is planted in SD units. It answers "what
paired difference in per-assay Spearman can this design detect?" Every outcome in it is simulated.

## scripts/route_a_common.py

Shared loading, alignment-based mapping and flank geometry for Route A.

Route A asks whether the one-parameter correction helps on ProteinGym deep mutational scans.
Every ProteinGym assay has to be placed on the topology census before any score is looked up.
This module does that placement and nothing else: it never reads a DMS_score value.

Mapping is an alignment, not a global offset. The SERT cassette (scripts/load_sert.py) is the
reason: a best-fit offset of -24 mismaps 215 residues, while the two-segment alignment is exact.

Geometry (side and distance to the nearest TMD edge) comes from read_topology() in run_gate1.py
followed by build_flank_map() in extract_am_flank_scores.py, which is what produced
data/processed/gate0/am_full_flank_scores.npz. Using the same path keeps Route A on the geometry
that every earlier analysis used.

## scripts/route_a_fetch_proteingym.py

Route A, step A1: fetch and verify ProteinGym v1.3 substitutions, and write the download manifest.

Each file goes into its own directory under data/raw/proteingym/, which must be empty before the
first download. Downloads are checked against the size and SHA-256 recorded here, so a changed
upstream file stops the run instead of flowing into the analysis. Archive members are checked for
unsafe paths before extraction, and nothing in the archive is executed.

Sources, both from the ProteinGym repository (https://github.com/OATML-Markslab/ProteinGym):
  * DMS_ProteinGym_substitutions.zip, release v1.3, as listed in the repository README.
  * reference_files/DMS_substitutions.csv, the reference table on the main branch (blob 8d1ea9a1).

## scripts/route_a_resolve_uniprot.py

Resolve ProteinGym UniProt entry names to accessions (Route A, step A2).

ProteinGym's `UniProt_ID` column holds UniProt entry names such as `PAI1_HUMAN`, while the
topology census is keyed by accession such as `P05121`. This script asks the UniProt REST
search for every distinct entry name in the reference table and writes the answer, with the
query and the UniProt release, to a JSON file. Nothing here reads DMS scores.

Input:  data/raw/proteingym/reference_github_main/DMS_substitutions.csv
Output: data/raw/proteingym/uniprot_entry_map_v1.3.json (raw API answer, gitignored with data/raw/)

## scripts/route_b_chunk_plan.py

Route B, step B1a: build and verify the ProtT5 chunk plan before any GPU time is spent.

Why this runs first
-------------------
ProtT5 attention is quadratic and the census holds proteins up to 14,507 residues, so sequences are
embedded in overlapping chunks. Only chunks containing panel sites need to be run: 726,444 sites across
4,837 proteins, which is roughly 5,100 chunks rather than 4,863 whole proteins.

The failure mode is silent. An off-by-one in the chunk offset misaligns every embedding to the wrong
residue, the model still trains, and nothing downstream complains. SERT is the precedent: a best-fit
global offset mismapped 215 residues at 66.8% agreement and only an explicit reference-residue check
caught it.

The check here is free and exact. The panel stores `ref_aa` for every site, so for each site the residue
the census sequence carries at that position must equal the residue the panel says is native. Any
mismatch is a mapping bug, and the plan is refused rather than written.

Output
------
`chunk_plan.npz` holds, per panel site, the chunk it belongs to and its offset inside that chunk, plus the
chunk table itself (accession, start, end). The embedding job reads the chunk table; the trainer reads the
per-site index.

## scripts/route_b_embed_modal.py

Route B, step B1b: ProtT5-XL per-residue embeddings for the panel sites, on Modal.

Reads `results/route_b_20261012/chunk_plan.npz`, which `route_b_chunk_plan.py` wrote after verifying every
one of the 726,444 panel sites against the census sequence. Embeds the 5,043 chunks, gathers the
per-residue vector for each panel site, and writes one fp16 array per batch to a Modal Volume.

The gather happens inside the GPU function on purpose. Returning whole chunk embeddings would move
2,605,296 x 1024 x 2 bytes, about 5.3 GB; gathering first reduces it to 726,444 x 1024 x 2, about 1.49 GB.

Frontend choice. ProtT5-XL-UniRef50, not ESM. ESM1b already satisfies the K/R-gain rule
(comparator-robust at -0.073, 2026-09-29 entry), so building on ESM representations risks rule compliance
arriving from the frontend rather than from the loss term, which would make the central claim
unfalsifiable. ESM2 is the same family and the same risk. ProtT5 is a different lineage, and the plain
student control measures its baseline contrast rather than assuming anything.

Run
---
    pip install modal && modal setup          # one time, authenticates to your account
    modal run scripts/route_b_embed_modal.py

Then fetch the output and assemble it:

    modal volume get route-b-embeddings / ./results/route_b_20261012/embeddings
    python scripts/route_b_embed_modal.py --assemble

Cost and GPU choice
-------------------
The weights are baked into the image by `_bake_model`, so the ~11 GB pull happens once at image build and
no container ever downloads them. Chunks are length-bucketed into padded forwards under TOKEN_BUDGET
tokens, which turns 5,043 single-sequence forwards into a few hundred batched ones; without that the job
is launch-overhead-bound and a faster GPU buys nothing.

With both fixes, 2,605,296 tokens is roughly 1-2 minutes of compute plus about a minute of model load:

    A100 40GB  $2.10/hr   ~4 min container   ~$0.14
    H100 SXM5  $3.95/hr   ~3 min container   ~$0.20
    L4 24GB    $0.80/hr   ~9 min container   ~$0.12

The fixed overhead (cold start, model load to GPU) does not shrink with a faster card, so H100 saves about
a minute over A100 and costs a few cents more. Any of the three is fine; pick whichever has capacity.
Modal bills per second and the batch files are idempotent, so a failed or partial run is cheap and a rerun
only redoes missing batches.

### Retained implementation rationale for scripts/route_b_embed_modal.py

One task for the whole job. The model is loaded inside the function, so each task costs another load (about a minute); with only ~2 minutes of compute in total, splitting the work would cost more in repeated loads than it saves in parallelism. 5,043 chunks pickle to ~8 MB, which Modal moves as a blob.

## scripts/run_composition_profile.py

Positive control: are the charge-topology rules visible in the census sequences?

Label-free and model-free. This reads only the amino-acid sequences and the
topology annotations, and measures how often K/R and D/E actually occur at each
distance from a transmembrane boundary, on each side.

If the positive-inside rule is real in this proteome, cytosolic flanks should be
K/R-enriched relative to outer flanks, and the enrichment should fade with
distance. That gives an evolutionary decay length to compare against the decay
length in AlphaMissense's scores: a predictor that learned the rule from sequence
should show a similar profile.

Nothing here depends on a predictor or a clinical label, so it is the reference
both of those are measured against.

## scripts/run_form_a.py

Train Form A student models on planted rules and audit held-out predictions.

The pilot stage uses separate seeds to verify the pipeline and model capacity.
Synthetic targets use observed membrane geometry without
reading AlphaMissense scores. `summarize_form_a.py` reads the outputs.

## scripts/run_form_b.py

Run the Form B correction and evaluate it on held-out clusters.

The pilot stage fits on training clusters and evaluates on validation clusters.
The full stage evaluates on test clusters. Corrected scores are temporary audit
inputs and are deleted after use. No clinical labels are read.

## scripts/run_gate0.py

Gate 0: does AlphaMissense itself encode the positive-inside rule?

Label-free. This script never reads ClinVar. It asks whether AlphaMissense scores
the *same* amino-acid substitution in the *same* protein differently depending on
which side of the membrane the flank faces.

Side difference (per substitution class, per scope):
  Stratify by (protein, exact substitution). Keep strata with at least one scored
  variant on each side. Within a stratum take the cytosolic-minus-outer difference
  in AlphaMissense pathogenicity, then combine strata with efficient weights
  n_cyto*n_outer/(n_cyto+n_outer). Also report the stratified rank statistic
  P(cytosolic score > outer score), which assumes no scale.

  Positive = AlphaMissense calls the substitution more damaging on the cytosolic side.

A raw side difference is not evidence for the positive-inside rule on its own: if
cytosolic flanks are simply more constrained, every substitution class shifts the
same way. The rule makes a *directional* prediction, so the primary quantities are
differences in differences:

  positive_loss_minus_gain     losing K/R should hurt more inside while gaining K/R
                               should hurt more outside, so the contrast is signed
                               and cancels any side baseline shared by all classes.
  positive_loss_minus_placebo  excess over substitutions that change no charge.

Uncertainty is a protein-clustered bootstrap. Every class in a scope is resampled
under the same protein draw, so contrasts keep the correlation between classes.

## scripts/run_gate0_within_site.py

Gate 0 robustness check: within-site charge contrasts in AlphaMissense.

Label-free. Gate 0 stratified by protein and exact substitution, which controls
protein identity and substitution identity but *not* how constrained the
individual residue is. If cytosolic K/R residues sit in more constrained sites
than outer ones, Gate 0's side difference is still confounded.

This removes that. At each K or R site, AlphaMissense scores both the
charge-preserving substitution (K<->R) and charge-losing ones. Their difference
is a charge-loss penalty computed entirely inside one site, so any constraint
shared by the whole residue cancels:

    penalty(site) = mean(score | charge lost) - mean(score | charge preserved)

The side comparison is then a second difference, taken across sites within the
same protein and the same reference residue. The same construction runs at D/E
sites using D<->E as the preserving substitution.

If the positive-inside rule is what drives Gate 0, the penalty should be larger
on the cytosolic side and should shrink with distance from the membrane.

## scripts/run_gate1.py

Compute outcome-blind feasibility counts for the membrane-topology study.

Clinical labels select eligible variants, but cell-level outputs combine benign
and pathogenic events so the analysis remains blinded before preregistration.

## scripts/run_method_label_controls.py

Side-label destruction controls and declared relative-score candidate counts.

The original observational contrast receives no randomization p-value. The
negative control tests the pipeline after its topology labels are destroyed.

## scripts/run_phaseA_rule_audit.py

Phase A: does AlphaMissense already encode the positive-inside rule?

Label-free audit. Reads only the Oxford-MRC census topology and AlphaMissense
pathogenicity scores. No ClinVar classification is opened, so this cannot
unblind the Gate 1 side-by-charge-change cells.

The audit asks how AlphaMissense's score varies with (a) which side of the
membrane an aqueous residue sits on, (b) the signed change in formal charge the
substitution makes, and (c) how far the residue is from the nearest TMD edge.
If the predictor already encodes the rule, cytosolic positive-charge losses
should score higher than the same losses on the outer side, and the gap should
decay with distance from the membrane.

## scripts/run_phaseB_dms_inventory.py

Phase B: inventory deep mutational scans that map onto the census topology.

The gate asks whether at least four census membrane proteins have an abundance,
surface-expression or trafficking scan with enough charge changes in the aqueous
flanks to fit the Phase C model. Function readouts (currents, signalling, binding,
catalysis) are reported but flagged, because mixing them with expression readouts
would muddy the model.

Numbering is resolved against the target sequence that MaveDB stores with each
score set, not against an accession, so an isoform offset cannot silently shift
every residue. A score set is only counted when its target sequence can be
located exactly inside the census sequence.

Label-free with respect to ClinVar: no clinical classification is read here.

### Retained implementation rationale for scripts/run_phaseB_dms_inventory.py

Most MaveDB targets store the coding DNA, so translate it. The deposited region is not always in frame 0, so every frame is tried and the one that lands inside a census sequence wins.

Locating the target sequence gives a candidate offset, but some deposits omit the start codon while still numbering hgvs_pro in full-protein coordinates. Choose the offset that actually reconciles the reference residues, and keep the agreement rate as the evidence that it is right.

## scripts/run_phase_c.py

Phase C, run exactly as preregistered.

Frozen protocol: docs/PHASE_C_PREREGISTRATION.md
SHA-256 754a9a713a2afc7fd2e3e43be9af55bb762805feeb3d3991be0667c207635d79
Frozen 2026-09-28 23:16:29 PDT, verified before this script was written.

Primary: the within-site K/R gain side difference over the <=10 residue window,
cytosolic minus outer, with a restricted wild cluster bootstrap-t clustered on
gene, reported under both Webb and Rademacher weights.

Every residue set, sensitivity fit and decision rule here is transcribed from the
frozen text. Deviations are recorded in the output under `deviations` and belong
in the decision log, never in the frozen document.

### Retained implementation rationale for scripts/run_phase_c.py

Frozen rule: "a score set whose anchors do not separate is excluded". Separation is the rank statistic; 0.5 is no separation at all. The observed values fall in two clean groups, 0.501 and then 0.615 upward, so this threshold does not sit near any set's value.

## scripts/run_phase_c_corrected.py

Run Phase C under the frozen correction plan.

Correction plan: docs/PHASE_C_CORRECTION_PLAN.md
SHA-256 c74e03ec741bcd08c4062e6a282db972e145d486b4e54540c4295d0a331dca41
frozen before any corrected number was computed.

Protocol unchanged: docs/PHASE_C_PREREGISTRATION.md, SHA-256 754a9a71...
Implements matched substitutions, consistent inference weights, a shared analysis
and simulation estimator, test-inversion intervals, and contributing-site counts.

## scripts/run_phase_c_matched.py

Compare lab and predictor side effects on identical variants.

Implements section 7 of the frozen preregistration. The same within-site estimator
is applied to lab scores and each predictor on a shared set of protein, position,
and substitution tuples. Each bootstrap draw resamples both measurements together
to preserve their correlation.

Each source is standardized by its within-site penalty spread before differencing.
ESM1b is the positive control.

## scripts/run_rule_audit.py

Demonstrate the reusable tool on two models using the earlier raw settings.

Only the original raw-score, all-uncharged K/R-gain comparison is run.
No comparator, scale, subgroup or simulation robustness analyses are invoked.

## scripts/run_step2_arms.py

Estimate slope uncertainty and decompose the gain penalty by arm.

The five folds share training data, so their slopes are not independent estimates.
This script computes a protein-clustered bootstrap interval for the slope.

The gain penalty is score(X -> K/R) - score(X -> uncharged). The analysis regresses
the K/R and comparator arms separately against their B1 counterparts, then repeats
the comparator analysis for polar-only and hydrophobic-only residue sets. This
distinguishes a charge-specific signal from constraints in the comparator arm.

Label-free throughout.

## scripts/run_step2_calibrated.py

Compare each predictor with the compositional footprint after calibration.

The preregistered self-normalized metric is unsuitable here because B1 is a
deterministic function of side, distance bin, reference residue, and alternate
residue. Within a stratum it has no measurement noise, so comparisons with noisy
predictors primarily reflect noise rather than footprint use. The deviation is
recorded in the decision log.

Primary analysis: fit the predictor's within-site penalties on B0 and B1 using
training families. On held-out families, compare the observed side difference
with beta times the B1 side difference. The residual measures side-specific over-
or under-use of the footprint in the predictor's own units.

Secondary analysis: average penalties within reference, alternate, and distance
strata, then correlate the predictor and B1 side differences.

Label-free throughout.

### Retained implementation rationale for scripts/run_step2_calibrated.py

The older Gate 0 helper defines loss as every alternate outside the reference charge group, including reversals and histidine. Step 2 must use the clean, uncharged comparator used by the project's later rules.

## scripts/run_step2_composition_control.py

Step 2: does a predictor follow the compositional footprint, fall short, or go beyond it?

Frozen design: docs/STEP2_COMPOSITION_CONTROL_PREREG.md
SHA-256 06cdc2e8449ab03ced939b6cf1599d7d754a5fcec36c4a2ca992c5e9d1de5385

Composition by side and distance is the evolutionary footprint of the topology
rules, not a rival explanation to them. So the baselines here are not controls to
be subtracted; they are alternative "predictors" scored through the same
within-site estimator, and the question is whether a real predictor reproduces
them, falls short, or exceeds them.

Label-free: no clinical label and no deep-mutational-scanning effect is read.

Within-site penalties use the audit's residue definitions (run_gate0_within_site),
so the predictor numbers here reproduce the ones disclosed in section 2 of the
frozen design rather than a differently-defined variant of them.

## scripts/run_step3_physics_comparison.py

Compare cytosolic K/R gain profiles with ER insertion measurements.

Frozen design: docs/STEP3_PHYSICAL_ANCHOR_PREREG.md
SHA-256 40a411e3968bb1c0f0ad1365bc4cf47878411a82894371c1814e75aa701021af

Lerch-Bader et al. (2008) measured the free-energy effect of adding a positive
charge at different distances on the cytosolic side. The corresponding model
quantity is the cytosolic within-site gain penalty by distance.

Negative insertion ddG and negative model penalties both indicate favorability.
Both values are negated for comparisons on a common favorability axis. These
descriptive comparisons do not establish model failure.

Label-free.

## scripts/run_within_site_controls.py

Controls and scale fixes for the within-site charge audit.

Five things the headline result needs before it can be frozen.

**Standardised units.** AlphaMissense scores run 0 to 1; ESM1b is an unbounded
log-likelihood ratio. "+1.11 versus +0.051" compares nothing. Every side
difference here is also divided by the standard deviation of that model's own
within-site penalties, and reported as a rank statistic P(cytosolic > outer),
which assumes no scale at all.

**Glycosylation control.** On the outer side, S/T to K/R destroys an N-X-S/T
sequon, and N to K/R destroys it too. Losing glycosylation is damaging for
reasons that have nothing to do with the positive-inside rule, and it only
happens on the outer side, so it biases the gain contrast toward the hypothesis.
Sites participating in a sequon are excluded here.

**Lysine versus arginine.** Lysine is a ubiquitination target, so losing one can
look benign and gaining one damaging for degradation reasons. That biases against
the hypothesis. Arginine-only is the cleaner test.

**Floor effects.** AlphaMissense is bounded, so a side difference could be
compression near an endpoint. Repeated on the logit scale and on ranks.

**Ion channels.** The census carries no pore or re-entrant loop annotation, so
pore loops cannot be excluded directly. Excluding ion channels entirely is the
available proxy, and it matters because Kir2.1 has a pore loop.

## scripts/search_ortholog_score_sets.py

Ortholog-aware rerun of the Phase B search.

Phase B matched MaveDB deposits to the census by human gene symbol. That hid
`mKir2.1`, a mouse deposit of a protein whose human ortholog KCNJ2 is in the
census, and the same bug could hide other rodent scans.

This enumerates every MaveDB score set, aligns each target sequence to the census
by global alignment, and applies the ortholog rule frozen in
PHASE_B_INCLUSION_RULES_FROZEN.md section 4a:

* at least 95% identity over the mapped region
* topology conserved, checked here as: the alignment covers every transmembrane
  segment and its +/-10 flank with no gaps, so census topology transfers position
  by position
* only positions where the two reference residues are identical are usable
* cross-species sets flagged

Counts only. No measured effect is read.

## scripts/simulate_phase_c_design.py

Evaluate inference methods on the planned Phase C design.

Phase C has 16 genes, or 12 after excluding GPCRs. The simulation uses these
cluster counts and the blinded gene-level cell sizes.

The following methods are compared:

1. Percentile bootstrap over genes.
2. Restricted wild cluster bootstrap-t with Rademacher weights.
3. Restricted wild cluster bootstrap-t with Webb weights.

Rademacher weights admit only 2^G distinct draws. At G = 12 that is 4,096, fewer
than the 9,999 replicates planned, so the bootstrap distribution is coarse and the
attainable p-values are granular. Webb's 6-point distribution is the standard fix
below about 12 clusters.

No measured outcomes are used.

## scripts/simulate_realized_design.py

Simulate Phase C using the realized design and shared analysis engine.

Uses the observed per-(gene, reference residue) cell sizes from the corrected run.
The analysis and simulation share the estimator implementation so the reported
size and power apply to the fitted procedure.

## scripts/sweep_ortholog_names.py

Complete the ortholog sweep by querying every MaveDB target name.

The score-set search endpoint returns at most 100 results and exposes no paging
parameter, so the first sweep screened only 100 deposits and found Kir2.1 only by
direct URN lookup. This queries each of the ~1,157 target gene names in turn,
caching every response, so the screen covers the whole catalogue.

Applies the ortholog rule frozen in PHASE_B_INCLUSION_RULES_FROZEN.md section 4a.
Counts only; no measured effect is read.

## scripts/validate_method_benchmark.py

Known-answer benchmark of the retained K/R-gain estimator and explicit baselines.

The synthetic design is saved before execution in METHOD_COMPLETION_PLAN_20261004.
All scenarios survive reporting. No real-data comparator or scale sweep is run.

## scripts/verify_estimator.py

Verify the primary estimator with an independent implementation and simulation.

The independent implementation uses fixed-effects least squares rather than the
production stratified weighted mean. The two algebraically equivalent estimates
must agree without sharing estimator code.

The simulation preserves the observed sites, sides, strata, proteins, and
distances while replacing outcomes with values generated from a known effect. It
checks estimator bias and protein-clustered interval coverage. No measured outcome
is used.

## scripts/verify_wild_bootstrap.py

Compare inference methods on the real design: percentile bootstrap vs wild cluster bootstrap-t.

The audit currently uses a protein-clustered percentile bootstrap, which the
simulation showed rejects at 6.0% against a nominal 5%. That is typical with a
few dozen clusters, but the restricted wild cluster bootstrap-t usually controls
size better, because it imposes the null when generating replicates and
studentises the statistic.

Implementation, for the fixed-effects form of the estimator:

    y~ and x~ are stratum-demeaned outcome and side indicator
    beta = sum(x~ y~) / sum(x~^2)
    cluster-robust se from per-protein score sums

Restricted (null-imposed) replicates take the residual under H0: beta = 0, which
is simply y~, and flip it by a protein-level Rademacher weight. The bootstrap t
distribution is then compared with the observed t.

Uses no outcome data: the outcome is simulated on the real design.

## tests/test_form_ab.py

Test the Form A/B generator, oracle, students, split, guards, and correction.

Run with: .venv_ml/Scripts/python.exe -m unittest discover -s tests -p test_form_ab.py -v

## tests/test_rule_audit_placebo.py

Checks for the placebo rules, disjoint arms and shared-draw contrasts.

Run with: python -m unittest discover -s tests -p test_rule_audit_placebo.py -v

## tools/membrane_rule_audit/audit.py

Audit mutation-score contrasts near membrane boundaries.

For each charge and neutral-placebo rule, the audit reports:

1. A within-site side difference that removes additive position effects.
2. A descriptive, cross-fitted composition baseline.
3. A distance profile relative to the nearest membrane boundary.

Each rule has a declared direction; placebo rules expect zero. Differences against
a placebo use shared protein resamples to preserve correlation. The composition
calibration is descriptive and does not identify an internal model mechanism.

Usage
-----
    python -m membrane_rule_audit --topology census.html --scores scores.tsv \
        --score-direction higher-is-damaging --out results/

See README.md for the score-file formats and for how to add a new predictor.

### Retained implementation rationale for tools/membrane_rule_audit/audit.py

Uncharged, conservative, and with no side asymmetry of their own: whatever these pick up is a side baseline the charge rules also inherit. Two sets, so a placebo reading cannot rest on one choice of residues.

Polar probes rather than placebos. If a side effect tracks polarity instead of charge, these carry it and the hydrophobic pair above does not. N is the N-glycosylation residue, so gaining it on the outer side can create a sequon; Q is its non-glycosylatable isostere and isolates polarity on its own.

What the positive-inside rule predicts for each cytosolic-minus-outer contrast. Losing a cytosolic K/R should hurt more inside; gaining a K/R should hurt more outside. None means this tool declares no direction, so no verdict is issued.

Every gain rule's group held out at once, so one comparator serves them all and a contrast between two gain rules differs only in its numerator. Narrower than the per-rule comparator, so it is reported beside the 14-residue definition, not instead.

## tools/membrane_rule_audit/disulfide.py

Evidence-aware, matched-site audit of structural disulfide loss.

The estimand compares an IDENTICAL Cys replacement at bonded and experimentally
verified free sites. It never subtracts the other replacements at that position:
every non-Cys replacement removes the original bond. Unknown is not free.

Only NumPy and the standard library are required. Score files are not opened
until annotation, matching, exposure and independent-family gates have passed.

## tools/membrane_rule_audit/glycosylation_constraints.py

A label-free lower bound for the existing glycosylation-loss coefficient.

The inputs are matched score contrasts, never experimental outcome labels.
Each constraint is affine in the coefficient: delta(beta) = delta(0) +
beta * feature_delta.  A expects a nonnegative delta; B a nonpositive delta.
The smallest shared nonnegative beta satisfying the calibration constraints
is the maximum of their individual lower bounds.  Application to a new
protein is not guaranteed to satisfy that protein's constraints.

## tools/membrane_rule_audit/panel.py

Reusable raw-score panel analysis using the earlier uncharged comparator.

This module packages the existing within-site estimator. It does not run
alternative comparators, transformed scales, subgroups, or calibration studies.

## tools/membrane_rule_audit/proline.py

Conditional Pro-gain audit in independently defined native structural contexts.

Native Pro loss is excluded. A TM-context contrast is not by itself proof of a
membrane-specific physical mechanism or clinical pathogenicity. Geometry and
matching are fixed before score access; comparator/scale sensitivity is retained.

## tools/membrane_rule_audit/refinement.py

Annotation-driven additions to a predictor that consumes no predictor scores.

These are directional correction features, not clinical probabilities or a claim
that every mutation following a chemical pattern has the same functional effect.
Real coefficients must be calibrated and the complete predictions evaluated
later. In particular, the control-free disulfide arm tests predictive utility;
it does not identify a causal disulfide-specific penalty.

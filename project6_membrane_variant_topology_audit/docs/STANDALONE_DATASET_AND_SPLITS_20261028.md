# Standalone dataset, hold-outs and remaining data

28 October 2026. This is an accounting of saved inputs and evaluation records, not a new training or benchmarking run. No new experimental target table or reserved outcome was opened.

## 1. What the model learns from

The supervised dataset is the existing processed ProteinGym DMS subset: **127,811 variant–assay measurements, 44 assays, 39 human proteins, 34 dependence groups**. It contains **43,991 measurements across 15 assays in 13 membrane proteins**, plus **83,820 measurements across 29 assays in 26 other proteins**. Counts are measurements, not necessarily unique substitutions across repeated assays.

The target is within-assay rank-normalized experimental damage: `rankdata(-DMS_score)/n`. It is not a clinical pathogenicity probability. The frozen ProtT5 encoder supplies representations and is not fine-tuned here. Training fits the PCA projection and regression head; correction calibration fits coefficients against inner out-of-fold DMS prediction errors. The standalone feature set has 18 substitution/topology features plus 64 ProtT5 PCs. AlphaMissense/ESM1b prediction scores are not inputs or targets.

The existing saved full-fit head uses the whole 39-protein development pool for scoring/behavior characterization. The quoted accuracy figures, including 0.4626, instead come from cross-validated refits; they are not in-sample accuracy from that saved full-fit head.

## 2. What was held out in the main benchmark

**The 13 membrane proteins were not a permanent test set.** All 39 proteins rotated through guarded group-out evaluation. For each of 34 outer folds, one entire dependence group was excluded from PCA fitting, head training and correction calibration. Its proteins were scored only after those fits. This left **36–38 proteins in training and 1–3 in testing**, depending on group size. All outcomes for one protein, across all its assays, stay together.

For example, OPSD (rhodopsin) and CCR5 belong to one dependence group: both were held out together and the other 37 proteins supplied the outer training pool. Inner three-fold group-out predictions within that training pool supplied correction-calibration residuals. The outer test proteins never supplied those calibration targets.

Every one of the 39 proteins received held-out predictions once under this scheme. The headline membrane score averages the 13 membrane proteins, with assays first averaged within protein. The 26 other proteins form the secondary evaluation. Historical project-wide exposure remains separate from this within-run separation: this is cross-validation on previously used development data, not a newly untouched confirmation set.

### The 13 membrane proteins in the headline table

| Dataset label | Accession | Measurements | Assays | Group |
|---|---|---:|---:|---|
| GLPA | P02724 | 245 | 1 | DG_P02724 |
| ERBB2 | P04626 | 326 | 1 | DG_P04626 |
| OPSD | P08100 | 165 | 1 | DG_P08100 |
| MET | P08581 | 5,393 | 1 | DG_P08581 |
| CP2C9 | P11712 | 12,512 | 2 | DG_P11712 |
| LYAM1 | P14151 | 359 | 1 | DG_P14151 |
| EPHB2 | P29323 | 1,239 | 1 | DG_P29323 |
| SC6A4 | P31645 | 11,576 | 1 | DG_P31645 |
| CCR5 | P51681 | 6,137 | 1 | DG_P08100 |
| KCNH2 | Q12809 | 200 | 1 | DG_Q12809 |
| SCN5A | Q14524 | 224 | 1 | DG_Q14524 |
| VKOR1 | Q9BQB6 | 3,392 | 2 | DG_Q9BQB6 |
| ACE2 | Q9BYF1 | 2,223 | 1 | DG_Q9BYF1 |

## 3. What was outside training and what is left

The historical 89-protein inventory can be accounted for without treating every unfitted protein as untouched:

| Status in the historical inventory | Count | Meaning now |
|---|---:|---|
| Current development/cross-validation pool | 39 | Already used in the standalone results. |
| Previously reserved, evaluated with an older model | 27 | Seven membrane and 20 other proteins. Not trained into the current 39-protein pool, but outcomes were already used in the project. Available for an explicitly exploratory standalone expansion. |
| Previously reserved, skipped by the older evaluation | 4 | Sequence/embedding or mapping failures; not four ready benchmark additions. |
| Pending confirmation candidate role | 6 | Not in the current 39-protein pool. Role name does not certify usable data or untouched outcomes. |
| Other historical exclusions | 13 | Mapping, coverage, relatedness or other curation barriers; require reassessment before admission. |
| Total historical inventory | 89 | This is not an exhaustive census of all public DMS data. |

The seven earlier membrane evaluation proteins are **NPC1, SLC22A1 (S22A1), HMGCR (HMDH), APP (A4), CD19, ADRB2 and KCNQ4**. The current standalone benchmark has not evaluated those seven. The 20 usable older non-membrane evaluation proteins likewise remain outside the current standalone benchmark. The four skipped accessions are HECD1, Q53Z42, KRAS (RASK) and alpha-synuclein (SYUA). None of this transfers the older model's performance or preregistered status to the standalone model.

The six pending-role proteins are **MC4R, P2RY8, SLC13A5, TECR, LARGE1 and FKRP**. The older local curation has no usable DMS score files for the first three; that is not a claim that no newer public assay could exist. LARGE1/FKRP have recovered and mapped assays with prior QC inspection, but are not integrated into the standalone comparison. TECR has a recovered table but unresolved assay-provenance/quality concerns. These six are therefore not six ready, untouched test datasets.

The other 13 exclusions are: AMFR, AVPR2, ILF3, KCNE1, KCNQ1, NKX31, PITX2, PKN1, RBP1, SRBS1, TNKS2, TPOR, YAP1.

**Important exposure-record conflict:** a later option-C protocol calls LDLR, SGCB and SGCA unopened, but the earlier exposure register explicitly lists all three as already exposed in Phase C. The disulfide-status review also records LDLR's historical exposure. This accounting does not certify those three as fresh confirmation. They may still be useful for an expanded, transparently exploratory comparison. A genuinely untouched test requires a reconciled exposure ledger and a fixed evaluation plan, not just a `pending` or `confirmation` label.

## 4. What the five-fold sensitivity check did

It repeated the same prediction/correction procedure on the same data using **five larger held-out blocks**, with whole dependence groups kept together. For each round, the other four blocks supplied training data; training-fold PCA, the regression head and inner three-fold correction calibration were fitted again. The held-out block was then scored. After five rounds, every measurement had one outer held-out prediction. This was not a random variant split, new data, a test of the full-fit model, or five retrainings on all 39 proteins.

| Fold | Training proteins | Held-out proteins | Held-out membrane proteins | Training measurements | Held-out measurements |
|---|---:|---:|---:|---:|---:|
| 1 | 30 | 9 | 3 | 102,200 | 25,611 |
| 2 | 32 | 7 | 3 | 102,221 | 25,590 |
| 3 | 33 | 6 | 3 | 102,190 | 25,621 |
| 4 | 31 | 8 | 3 | 102,353 | 25,458 |
| 5 | 30 | 9 | 1 | 102,280 | 25,531 |

This asks whether conclusions depend strongly on holding out one group at a time versus holding out larger blocks and having less training data. The settings and data remain fixed; only the outer grouping scheme changes, with corresponding refits.

| Model | Guarded group-out membrane Spearman | Five-fold membrane Spearman |
|---|---:|---:|
| Raw standalone | 0.462330 | 0.445566 |
| PI + disulfide | 0.462596 | 0.445561 |
| All four corrections | 0.459003 | 0.442692 |
| AlphaMissense | 0.455125 | 0.455125 |
| VespaG | 0.450895 | 0.450895 |

The tiny PI-plus-disulfide advantage over raw does not persist under five-fold evaluation, and the ordering relative to AlphaMissense/VespaG reverses. External baseline score files remain unchanged; the standalone fits change. The check therefore limits a broad superiority claim.

**Disulfide coverage caveat:** all measured disulfide-supported examples—OPSD, CCR5 and SC6A4—land in five-fold test fold 3. The remaining training/calibration data have no eligible disulfide-feature variation, so that fold's learned disulfide coefficient is zero. Other folds have no eligible disulfide examples to change at test time. Thus the direct disulfide term changes zero measured held-out rows in this five-fold run. It is not a direct test of a nonzero disulfide correction's predictive benefit. Guarded group-out separates its two supported dependence groups and yields nonzero held-out correction terms.

## 5. The rule audits are a separate evaluation

The label-free standalone rule audits score all 193,515 possible missense substitutions at 10,185 positions in the existing 13 membrane proteins. These are model predictions, not 193,515 additional experimental training labels. The positive-inside primary contrast has 596 contributing positions, 11 proteins and 10 dependence groups. Neither those audits nor the annotation-defined disulfide check create a fresh held-out DMS test set.

## Complete inventories and sources

The accompanying folder `results/standalone_dataset_accounting_20261028/` contains the exact 39-protein list, 44-assay list, five-fold assignments, all 89 historical role statuses, three additional exposure conflicts, and the scheme comparison values. Gene labels are retained exactly as stored in the dataset.

Source records include the current benchmark plan/results, `features/l1_table.npz` metadata only, the frozen historical role manifest, the old evaluation accession list, the exposure-history register and the option-C protocol. No new model was fitted, new performance benchmark run, new experimental target read or reserved dataset fetched.

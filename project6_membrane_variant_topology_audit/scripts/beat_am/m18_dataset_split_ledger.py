"""Reconcile the existing standalone model's datasets and split assignments."""
from collections import Counter
from datetime import datetime,timezone
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'results/standalone_multirule_20261025/benchmark_20261025'
OUT=ROOT/'results/standalone_dataset_accounting_20261028'
SOURCES={}


def source(path):
    path=ROOT/path
    SOURCES[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
    return path


def json_read(path):return json.loads(source(path).read_text(encoding='utf-8'))


def csv_read(path):
    with source(path).open(encoding='utf-8-sig',newline='') as stream:return list(csv.DictReader(stream))


def csv_write(name,rows):
    with (OUT/name).open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    plan=json_read('results/standalone_multirule_20261025/benchmark_20261025/BENCHMARK_PLAN.json')
    table_path=source('results/beat_am_20261016/features/l1_table.npz')
    if SOURCES[str(table_path.relative_to(ROOT))]!=plan['table_sha256']:raise AssertionError('Dataset differs from benchmark')
    # Do not access target, AlphaMissense or ESM1b arrays during dataset accounting.
    with np.load(table_path,allow_pickle=False) as table:
        accession=table['accession'].astype(str);gene=table['gene'].astype(str)
        membrane=table['is_membrane'].astype(int);assay=table['dms_id'].astype(str)
    assert len(accession)==127811 and len(set(accession))==39 and len(set(assay))==44
    historical=csv_read('results/beat_am_20261016/inventory_frozen/group_manifest_frozen.csv')
    roles={r['accession']:r for r in historical}
    old_evaluation=csv_read('results/beat_am_20261016/confirmation/confirmation_per_protein.csv')
    evaluated_before={r['accession'] for r in old_evaluation}
    assert len(evaluated_before)==27 and not evaluated_before.intersection(accession)
    group_for={a:g for g,members in plan['dependence_members'].items() for a in members}
    fold_for={mode:{a:f['fold'] for f in folds for a in f['test_accessions']} for mode,folds in plan['folds'].items()}
    proteins=[];assays=[]
    for a in sorted(set(accession)):
        selected=accession==a
        genes=sorted(set(gene[selected]));assert len(genes)==1
        row=dict(accession=a,gene=genes[0],membrane=int(membrane[selected][0]),
                 measurements=int(selected.sum()),assays=len(set(assay[selected])),dependence_group=group_for[a],
                 guarded_test_fold=fold_for['guarded_loocv'][a],five_fold_test_fold=fold_for['grouped_cv'][a],
                 historical_role=roles[a]['role'])
        proteins.append(row)
        for name in sorted(set(assay[selected])):
            assays.append(dict(accession=a,gene=genes[0],assay=name,membrane=row['membrane'],
                               measurements=int((selected & (assay==name)).sum()),dependence_group=group_for[a]))
    assert len(assays)==44
    used=set(accession)
    historical_rows=[]
    for row in historical:
        a=row['accession'];role=row['role']
        if a in used:status='current_39_training_and_cross_validation_pool'
        elif role.startswith('confirmation_'):
            status='outside_current_pool_old_model_evaluated' if a in evaluated_before else 'outside_current_pool_old_evaluation_skipped'
        elif role=='pending_confirmation_candidate':status='pending_candidate_not_in_current_pool'
        else:status='other_historical_exclusion_not_in_current_pool'
        historical_rows.append(dict(**row,current_status=status,current_standalone_benchmark=(a in used),
                                    evaluated_in_older_model=(a in evaluated_before),
                                    fresh_untouched_status='not_certified_by_manifest_role'))
    counts=Counter(r['current_status'] for r in historical_rows)
    assert sorted(counts.values())==[4,6,13,27,39]
    folds=[]
    for f in plan['folds']['grouped_cv']:
        test=np.isin(accession,f['test_accessions'])
        folds.append(dict(fold=f['fold'],train_measurements=int((~test).sum()),test_measurements=int(test.sum()),
                          train_proteins=len(set(accession[~test])),test_proteins=len(set(accession[test])),
                          train_groups=len(f['train_groups']),test_groups=len(f['test_groups']),
                          test_membrane_proteins=len(set(accession[test & (membrane==1)])),
                          held_out_genes=';'.join(sorted(set(gene[test])))))
    assert sum(f['test_measurements'] for f in folds)==127811
    csv_write('current_39_proteins.csv',proteins);csv_write('current_44_assays.csv',assays)
    csv_write('historical_89_status.csv',historical_rows);csv_write('five_fold_assignments.csv',folds)
    exposure=source('results/beat_am_20261016/inventory/exposure_history_register.md').read_text(encoding='utf-8')
    protocol=source('docs/OPTION_C_CALIBRATION_PREREG_20261025.md').read_text(encoding='utf-8')
    exposure_conflicts=[]
    for name,a in [('LDLR','P01130'),('SGCB','Q16585'),('SGCA','Q16586')]:
        assert a in exposure and a in protocol
        exposure_conflicts.append(dict(gene=name,accession=a,
            protocol_label='unopened; protocol not executed',
            prior_exposure_register='Phase C measured-outcome exposure recorded',
            accounting='historically exposed candidate, not certified fresh confirmation'))
    csv_write('additional_candidate_exposure_conflicts.csv',exposure_conflicts)
    ext=json_read('results/external_predictor_comparison_20261026/results.json')
    comparison=[]
    labels={'raw':'Raw standalone','positive_inside+disulfide_loss':'PI + disulfide',
            'positive_inside+charge_burial+glycosylation_loss+disulfide_loss':'All four corrections',
            'AlphaMissense':'AlphaMissense','VespaG':'VespaG'}
    for model,label in labels.items():
        row=dict(model=label)
        for mode in ['guarded_loocv','grouped_cv']:
            block=ext[mode]['paired_to_VespaG' if model=='VespaG' else 'requested_common']
            values=[r for r in block['summaries'] if r['scope']=='membrane' and r['assay_type']=='all' and r['model']==model]
            assert len(values)==1
            row[mode]=values[0]['mean_spearman']
        comparison.append(row)
    csv_write('evaluation_scheme_means.csv',comparison)
    membrane_proteins=[r for r in proteins if r['membrane']==1]
    lines=['# Standalone dataset, hold-outs and remaining data','',
        '28 October 2026. This is an accounting of saved inputs and evaluation records, not a new training or benchmarking run. No new experimental target table or reserved outcome was opened.','',
        '## 1. What the model learns from','',
        'The supervised dataset is the existing processed ProteinGym DMS subset: **127,811 variant–assay measurements, 44 assays, 39 human proteins, 34 dependence groups**. It contains **43,991 measurements across 15 assays in 13 membrane proteins**, plus **83,820 measurements across 29 assays in 26 other proteins**. Counts are measurements, not necessarily unique substitutions across repeated assays.','',
        'The target is within-assay rank-normalized experimental damage: `rankdata(-DMS_score)/n`. It is not a clinical pathogenicity probability. The frozen ProtT5 encoder supplies representations and is not fine-tuned here. Training fits the PCA projection and regression head; correction calibration fits coefficients against inner out-of-fold DMS prediction errors. The standalone feature set has 18 substitution/topology features plus 64 ProtT5 PCs. AlphaMissense/ESM1b prediction scores are not inputs or targets.','',
        'The existing saved full-fit head uses the whole 39-protein development pool for scoring/behavior characterization. The quoted accuracy figures, including 0.4626, instead come from cross-validated refits; they are not in-sample accuracy from that saved full-fit head.','',
        '## 2. What was held out in the main benchmark','',
        '**The 13 membrane proteins were not a permanent test set.** All 39 proteins rotated through guarded group-out evaluation. For each of 34 outer folds, one entire dependence group was excluded from PCA fitting, head training and correction calibration. Its proteins were scored only after those fits. This left **36–38 proteins in training and 1–3 in testing**, depending on group size. All outcomes for one protein, across all its assays, stay together.','',
        'For example, OPSD (rhodopsin) and CCR5 belong to one dependence group: both were held out together and the other 37 proteins supplied the outer training pool. Inner three-fold group-out predictions within that training pool supplied correction-calibration residuals. The outer test proteins never supplied those calibration targets.','',
        'Every one of the 39 proteins received held-out predictions once under this scheme. The headline membrane score averages the 13 membrane proteins, with assays first averaged within protein. The 26 other proteins form the secondary evaluation. Historical project-wide exposure remains separate from this within-run separation: this is cross-validation on previously used development data, not a newly untouched confirmation set.','',
        '### The 13 membrane proteins in the headline table','',
        '| Dataset label | Accession | Measurements | Assays | Group |',
        '|---|---|---:|---:|---|']
    for r in membrane_proteins:lines.append(f"| {r['gene']} | {r['accession']} | {r['measurements']:,} | {r['assays']} | {r['dependence_group']} |")
    lines.extend(['','## 3. What was outside training and what is left','',
        'The historical 89-protein inventory can be accounted for without treating every unfitted protein as untouched:', '',
        '| Status in the historical inventory | Count | Meaning now |',
        '|---|---:|---|',
        '| Current development/cross-validation pool | 39 | Already used in the standalone results. |',
        '| Previously reserved, evaluated with an older model | 27 | Seven membrane and 20 other proteins. Not trained into the current 39-protein pool, but outcomes were already used in the project. Available for an explicitly exploratory standalone expansion. |',
        '| Previously reserved, skipped by the older evaluation | 4 | Sequence/embedding or mapping failures; not four ready benchmark additions. |',
        '| Pending confirmation candidate role | 6 | Not in the current 39-protein pool. Role name does not certify usable data or untouched outcomes. |',
        '| Other historical exclusions | 13 | Mapping, coverage, relatedness or other curation barriers; require reassessment before admission. |',
        '| Total historical inventory | 89 | This is not an exhaustive census of all public DMS data. |','',
        'The seven earlier membrane evaluation proteins are **NPC1, SLC22A1 (S22A1), HMGCR (HMDH), APP (A4), CD19, ADRB2 and KCNQ4**. The current standalone benchmark has not evaluated those seven. The 20 usable older non-membrane evaluation proteins likewise remain outside the current standalone benchmark. The four skipped accessions are HECD1, Q53Z42, KRAS (RASK) and alpha-synuclein (SYUA). None of this transfers the older model\'s performance or preregistered status to the standalone model.','',
        'The six pending-role proteins are **MC4R, P2RY8, SLC13A5, TECR, LARGE1 and FKRP**. The older local curation has no usable DMS score files for the first three; that is not a claim that no newer public assay could exist. LARGE1/FKRP have recovered and mapped assays with prior QC inspection, but are not integrated into the standalone comparison. TECR has a recovered table but unresolved assay-provenance/quality concerns. These six are therefore not six ready, untouched test datasets.','',
        'The other 13 exclusions are: '+', '.join(sorted(r['gene'] for r in historical_rows if r['current_status']=='other_historical_exclusion_not_in_current_pool'))+'.','',
        '**Important exposure-record conflict:** a later option-C protocol calls LDLR, SGCB and SGCA unopened, but the earlier exposure register explicitly lists all three as already exposed in Phase C. The disulfide-status review also records LDLR\'s historical exposure. This accounting does not certify those three as fresh confirmation. They may still be useful for an expanded, transparently exploratory comparison. A genuinely untouched test requires a reconciled exposure ledger and a fixed evaluation plan, not just a `pending` or `confirmation` label.','',
        '## 4. What the five-fold sensitivity check did','',
        'It repeated the same prediction/correction procedure on the same data using **five larger held-out blocks**, with whole dependence groups kept together. For each round, the other four blocks supplied training data; training-fold PCA, the regression head and inner three-fold correction calibration were fitted again. The held-out block was then scored. After five rounds, every measurement had one outer held-out prediction. This was not a random variant split, new data, a test of the full-fit model, or five retrainings on all 39 proteins.','',
        '| Fold | Training proteins | Held-out proteins | Held-out membrane proteins | Training measurements | Held-out measurements |',
        '|---|---:|---:|---:|---:|---:|'] )
    for f in folds:lines.append(f"| {f['fold']} | {f['train_proteins']} | {f['test_proteins']} | {f['test_membrane_proteins']} | {f['train_measurements']:,} | {f['test_measurements']:,} |")
    lines.extend(['','This asks whether conclusions depend strongly on holding out one group at a time versus holding out larger blocks and having less training data. The settings and data remain fixed; only the outer grouping scheme changes, with corresponding refits.','',
        '| Model | Guarded group-out membrane Spearman | Five-fold membrane Spearman |',
        '|---|---:|---:|'])
    for r in comparison:lines.append(f"| {r['model']} | {r['guarded_loocv']:.6f} | {r['grouped_cv']:.6f} |")
    lines.extend(['','The tiny PI-plus-disulfide advantage over raw does not persist under five-fold evaluation, and the ordering relative to AlphaMissense/VespaG reverses. External baseline score files remain unchanged; the standalone fits change. The check therefore limits a broad superiority claim.','',
        '**Disulfide coverage caveat:** all measured disulfide-supported examples—OPSD, CCR5 and SC6A4—land in five-fold test fold 3. The remaining training/calibration data have no eligible disulfide-feature variation, so that fold\'s learned disulfide coefficient is zero. Other folds have no eligible disulfide examples to change at test time. Thus the direct disulfide term changes zero measured held-out rows in this five-fold run. It is not a direct test of a nonzero disulfide correction\'s predictive benefit. Guarded group-out separates its two supported dependence groups and yields nonzero held-out correction terms.','',
        '## 5. The rule audits are a separate evaluation','',
        'The label-free standalone rule audits score all 193,515 possible missense substitutions at 10,185 positions in the existing 13 membrane proteins. These are model predictions, not 193,515 additional experimental training labels. The positive-inside primary contrast has 596 contributing positions, 11 proteins and 10 dependence groups. Neither those audits nor the annotation-defined disulfide check create a fresh held-out DMS test set.','',
        '## Complete inventories and sources','',
        'The accompanying folder `results/standalone_dataset_accounting_20261028/` contains the exact 39-protein list, 44-assay list, five-fold assignments, all 89 historical role statuses, three additional exposure conflicts, and the scheme comparison values. Gene labels are retained exactly as stored in the dataset.','',
        'Source records include the current benchmark plan/results, `features/l1_table.npz` metadata only, the frozen historical role manifest, the old evaluation accession list, the exposure-history register and the option-C protocol. No new model was fitted, new performance benchmark run, new experimental target read or reserved dataset fetched.'])
    text='\n'.join(lines)+'\n';(OUT/'REPORT.md').write_text(text,encoding='utf-8')
    for filename,digest in SOURCES.items():
        if hashlib.sha256((ROOT/filename).read_bytes()).hexdigest()!=digest:raise AssertionError('Source changed')
    (OUT/'accounting.json').write_text(json.dumps(dict(created_utc=datetime.now(timezone.utc).isoformat(),
        current_measurements=len(accession),current_assays=len(assays),current_proteins=len(proteins),
        membrane_measurements=int((membrane==1).sum()),non_membrane_measurements=int((membrane==0).sum()),
        historical_status_counts=dict(counts),five_fold=folds,evaluation_scheme_means=comparison,
        no_new_training=True,no_new_benchmark=True,no_new_experimental_targets_read=True,
        additional_exposure_conflicts=exposure_conflicts,source_sha256=SOURCES),indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(report=str(DOC),current_proteins=len(proteins),current_assays=len(assays),
                         historical_status_counts=dict(counts),five_fold=folds),indent=2))


if __name__=='__main__':main()

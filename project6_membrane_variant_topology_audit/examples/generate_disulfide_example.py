"""Generate synthetic matched-site inputs for demonstrating the disulfide audit."""
import csv
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'examples/disulfide_toy'


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    sequences={}
    sites=[]
    scores=[]
    for i in range(6):
        accession=f'SYNTHETIC_{i}'
        sequence='ACACACACAC'
        sequences[accession]=sequence
        for pos,state,rsa in [(2,'bonded',0.2),(4,'bonded',0.3),(6,'free',0.22),(8,'free',0.32)]:
            sites.append({'accession':accession,'position':pos,'native':'C',
                'dataset_kind':'synthetic',
                'sequence_sha256':hashlib.sha256(sequence.encode()).hexdigest(),
                'sequence_version':'synthetic v1','protein_family':f'SYNTHETIC_FAMILY_{i}',
                'role':'confirmation','exposure_status':'unseen_and_homology_checked',
                'bond_state':state,'bond_partner':4 if pos==2 else 2 if pos==4 else '',
                'bond_role':'structural','evidence_type':'experimental',
                'evidence_source':'SYNTHETIC known-answer fixture; not an experiment',
                'native_context_verified':True,'mapping_verified':True,
                'control_reviewed':True,'metal_ligand':False,'signal_peptide':False,
                'other_functional_cysteine':False,'redox_or_catalytic':False,
                'domain':'synthetic_domain','compartment':'extracellular',
                'secondary_structure':'C','rsa':rsa,'context_source':'SYNTHETIC geometry'})
            for alternate in ('S','A'):
                scores.append({'accession':accession,'position':pos,'reference':'C',
                    'alternate':alternate,'score':i*10+(0.3 if state=='bonded' else 0)})
    for name,rows in [('sites.csv',sites),('scores.csv',scores)]:
        with (OUT/name).open('w',encoding='utf-8',newline='') as fh:
            writer=csv.DictWriter(fh,fieldnames=list(rows[0]))
            writer.writeheader();writer.writerows(rows)
    (OUT/'sequences.json').write_text(json.dumps(sequences,indent=2),encoding='utf-8')
    (OUT/'expected.json').write_text(json.dumps({'data_kind':'synthetic',
        'real_protein_or_predictor_data':False,'families':6,'pairs':12,
        'expected_C>S':0.3,'expected_C>A':0.3},indent=2),encoding='utf-8')
    print('Wrote synthetic disulfide fixture:',OUT)


if __name__=='__main__': main()

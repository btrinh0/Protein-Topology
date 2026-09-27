"""Independent pandas reconstruction of the expanded positive-inside results."""
from __future__ import annotations

import argparse
import io
import json
from pathlib import Path
import re
import zipfile

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
ALPHABET='ACDEFGHIKLMNPQRSTVWY'
NEUTRAL=[ALPHABET.index(a) for a in 'ACFGILMNPQSTVWY']
CHARGED=[ALPHABET.index(a) for a in 'KR']
SEED,DRAWS=20261003,1999


def source_native_check(directory):
    """Check the source sequences and matrix reference diagonals independently."""
    from Bio import SeqIO
    data=ROOT/'data/external/gemme_membrane_census_20261031'
    with zipfile.ZipFile(data/'vespag_fasta_training_datasets.zip') as archive:
        text=archive.read('vespag_fasta_training_datasets/seq_human5k.fasta').decode()
    source={r.id:str(r.seq) for r in SeqIO.parse(io.StringIO(text),'fasta')}
    topology=(ROOT/'data/raw/topology_viewer.html').read_text(encoding='utf-8-sig')
    match=re.search(r'const PROTEINS_INLINE = (\[.*?\]);',topology,flags=re.S)
    canonical={r['uniprot']:r['sequence'] for r in json.loads(match[1])}
    protocol=json.loads((directory/'PROTOCOL.json').read_text())
    files={r['accession']:r for r in json.loads((data/'download_manifest.json').read_text())['entries']}
    missing=zero=nonzero=positions=0
    for accession in protocol['included_accessions']:
        assert source[accession]==canonical[accession]
        frame=pd.read_csv(ROOT/files[accession]['path'],sep=r'\s+',index_col=0)
        assert frame.shape==(20,len(source[accession]))
        matrix=frame.loc[list(ALPHABET.lower())].to_numpy().T
        ids=np.array([ALPHABET.find(a) for a in source[accession]])
        valid=ids>=0
        diagonal=matrix[np.flatnonzero(valid),ids[valid]]
        missing+=int(np.isnan(diagonal).sum())
        zero+=int((np.isfinite(diagonal)&(np.abs(diagonal)<=1e-12)).sum())
        nonzero+=int((np.isfinite(diagonal)&(np.abs(diagonal)>1e-12)).sum())
        positions+=len(diagonal)
    assert nonzero==0,'Source matrix reference residues disagree with the source sequences'
    report=dict(proteins_checked=len(protocol['included_accessions']),canonical_sequences_exact=True,
        positions_checked=positions,source_native_missing=missing,source_native_zero=zero,
        source_native_nonzero=nonzero,reference_diagonal_check_passed=True)
    (directory/'SOURCE_NATIVE_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


def verify(directory):
    result=pd.read_csv(directory/'results.csv')
    checks=[];per_protein=[];saved_meta={}
    for row in result.to_dict('records'):
        source=row.get('input_npz')
        path=ROOT/source if isinstance(source,str) else directory/f"{row['model']}__{row['scope']}__{row['coverage']}.npz"
        with np.load(path,allow_pickle=False) as z:
            assert ''.join(z['aa_alphabet'])==ALPHABET
            names=['protein','position','side','ref_aa','alt_aa','distance']
            fields={k:z[k].copy() for k in names}
            fields['score']=z['score'].copy()
            accessions=z['accessions'].copy()
        f=pd.DataFrame(fields)
        assert (f.ref_aa!=f.alt_aa).all() and f.score.notna().all()
        assert f.side.isin([0,1]).all() and f.distance.between(1,10).all()
        assert not f.duplicated(['protein','position','alt_aa']).any()
        shared_key=(row['scope'],row['coverage'])
        metadata=f.drop(columns='score')
        if shared_key in saved_meta:
            pd.testing.assert_frame_equal(saved_meta[shared_key],metadata)
        else:
            saved_meta[shared_key]=metadata
        f=f[f.ref_aa.isin(NEUTRAL)]
        keys=['protein','position','ref_aa','side']
        kr=f[f.alt_aa.isin(CHARGED)].groupby(keys).score.mean()
        neutral=f[f.alt_aa.isin(NEUTRAL)].groupby(keys).score.mean()
        sites=(kr-neutral).dropna().reset_index(name='penalty')
        assert len(sites)==row['n_eligible_sites']
        grouped=sites.groupby(['protein','ref_aa','side']).penalty.agg(['mean','size']).unstack('side').dropna()
        weight=grouped['size'][0]*grouped['size'][1]/(grouped['size'][0]+grouped['size'][1])
        differences=grouped['mean'][0]-grouped['mean'][1]
        parts=pd.DataFrame({'num':weight*differences,'den':weight}).groupby('protein').sum()
        point=float(parts.num.sum()/parts.den.sum())
        assert abs(point-row['estimate'])<1e-12
        assert len(grouped)==row['n_strata']
        assert len(parts)==row['n_proteins']
        assert int(grouped['size'].sum().sum())==row['n_contributing_sites']
        rng=np.random.default_rng(SEED);samples=[];n=len(parts)
        for start in range(0,DRAWS,64):
            counts=rng.multinomial(n,np.full(n,1/n),size=min(64,DRAWS-start))
            samples.extend((counts@parts.num.to_numpy()/(counts@parts.den.to_numpy())).tolist())
        ci=np.percentile(samples,[2.5,97.5])
        assert np.allclose(ci,[row['ci_low'],row['ci_high']],atol=1e-12,rtol=0)
        for pid,p in parts.iterrows():
            per_protein.append(dict(model=row['model'],scope=row['scope'],coverage=row['coverage'],
                accession=str(accessions[int(pid)]),estimate=float(p.num/p.den),
                matching_weight=float(p.den),role='descriptive protein point estimate; no protein-specific interval'))
        checks.append(dict(model=row['model'],scope=row['scope'],coverage=row['coverage'],
                           point_error=abs(point-row['estimate']),interval_error=float(np.max(np.abs(ci-[row['ci_low'],row['ci_high']]))),
                           contribution_counts_match=True))
    report=dict(method='independent pandas grouping; original estimator functions not called',
                rows_verified=len(checks),max_point_error=max(c['point_error'] for c in checks),
                all_estimates_intervals_and_counts_match=True,matched_model_metadata_identical=True,
                checks=checks)
    (directory/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    pd.DataFrame(per_protein).to_csv(directory/'independent_per_protein_results.csv',index=False)
    print(json.dumps({k:v for k,v in report.items() if k!='checks'},indent=2),flush=True)
    if directory.name=='gemme_membrane_census_20261031':
        source_native_check(directory)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('directory',type=Path)
    args=parser.parse_args()
    verify(args.directory)


if __name__=='__main__':
    main()

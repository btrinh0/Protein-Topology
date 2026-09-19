"""Audit positive-inside contrasts in the authors' public full-protein GEMME scores."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
from pathlib import Path
import struct
import sys
import zipfile
import zlib

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
sys.path.insert(0,str(ROOT/'scripts/beat_am'))
sys.path.insert(0,str(ROOT/'tools'))
from route_a_common import load_census, census_flanks
import route_a_common
from membrane_rule_audit import panel
from m10_fetch_external_baselines import request

DATA=ROOT/'data/external/gemme_membrane_census_20261031'
OUT=ROOT/'results/gemme_membrane_census_20261031'
URL='https://zenodo.org/records/11085958/files/gemme_rawpreds_by_training_dataset.zip?download=1'
AM=ROOT/'data/processed/gate0/am_full_flank_scores.npz'
IMPLEMENTATION_SHA='c04805510465606b5590275f646d4a5f387e3b5521f4e848157ff26eccf71117'
SEED,DRAWS,WINDOW=20261003,1999,10


def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def sequences():
    with zipfile.ZipFile(DATA/'vespag_fasta_training_datasets.zip') as z:
        text=z.read('vespag_fasta_training_datasets/seq_human5k.fasta').decode()
    records={};key=None
    for line in text.splitlines():
        if line.startswith('>'):
            key=line[1:].split()[0]
            assert key not in records
            records[key]=''
        elif key:
            records[key]+=line.strip()
    return records


def prepare():
    OUT.mkdir(parents=True,exist_ok=True)
    assert panel.sha256(Path(panel.__file__))==IMPLEMENTATION_SHA
    census=load_census();seq=sequences()
    members=json.loads((DATA/'candidate_members.json').read_text())
    eligible=[];excluded=[]
    for row in sorted(members,key=lambda r:r['accession']):
        a=row['accession']
        if a not in seq:
            excluded.append(dict(accession=a,reason='missing source sequence'))
        elif seq[a]!=census[a]['sequence']:
            excluded.append(dict(accession=a,reason='source and topology sequences differ'))
        else:
            eligible.append(row)
    assert len({r['accession'] for r in eligible})==len(eligible)
    sources=[Path(__file__),Path(panel.__file__),Path(route_a_common.__file__),
             DATA/'candidate_members.json',DATA/'vespag_fasta_training_datasets.zip',
             route_a_common.CENSUS,AM,ROOT/'scripts/extract_am_flank_scores.py',
             ROOT/'scripts/run_gate1.py']
    protocol=dict(study='label-free positive-inside audit of public full-protein human GEMME predictions',
        source_url=URL,source_dataset='human/ normalized GEMME evolCombi predictions, Zenodo 11085958',
        cohort='every source protein matching the membrane census sequence exactly',
        candidate_proteins=len(members),exact_sequence_proteins=len(eligible),
        included_accessions=[r['accession'] for r in eligible],
        selection_uses_scores=False,excluded_sequences=excluded,
        raw_score_direction='source higher=fitness; negate once for higher=damage',
        primary=dict(window=WINDOW,comparator=panel.NEUTRAL,availability='>=1 K/R and >=1 neutral alternative',
                     matching='protein/native, both membrane sides',expected_sign='negative',
                     bootstrap='whole protein',draws=DRAWS,seed=SEED),
        sensitivity='complete 19-substitution positions',
        comparison='AlphaMissense and GEMME use identical finite mutation cells on exact same positions',
        fits=0,experimental_outcomes_read=False,
        implementation_sha256=IMPLEMENTATION_SHA,
        source_sha256={str(p.relative_to(ROOT)):panel.sha256(p) for p in sources})
    path=OUT/'PROTOCOL.json'
    if path.exists():
        assert json.loads(path.read_text())==protocol,'Protocol or input changed'
    else:
        save(path,protocol)
    pd.DataFrame(excluded).to_csv(OUT/'sequence_exclusions.csv',index=False)
    print(json.dumps(dict(candidates=len(members),exact_sequence_matches=len(eligible),excluded=len(excluded))),flush=True)
    return census,eligible,protocol


def read_range(start,end):
    with request(URL,headers={'Range':f'bytes={start}-{end}'}) as response:
        if response.status!=206:
            raise ValueError('Byte range not honored; full download refused')
        content_range=response.headers.get('Content-Range','')
        if not content_range.startswith(f'bytes {start}-{end}/'):
            raise ValueError('Wrong byte range returned')
        value=response.read(end-start+2)
    assert len(value)==end-start+1
    return value


def blocks(members):
    pending=[];start=end=None
    for row in sorted(members,key=lambda r:r['offset']):
        lo=row['offset'];hi=lo+30+len(row['member'].encode())+row['compressed']+512
        if pending and (lo>end+131072 or hi-start>4*1024*1024):
            yield start,end,pending
            pending=[];start=end=None
        if not pending:
            start=lo
        pending.append(row);end=hi
    if pending:
        yield start,end,pending


def fetch_block(block):
    start,end,members=block
    blob=read_range(start,end)
    result=[]
    for row in members:
        relative=row['offset']-start
        header=struct.unpack_from('<4s5H3I2H',blob,relative)
        assert header[0]==b'PK\x03\x04'
        filename_length,extra_length=header[-2:]
        data_start=relative+30+filename_length+extra_length
        filename=blob[relative+30:relative+30+filename_length].decode()
        assert filename==row['member']
        raw=blob[data_start:data_start+row['compressed']]
        if len(raw)!=row['compressed']:
            raw=read_range(row['offset']+30+filename_length+extra_length,
                           row['offset']+30+filename_length+extra_length+row['compressed']-1)
        data=zlib.decompress(raw,-15) if row['method']==8 else raw
        assert len(data)==row['bytes'] and zlib.crc32(data)==row['crc32']
        path=DATA/'matrices'/(row['accession']+'.txt.gz')
        path.parent.mkdir(exist_ok=True,parents=True)
        with gzip.open(path,'wb') as stream:
            stream.write(data)
        result.append(dict(accession=row['accession'],path=str(path.relative_to(ROOT)),
                           sha256=panel.sha256(path),source_sha256=hashlib.sha256(data).hexdigest(),
                           source_member=row['member'],source_crc32=row['crc32'],source_bytes=len(data)))
    return result,len(blob)


def fetch(members):
    path=DATA/'download_manifest.json'
    manifest=json.loads(path.read_text()) if path.exists() else dict(entries=[],errors=[])
    entries={r['accession']:r for r in manifest['entries']}
    needed=[]
    for row in members:
        a=row['accession']
        if a in entries and (ROOT/entries[a]['path']).exists():
            assert panel.sha256(ROOT/entries[a]['path'])==entries[a]['sha256']
        else:
            needed.append(row)
    groups=list(blocks(needed))
    print(f'GEMME source: {len(needed)} uncached proteins, {len(groups)} byte-range blocks, '
          f'{sum(end-start+1 for start,end,_ in groups):,} transfer bytes',flush=True)
    manifest['errors']=[];total_bytes=0
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures={pool.submit(fetch_block,group):group for group in groups}
        for n,future in enumerate(as_completed(futures),1):
            group=futures[future]
            try:
                new,transfer=future.result();total_bytes+=transfer
                entries.update({r['accession']:r for r in new})
            except Exception as exc:
                manifest['errors'].append(dict(accessions=[r['accession'] for r in group[2]],
                                               error=f'{type(exc).__name__}: {exc}'))
                print('Download block failed:',str(exc),flush=True)
            manifest['entries']=[entries[a] for a in sorted(entries)]
            manifest['source_url']=URL
            manifest['transfer_bytes_this_run']=total_bytes
            manifest['updated_utc']=datetime.now(timezone.utc).isoformat()
            save(path,manifest)
            if n%10==0 or n==len(groups):
                print(f'Cached GEMME blocks {n}/{len(groups)}; proteins {len(entries)}/{len(members)}',flush=True)
    if manifest['errors']:
        raise RuntimeError('Some source blocks failed; rerun to resume')
    assert {r['accession'] for r in members} <= set(entries)
    return manifest


def write_archive(loaded,matrix,path):
    rows,cols=np.where(np.isfinite(matrix))
    arrays={k:loaded[k][rows] for k in ['protein','position','side','distance','tmd_index','ref_aa']}
    np.savez_compressed(path,**arrays,alt_aa=cols,score=matrix[rows,cols],
                        accessions=loaded['accessions'],aa_alphabet=np.array(list(panel.ALPHABET)))


def audit(census,members,protocol):
    manifest=json.loads((DATA/'download_manifest.json').read_text())
    files={r['accession']:r for r in manifest['entries']}
    accessions=np.array([r['accession'] for r in members])
    arrays={k:[] for k in ['protein','position','side','distance','tmd_index','ref_aa','alt_aa','score']}
    inventory=[]
    for pid,a in enumerate(accessions):
        source=files[str(a)];path=ROOT/source['path']
        assert panel.sha256(path)==source['sha256']
        frame=pd.read_csv(path,sep=r'\s+',index_col=0)
        assert set(frame.index)==set(panel.ALPHABET.lower())
        assert list(frame.columns)==['V'+str(i) for i in range(1,len(census[str(a)]['sequence'])+1)]
        score=-frame.loc[list(panel.ALPHABET.lower())].to_numpy(dtype=float).T
        seq=census[str(a)]['sequence']
        native=np.array([panel.ALPHABET.index(n) if n in panel.ALPHABET else -1 for n in seq])
        valid=native>=0
        score[np.flatnonzero(valid),native[valid]]=np.nan
        geometry=census_flanks(census[str(a)])
        positions=np.array(sorted(p for p,g in geometry.items() if 1<=g[1]<=WINDOW and valid[p-1]),dtype=int)
        if not len(positions):
            inventory.append(dict(accession=str(a),gene=census[str(a)]['gene'],n_scored_flank_sites=0,
                                   n_eligible_sites=0,n_contributing_sites=0,n_strata=0))
            continue
        subset=score[positions-1]
        rows,cols=np.where(np.isfinite(subset))
        selected=positions[rows]
        geo=np.array([geometry[int(p)] for p in selected])
        current=dict(protein=np.full(len(rows),pid),position=selected,side=geo[:,0],distance=geo[:,1],
                     tmd_index=geo[:,2],ref_aa=native[selected-1],alt_aa=cols,score=subset[rows,cols])
        for k in arrays:
            arrays[k].append(current[k])
        # Compute eligibility with the frozen estimator and no new fit.
        local=dict(protein=np.full(len(positions),pid),position=positions,
                   side=np.array([geometry[int(p)][0] for p in positions]),
                   ref_aa=native[positions-1],scores=subset)
        contrast=panel.gain_contrast(local)
        parts=panel.stratum_parts(local,contrast['penalty'],contrast['eligible'])
        inventory.append(dict(accession=str(a),gene=census[str(a)]['gene'],
                               n_scored_flank_sites=int(np.isfinite(subset).any(axis=1).sum()),
                               n_eligible_sites=parts['n_sites'],n_contributing_sites=parts['n_contributing_sites'],
                               n_strata=parts['n_strata']))
    full_path=OUT/'GEMME__public_human_census__available_cells.npz'
    np.savez_compressed(full_path,**{k:np.concatenate(v) for k,v in arrays.items()},
                        accessions=accessions,aa_alphabet=np.array(list(panel.ALPHABET)))
    gemme=panel.load_panel(full_path,WINDOW)
    results=[]
    def run(model,scope,coverage,loaded,matrix):
        path=OUT/f'{model}__{scope}__{coverage}.npz'
        if path!=full_path:
            write_archive(loaded,matrix,path)
        report=panel.run_panel_audit(path,OUT/path.stem,model,'higher-is-damaging',WINDOW,DRAWS,SEED)
        r=report['result'];ci=r['ci95_raw']
        verdict=('Expected direction' if ci and ci[1]<0 else
                 'Opposite direction' if ci and ci[0]>0 else 'Inconclusive')
        results.append(dict(model=model,scope=scope,coverage=coverage,input_npz=str(path.relative_to(ROOT)),
            estimate=r['estimate_raw'],ci_low=ci[0] if ci else None,ci_high=ci[1] if ci else None,
            verdict=verdict,n_proteins=r['n_proteins'],n_eligible_sites=r['n_sites'],
            n_contributing_sites=r['n_contributing_sites'],n_strata=r['n_strata']))
    run('GEMME','public_human_census','available_cells',gemme,gemme['scores'])
    complete=np.isfinite(gemme['scores']).sum(axis=1)==19
    run('GEMME','public_human_census','complete19_positions',gemme,
        np.where(complete[:,None],gemme['scores'],np.nan))
    print('Full-protein GEMME audit complete; aligning cached AlphaMissense cells',flush=True)
    am=panel.load_panel(AM,WINDOW)
    index={(str(accessions[int(p)]),int(q)):i for i,(p,q) in enumerate(zip(gemme['protein'],gemme['position']))}
    am_matrix=np.full_like(gemme['scores'],np.nan);matched_positions=0;metadata_mismatches=0
    for i,(pid,pos) in enumerate(zip(am['protein'],am['position'])):
        accession=str(am['accessions'][int(pid)])
        j=index.get((accession,int(pos)))
        if j is None:
            continue
        if any(am[k][i]!=gemme[k][j] for k in ['ref_aa','side','distance','tmd_index']):
            metadata_mismatches+=1
            continue
        am_matrix[j]=am['scores'][i];matched_positions+=1
    common=np.isfinite(am_matrix)&np.isfinite(gemme['scores'])
    g_common=np.where(common,gemme['scores'],np.nan)
    am_matrix=np.where(common,am_matrix,np.nan)
    complete=common.sum(axis=1)==19
    for model,matrix in [('GEMME',g_common),('AlphaMissense',am_matrix)]:
        for coverage in ['available_cells','complete19_positions']:
            selected=matrix if coverage=='available_cells' else np.where(complete[:,None],matrix,np.nan)
            run(model,'identical_GEMME_AM_cells',coverage,gemme,selected)
    write_json=save
    write_json(OUT/'results.json',results)
    pd.DataFrame(results).to_csv(OUT/'results.csv',index=False)
    pd.DataFrame(inventory).to_csv(OUT/'protein_eligibility.csv',index=False)
    assert all(panel.sha256(ROOT/p)==digest for p,digest in protocol['source_sha256'].items())
    save(OUT/'VERIFICATION.json',dict(original_audit_unchanged=True,all_sources_unchanged=True,
        source_matrices=len(members),source_sequences_exact=True,matched_am_positions=matched_positions,
        am_metadata_mismatches=metadata_mismatches,shared_mutation_cells=int(common.sum()),
        identity_substitutions_excluded=True,no_outcome_labels_or_fitting=True,
        protocol_sha256=panel.sha256(OUT/'PROTOCOL.json')))
    lines=['# Positive-inside audit of full-protein public GEMME predictions','',
           '31 October 2026. The authors\' public human prediction archive supplied '
           f'{len(members)} proteins matching the topology sequence exactly. '
           'The unchanged original positive-inside audit was applied with 10-residue flanks, '
           'the uncharged comparator and 1,999 protein-bootstrap replicates (seed 20261003). '
           'Negative is the expected direction.','',
           '| Model | Cohort | Coverage | Delta | 95% protein interval | Contributing proteins | Contributing sites | Result |',
           '|---|---|---|---:|---|---:|---:|---|']
    for r in results:
        ci=f"[{r['ci_low']:+.6f}, {r['ci_high']:+.6f}]" if r['ci_low'] is not None else 'Not estimable'
        lines.append(f"| {r['model']} | {r['scope']} | {r['coverage']} | {r['estimate']:+.6f} | "
                     f"{ci} | {r['n_proteins']} | {r['n_contributing_sites']} | {r['verdict']} |")
    lines+=['','The GEMME archive contains normalized evolCombi predictions used by the VespaG authors; '
            'this is a separate score source and protein cohort from the ProteinGym table analysis. '
            'GEMME and AlphaMissense matched rows use the same mutation cells, native residues and '
            'membrane geometry. Cross-model raw magnitudes are not an effect-size comparison. '
            'VespaG was not substituted with GEMME teacher scores and was not extrapolated to this cohort.','',
            'Sources: [authors\' data](https://zenodo.org/records/11085958), '
            '[GEMME](https://doi.org/10.1093/molbev/msz179).','']
    text='\n'.join(lines)
    (OUT/'REPORT.md').write_text(text,encoding='utf-8')
    (ROOT/'docs/GEMME_MEMBRANE_CENSUS_20261031.md').write_text(text,encoding='utf-8')
    print(pd.DataFrame(results).to_string(index=False),flush=True)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--fetch',action='store_true')
    args=parser.parse_args()
    census,members,protocol=prepare()
    if args.prepare_only:
        return
    if args.fetch:
        fetch(members)
    audit(census,members,protocol)


if __name__=='__main__':
    main()

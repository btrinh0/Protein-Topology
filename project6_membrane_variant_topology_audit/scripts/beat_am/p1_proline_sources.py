"""Fetch experimental structural annotations for the proline audit without reading variant
scores.
"""
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
from urllib import request

ROOT=Path(__file__).resolve().parents[2]
RAW=ROOT/'data/external/proline_rule_20261025'
OUT=ROOT/'results/proline_rule_20261025'
PROTOCOL=ROOT/'docs/PROLINE_TM_PROTOCOL_20261025.json'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def download(url,path):
    if not path.exists():
        path.parent.mkdir(parents=True,exist_ok=True)
        req=request.Request(url,headers={'User-Agent':'ProlineContextAudit/1.0'})
        with request.urlopen(req,timeout=50) as response: data=response.read()
        temporary=path.with_suffix(path.suffix+'.partial')
        temporary.write_bytes(data);temporary.replace(path)
    return {'url':url,'path':str(path.resolve()),'sha256':sha(path),'bytes':path.stat().st_size}


def candidates():
    text=(ROOT/'data/raw/topology_viewer.html').read_text(encoding='utf-8-sig')
    census={r['uniprot']:r for r in json.loads(re.search(r'const PROTEINS_INLINE = (\[.*?\]);',text,re.S)[1])}
    registry=(ROOT/'scripts/run_phase_c.py').read_text(encoding='utf-8').split('SCORE_SETS = {',1)[1].split('\n}',1)[0]
    result={acc:{'accession':acc,'gene':gene,'source':'previous Phase C registry'}
            for gene,urn,acc in re.findall(r'"([A-Z0-9]+)":\s*\("([^"]+)",\s*"([A-Z0-9]+)"',registry)}
    manifest=list(csv.DictReader((ROOT/'results/beat_am_20261016/inventory_frozen/group_manifest_frozen.csv').open(encoding='utf-8')))
    for row in manifest:
        if row['was_fitted_in_l1']=='True' and row['accession'] in census:
            result.setdefault(row['accession'],{'accession':row['accession'],'gene':row['gene'],'source':'fitted membrane pool'})
    # Amendment 02 expands the structural audit through deterministic, outcome-independent protein selection.
    links=RAW/'human_reviewed_pdb_links.tsv'
    if links.exists():
        with links.open(encoding='utf-8',newline='') as fh:
            linked={r['Entry'] for r in csv.DictReader(fh,delimiter='\t') if r.get('PDB','').strip()}
        eligible=[acc for acc,rec in census.items() if acc in linked and acc not in result and len(rec['tmds'])>=4]
        eligible.sort(key=lambda acc:hashlib.sha256(f'20261025:{acc}'.encode()).hexdigest())
        for acc in eligible[:24]:
            result[acc]={'accession':acc,'gene':census[acc]['gene'],
                         'source':'amendment_02_deterministic_label_free_panel'}
    for acc,info in result.items():
        info['sequence']=census[acc]['sequence']
        info['tmds']=census[acc]['tmds']
        info['loops']=census[acc]['loops']
        info['frozen_group']=next((r['frozen_group'] for r in manifest if r['accession']==acc),None)
    return result


def prepare_one(info):
    acc=info['accession']
    record=download(f'https://rest.uniprot.org/uniprotkb/{acc}.json',RAW/'uniprot'/f'{acc}.json')
    data=json.loads(Path(record['path']).read_text(encoding='utf-8'))
    # Pin the sequence used by the prediction cache and report current UniProt disagreements.
    agreement=data['sequence']['value']==info['sequence']
    ranked=[]
    for ref in data.get('uniProtKBCrossReferences',[]):
        if ref['database']!='PDB': continue
        prop={p['key']:p['value'] for p in ref.get('properties',[])}
        if prop.get('Method') not in {'X-ray','EM'}: continue
        try: resolution=float(prop.get('Resolution','').split()[0])
        except (ValueError,IndexError): continue
        if resolution>3.5: continue
        for chains,start,end in re.findall(r'([A-Za-z0-9/]+)=(\d+)-(\d+)',prop.get('Chains','')):
            for chain in chains.split('/'):
                ranked.append({'pdb_id':ref['id'].lower(),'chain':chain,'method':prop['Method'],
                    'resolution':resolution,'range':[int(start),int(end)],
                    'coverage':(int(end)-int(start)+1)/len(info['sequence'])})
    ranked.sort(key=lambda r:(-r['coverage'],r['resolution'],r['pdb_id'],r['chain']))
    selected=[];seen=set()
    for row in ranked:
        if row['pdb_id'] in seen: continue
        seen.add(row['pdb_id']);selected.append(row)
        if len(selected)==3: break
    attempts=[]
    for row in selected:
        pdb=row['pdb_id']
        local=ROOT.parent/'data/raw/structures'/f'{pdb}.cif.gz'
        try:
            if local.exists():
                source={'url':f'https://files.rcsb.org/download/{pdb.upper()}.cif.gz',
                        'path':str(local.resolve()),'sha256':sha(local),'bytes':local.stat().st_size,
                        'existing_workspace_structure':True}
            else:
                source=download(f'https://files.rcsb.org/download/{pdb.upper()}.cif.gz',RAW/'structures'/f'{pdb}.cif.gz')
            # Reject HTML error responses returned in place of structural files.
            with gzip.open(source['path'],'rt',encoding='utf-8') as fh:
                first=fh.readline()
            if not first.startswith('data_'): raise ValueError('Not an mmCIF data block')
            attempts.append({**row,'source':source,'status':'downloaded'})
        except Exception as exc:
            attempts.append({**row,'status':'fetch_failed','error':str(exc)})
    return {**info,'canonical_sequence_matches_current_uniprot':agreement,
            'uniprot_source':record,'ranked_structure_candidates':len(ranked),'attempts':attempts}


def main():
    RAW.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
    frozen=OUT/'protocol_freeze.json'
    if frozen.exists():
        assert json.loads(frozen.read_text())['protocol_sha256']==sha(PROTOCOL),'Frozen protocol changed'
    else:
        frozen.write_text(json.dumps({'protocol_sha256':sha(PROTOCOL),
            'frozen_utc':datetime.now(timezone.utc).isoformat(),'proline_outcome_contrasts_computed':False},indent=2))
    original=RAW/'source_manifest.json'
    archive=RAW/'source_manifest_initial24.json'
    if original.exists() and not archive.exists(): archive.write_bytes(original.read_bytes())
    download('https://rest.uniprot.org/uniprotkb/stream?query=%28organism_id%3A9606%29%20AND%20%28reviewed%3Atrue%29&format=tsv&fields=accession%2Cxref_pdb',RAW/'human_reviewed_pdb_links.tsv')
    pool=candidates();results=[]
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures={executor.submit(prepare_one,row):acc for acc,row in pool.items()}
        for future in as_completed(futures):
            acc=futures[future]
            try:
                result=future.result();results.append(result)
                print(acc,result['gene'],'structure attempts',len(result['attempts']),flush=True)
            except Exception as exc:
                results.append({**pool[acc],'status':'source_failed','error':str(exc)})
                print(acc,'SOURCE FAILED',str(exc),flush=True)
            (RAW/'source_manifest.json').write_text(json.dumps({'protocol_sha256':sha(PROTOCOL),
                'amendment_02_sha256':sha(ROOT/'docs/PROLINE_TM_AMENDMENT_02_20261025.md'),
                'variant_score_files_opened':False,'proteins':sorted(results,key=lambda r:r['accession'])},indent=2),encoding='utf-8')
    print(json.dumps({'candidate_proteins':len(pool),'source_records':len(results)},indent=2))


if __name__=='__main__': main()

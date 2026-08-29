"""Fetch public annotation/identifier data for the frozen disulfide protocol."""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
from urllib import parse, request

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / 'data/external/disulfide_rule_20261025'
ACCESSIONS = ('P00740', 'P06213', 'P01130', 'P02768', 'P08100', 'Q9BYF1')
FIELDS = ('accession','gene_primary','protein_name','sequence','ft_disulfid',
          'ft_domain','ft_signal','ft_topo_dom','ft_transmem','ft_binding',
          'ft_act_site','ft_mod_res','cc_ptm')
MAVEDB_METADATA = ('urn', 'title', 'shortDescription', 'abstractText', 'methodText',
    'numVariants', 'publishedDate', 'primaryPublicationIdentifiers',
    'targetGenes', 'targetGene', 'experiment', 'doiIdentifiers')
BANNED_IDENTIFIERS = ('score','clinvar','classif','pathogen','acmg','clnsig',
                      'gnomad','hgmd','dbsnp','rsid')


def digest(data):
    return hashlib.sha256(data).hexdigest()


def fetch(url, path, *, transform=None, validator=None):
    path = Path(path)
    if path.exists():
        data = path.read_bytes()
        cached = True
    else:
        print('Fetching',url,flush=True)
        req = request.Request(url,headers={'User-Agent':'DisulfideRuleAudit/1.0'})
        with request.urlopen(req,timeout=55) as response:
            data = response.read()
        if transform:
            data = transform(data)
        if validator:
            validator(data)
        path.parent.mkdir(parents=True,exist_ok=True)
        temp = path.with_suffix(path.suffix+'.partial')
        temp.write_bytes(data)
        temp.replace(path)
        cached = False
    if validator:
        validator(data)
    return {'url':url,'path':str(path.resolve()),'sha256':digest(data),
            'bytes':len(data),'cached':cached,
            'retrieved_utc':datetime.now(timezone.utc).isoformat()}


def validate_identifiers(raw):
    header = next(csv.reader(io.StringIO(raw.decode('utf-8-sig'))))
    bad = [c for c in header if any(x in c.lower() for x in BANNED_IDENTIFIERS)]
    if bad or 'hgvs_pro' not in header:
        raise ValueError(f'Refusing identifier response with unexpected columns: {header}')


def metadata_only(raw):
    d=json.loads(raw)
    # The acquisition allowlist excludes summary statistics and score distributions.
    return json.dumps({k:d[k] for k in MAVEDB_METADATA if k in d},
                      ensure_ascii=False,indent=2).encode('utf-8')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,default=DEFAULT_OUT)
    p.add_argument('--skip-bulk',action='store_true')
    p.add_argument('--f9-identifiers',action='store_true')
    args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    records=[]
    failures=[]
    jobs=[]
    if not args.skip_bulk:
        query=parse.urlencode({'query':'(organism_id:9606) AND (reviewed:true)',
                             'format':'tsv','fields':','.join(FIELDS)})
        jobs.append(('https://rest.uniprot.org/uniprotkb/stream?'+query,
                     args.out/'human_reviewed_features.tsv',{}))
    for acc in ACCESSIONS:
        jobs.append((f'https://rest.uniprot.org/uniprotkb/{acc}.json',
                     args.out/'uniprot'/f'{acc}.json',{}))
    if args.f9_identifiers:
        # Request variant identifiers in the MaveDB namespace for the four F9 endpoints.
        for suffix in ('a-1','b-1','c-1','d-1'):
            urn='urn:mavedb:00001200-'+suffix
            api=f'https://api.mavedb.org/api/v1/score-sets/{urn}'
            jobs.append((api,args.out/'mavedb'/f'{urn.replace(":","_")}.metadata.json',
                         {'transform':metadata_only}))
            jobs.append((api+'/variants/data?namespaces=mavedb',
                         args.out/'mavedb'/f'{urn.replace(":","_")}.identifiers.csv',
                         {'validator':validate_identifiers}))
    for url,path,kwargs in jobs:
        try:
            record=fetch(url,path,**kwargs)
            records.append(record)
            print('Saved',path.name,record['bytes'],'bytes',flush=True)
        except Exception as exc:
            failures.append({'url':url,'error':str(exc)})
            print('FAILED',url,str(exc),flush=True)
        (args.out/'source_manifest.json').write_text(json.dumps({
            'outcome_scores_requested':False,'records':records,'failures':failures},
            indent=2),encoding='utf-8')
    print(json.dumps({'sources':len(records),'failures':failures},indent=2))
    return 0 if not failures else 1


if __name__=='__main__':
    raise SystemExit(main())

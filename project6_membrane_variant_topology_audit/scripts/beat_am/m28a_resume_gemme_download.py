"""Resume the frozen GEMME audit with fewer requests and explicit backoff."""
from __future__ import annotations

from datetime import datetime,timezone
import json
from pathlib import Path
import time
import urllib.error
import urllib.request

import m28_gemme_membrane_census as study


def respectful_range(start,end):
    for attempt in range(4):
        req=urllib.request.Request(study.URL,headers={
            'User-Agent':'protein-rule-audit-comparison/1.0',
            'Accept-Encoding':'identity','Range':f'bytes={start}-{end}'})
        try:
            with urllib.request.urlopen(req,timeout=90) as response:
                assert response.status==206
                assert response.headers.get('Content-Range','').startswith(f'bytes {start}-{end}/')
                value=response.read(end-start+2)
            assert len(value)==end-start+1
            return value
        except urllib.error.HTTPError as exc:
            if exc.code not in [429,500,502,503,504] or attempt==3:
                raise
            header=exc.headers.get('Retry-After','')
            delay=max(60,int(header)) if header.isdigit() else 60*(attempt+1)
            print(f'Source returned {exc.code}; respecting a {delay}-second backoff',flush=True)
            deadline=time.monotonic()+delay
            while time.monotonic()<deadline:
                time.sleep(min(30,deadline-time.monotonic()))
                if time.monotonic()<deadline:
                    print('Waiting for source rate-limit reset',flush=True)
    raise RuntimeError('HTTP retry limit reached')


def coalesced_blocks(members):
    pending=[];start=end=None
    for row in sorted(members,key=lambda r:r['offset']):
        lo=row['offset'];hi=lo+30+len(row['member'].encode())+row['compressed']+512
        if pending and (lo>end+2*1024**2 or hi-start>16*1024**2):
            yield start,end,pending
            pending=[];start=end=None
        if not pending:
            start=lo
        pending.append(row);end=hi
    if pending:
        yield start,end,pending


def main():
    census,members,protocol=study.prepare()
    path=study.DATA/'download_manifest.json'
    manifest=json.loads(path.read_text())
    entries={r['accession']:r for r in manifest['entries']}
    pending=[]
    for r in members:
        a=r['accession']
        if a in entries and (study.ROOT/entries[a]['path']).exists():
            assert study.panel.sha256(study.ROOT/entries[a]['path'])==entries[a]['sha256']
        else:
            pending.append(r)
    groups=list(coalesced_blocks(pending))
    study.save(study.OUT/'TRANSFER_STRATEGY.json',dict(
        reason='source rate limit on many small byte-range requests',
        change='coalesce requests into <=16 MiB blocks; sequential requests; respect Retry-After',
        analysis_protocol_unchanged=True,source_crc_checks_unchanged=True,
        resume_script_sha256=study.panel.sha256(Path(__file__)),
        remaining_proteins=len(pending),requests=len(groups),
        requested_bytes=sum(end-start+1 for start,end,_ in groups)))
    print(f'Resuming {len(pending)} proteins with {len(groups)} sequential requests',flush=True)
    study.read_range=respectful_range
    manifest['errors']=[];transfer=0
    for i,group in enumerate(groups,1):
        new,size=study.fetch_block(group);transfer+=size
        entries.update({r['accession']:r for r in new})
        manifest['entries']=[entries[a] for a in sorted(entries)]
        manifest['resume_transfer_bytes']=transfer
        manifest['updated_utc']=datetime.now(timezone.utc).isoformat()
        study.save(path,manifest)
        print(f'GEMME transfer {i}/{len(groups)}; cached {len(entries)}/{len(members)} proteins',flush=True)
    assert {r['accession'] for r in members}<=set(entries)
    study.audit(census,members,protocol)


if __name__=='__main__':
    main()

"""Finish the same public GEMME transfer in four bounded concurrent batches."""
from pathlib import Path

import m28_gemme_membrane_census as study
from m28a_resume_gemme_download import respectful_range,coalesced_blocks


def main():
    census,members,protocol=study.prepare()
    study.save(study.OUT/'TRANSFER_STRATEGY_CONCURRENT.json',dict(
        reason='large sequential transfers were slow; four bounded downloads use the existing thread pool',
        change='<=16 MiB coalesced ranges; at most four concurrent requests; explicit Retry-After backoff',
        analysis_protocol_unchanged=True,source_crc_checks_unchanged=True,
        resume_script_sha256=study.panel.sha256(Path(__file__)),
        backoff_script_sha256=study.panel.sha256(Path(__file__).with_name('m28a_resume_gemme_download.py'))))
    study.blocks=coalesced_blocks
    study.read_range=respectful_range
    study.fetch(members)
    study.audit(census,members,protocol)


if __name__=='__main__':
    main()

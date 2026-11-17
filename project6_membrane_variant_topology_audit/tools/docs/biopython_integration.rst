Biopython integration for the positive-inside audit
===================================================

This optional interface reads annotated protein records and existing variant
scores. It measures whether introducing lysine or arginine has a lower
relative predicted penalty on the cytosolic side of a membrane. The interface
uses the package's existing positive-inside estimator.

Installation and offline example
--------------------------------

Use Python 3.10 or later. From the package directory (the source archive root, or tools/ in a checkout)::

    python -m pip install ".[biopython]"
    python examples/biopython_positive_inside.py

The example supplies two synthetic proteins, topology annotations, and a
tab-separated score table. Its expected raw contrast is -0.30, with 16
contributing positions, 8 matched strata, and 2 contributing proteins.
Both proteins have the same planted effect, so the bootstrap interval is
degenerate in this example. All sequences, annotations, and scores are
artificial software fixtures.

The core package requires NumPy. The optional extra adds Biopython 1.85 or
later. Supported minimum and current versions are checked in the submission
verification report. No model, GPU, remote API, or DMS labels are required
to audit a supplied score table.

Annotation contract
-------------------

Records must be Bio.SeqRecord.SeqRecord objects with unique explicit IDs
and known protein sequences. If molecule_type is supplied, it must identify
a protein. The accepted feature types are "transmembrane region" and
"topological domain", as produced by Bio.SeqIO's UniProt XML parser.
FASTA alone supplies no membrane topology.

Relevant features must have exact local contiguous locations and strand
None or +1. Fuzzy, uncertain, compound, remote, out-of-bounds, or overlapping
relevant annotations raise ValueError. Unrelated feature types are ignored.
Biopython's [5:15] location becomes residues 6 through 15 in the audit's
one-based inclusive coordinates.

Topological-domain description qualifiers use the following vocabulary.
Matching ignores case and surrounding whitespace.

* Cytoplasmic, Cytosolic, Inside: cytosolic side.
* Extracellular, Lumenal, Luminal, Periplasmic, Outside: noncytosolic side.
* Other or absent descriptions: unknown, counted and excluded.

A description may be a string or a singleton list. Side is never inferred
by assuming alternating membrane crossings. Only a topological domain
directly adjoining a transmembrane interval contributes sites. Sites within
the selected window of the nearest adjacent membrane edge are retained;
the default window is 10 residues. Transmembrane-core sites are excluded.

Score contract
--------------

Each score row has five fields: record ID, one-based position, reference
residue, alternate residue, and numerical score. IDs must match exactly.
The supplied reference must match the record sequence. Residues must be
uppercase standard amino acids, and reference and alternate must differ.
Sequence isoforms and score numbering must be reconciled before calling
the API.

Scores may be finite or NaN; infinity is rejected. Missing substitutions
can also be omitted. Identical retained duplicate substitutions do not
increase sample size. A finite duplicate supersedes a missing value.
Conflicting finite duplicates within retained flanks raise ValueError
when their absolute difference exceeds 1e-12.

Declare higher-is-damaging or higher-is-tolerated explicitly. The latter
is negated, giving a common higher-is-damaging orientation. Raw contrast
magnitudes inherit the model's score scale.

Calculation and output
----------------------

At each eligible site, the estimator subtracts the mean score of available
uncharged substitutions (ACFGILMNPQSTVWY) from the mean of available K/R
substitutions. Eligibility requires an uncharged reference residue, at
least one finite K/R score, and at least one finite comparator score.

The calculation then subtracts the noncytosolic mean penalty from the
cytosolic mean within each protein/reference-residue stratum. Strata need
eligible sites on both sides. Each difference receives the weight
n_inside*n_outside/(n_inside+n_outside), where n counts eligible sites.
The final estimate is the weighted average across contributing strata.

A negative estimate agrees with the specified positive-inside expectation.
The report returns the estimate, comparator, window, score direction, and
site, stratum, and protein counts. It reports no automatic pass/fail label.
No matched strata yields estimate_raw=None and ci95_raw=None.

The interval resamples whole contributing proteins with replacement,
recalculates the weighted estimate, and returns the 2.5th and 97.5th
percentiles. It requires at least two contributing proteins and two
bootstrap draws. Set draws=0 to disable it; seed controls reproducibility.
The interval describes stability across the supplied proteins and does not
resolve annotation or comparator sensitivity.

Run the integration tests from the package directory after installation::

    python -W error::ResourceWarning -m unittest discover -s tests -p test_biopython_integration.py -v

Scope
-----

The integration covers positive-inside. The package's other audits retain
their existing interfaces and annotation requirements.

Copyright (c) 2026, btrinh0. This documentation is dual licensed under the
Biopython License Agreement or BSD 3-Clause; see ../LICENSE.

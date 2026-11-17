# Synthetic Biopython example

These files contain two artificial proteins and artificial substitution
scores generated from the project's existing positive-inside toy example.
They contain no human participants, patient measurements, clinical labels,
or downloaded predictor outputs.

- synthetic_proteins.xml: minimal UniProt-style XML for SeqIO, with explicit
  transmembrane and topological-domain annotations.
- toy_scores.tsv: record ID, one-based protein substitution, and score.

The expected contrast is -0.30, with 16 contributing sites, 8 strata,
and 2 proteins. The two proteins have identical planted contrasts.
The input fixture predates this integration.

After installing the optional extra, run from the package directory:

    python examples/biopython_positive_inside.py

Copyright (c) 2026, btrinh0. These synthetic fixtures and notes are dual
licensed under the Biopython License Agreement or BSD 3-Clause; see
../../LICENSE.


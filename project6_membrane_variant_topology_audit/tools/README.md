# membrane-rule-audit

Audit biological patterns in protein variant-effect scores using supplied
sequence, membrane-topology, and structural annotations.

## Installation

Requires Python 3.10+ and NumPy. From the package directory, run:

    python -m pip install .

The package directory is project6_membrane_variant_topology_audit/tools in
the GitHub release archive, and the root of the separate source distribution.

Install the optional Biopython interface with:

    python -m pip install ".[biopython]"

## Commands

| Command | Purpose |
|---|---|
| membrane-rule-audit | Charge gain, charge loss, and comparator checks near membrane boundaries |
| membrane-rule-panel | Positive-inside contrast with the uncharged comparator |
| disulfide-rule-audit | Disulfide-loss comparison with matched, verified free-cysteine controls |
| proline-rule-audit | Proline introduction in independently annotated membrane-helix contexts |
| multi-rule-refine | Fit or apply additive score adjustments using supplied rule features |

Each command supports --help. Audits require an explicit score direction:
higher-is-damaging or higher-is-tolerated.

## Offline positive-inside example

From the package directory:

    python examples/generate_toy_example.py --out toy
    membrane-rule-panel --scores toy/toy_mapped_scores.npz --model SyntheticExample --score-direction higher-is-damaging --out toy_panel

The expected contrast is -0.30. Read toy_panel/membrane_rule_panel.json
for the estimate, counts, and interval.

For CSV or TSV input, the general command accepts uniprot_id, protein_variant,
and score columns, with substitutions such as R123W:

    membrane-rule-audit --topology toy/toy_topology.json --scores toy/toy_scores.tsv --score-direction higher-is-damaging --out toy_general

The example generator documents the JSON and mapped NPZ input layouts.
Residue positions are one-based in the supplied protein sequence.

## Biopython records

The optional interface accepts annotated SeqRecords and precomputed missense
scores. Run its offline UniProt XML example after installing the extra:

    python examples/biopython_positive_inside.py

In Python:

    from Bio import SeqIO
    from membrane_rule_audit.biopython import audit_seqrecords

    with open("proteins.xml", "rb") as handle:
        records = list(SeqIO.parse(handle, "uniprot-xml"))

    # Each row is: record ID, one-based position, reference, alternate, score.
    report = audit_seqrecords(
        records, score_rows, score_direction="higher-is-damaging"
    )

Records require exact transmembrane and topological-domain annotations,
with score IDs and reference residues matching the sequence. Unknown sides
are counted and excluded. The SeqRecord interface covers positive-inside.
See [the annotation contract](docs/biopython_integration.rst).

## Disulfide example

From the package directory:

    python examples/generate_disulfide_example.py
    disulfide-rule-audit --sites examples/disulfide_toy/sites.csv --sequences examples/disulfide_toy/sequences.json --scores examples/disulfide_toy/scores.csv --score-direction higher-is-damaging --out toy_disulfide

The synthetic disulfide comparison has an expected contrast of +0.30.
The matched analysis requires evidence that control cysteines are free.

## Tests

With the Biopython extra installed, run:

    python -W error::ResourceWarning -m unittest discover -s tests -p test_biopython_integration.py -v

## Citation and licence

See [CITATION.cff](CITATION.cff) for the software citation.

Copyright (c) 2026, btrinh0. The package and accompanying integration
materials are dual licensed under the Biopython License Agreement or
BSD 3-Clause. See [LICENSE](LICENSE) for the terms and scope.


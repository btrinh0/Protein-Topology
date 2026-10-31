# Protein Topology

Membrane Rule Audit measures biological-rule contrasts in protein variant scores
and evaluates annotation-derived corrections on a standalone predictor.

## Install

Use Python 3.10 or newer and NumPy. Installation directly from GitHub also needs Git.

```sh
python -m pip install "membrane-rule-audit @ git+https://github.com/btrinh0/Protein-Topology.git@main#subdirectory=project6_membrane_variant_topology_audit/tools"
```

From a clone, install the package with:

```sh
python -m pip install ./project6_membrane_variant_topology_audit/tools
```

## Commands

| Command | Purpose |
| --- | --- |
| `membrane-rule-audit` | Charge contrasts near annotated membrane boundaries |
| `membrane-rule-panel` | Within-site positive-inside contrasts from mapped score arrays |
| `disulfide-rule-audit` | Bonded cysteine substitutions against matched, verified free cysteines |
| `proline-rule-audit` | Proline introduction in regular membrane helices against matched controls |
| `multi-rule-refine` | Fit or apply additive rule adjustments to supplied predictions |

Each command provides `--help`. Audits require an explicit score direction:
`higher-is-damaging` or `higher-is-tolerated`.

Run the synthetic positive-inside example from the repository root:

```sh
python project6_membrane_variant_topology_audit/examples/generate_toy_example.py --out toy
membrane-rule-panel --scores toy/toy_mapped_scores.npz --model SyntheticExample --score-direction higher-is-damaging --out toy_panel
```

The expected raw contrast in `toy_panel/membrane_rule_panel.json` is -0.30.
The general audit accepts CSV or TSV columns `uniprot_id`, `protein_variant` and
`score`; variants such as `R123W` use one-based sequence positions. The
[example generator](project6_membrane_variant_topology_audit/examples/generate_toy_example.py)
shows the topology records and mapped NPZ fields. Mapped `side` values are 0 for
cytosolic and 1 for noncytosolic, with alphabet `ACDEFGHIKLMNPQRSTVWY`.

# Membrane Rule Audit

Analyze existing protein variant scores using membrane topology or curated
annotations. The audit reports contrasts, uncertainty intervals when enough
groups contribute, and counts for interpreting each comparison.

## Install

Use Python 3.10 or newer. The GitHub install also requires Git.

```sh
python -m pip install "membrane-rule-audit @ git+https://github.com/btrinh0/Protein-Topology.git@main#subdirectory=project6_membrane_variant_topology_audit/tools"
```

From a clone of the repository, run the following from its root.

```sh
python -m pip install ./project6_membrane_variant_topology_audit/tools
```

## Commands

| Command | Purpose | Inputs |
|---|---|---|
| `membrane-rule-audit` | Charge gain, charge loss and comparator checks near membrane boundaries | Topology JSON or viewer HTML and a variant-score table or NPZ |
| `membrane-rule-panel` | Positive-inside contrast with the all-uncharged comparator | A mapped score NPZ with membrane-side and distance fields |
| `disulfide-rule-audit` | Compare substitutions at annotated bonds with matched, verified free cysteines | Canonical sequences, curated site evidence and variant scores |
| `proline-rule-audit` | Compare proline introduction in regular membrane helices with matched controls | Structural site JSON and a position-by-substitution score NPZ |
| `multi-rule-refine` | Fit or apply additive adjustments to supplied predictions | Prepared calibration data or scores, rule features and coefficients |

Use `--help` on each command for its arguments. For the audits, declare the
input score direction as `higher-is-damaging` or `higher-is-tolerated`.

## Try the positive-inside audit

From the repository root, generate a small set of artificial inputs with a
known contrast of −0.30. The two examples below use different input formats
for the same synthetic scores.

```sh
python project6_membrane_variant_topology_audit/examples/generate_toy_example.py --out toy
membrane-rule-audit --topology toy/toy_topology.json --scores toy/toy_scores.tsv --score-direction higher-is-damaging --out toy_general
membrane-rule-panel --scores toy/toy_mapped_scores.npz --model SyntheticExample --score-direction higher-is-damaging --out toy_panel
```

Read `toy_panel/membrane_rule_panel.json` for the estimate and counts. A negative
contrast means introducing lysine or arginine has a lower relative penalty on
the cytosolic side.

For your own data, the general command accepts CSV or TSV with columns
`uniprot_id`, `protein_variant` and `score`. Write a variant as `R123W`, with a
one-based position in the sequence supplied with the topology. The
[example generator](../examples/generate_toy_example.py) shows the topology
records and mapped NPZ fields. In the mapped NPZ, `side` is 0 for cytosolic and
1 for outer, `distance` is in residues, and the alphabet is
`ACDEFGHIKLMNPQRSTVWY`.

## Try the disulfide audit

The [disulfide example generator](../examples/generate_disulfide_example.py)
creates a synthetic sequence map, site evidence and score table. From the
repository root, run the following.

```sh
python project6_membrane_variant_topology_audit/examples/generate_disulfide_example.py
disulfide-rule-audit --sites project6_membrane_variant_topology_audit/examples/disulfide_toy/sites.csv --sequences project6_membrane_variant_topology_audit/examples/disulfide_toy/sequences.json --scores project6_membrane_variant_topology_audit/examples/disulfide_toy/scores.csv --score-direction higher-is-damaging --out toy_disulfide
```

The expected contrast is +0.30 for each tested substitution. For real data,
the matched analysis requires evidence that the control cysteines are free.

## Python API

```python
from pathlib import Path

from membrane_rule_audit import run_panel_audit

summary = run_panel_audit(
    Path("toy/toy_mapped_scores.npz"),
    Path("toy_panel"),
    model="SyntheticExample",
    score_direction="higher-is-damaging",
)
print(summary["result"]["estimate_raw"])
```

The package takes variant scores and annotations as inputs. Audit results
describe patterns in those scores on the supplied sites and comparators.
The refinement command applies supplied coefficients or fits them from a
calibration table.

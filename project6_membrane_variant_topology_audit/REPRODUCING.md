# Reproducing the package and paper checks

Run commands from `project6_membrane_variant_topology_audit/`.

## Environments

The current package requires Python 3.10 or newer and NumPy. The recorded
analysis environment uses Python 3.12.13 and the observed versions in
`requirements-analysis.txt`. The embedding requirements are separate; exact
versions of the original GPU environment were not recovered.

```sh
python -m pip install -r requirements-analysis.txt
```

Use the current package for its public CLI:

```sh
python -m pip install ./tools
membrane-rule-audit --help
membrane-rule-panel --help
```

Form A/B intentionally loads the earlier validated panel from the installed
distribution. Use a separate environment for those experiments and the full
research test suite. The original wheel and release manifest are included under
`provenance/validated_release/` without modification.

Restore their historical paths with PowerShell:

```powershell
$releasePath = 'dist/membrane-rule-audit-v0.1.0-validated-20261004'
New-Item -ItemType Directory -Force -Path "$releasePath/wheel" | Out-Null
Copy-Item -LiteralPath 'provenance/validated_release/membrane_rule_audit-0.1.0-py3-none-any.whl' -Destination "$releasePath/wheel/"
Copy-Item -LiteralPath 'provenance/validated_release/RELEASE_MANIFEST.json' -Destination "$releasePath/"
python -m pip install --no-deps --force-reinstall provenance/validated_release/membrane_rule_audit-0.1.0-py3-none-any.whl
```

On other platforms, create the same directory and copy the same two files to
those locations before installing the wheel. The manifest pins the installed
validated modules; the tests import current working-tree APIs separately.
The release loader excludes working-tree build metadata when locating the
installed validated panel.

## Software and saved-result checks

```sh
python -m unittest discover -s tests -p "test_*.py" -v
python scripts/paper_support/recheck_final_paper.py --out results/reproduction_saved_checks.json
```

The three published-model adapter tests require the assets listed in
`REQUIRED_INPUTS.csv`; they are skipped when the cached ESM2 input is absent.
Tests use synthetic data except for explicitly asset-dependent checks.

For the added score-array checks, restore the original data and provide the
project root containing those inputs:

```sh
python scripts/paper_support/recheck_final_paper.py --project-root PATH_TO_RESTORED_PROJECT --score-audits --out results/reproduction_score_checks.json
```

The worked example reads the restored standalone position-score archive and
annotation bank, and writes to a fresh destination:

```sh
python scripts/paper_support/export_worked_example.py --out results/worked_example
```

## Source downloads and frozen settings

`REQUIRED_INPUTS.csv` identifies each missing input, its producing script or
source, byte counts and available hashes. Download manifests under `data/` pin
source URLs and identifiers. Keep the included group manifests, assay lists,
model settings and benchmark plans when checking the paper's exact cohorts.

The original acquisition entry points are:

- `route_a_fetch_proteingym.py` and `route_a_resolve_uniprot.py` for the pinned assay release and accession mapping.
- `beat_am/c1_disulfide_sources.py` and `beat_am/p1_proline_sources.py` for annotation sources.
- `beat_am/m10_fetch_external_baselines.py` for matched public predictor scores.
- `beat_am/m12_prepare_speed_assets.py` and `beat_am/m14_fetch_cached_esm2.py` for published-model timing assets.
- `beat_am/m28_gemme_membrane_census.py` and its recovery helpers for the frozen GEMME candidate set.

These paths are relative to `scripts/`. Some historical runners have fixed Linux
paths, protected plans or fixed output destinations; follow `RUN_GUIDE.md` before
running them. Restoring a newer source release is not equivalent to restoring a
pinned snapshot. Respect each source distribution's license and access terms.

## Timing provenance

`provenance/runner_at_measurement_20261026.py` is an immutable copy of the actual
measured runner, with SHA-256
`6c1e5f4a42afa5008b5b9465262ff800636f558fc43a5d13fe6120b0a6e4c1ac`.
It assumes its original location under `results/speed_comparison_20261026/` and
is retained for historical reproduction rather than as a current entry point.

The maintained timing runner has a shorter header. The report finalizer verifies
the retained measured source directly, preserves timing samples and historical
test evidence, and records the publication sources separately. Regenerate reports
only against copied matching measurement outputs. New timings require a fresh
measurement record and their actual environment information.

The published-model adapters explicitly construct float32 heads so callers that
change PyTorch's default dtype do not create layers incompatible with the cached
float32 inputs. The published-weight checks cover this interaction; the original
timing samples remain unchanged.

## Known gaps

The final Figure 2 TikZ source and original external comparator-sensitivity
resampling runner/settings were not recovered. The saved comparator intervals
are retained, while the added checker reconstructs point estimates. The original
GPU and typesetting environments are incompletely recorded.

Historical records may contain source hashes from earlier code states. Current
publication hashes are recorded separately; old plans and results retain their
original hashes. The [archive record](docs/ARCHIVED_EXPERIMENTS.md) describes the
completed hold-out and topology-feature experiments and their retained exposure
history.

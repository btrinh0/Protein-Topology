# CPU speed benchmark with cached inputs

Date: 26 October 2026 (America/Los_Angeles).

Scope selected by the user: prediction with cached inputs, plus lookup for every model. No foundation-model encoder or MSA-search runtime is included. No new training or coefficient fitting was performed.

The main workload is the same 43,991 membrane-panel measurement rows used in the accuracy comparison (13 proteins, 15 assays). Duplicated variant requests across assays remain duplicated for every method. Smaller batches use the same fixed-seed requests.

Model initialization, disk loading, imports and tensor construction are excluded from warm timings. Single-query lookups are repeated within each timed sample to reduce clock noise. Methods are warmed and interleaved in random fixed-seed order; results are medians of seven measurements with 10th/90th percentiles retained in JSON.

## Hardware

```json
{
  "platform": "Windows-11-10.0.26100-SP0",
  "processor": "Intel64 Family 6 Model 151 Stepping 2, GenuineIntel",
  "python": "3.12.13 (main, Mar 10 2026, 18:15:39) [MSC v.1944 64 bit (AMD64)]",
  "logical_cpus": 24,
  "torch": "2.14.0+cpu",
  "numpy": "2.5.3",
  "sklearn": "1.9.1",
  "torch_cuda_available": false,
  "cpu_name": "12th Gen Intel(R) Core(TM) i9-12900KF",
  "physical_ram_gib": 47.795555114746094
}
```

## Prediction after cached embeddings

This includes each measured prediction head and its output selection. For the standalone model, PCA and supplied biochemical/topology features are included; rebuilding those features is also timed separately. VESPA-light includes its published conservation CNN plus its ten logistic-regression models. VespaG uses its published v1 ESM2 head. Cached embeddings have different dimensions and were produced by each method's own encoder.

The existing ProtT5 cache uses canonical full-protein context; the published ESM2 cache uses ProteinGym assay constructs. Queried substitutions are matched through the verified construct alignment, but the encoder contexts are not asserted identical. This is a post-embedding component comparison, not a controlled end-to-end sequence benchmark.

| Threads | Stage / model | Batch size | Median ms | Variants / second |
|---:|---|---:|---:|---:|
| 1 | head_only/Standalone_raw | 43,991 | 853.016 | 51,571 |
| 1 | cached_embedding/Standalone_raw | 43,991 | 846.833 | 51,948 |
| 1 | cached_embedding_rebuild_features/Standalone_raw | 43,991 | 994.224 | 44,247 |
| 1 | cached_embedding/Standalone_PI_disulfide | 43,991 | 864.492 | 50,887 |
| 1 | correction_only/Standalone_PI_disulfide | 43,991 | 0.742 | 59,271,086 |
| 1 | cached_embedding/Standalone_PI_glyco_disulfide | 43,991 | 860.110 | 51,146 |
| 1 | correction_only/Standalone_PI_glyco_disulfide | 43,991 | 0.701 | 62,763,590 |
| 1 | cached_embedding/VESPAl | 43,991 | 190.931 | 230,403 |
| 1 | cached_embedding/VespaG | 43,991 | 29.165 | 1,508,359 |
| 8 | head_only/Standalone_raw | 43,991 | 193.765 | 227,033 |
| 8 | cached_embedding/Standalone_raw | 43,991 | 170.105 | 258,610 |
| 8 | cached_embedding_rebuild_features/Standalone_raw | 43,991 | 307.091 | 143,251 |
| 8 | cached_embedding/Standalone_PI_disulfide | 43,991 | 196.872 | 223,449 |
| 8 | correction_only/Standalone_PI_disulfide | 43,991 | 0.894 | 49,228,961 |
| 8 | cached_embedding/Standalone_PI_glyco_disulfide | 43,991 | 170.953 | 257,327 |
| 8 | correction_only/Standalone_PI_glyco_disulfide | 43,991 | 0.790 | 55,705,964 |
| 8 | cached_embedding/VESPAl | 43,991 | 74.387 | 591,379 |
| 8 | cached_embedding/VespaG | 43,991 | 11.145 | 3,947,045 |

## Lookup for every model

All models use the same in-memory indexed-array retrieval implementation on the same requests. These timings measure cached lookup, not original model inference; differences here are memory/access noise and do not establish a model speed advantage.

| Threads | Model | Requests | Median ms |
|---:|---|---:|---:|
| 1 | Standalone_raw | 43,991 | 0.127 |
| 1 | Standalone_PI_disulfide | 43,991 | 0.120 |
| 1 | Standalone_PI_glyco_disulfide | 43,991 | 0.128 |
| 1 | AlphaMissense | 43,991 | 0.111 |
| 1 | ESM1b | 43,991 | 0.106 |
| 1 | GEMME | 43,991 | 0.113 |
| 1 | VESPA | 43,991 | 0.121 |
| 1 | VESPAl | 43,991 | 0.125 |
| 1 | VespaG | 43,991 | 0.128 |
| 8 | Standalone_raw | 43,991 | 0.129 |
| 8 | Standalone_PI_disulfide | 43,991 | 0.110 |
| 8 | Standalone_PI_glyco_disulfide | 43,991 | 0.118 |
| 8 | AlphaMissense | 43,991 | 0.109 |
| 8 | ESM1b | 43,991 | 0.133 |
| 8 | GEMME | 43,991 | 0.101 |
| 8 | VESPA | 43,991 | 0.099 |
| 8 | VESPAl | 43,991 | 0.124 |
| 8 | VespaG | 43,991 | 0.102 |

## Stages that could not be measured

- **AlphaMissense fresh inference:** Official trained weights are not released; lookup was timed.
- **ESM1b fresh inference:** Masked language-model forward passes are outside the selected cached-input scope; lookup was timed.
- **GEMME cached-alignment prediction:** The official executable/container is not configured locally and the official download host was unreachable. Lookup was timed; no inferred or literature runtime is substituted.
- **VESPA full cached-input prediction:** No cached ProtT5 mutation log-odds inputs are available. Its complete pipeline is not replaced by VESPA-light or by synthetic feature values; lookup was timed.

## Interpretation

A fast prediction head does not establish a fast complete predictor. ProtT5 generation, ESM2 generation, ESM1b masked-sequence inference, MSA search, annotation acquisition, and loading an entire public proteome archive are outside this selected benchmark. Lookup timings must not be compared to another method's sequence-to-score timing.

The saved standalone head is the existing full-development 18+64-feature bundle, with the same fixed architecture used for the accuracy comparison. Runtime is not measured by refitting the head or mixing in AlphaMissense features. The additional layer's measured coefficients and features come from the existing guarded folds; their effects reproduce the saved refined-minus-raw predictions.

Published model classes are loaded from pinned upstream source definitions with strict state-dictionary matching. No synthetic embeddings or random model weights are timed as model predictions. The LR compatibility check is recorded in verification. Numeric outputs from these runtime adapters are not used to replace the published scores in the accuracy comparison.

## Reproduce

```powershell
.\.venv_ml\Scripts\python.exe scripts\beat_am\m13_local_speed_benchmark.py
```

Sources: [AlphaMissense availability](https://github.com/google-deepmind/alphamissense), [VESPA](https://github.com/Rostlab/VESPA), [VespaG](https://github.com/JSchlensok/VespaG), [published cached ESM2 inputs](https://zenodo.org/records/11085958).

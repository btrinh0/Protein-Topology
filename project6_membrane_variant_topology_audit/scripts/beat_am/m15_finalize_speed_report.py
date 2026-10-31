"""Finalize timing reports and provenance while preserving the measured samples."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(Path(__file__).resolve().parent))
import m13_local_speed_benchmark as speed


def timing_snapshot(expected_sha256):
    """Read the retained timing source and verify its recorded measurement hash."""
    path = ROOT / "tests/fixtures/runner_at_measurement_20261026.py"
    snapshot = path.read_bytes()
    if hashlib.sha256(snapshot).hexdigest() != expected_sha256:
        raise ValueError("Retained timing source differs from the recorded measurement hash")
    return snapshot


def main():
    out = speed.OUT
    result = json.loads((out/"results.json").read_text(encoding="utf-8"))
    source_path = Path(speed.__file__)
    new = 'The official executable/container is not configured locally and the official download host was unreachable. Lookup was timed; no inferred or literature runtime is substituted.'

    expected = result["source_hashes"]["source_code"]
    snapshot = timing_snapshot(expected)
    (out/"runner_at_measurement.py").write_bytes(snapshot)
    result["unavailable_stages"]["GEMME cached-alignment prediction"] = new
    result["presentation_amendment"] = dict(
        changes="Clarified GEMME runtime availability and differing cached encoder contexts; timing samples unchanged.",
        timing_runner_sha256=expected,report_renderer_sha256=speed.sha(source_path))
    speed.dump(out/"results.json",result)
    speed.write_report(result)
    verification = json.loads((out/"verification.json").read_text(encoding="utf-8"))
    verification["publication_source_hashes"] = dict(
        timing_runner=speed.sha(source_path),
        report_finalizer=speed.sha(Path(__file__)),
        timing_tests=speed.sha(ROOT/"tests/test_speed_comparison.py"))
    verification["timing_runner_snapshot_sha256"] = speed.sha(out/"runner_at_measurement.py")
    verification["cached_input_scope_selected_by_user"] = True
    if "thread_libraries" in verification:
        verification["thread_libraries_snapshot_after_timing"] = verification.pop("thread_libraries")
    profiles = {}
    for threads in result["threads"]:
        speed.torch.set_num_threads(threads)
        with speed.threadpool_limits(limits=threads):
            details = speed.threadpool_info()
            if any(item["num_threads"] > threads for item in details):
                raise ValueError("Thread-limit context did not constrain a loaded numerical library")
            profiles[str(threads)] = details
    verification["thread_limit_context_verified_after_run"] = profiles
    speed.dump(out/"verification.json",verification)
    measurements = {(r["threads"],r["method"]):r for r in result["measurements"] if r["rows"] == result["primary_rows"]}
    raw = measurements[1,"cached_embedding/Standalone_raw"]["median_seconds"]
    delta = measurements[1,"correction_only/Standalone_PI_disulfide"]["median_seconds"]
    text = f"""# Current paper contribution and speed interpretation

Date: 26 October 2026 (America/Los_Angeles).

The strongest candidate contribution is an explicit, reproducible biological-rule audit and selective correction framework. It tests defined, context-dependent expectations on a predictor's outputs, reports matched comparisons and uncertainty, and makes correction terms inspectable. The generic use of frozen protein-language-model representations and a small prediction head already exists in VESPA/VespaG; this architecture alone is not the differentiating claim.

The paper can demonstrate the workflow on a standalone DMS-supervised ProtT5 predictor using explicit substitution/topology features. It can report where an expected contrast is present, opposite, inconclusive, or shifted by a correction. Those categories must remain distinct from better experimental prediction accuracy. Direct superiority in biological consistency over GEMME/VespaG has not been established by the accuracy comparison; their score correlations alone do not answer that question.

The current best observed membrane correlation is 0.462596 for inside-out plus disulfide versus 0.462330 raw. This tiny difference has an interval crossing zero. The targeted disulfide MSE reduction of 10.54% versus inside-out alone is exploratory, based on three proteins and two dependence groups. Glycosylation remains a conditional audit/calibration finding; a universal glycan-loss penalty conflicts with some assay objectives, including ACE2 binding. The present evidence does not support a general accuracy-breakthrough claim or three independently validated repairs.

The new speed experiment supports inexpensive correction application: the cached-feature inside-out/disulfide layer takes {delta*1000:.3f} ms for 43,991 requests on one CPU thread, compared with {raw*1000:.3f} ms for the base cached-input prediction pipeline (about {delta/raw*100:.3f}% as an isolated component ratio). This is not the noisy difference between two complete pipeline medians. Annotation-feature construction and the sequence encoder are excluded from that isolated overhead.

VespaG and VESPA-light were faster than the standalone predictor in the measured post-embedding stages. Common cached-score lookup is inexpensive for every model. The paper should therefore emphasize biological diagnostics, controlled refinement, and transparent limits, rather than claim a generally faster predictor. These measurements do not establish full sequence-to-score runtime.

Novelty relative to all prior audit/constraint methods has not been established by these baseline papers alone. A final manuscript should compare the audit definition, controls, context handling, and repair/evaluation design with the closest prior methods before making a first-of-its-kind claim.

Sources: [VESPA](https://doi.org/10.1007/s00439-021-02411-y), [VespaG](https://doi.org/10.1093/bioinformatics/btae621), `results/external_predictor_comparison_20261026/`, `results/standalone_multirule_20261025/label_free_audit_20261025/`, and `results/speed_comparison_20261026/`.
"""
    (out/"CONTRIBUTION.md").write_text(text,encoding="utf-8")
    print("Timing samples unchanged; original source verified; historical test record retained; reports finalized.")


if __name__ == "__main__":
    main()

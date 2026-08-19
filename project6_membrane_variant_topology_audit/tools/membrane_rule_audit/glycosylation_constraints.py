"""Calibrate a label-free lower bound for the glycosylation-loss coefficient."""
from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from .refinement import RULES


def fit_glycosylation_floor(
    constraints: Sequence[Mapping], *, calibration_groups: Sequence[str],
    excluded_groups: Sequence[str] = (), margin: float = 0.0,
) -> dict:
    """Solve a one-dimensional constrained problem on declared groups only.

    ``baseline_delta`` is the protein's matched A or B contrast with the
    glycosylation coefficient set to zero. ``feature_delta`` is the same
    linear contrast applied to the binary loss feature. Other rule terms
    may already be included in baseline_delta and are held fixed.
    """
    allowed, excluded = set(calibration_groups), set(excluded_groups)
    if not allowed or allowed & excluded:
        raise ValueError("Calibration and excluded groups must be nonempty/disjoint")
    if not math.isfinite(margin) or margin < 0:
        raise ValueError("A finite nonnegative directional margin is required")
    if not constraints:
        raise ValueError("No supported glycosylation constraints")
    rows, beta = [], 0.0
    seen = set()
    for source in constraints:
        row = dict(source)
        group, contrast = str(row["group"]), row["contrast"]
        identity = (str(row["accession"]), contrast)
        if group not in allowed or group in excluded:
            raise ValueError("Constraint comes from an undeclared or excluded group")
        if contrast not in ("A", "B") or identity in seen:
            raise ValueError("Unique protein/contrast A or B constraints required")
        seen.add(identity)
        delta, slope = float(row["baseline_delta"]), float(row["feature_delta"])
        if not math.isfinite(delta) or not math.isfinite(slope):
            raise ValueError("Finite score and feature contrasts required")
        direction = 1.0 if contrast == "A" else -1.0
        value, gain = direction * delta, direction * slope
        if gain < -1e-12:
            raise ValueError("Glycosylation feature moves a constraint the wrong way")
        if gain <= 1e-12:
            if value < margin - 1e-12:
                raise ValueError("The feature cannot repair an unsupported constraint")
            lower = 0.0
        else:
            lower = max(0.0, (margin - value) / gain)
        beta = max(beta, lower)
        rows.append({**row, "minimum_beta": lower})
    residuals = [
        (1.0 if r["contrast"] == "A" else -1.0)
        * (r["baseline_delta"] + beta * r["feature_delta"]) - margin
        for r in rows
    ]
    if min(residuals) < -1e-10:
        raise RuntimeError("The solved coefficient does not satisfy calibration constraints")
    return dict(
        version=1, method="minimum_nonnegative_affine_constraint_floor",
        rule="glycosylation_loss", beta=beta, margin=margin,
        calibration_groups=sorted(allowed), excluded_groups=sorted(excluded),
        constraint_count=len(rows), constraints=rows,
        binding_constraints=[dict(accession=r["accession"], contrast=r["contrast"])
                             for r in rows if abs(r["minimum_beta"] - beta) < 1e-10],
        minimum_signed_calibration_slack=min(residuals),
        experimental_labels_used=False,
        guarantee="calibration protein-level matched contrasts only",
    )


def coefficients_with_glycosylation_floor(coefficients: Mapping, fitted: Mapping) -> dict:
    """Return coefficients accepted by the existing standalone refinement API.

    The floor was derived with the other three coefficients held fixed.
    Callers must use the same other coefficients used during calibration.
    """
    if set(coefficients) != set(RULES):
        raise ValueError("Every existing rule coefficient must be supplied")
    if fitted.get("method") != "minimum_nonnegative_affine_constraint_floor":
        raise ValueError("Unexpected glycosylation floor artifact")
    beta = float(fitted["beta"])
    result = {k: float(v) for k, v in coefficients.items()}
    if any(not math.isfinite(v) or v < 0 for v in [beta, *result.values()]):
        raise ValueError("Finite nonnegative rule coefficients required")
    result["glycosylation_loss"] = max(result["glycosylation_loss"], beta)
    return result

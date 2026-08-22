"""Export the method-validation figure and label-control supplement."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results/method_completion_20261004"
OUT = ROOT / "docs/figures/method_completion_20261004"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run() -> dict:
    benchmark = json.loads((RESULT / "benchmark_results.json").read_text(encoding="utf-8"))
    controls = json.loads((RESULT / "label_controls_and_counts.json").read_text(encoding="utf-8"))
    demo_path = ROOT / "results/tool_demo_20261003/two_model_demo.json"
    demo = json.loads(demo_path.read_text(encoding="utf-8"))
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
        "axes.spines.top": False, "axes.spines.right": False,
        "svg.fonttype": "none", "savefig.facecolor": "white"})
    scenarios = [
        ("null", "Balanced null"),
        ("shared_position_confound", "Shared position constraint"),
        ("reference_composition_confound", "Reference imbalance"),
        ("true_outer_preference", "Planted outer preference"),
        ("opposite_preference", "Planted opposite preference"),
        ("distance_decay", "Planted distance decay"),
        ("uneven_clustered_null", "Unequal protein sizes"),
    ]
    methods = [
        ("naive_charge_arm_only", "Charged arm alone", "#8b5b81", "^", -0.18),
        ("naive_unstratified_paired", "Paired, no matching", "#be8322", "s", 0.0),
        ("packaged_within_site_stratified", "Packaged paired + matched", "#126c8a", "o", 0.18),
    ]
    index = {(c["scenario"], c["method"]): c for c in benchmark["cells"]}
    fig, axes = plt.subplots(1, 2, figsize=(15.3, 7.4), gridspec_kw={"width_ratios": [1.65, 1]})
    fig.subplots_adjust(left=0.19, right=0.975, top=0.82, bottom=0.31, wspace=0.35)
    left, right = axes
    for method, label, color, marker, offset in methods:
        cells = [index[(s, method)] for s, _ in scenarios]
        point = np.array([c["bias"] for c in cells])
        lo = np.array([c["bias_mc95"][0] for c in cells])
        hi = np.array([c["bias_mc95"][1] for c in cells])
        left.errorbar(point, np.arange(len(cells)) + offset, xerr=[point - lo, hi - point],
                      fmt=marker, color=color, ms=6.5, capsize=2.5, lw=1.2, label=label)
    left.axvline(0, color="#2c3c46", lw=1, ls="--", alpha=0.65)
    left.set_yticks(np.arange(len(scenarios)), [label for _, label in scenarios])
    left.invert_yaxis()
    left.set_xlim(-0.35, 0.67)
    left.set_xlabel("Bias relative to known charge target\n(synthetic score units)", labelpad=10)
    left.set_title("A   Known-answer validation", loc="left", fontweight="bold", pad=15)
    left.grid(axis="x", alpha=0.12)
    handles, labels = left.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.19, 0.17),
               ncol=3, frameon=False, fontsize=10)
    points = []
    for y, name, label, color in [(1, "alphamissense", "AlphaMissense", "#ac4538"),
                                   (0, "esm1b", "ESM1b", "#276b58")]:
        result = demo["results"][name]["result"]
        sd = result["penalty_sd"]
        estimate = result["estimate_raw"] / sd
        lo, hi = np.array(result["ci95_raw"]) / sd
        right.errorbar(estimate, y, xerr=[[estimate - lo], [hi - estimate]], fmt="o", color=color,
                       capsize=5, ms=8, lw=2)
        right.text(estimate, y + 0.19, f"{estimate:+.3f} SD", color=color, ha="center", fontweight="bold")
        points.append({"model": name, "estimate_raw": result["estimate_raw"],
            "normalizing_penalty_sd": sd, "estimate_sd": estimate, "ci95_sd": [float(lo), float(hi)],
            "contributing_sites": result["n_contributing_sites"], "contributing_proteins": result["n_proteins"]})
    right.axvline(0, color="#2c3c46", lw=1, ls="--", alpha=0.65)
    right.set_yticks([0, 1], ["ESM1b", "AlphaMissense"])
    right.set_ylim(-0.55, 1.55)
    right.set_xlim(-0.12, 0.12)
    right.set_xticks([-0.10, -0.05, 0, 0.05, 0.10])
    right.set_xlabel("Cytosolic − outer charge penalty\n(each model's own penalty SD)", labelpad=10)
    right.set_title("B   Opposite score preferences", loc="left", fontweight="bold", pad=15)
    right.grid(axis="x", alpha=0.12)
    right.text(0.5, 0.025, "− Greater relative outer penalty\n+ Greater relative cytosolic penalty",
               transform=right.transAxes, ha="center", va="bottom", fontsize=9, color="#364750")
    fig.suptitle("A validated membrane charge-score comparison", x=0.54, y=0.965,
                 fontsize=20, fontweight="bold", color="#182e3c")
    fig.text(0.54, 0.895, "Synthetic rules test the calculation; released scores reveal conditional behavior.",
             ha="center", fontsize=12, color="#455965")
    fig.text(0.19, 0.097, "A: 1,000 replicates per scenario; bars show Monte Carlo uncertainty of mean bias.\n"
             "Standard fixed effects coincides with the packaged estimator (maximum difference 1.4 × 10⁻¹⁶).",
             fontsize=10, color="#455965", va="top")
    fig.text(0.19, 0.018, "B: raw K/R scores minus 14 other uncharged alternatives at each site; distances 1–10.\n"
             "Bars are protein-resampling stability intervals. This is score behavior, not insertion or clinical accuracy.",
             fontsize=10, color="#455965", va="bottom")
    main_paths = []
    for suffix in ("png", "svg", "pdf"):
        path = OUT / f"method_plus_discovery.{suffix}"
        fig.savefig(path, dpi=220, bbox_inches="tight")
        main_paths.append(path)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.2))
    fig.subplots_adjust(top=0.72, bottom=0.2, wspace=0.28)
    for ax, name, label, color in zip(axes, ("alphamissense", "esm1b"),
                                    ("AlphaMissense", "ESM1b"), ("#ac4538", "#276b58")):
        result = controls["models"][name]
        values = np.array(result["shuffled_estimates_raw"])
        ax.hist(values, bins=35, color=color, alpha=0.75, edgecolor="white", linewidth=0.5)
        ax.axvline(0, color="#293d47", ls="--", lw=1.4)
        ax.set_title(label, fontweight="bold")
        ax.set_xlabel("Shuffled cytosolic − outer contrast (raw units)")
        ax.set_ylabel("Shuffles")
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.text(0.02, 0.98, f"Mean {result['shuffled_mean_raw']:+.2g}\n"
                f"Variance / exact oracle {result['empirical_to_exact_variance_ratio']:.3f}",
                transform=ax.transAxes, va="top", fontsize=10)
    fig.suptitle("Side-label destruction removes the score association", y=0.97, fontsize=17, fontweight="bold")
    fig.text(0.5, 0.855, "999 shuffles per model; protein/reference strata and their side counts preserved.",
             ha="center", fontsize=11)
    fig.text(0.5, 0.065, "These are pipeline negative controls. No causal randomization p-value is assigned to the original scores.",
             ha="center", fontsize=10)
    supplement_paths = []
    for suffix in ("png", "svg", "pdf"):
        path = OUT / f"shuffled_label_controls.{suffix}"
        fig.savefig(path, dpi=220, bbox_inches="tight")
        supplement_paths.append(path)
    plt.close(fig)
    summary = {"date": "2026-10-04", "main_figure_models": points,
        "main_figure_scenarios": [s for s, _ in scenarios],
        "input_hashes": {p.relative_to(ROOT).as_posix(): sha(p) for p in
                         (RESULT / "benchmark_results.json", RESULT / "label_controls_and_counts.json", demo_path)},
        "figures": [{"path": p.relative_to(ROOT).as_posix(), "sha256": sha(p)}
                    for p in main_paths + supplement_paths], "script_sha256": sha(Path(__file__))}
    (RESULT / "figure_manifest.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"main_figure": str(main_paths[0]), "supplement": str(supplement_paths[0])}, indent=2))
    return summary


if __name__ == "__main__":
    run()

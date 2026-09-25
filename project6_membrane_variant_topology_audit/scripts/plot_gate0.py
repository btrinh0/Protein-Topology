"""Figure for Gate 0: side asymmetry in AlphaMissense scores versus distance to the TMD.

Reads results/gate0_*/gate0_summary.json and writes gate0_side_asymmetry.png.
The matching table view is gate0_estimates.csv in the same directory.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RESULTS = ROOT / "results" / "gate0_20260924"

SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
TEXT_MUTED = "#8a8980"
GRID = "#e6e5e1"

# Colour follows the role (loss / gain / placebo), so it means the same thing in
# every panel. Validated: worst all-pairs CVD dE 9.2, normal-vision dE 24.0.
ROLE_COLOUR = {"loss": "#2a78d6", "gain": "#eb6834", "placebo": "#1baf7a"}

BINS = ["1-5", "6-10", "11-15", "16-30", "31-60"]


def series(results: dict, class_name: str) -> tuple[list[float], list[float], list[float]]:
    point, low, high = [], [], []
    for key in BINS:
        entry = results["distance_dose_response"][key]["classes"][class_name]
        point.append(entry["mean_difference_cytosolic_minus_outer"])
        low.append(entry["mean_difference_ci95"][0])
        high.append(entry["mean_difference_ci95"][1])
    return point, low, high


def contrast_series(results: dict, name: str) -> tuple[list[float], list[float], list[float]]:
    point, low, high = [], [], []
    for key in BINS:
        entry = results["distance_dose_response"][key]["contrasts"][name]
        point.append(entry["difference_in_differences"])
        low.append(entry["ci95"][0])
        high.append(entry["ci95"][1])
    return point, low, high


def style_axis(ax, title: str, subtitle: str, ylabel: str | None) -> None:
    ax.set_facecolor(SURFACE)
    ax.set_title(title, color=TEXT_PRIMARY, fontsize=11, fontweight="600", loc="left", pad=18)
    ax.text(0, 1.02, subtitle, transform=ax.transAxes, color=TEXT_SECONDARY,
            fontsize=8.5, va="bottom", ha="left")
    if ylabel:
        ax.set_ylabel(ylabel, color=TEXT_SECONDARY, fontsize=9)
    ax.axhline(0, color=TEXT_MUTED, linewidth=1, linestyle=(0, (4, 3)), zorder=1)
    ax.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(GRID)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=8.5, length=0)


def draw_lines(ax, entries: list[tuple[str, str, tuple]], label_offsets: dict[str, float]) -> None:
    x = range(len(BINS))
    for label, role, (point, low, high) in entries:
        colour = ROLE_COLOUR[role]
        ax.fill_between(x, low, high, color=colour, alpha=0.16, linewidth=0, zorder=2)
        ax.plot(x, point, color=colour, linewidth=2, zorder=3, solid_capstyle="round")
        ax.plot(x, point, "o", color=colour, markersize=5, markeredgecolor=SURFACE,
                markeredgewidth=2, zorder=4)
        ax.annotate(label, (len(BINS) - 1, point[-1]),
                    textcoords="offset points", xytext=(9, label_offsets.get(label, 0)),
                    color=colour, fontsize=8.5, fontweight="600", va="center", zorder=5)
    ax.set_xticks(list(x))
    ax.set_xticklabels(BINS)
    ax.set_xlim(-0.25, len(BINS) - 0.35)


def run(args: argparse.Namespace) -> Path:
    results_dir = Path(args.results)
    summary = json.loads((results_dir / "gate0_summary.json").read_text(encoding="utf-8"))
    results = summary["results"]

    figure, axes = plt.subplots(2, 2, figsize=(12.5, 9.2), facecolor=SURFACE)
    figure.subplots_adjust(hspace=0.42, wspace=0.28, left=0.07, right=0.88, top=0.87, bottom=0.11)

    style_axis(axes[0][0], "Positive charge: the placebo is not flat",
               "Cytosolic minus outer AlphaMissense score, within protein and substitution",
               "more damaging on the cytosolic side  →")
    draw_lines(axes[0][0], [
        ("K/R loss", "loss", series(results, "positive_charge_loss")),
        ("K/R gain", "gain", series(results, "positive_charge_gain")),
        ("no charge change", "placebo", series(results, "charge_neutral_placebo")),
    ], {"K/R gain": 6, "no charge change": -6})

    style_axis(axes[0][1], "Negative charge behaves the same way",
               "D/E substitutions, same estimator and same placebo", None)
    draw_lines(axes[0][1], [
        ("D/E loss", "loss", series(results, "negative_charge_loss")),
        ("D/E gain", "gain", series(results, "negative_charge_gain")),
        ("no charge change", "placebo", series(results, "charge_neutral_placebo")),
    ], {"D/E gain": 9, "no charge change": -9})

    style_axis(axes[1][0], "The charge-specific part does not decay",
               "Difference in differences against gain and against the placebo",
               "charge-specific side asymmetry  →")
    draw_lines(axes[1][0], [
        ("loss − gain", "loss", contrast_series(results, "positive_loss_minus_gain")),
        ("loss − placebo", "gain", contrast_series(results, "positive_loss_minus_placebo")),
    ], {"loss − gain": -7, "loss − placebo": 7})
    axes[1][0].set_xlabel("residues from the transmembrane boundary", color=TEXT_SECONDARY, fontsize=9)

    # Panel D: the subgroup where the confound disappears.
    ax = axes[1][1]
    style_axis(ax, "First transmembrane segment only",
               "Window ≤10 residues; faded = all segments, solid = first segment only", None)
    rows = [("K/R loss", "loss", "positive_charge_loss"),
            ("K/R gain", "gain", "positive_charge_gain"),
            ("no charge change", "placebo", "charge_neutral_placebo")]
    for index, (_, role, class_name) in enumerate(rows):
        colour = ROLE_COLOUR[role]
        for offset, scope, alpha in ((-0.16, results["by_window"]["10"], 0.35),
                                     (0.16, results["first_tmd_only"], 1.0)):
            entry = scope["classes"][class_name]
            value = entry["mean_difference_cytosolic_minus_outer"]
            low, high = entry["mean_difference_ci95"]
            ax.plot([index + offset, index + offset], [low, high], color=colour,
                    linewidth=2, alpha=alpha, solid_capstyle="round", zorder=3)
            ax.plot([index + offset], [value], "o", color=colour, markersize=8, alpha=alpha,
                    markeredgecolor=SURFACE, markeredgewidth=2, zorder=4)
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels([label for label, _, _ in rows], fontsize=8.5, fontweight="600")
    for tick, (_, role, _) in zip(ax.get_xticklabels(), rows):
        tick.set_color(ROLE_COLOUR[role])
    ax.set_xlim(-0.5, len(rows) - 0.5)
    ax.margins(y=0.18)

    handles = [plt.Line2D([], [], color=ROLE_COLOUR[role], linewidth=2, marker="o",
                          markersize=5, markeredgecolor=SURFACE, markeredgewidth=1.5, label=label)
               for label, role in (("charge loss", "loss"), ("charge gain", "gain"),
                                   ("no charge change / reference", "placebo"))]
    figure.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
                  fontsize=9, labelcolor=TEXT_SECONDARY, bbox_to_anchor=(0.5, 0.015))

    figure.suptitle("AlphaMissense scores membrane flanks asymmetrically, but not only by charge",
                    color=TEXT_PRIMARY, fontsize=14, fontweight="600", x=0.07, ha="left", y=0.965)
    figure.text(0.07, 0.925,
                "Stratified by protein and exact substitution; bands and bars are 95% "
                "protein-clustered bootstrap intervals. No clinical labels used.",
                color=TEXT_SECONDARY, fontsize=9.5, ha="left")

    out_path = results_dir / "gate0_side_asymmetry.png"
    figure.savefig(out_path, dpi=200, facecolor=SURFACE)
    plt.close(figure)
    print(out_path)
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", default=str(DEFAULT_RESULTS))
    run(parser.parse_args())


if __name__ == "__main__":
    main()

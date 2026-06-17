"""Figure: within-site charge penalties by membrane side and distance.

Reads results/gate0_*/gate0_within_site.json and writes within_site_charge_rules.png.
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
LOSS, GAIN = "#2a78d6", "#eb6834"

BINS = ["1-5", "6-10", "11-15", "16-30", "31-60"]


def pull(block: dict, path: list[str]) -> tuple[list, list, list]:
    node = block
    for key in path:
        node = node[key]
    point, low, high = [], [], []
    for key in BINS:
        entry = node[key]
        point.append(entry["side_difference"])
        low.append(entry["ci95"][0])
        high.append(entry["ci95"][1])
    return point, low, high


def panel(ax, title: str, subtitle: str, entries, ylabel: str | None) -> None:
    ax.set_facecolor(SURFACE)
    ax.set_title(title, color=TEXT_PRIMARY, fontsize=11.5, fontweight="600", loc="left", pad=20)
    ax.text(0, 1.03, subtitle, transform=ax.transAxes, color=TEXT_SECONDARY,
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

    x = range(len(BINS))
    for label, colour, (point, low, high), offset in entries:
        ax.fill_between(x, low, high, color=colour, alpha=0.16, linewidth=0, zorder=2)
        ax.plot(x, point, color=colour, linewidth=2, zorder=3, solid_capstyle="round")
        ax.plot(x, point, "o", color=colour, markersize=5, markeredgecolor=SURFACE,
                markeredgewidth=2, zorder=4)
        ax.annotate(label, (len(BINS) - 1, point[-1]), textcoords="offset points",
                    xytext=(9, offset), color=colour, fontsize=8.5, fontweight="600",
                    va="center", zorder=5)
    ax.set_xticks(list(x))
    ax.set_xticklabels(BINS)
    ax.set_xlim(-0.25, len(BINS) - 0.35)


def run(args: argparse.Namespace) -> Path:
    results_dir = Path(args.results)
    data = json.loads((results_dir / "gate0_within_site.json").read_text(encoding="utf-8"))
    positive = data["results"]["positive"]
    negative = data["results"]["negative"]

    figure, axes = plt.subplots(1, 2, figsize=(12.5, 5.6), facecolor=SURFACE)
    figure.subplots_adjust(wspace=0.3, left=0.075, right=0.86, top=0.76, bottom=0.17)

    panel(axes[0], "Positive charge: loss is encoded, gain is not",
          "Rule predicts loss positive and decaying, gain negative and decaying",
          [("K/R loss  ✓", LOSS, pull(positive, ["distance_dose_response"]), 8),
           ("K/R gain  ✗", GAIN, pull(positive, ["gain_at_uncharged_sites",
                                                      "distance_dose_response"]), -8)],
          "penalty larger on the cytosolic side  →")
    axes[0].set_xlabel("residues from the transmembrane boundary",
                       color=TEXT_SECONDARY, fontsize=9)

    panel(axes[1], "Negative charge: gain is encoded, loss is not",
          "Negative-inside depletion predicts gain positive and decaying",
          [("D/E gain  ✓", GAIN, pull(negative, ["gain_at_uncharged_sites",
                                                      "distance_dose_response"]), 0),
           ("D/E loss  ✗", LOSS, pull(negative, ["distance_dose_response"]), 0)],
          None)
    axes[1].set_xlabel("residues from the transmembrane boundary",
                       color=TEXT_SECONDARY, fontsize=9)

    handles = [plt.Line2D([], [], color=colour, linewidth=2, marker="o", markersize=5,
                          markeredgecolor=SURFACE, markeredgewidth=1.5, label=label)
               for label, colour in (("charge removed", LOSS), ("charge introduced", GAIN))]
    figure.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
                  fontsize=9, labelcolor=TEXT_SECONDARY, bbox_to_anchor=(0.5, 0.005))

    figure.suptitle("AlphaMissense encodes two of the four charge-topology rules",
                    color=TEXT_PRIMARY, fontsize=14, fontweight="600", x=0.075, ha="left", y=0.965)
    figure.text(0.075, 0.895,
                "Charge penalty measured inside a single residue, so site constraint cancels; "
                "side difference stratified by protein and reference residue.\n"
                "95% protein-clustered bootstrap intervals. No clinical labels used.",
                color=TEXT_SECONDARY, fontsize=9.5, ha="left", va="top")

    out_path = results_dir / "within_site_charge_rules.png"
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

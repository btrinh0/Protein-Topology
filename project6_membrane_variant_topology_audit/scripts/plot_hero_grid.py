"""Hero figure: evolution's footprint, two predictors, and the lab data.

Rows are sources of evidence, columns are the two positive-charge rules. Every
predictor curve is divided by that model's own spread of within-site penalties,
because AlphaMissense is bounded 0-1 and ESM1b is an unbounded log-likelihood
ratio; raw numbers across the two mean nothing.

The composition row is a residue count, not a model output, so it keeps its own
units. In the gain column it is drawn sign-flipped, because the rule's prediction
flips: positive residues are depleted from outer flanks, so *adding* one there is
the costly direction.

The fourth row uses the finalized Phase C summary and reports its sample sizes.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
SURFACE = "#fcfcfb"
TEXT_PRIMARY, TEXT_SECONDARY, TEXT_MUTED = "#0b0b0b", "#52514e", "#8a8980"
GRID, PENDING = "#e6e5e1", "#b9b8b0"
LOSS, GAIN, COMPOSITION = "#2a78d6", "#eb6834", "#1baf7a"
BINS = ["1-5", "6-10", "11-15", "16-30", "31-60"]


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def curves(results_dir: Path, controls_dir: Path) -> dict[str, list[float]]:
    audit = load(results_dir / "gate0_within_site.json")["results"]["positive"]
    controls = load(controls_dir / "within_site_controls.json")["results"]
    loss_sd = controls["loss"]["all"]["penalty_sd"]
    gain_sd = controls["gain"]["all"]["penalty_sd"]
    loss = [audit["distance_dose_response"][b]["side_difference"] / loss_sd for b in BINS]
    gain = [audit["gain_at_uncharged_sites"]["distance_dose_response"][b]["side_difference"]
            / gain_sd for b in BINS]
    return {"loss": loss, "gain": gain}


def style(ax, ylabel: str | None, show_x: bool) -> None:
    ax.set_facecolor(SURFACE)
    ax.axhline(0, color=TEXT_MUTED, linewidth=1, linestyle=(0, (4, 3)), zorder=1)
    ax.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(GRID)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=8, length=0)
    ax.set_xticks(range(len(BINS)))
    ax.set_xticklabels(BINS if show_x else [""] * len(BINS))
    ax.set_xlim(-0.3, len(BINS) - 0.7)
    if ylabel:
        ax.set_ylabel(ylabel, color=TEXT_SECONDARY, fontsize=8, labelpad=2)


def draw(ax, values: list[float], colour: str, label: str) -> None:
    x = range(len(values))
    ax.plot(x, values, color=colour, linewidth=2.2, zorder=3, solid_capstyle="round")
    ax.plot(x, values, "o", color=colour, markersize=5.5, markeredgecolor=SURFACE,
            markeredgewidth=2, zorder=4)
    ax.annotate(label, (len(values) - 1, values[-1]), textcoords="offset points",
                xytext=(8, 0), color=colour, fontsize=8, fontweight="600", va="center")


def run(args: argparse.Namespace) -> Path:
    composition = load(Path(args.composition) / "composition_summary.json")
    kr = [composition["results"]["by_distance_bin"]["positive_KR"][b]["difference"]
          for b in BINS]
    am = curves(ROOT / "results" / "am_full_20260927", ROOT / "results" / "controls_am_full_20260927")
    esm = curves(ROOT / "results" / "esm1b_20260927", ROOT / "results" / "controls_esm1b_20260927")

    figure, axes = plt.subplots(4, 2, figsize=(11.5, 13.2), facecolor=SURFACE)
    figure.subplots_adjust(hspace=0.42, wspace=0.34, left=0.235, right=0.87,
                           top=0.845, bottom=0.07)

    rows = [
        ("Evolution", "residue frequency\ncytosolic − outer",
         COMPOSITION, kr, [-v for v in kr], "K/R frequency", "flipped to the\ngain prediction"),
        ("AlphaMissense", "standardised\nside difference",
         None, am["loss"], am["gain"], "K/R loss", "K/R gain"),
        ("ESM1b", "standardised\nside difference",
         None, esm["loss"], esm["gain"], "K/R loss", "K/R gain"),
    ]
    for row, (title, ylabel, fixed, left, right, left_label, right_label) in enumerate(rows):
        for column, (values, colour, label) in enumerate((
                (left, fixed or LOSS, left_label), (right, fixed or GAIN, right_label))):
            ax = axes[row][column]
            style(ax, ylabel if column == 0 else None, show_x=False)
            draw(ax, values, colour, label)
        box = axes[row][0].get_position()
        figure.text(0.012, (box.y0 + box.y1) / 2, title, color=TEXT_PRIMARY,
                    fontsize=11, fontweight="700", ha="left", va="center")

    lab = json.loads((ROOT / "results" / "phase_c_20260928"
                      / "lab_row_for_hero.json").read_text(encoding="utf-8"))
    lab_n = {"loss": [62, 44, 16, 29, 36], "gain": [292, 179, 112, 234, 305]}
    for column, key in enumerate(("loss", "gain")):
        ax = axes[3][column]
        style(ax, "standardised\nside difference" if column == 0 else None, show_x=True)
        draw(ax, lab[key], LOSS if key == "loss" else GAIN, "K/R " + key)
        for x, (value, n) in enumerate(zip(lab[key], lab_n[key])):
            ax.annotate(f"n={n}", (x, value), textcoords="offset points",
                        xytext=(0, -14), ha="center", color=TEXT_MUTED, fontsize=6.5)
        ax.set_xlabel("residues from the transmembrane boundary",
                      color=TEXT_SECONDARY, fontsize=8.5)
    axes[3][1].text(0.03, 0.04,
                    "UNINFORMATIVE\nprimary ≤10: −0.085, 95% CI [−0.41, +0.17], p = 0.61",
                    transform=axes[3][1].transAxes, color="#b03a2e", fontsize=7.5,
                    fontweight="700", va="bottom")
    box = axes[3][0].get_position()
    figure.text(0.012, (box.y0 + box.y1) / 2, "Lab data\n(10 proteins)", color=TEXT_PRIMARY,
                fontsize=11, fontweight="700", ha="left", va="center")

    for column, heading, prediction in (
            (0, "Losing a positive charge", "rule: positive, decaying"),
            (1, "Gaining a positive charge", "rule: negative, decaying")):
        axes[0][column].set_title(heading, color=TEXT_PRIMARY, fontsize=12.5,
                                  fontweight="700", loc="left", pad=26)
        axes[0][column].text(0, 1.055, prediction, transform=axes[0][column].transAxes,
                             color=TEXT_SECONDARY, fontsize=8.5, ha="left")

    figure.suptitle("The gain column is where the models part company",
                    color=TEXT_PRIMARY, fontsize=15, fontweight="700", x=0.11, ha="left", y=0.975)
    figure.text(0.11, 0.935,
                "Predictor rows are divided by each model's own spread of within-site penalties, so the two are comparable. "
                "The composition row is a\nresidue count in its own units. Evolution and ESM1b both predict a cost to adding "
                "a positive charge on the outer side; AlphaMissense does not.",
                color=TEXT_SECONDARY, fontsize=9.5, ha="left", va="top")

    out_path = Path(args.outdir) / "hero_grid_charge_rules.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(out_path, dpi=200, facecolor=SURFACE)
    plt.close(figure)
    print(out_path)
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--composition", default=str(ROOT / "results" / "composition_20260927"))
    parser.add_argument("--outdir", default=str(ROOT / "results" / "figures_20260927"))
    run(parser.parse_args())


if __name__ == "__main__":
    main()

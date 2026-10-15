"""Render the retained trial layouts for the final paper figures."""

import json
import os
import sys

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

import figures_sts as cur

OUT = os.path.join(cur.OUT, "trial_figures4papers")

PALETTE = {
    "blue_main": "#0F4D92", "blue_secondary": "#3775BA",
    "green_1": "#DDF3DE", "green_2": "#AADCA9", "green_3": "#8BCF8B",
    "red_1": "#F6CFCB", "red_2": "#E9A6A1", "red_strong": "#B64342",
    "neutral": "#CFCECE", "highlight": "#FFD700",
    "teal": "#42949E", "violet": "#9A4D8E",
}
COMPARATOR_REDS = ["#D4685F", "#E5A19B"]
ERR = {"elinewidth": 2, "capthick": 2, "capsize": 15}


def style(font_size=24, axes_linewidth=3):
    mpl.rcdefaults()
    mpl.rcParams.update({
        "font.family": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": font_size, "axes.linewidth": axes_linewidth,
        "axes.spines.right": False, "axes.spines.top": False,
        "legend.frameon": False, "svg.fonttype": "none", "pdf.fonttype": 42,
        "mathtext.fontset": "dejavusans",
    })


def finalize(fig, name, pad=2):
    fig.tight_layout(pad=pad)
    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=300, facecolor="white")
    plt.close(fig)
    print(f"  wrote trial_figures4papers/{name}.png/.pdf")


def annotate(ax, x, top, value, fmt, size=32):
    ax.text(x, top, fmt.format(value), ha="center", va="bottom", fontsize=size)


def legend_panel(ax, handles, ncol=1, size=32):
    ax.set_axis_off()
    ax.legend(handles=handles, loc="center", fontsize=size, ncol=ncol, frameon=False)



def trial_accuracy():
    style(24, 3)
    per = pd.read_csv(os.path.join(cur.PVALS, "per_protein.csv"))
    main = pd.read_csv(os.path.join(cur.PVALS, "highlighted_configuration_main_comparisons.csv"))
    per = per[per.scheme == "guarded_loocv"]
    wide = per.pivot(index="gene", columns="model", values="spearman")
    g2 = pd.read_csv(os.path.join(cur.GLYCO, "accuracy_per_protein.csv"))
    g2 = g2[g2.membrane == 1].pivot(index="accession", columns="model", values="spearman")
    gene_of = per.drop_duplicates("accession").set_index("accession").gene
    wide["PI_DS_constraint_glyco"] = g2.rename(index=gene_of).reindex(wide.index)["PI_DS_constraint_glyco"]

    methods = [("positive_inside+disulfide_loss", "Ours (PI+DS)", PALETTE["blue_main"]),
               ("PI_DS_constraint_glyco", "Ours (+ glyco)", PALETTE["blue_secondary"]),
               ("VespaG", "VespaG", COMPARATOR_REDS[0]),
               ("GEMME", "GEMME", COMPARATOR_REDS[1])]
    means = np.array([wide[m].mean() for m, _, _ in methods])
    half = np.array([cur.T_975[12] * wide[m].std(ddof=1) / np.sqrt(len(wide)) for m, _, _ in methods])

    fig = plt.figure(figsize=(30, 10))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.6, 0.55])
    ax = fig.add_subplot(gs[0])
    bars = ax.bar(np.arange(len(methods)), means, yerr=half, error_kw=ERR,
                  color=[c for _, _, c in methods], edgecolor="black", linewidth=2)
    for b, v, h in zip(bars, means, half):
        annotate(ax, b.get_x() + b.get_width() / 2, v + h + 0.004, v, "{:.3f}")

    lo = means.min() - means.std()
    ax.set_ylim(lo, (means + half).max() + 0.06)
    ax.set_xticks([])
    ax.set_ylabel("Mean Spearman", fontsize=40, labelpad=12)
    ax.tick_params(axis="y", labelsize=30, length=10, width=2)

    ax = fig.add_subplot(gs[1])
    diff = (wide["positive_inside+disulfide_loss"] - wide["VespaG"]).sort_values(ascending=False)
    cols = [PALETTE["green_3"] if d > 0 else COMPARATOR_REDS[0] for d in diff.values]
    bars = ax.bar(np.arange(len(diff)), diff.values, color=cols, edgecolor="black", linewidth=2)
    for b, v in zip(bars, diff.values):
        annotate(ax, b.get_x() + b.get_width() / 2, max(v, 0) + 0.004, v, "{:+.2f}", size=20)
    ax.axhline(0, color="black", lw=2)
    ax.set_xticks(np.arange(len(diff)))
    ax.set_xticklabels(diff.index, rotation=45, ha="right", fontsize=24)
    t = main[(main.scheme == "guarded_loocv") & (main.model == "positive_inside+disulfide_loss")
             & (main.reference == "VespaG")].iloc[0]
    ax.set_ylabel(r"$\Delta\rho$ (Ours $-$ VespaG)", fontsize=40, labelpad=12)
    ax.set_title(f"mean {t.delta_spearman:+.3f}, paired p = {t.p_two_sided:.2f}", fontsize=32)
    ax.tick_params(axis="y", labelsize=30, length=10, width=2)

    legend_panel(fig.add_subplot(gs[2]),
                 [Patch(facecolor=c, edgecolor="black", label=l) for _, l, c in methods]
                 + [Patch(facecolor=PALETTE["green_3"], edgecolor="black", label="Ours higher"),
                    Patch(facecolor=COMPARATOR_REDS[0], edgecolor="black", label="VespaG higher")])
    finalize(fig, "accuracy_f4p")



def trial_refinement():
    style(24, 3)
    rows = cur.refinement_rows()
    labels = ["Positive-\ninside", "Disulfide", "Glyco A\n(Pro at X)", "Glyco B\n(S" + r"$\leftrightarrow$" + "T)",
              "Charge burial\n(control)"]
    fig = plt.figure(figsize=(30, 10))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.22])
    ax = fig.add_subplot(gs[0])
    w = 0.38
    tops = []
    for k, (_, before, after, sign) in enumerate(rows):
        for vals, dx, face, hatch in ((before, -w / 2, PALETTE["neutral"], "//"),
                                      (after, w / 2, PALETTE["blue_main"], None)):
            if vals is None:
                continue
            est, lo, hi, sd = vals
            v = sign * est / sd
            o = sorted((sign * lo / sd, sign * hi / sd))
            x = k + (dx if before is not None else 0)
            ax.bar(x, v, w, color=face, hatch=hatch, edgecolor="black", linewidth=2,
                   yerr=[[v - o[0]], [o[1] - v]], error_kw=ERR)
            annotate(ax, x, (o[1] if v >= 0 else o[0]) + (0.03 if v >= 0 else -0.16), v,
                     "{:+.2f}", size=26)
            tops += o
    ax.axhline(0, color="black", lw=2)
    ax.set_xticks(np.arange(len(rows)))
    ax.set_xticklabels(labels, fontsize=28)
    ax.set_ylabel("Standardized contrast\n(+ = agrees with rule)", fontsize=36, labelpad=12)
    ax.set_ylim(min(tops) - 0.3, max(tops) + 0.25)
    ax.tick_params(axis="y", labelsize=30, length=10, width=2)
    legend_panel(fig.add_subplot(gs[1]),
                 [Patch(facecolor=PALETTE["neutral"], hatch="//", edgecolor="black", label="Before"),
                  Patch(facecolor=PALETTE["blue_main"], edgecolor="black", label="After")])
    finalize(fig, "refinement_f4p")



def trial_validation():
    style(16, 2)
    est = pd.read_csv(os.path.join(cur.FORM_A, "estimates.csv"))
    crit = cur.load_json(os.path.join(cur.FORM_A, "criteria.json"))
    bench = cur.load_json(os.path.join(cur.METHOD, "benchmark_results.json"))
    planted = ["P_neg05", "P_neg10", "P_pos10", "P_decay"]
    models = [("S1", "Full model", PALETTE["blue_main"]), ("S0", "Noise-free generator", PALETTE["green_3"]),
              ("S0n", "Noisy generator", PALETTE["teal"]), ("S3", "Side-free chemistry", PALETTE["red_2"])]

    fig = plt.figure(figsize=(28, 7))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 0.8, 1.2, 0.6])

    ax = fig.add_subplot(gs[0])
    ax.plot([-0.13, 0.13], [-0.13, 0.13], color=PALETTE["neutral"], lw=2.5, zorder=0)
    for m, _, c in models:
        sub = est[(est.student == m) & est.scenario.isin(planted)]
        g = sub.groupby("scenario").agg(x=("oracle_raw", "first"), y=("estimate_raw", "mean"))
        ax.scatter(g.x, g.y, s=220, color=c, edgecolor="black", linewidth=1.5, alpha=0.9,
                   zorder=3 if m == "S1" else 2)
    ax.text(0.05, 0.95, f"slope {crit['S1']['slope_through_origin']:.3f}", transform=ax.transAxes,
            va="top", fontsize=24)
    ax.set_xlabel("Planted magnitude")
    ax.set_ylabel(r"Recovered $\hat{\Delta}$")
    ax.set_title("Recovery", fontsize=22)

    ax = fig.add_subplot(gs[1])
    s1 = np.abs(est[(est.student == "S1") & est.scenario.isin(planted)].estimate_raw.values)
    vals = [crit["S2_identity"]["max_abs_estimate"], s1.min()]
    bars = ax.bar([0, 1], vals, color=[PALETTE["red_2"], PALETTE["blue_main"]],
                  edgecolor="black", linewidth=2)
    for b, v in zip(bars, vals):
        annotate(ax, b.get_x() + b.get_width() / 2, v * 3, v, "{:.1e}", size=22)
    ax.set_yscale("log")
    ax.set_ylim(1e-20, 1)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["side-blind\n(max)", "side-aware\n(min)"])
    ax.set_ylabel(r"$|\hat{\Delta}|$")
    ax.set_title("Refusal", fontsize=22)

    ax = fig.add_subplot(gs[2])
    cells = {(c["scenario"], c["method"]): c for c in bench["cells"]}
    scen = [("shared_position_confound", "Position effects"), ("reference_composition_confound", "Composition")]
    meths = [("packaged_within_site_stratified", "Full method", PALETTE["blue_main"]),
             ("naive_charge_arm_only", "Without Stage 2", COMPARATOR_REDS[0]),
             ("naive_unstratified_paired", "Without Stage 3", COMPARATOR_REDS[1])]
    w = 0.26
    for i, (s, _) in enumerate(scen):
        for j, (m, _, c) in enumerate(meths):
            cell = cells[(s, m)]
            rate, (lo, hi) = cell["zero_rejection"], cell["zero_rejection_mc95"]
            x = i + (j - 1) * w
            ax.bar(x, rate, w, color=c, edgecolor="black", linewidth=2,
                   yerr=[[rate - lo], [hi - rate]], error_kw={**ERR, "capsize": 8})
            annotate(ax, x, hi + 0.02, rate, "{:.2f}", size=20)
    ax.axhline(0.05, color="black", lw=1.5, ls="--")
    ax.set_xticks([0, 1])
    ax.set_xticklabels([l for _, l in scen])
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("False-rejection rate")
    ax.set_title("Control necessity", fontsize=22)

    legend_panel(fig.add_subplot(gs[3]),
                 [Patch(facecolor=c, edgecolor="black", label=l) for _, l, c in models]
                 + [Patch(facecolor=c, edgecolor="black", label=l) for _, l, c in meths[1:]],
                 size=20)
    finalize(fig, "validation_f4p", pad=1)



def trial_forest():
    style(24, 3)
    sens = cur.load_json(os.path.join(cur.RESULTS, "comparator_sensitivity_20261001",
                                      "comparator_sensitivity.json"))
    comps = [("all_uncharged", "All uncharged\n(declared)"), ("polar_neutral_QNST", "Polar"),
             ("hydrophobic_AVLIMF", "Hydrophobic"), ("small_AGS", "Small"), ("aromatic_FWY", "Aromatic")]
    fig = plt.figure(figsize=(32, 9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.3])
    for k, (key, label, c) in enumerate((("alphamissense", "AlphaMissense", PALETTE["teal"]),
                                         ("esm1b", "ESM1b", PALETTE["violet"]))):
        ax = fig.add_subplot(gs[k])
        res = sens["results"][key]
        v = np.array([res[cp]["estimate_raw"] for cp, _ in comps])
        lo = np.array([res[cp]["ci95_raw"][0] for cp, _ in comps])
        hi = np.array([res[cp]["ci95_raw"][1] for cp, _ in comps])
        bars = ax.bar(np.arange(len(comps)), v, color=c, edgecolor="black", linewidth=2,
                      yerr=[v - lo, hi - v], error_kw=ERR)

        for b, a in zip(bars, [1.0, 0.75, 0.6, 0.45, 0.3]):
            b.set_facecolor(mpl.colors.to_rgba(c, a))
        span = max(abs(lo).max(), abs(hi).max())
        for x, val, l_, h_ in zip(np.arange(len(comps)), v, lo, hi):
            annotate(ax, x, (h_ if val >= 0 else l_) + (0.06 if val >= 0 else -0.2) * span, val,
                     "{:+.3f}", size=24)
        ax.axhline(0, color="black", lw=2)
        ax.set_ylim(-span * 1.35, span * 1.35)
        ax.set_xticks(np.arange(len(comps)))
        ax.set_xticklabels([l for _, l in comps], fontsize=24)
        ax.set_title(label, fontsize=36, color=c)
        ax.set_ylabel(r"$\Delta$ ($-$ = agrees with rule)" if k == 0 else "", fontsize=32)
        ax.tick_params(axis="y", labelsize=26, length=10, width=2)
    legend_panel(fig.add_subplot(gs[2]),
                 [Patch(facecolor=PALETTE["teal"], edgecolor="black", label="AlphaMissense"),
                  Patch(facecolor=PALETTE["violet"], edgecolor="black", label="ESM1b")], size=28)
    finalize(fig, "forest_f4p")



def trial_robustness():
    style(16, 2)
    audits = pd.read_csv(os.path.join(cur.PI_AUDIT, "all_audits.csv"))
    sub = audits[audits.model == "guarded_oof__positive_inside+disulfide_loss"]
    order = {"all_uncharged": "All uncharged", "polar_neutral": "Polar", "hydrophobic": "Hydrophobic",
             "small": "Small", "aromatic": "Aromatic"}
    grid = sub.pivot(index="comparator_name", columns="window", values="estimate_raw").loc[list(order)]
    fig, ax = plt.subplots(figsize=(12, 7))
    im = ax.imshow(grid.values, cmap="magma", aspect="auto")
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            v = grid.values[i, j]
            ax.text(j, i, f"{v:+.3f}", ha="center", va="center", fontsize=18,
                    color="white" if v < grid.values.mean() else "black")
    ax.set_xticks(range(grid.shape[1]))
    ax.set_xticklabels([f"≤{w}" for w in grid.columns])
    ax.set_yticks(range(grid.shape[0]))
    ax.set_yticklabels(list(order.values()))
    ax.set_xlabel("Flank window (residues)")
    cb = fig.colorbar(im, ax=ax)
    cb.set_label(r"$\Delta$ after refinement")
    finalize(fig, "robustness_f4p")



def trial_overview():
    saved = (cur.ACCENT, cur.RULE_SHADE, cur.save)
    cur.ACCENT, cur.RULE_SHADE = PALETTE["blue_main"], "#E3ECF6"

    def save_trial(fig, name):
        os.makedirs(OUT, exist_ok=True)
        for ext in ("png", "pdf"):
            fig.savefig(os.path.join(OUT, f"overview_f4p.{ext}"), dpi=300, facecolor="white",
                        bbox_inches="tight")
        plt.close(fig)
        print("  wrote trial_figures4papers/overview_f4p.png/.pdf")

    cur.save = save_trial
    mpl.rcParams.update({"font.family": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
                         "mathtext.fontset": "dejavusans"})
    try:
        cur.fig1_overview()
    finally:
        cur.ACCENT, cur.RULE_SHADE, cur.save = saved
        mpl.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm"})



def trial_overview_arch():
    from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle
    mpl.rcdefaults()
    mpl.rcParams.update({"font.family": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
                         "mathtext.fontset": "dejavusans", "font.size": 7, "pdf.fonttype": 42})
    blue, light, frozen, edge = PALETTE["blue_main"], "#E3ECF6", "#EDEDED", "#767676"
    fig, ax = plt.subplots(figsize=(6.5, 3.15))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 48.5)
    ax.set_aspect("equal")
    ax.axis("off")

    def box(x0, x1, yc, h, title, shape=None, note=None, kind="data"):
        face, ec = {"data": ("white", edge), "frozen": (frozen, edge), "fitted": (light, blue)}[kind]
        ax.add_patch(FancyBboxPatch((x0, yc - h / 2), x1 - x0, h, boxstyle="round,pad=0,rounding_size=1.2",
                                    fc=face, ec=ec, lw=0.8))
        lines = [(title, 7, "black")] + ([(shape, 7, "black")] if shape else []) + \
                ([(note, 5.8, "#4D4D4D")] if note else [])
        step = 2.6
        top = yc + step * (len(lines) - 1) / 2
        for k, (text, size, colour) in enumerate(lines):
            ax.text((x0 + x1) / 2, top - k * step, text, ha="center", va="center", fontsize=size,
                    color=colour)

    def arrow(p0, p1):
        ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=7, lw=0.8,
                                     color="#4D4D4D", shrinkA=0, shrinkB=0))


    y1, h = 39.5, 8.0
    box(0.5, 15.5, y1, h, "Protein sequence", r"+ substitution", "length $n$")
    box(19.5, 37.5, y1, h, "ProtT5-XL-UniRef50", r"$\mathbb{R}^{n\times 1024}$", "frozen, float16", "frozen")
    ax.add_patch(Rectangle((34.6, y1 + 1.6), 1.6, 1.3, fc=edge, lw=0))
    ax.add_patch(Circle((35.4, y1 + 2.9), 0.55, fc="none", ec=edge, lw=0.7))
    box(41.5, 56.5, y1, h, "Mutated position", r"$\mathbb{R}^{1024}$")
    box(60.5, 75.5, y1, h, "PCA", r"$\mathbb{R}^{1024}\rightarrow\mathbb{R}^{64}$", "fit on training folds", "fitted")
    box(79.5, 99.5, y1, h, "Tabular features", r"$\mathbb{R}^{18}$", "chemistry 9, identity 2, geometry 7")
    for x0, x1 in ((15.5, 19.5), (37.5, 41.5), (56.5, 60.5)):
        arrow((x0, y1), (x1, y1))


    y2 = 25.5
    box(74.5, 92.5, y2, h, "Concatenate", r"$\mathbb{R}^{64+18}=\mathbb{R}^{82}$")
    arrow((68.0, y1 - h / 2), (79.0, y2 + h / 2))
    arrow((89.5, y1 - h / 2), (86.0, y2 + h / 2))
    box(51.5, 70.5, y2, h, "Gradient-boosted trees", "500 iterations, depth 6", "learning rate 0.05", "fitted")
    box(39.5, 47.5, y2, h, "Raw", r"score $\in\mathbb{R}$")
    box(15.5, 35.5, y2, h, "Rule corrections", r"raw $+\sum_k b_k f_k$",
        "PI, DS: residual fit; glyco: constraint", "fitted")
    box(0.5, 11.5, y2, h, "Refined", r"score $\in\mathbb{R}$")
    for x0, x1 in ((74.5, 70.5), (51.5, 47.5), (39.5, 35.5), (15.5, 11.5)):
        arrow((x0, y2), (x1, y2))


    for k, (face, ec, label) in enumerate(((frozen, edge, "frozen"), (light, blue, "fitted on training groups only"))):
        x = 1 + k * 13
        ax.add_patch(FancyBboxPatch((x, 46.2), 2.4, 1.6, boxstyle="round,pad=0,rounding_size=0.4",
                                    fc=face, ec=ec, lw=0.7))
        ax.text(x + 3.2, 47.0, label, va="center", fontsize=6.2, color="#4D4D4D")


    ax.plot([1, 1, 76, 76], [14.2, 15.0, 15.0, 14.2], color=PALETTE["neutral"], lw=0.8)
    ax.text(38.5, 15.6, "Label-free audit", ha="center", va="bottom", fontsize=7.5)
    ax.add_patch(Rectangle((3, 5), 26, 5, fc="#f2f2f2", ec="none"))
    ax.plot([3, 29], [10, 10], color="#9a9a9a", lw=0.6)
    ax.plot([3, 29], [5, 5], color="#9a9a9a", lw=0.6)
    yy = np.linspace(1.6, 12.8, 200)
    ax.plot(16 + 1.2 * np.sin((yy - 1.6) * 2.4), yy, color="#4D4D4D", lw=1.0)
    for x in (10.5, 13, 19.5, 22):
        ax.text(x, 2.6, "+", ha="center", va="center", fontsize=8, color="#4D4D4D", weight="bold")
    ax.text(4, 11.6, "out", fontsize=7, va="center")
    ax.text(4, 2.6, "in", fontsize=7, va="center")
    arrow((31, 7.5), (37, 7.5))
    for x0, kr, side in ((40, 2.0, "in"), (52, 5.5, "out")):
        ax.add_patch(Rectangle((x0, 3.5), 1.6, kr, fc="#4D4D4D", lw=0))
        for i, hh in enumerate((3.4, 3.0, 3.7)):
            ax.add_patch(Rectangle((x0 + 2.2 + i * 1.8, 3.5), 1.4, hh, fc="#d9d9d9", lw=0))
        ax.plot([x0 - 0.5, x0 + 8.0], [3.5, 3.5], color="#9a9a9a", lw=0.5)
        ax.text(x0 + 3.7, 1.8, side, ha="center", va="center", fontsize=7)
    ax.text(40.8, 6.3, "K/R", ha="center", va="bottom", fontsize=6.5, color="#4D4D4D", rotation=90)
    arrow((62, 7.5), (68, 7.5))
    ax.text(71, 7.5, r"$\Delta$", ha="center", va="center", fontsize=11)
    ax.text(80, 7.5, "reads any predictor's scores;\nno experimental labels", va="center",
            fontsize=6.2, color="#4D4D4D")

    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"overview_arch_f4p.{ext}"), dpi=300, facecolor="white",
                    bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("  wrote trial_figures4papers/overview_arch_f4p.png/.pdf")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    if "--arch" in sys.argv:
        trial_overview_arch()
        sys.exit()
    print("trial figures (figures4papers style) ->", OUT)
    for build in (trial_overview, trial_validation, trial_accuracy, trial_refinement, trial_forest,
                  trial_robustness):
        print(build.__name__)
        build()

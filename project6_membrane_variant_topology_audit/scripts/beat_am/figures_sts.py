"""Generate the revised STS paper figures and verify their plotted values."""

import json
import os
import sys

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Arc, Circle, FancyArrowPatch, Polygon, Rectangle, RegularPolygon

from figures_paper import FORM_A, PROJ

RESULTS = os.path.join(PROJ, "results")
METHOD = os.path.join(RESULTS, "method_completion_20261004")
PVALS = os.path.join(RESULTS, "dms_comparator_pvalues_20261028")
PI_AUDIT = os.path.join(RESULTS, "positive_inside_standalone_20261028")
MULTI = os.path.join(RESULTS, "standalone_multirule_20261025")
GLYCO = os.path.join(RESULTS, "glycosylation_constraint_repair_20261031", "run2")
OUT = os.path.join(PROJ, "figures", "sts_20261031")

mpl.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm",
    "font.serif": ["cmr10", "DejaVu Serif"],
    "axes.formatter.use_mathtext": True,        # Use an ASCII minus because cmr10 lacks a Unicode minus glyph.
    "font.size": 7.5, "axes.labelsize": 7.5, "axes.titlesize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "lines.linewidth": 1.2, "lines.markersize": 3.5,
    "errorbar.capsize": 0, "legend.frameon": False,
    "pdf.fonttype": 42, "savefig.dpi": 300, "savefig.bbox": "tight",
})
FULL_W, HALF_W = 6.5, 3.2
ACCENT = "#009E73"
AM     = "#0072B2"
ESM    = "#D55E00"
GRAY_D, GRAY_M, GRAY_L = "#4d4d4d", "#9a9a9a", "#d9d9d9"
RULE_SHADE = "#eef6f0"
EBAR = 0.8
T_975 = {12: 2.178813, 19: 2.093024}

VALUES = {}
MISMATCH = {}
TEXT_ISSUES = [] # Record unsupported paper claims separately from failed figure checks.


def load_json(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def note(fig, line):
    VALUES.setdefault(fig, []).append(line)


def check(fig, label, file_value, paper_value, tol):
    """Print file value next to paper value; record a mismatch beyond rounding."""
    ok = abs(file_value - paper_value) <= tol + 1e-15
    line = (f"CHECK {label}: file {file_value:+.8g} | paper {paper_value:+.8g} "
            f"| {'match' if ok else 'MISMATCH'}")
    note(fig, line)
    if not ok:
        MISMATCH.setdefault(fig, []).append(line)


def panel_label(ax, letter):
    # Use math bold because cmr10 has no bold text face.
    ax.text(-0.18, 1.05, rf"$\mathbf{{{letter}}}$", transform=ax.transAxes, fontsize=9)


def save(fig, name):
    for line in VALUES.get(name, []):
        print("  " + line)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, f"VALUES_{name}.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(VALUES.get(name, [])) + "\n")
    if MISMATCH.get(name):
        print(f"  BLOCKED {name}: mismatch with paper text, figure not saved")
        plt.close(fig)
        return
    for ext in ("pdf", "png"):
        os.makedirs(os.path.join(OUT, ext), exist_ok=True)
        fig.savefig(os.path.join(OUT, ext, f"{name}.{ext}"), facecolor="white")
    plt.close(fig)
    print(f"  wrote pdf/{name}.pdf and png/{name}.png")



def fig1_overview():
    name = "fig1_overview"
    fig, ax = plt.subplots(figsize=(FULL_W, 2.2))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 34)
    ax.set_aspect("equal")
    ax.axis("off")

    def arrow(x0, x1, y):
        ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>",
                                     mutation_scale=6, lw=0.6, color=GRAY_D))

    def bracket(x0, x1, y, label):
        ax.plot([x0, x0, x1, x1], [y - 0.8, y, y, y - 0.8], color=GRAY_M, lw=0.6)
        ax.text((x0 + x1) / 2, y + 0.6, label, ha="center", va="bottom", fontsize=8)

    def number(x, y, n):
        ax.add_patch(Circle((x, y), 1.2, fc="white", ec=GRAY_D, lw=0.6))
        ax.text(x, y, str(n), ha="center", va="center", fontsize=7)


    bracket(1, 19, 31.5, "Input")
    bracket(23, 79, 31.5, "Model")
    bracket(83, 99, 31.5, "Output")

    seq_ref, seq_var, k = "MLKVAGTL", "MLKVKGTL", 4
    for row, (seq, y) in enumerate(((seq_ref, 26.5), (seq_var, 21.5))):
        for i, ch in enumerate(seq):
            x = 2.5 + i * 2.1
            if row == 1 and i == k:
                ax.add_patch(Rectangle((x - 1.0, y - 1.3), 2.0, 2.6, fc=GRAY_L, lw=0))
            ax.text(x, y, ch, ha="center", va="center", fontsize=7,
                    weight="bold" if (row == 1 and i == k) else "normal",
                    color="black" if (row == 1 and i == k) else GRAY_D)
    arrow(19.5, 23.5, 24)


    rng = np.random.default_rng(3)
    for i in range(6):
        for j in range(4):
            shade = str(0.55 + 0.4 * rng.random())
            ax.add_patch(Rectangle((26 + i * 1.5, 19.5 + j * 1.5), 1.5, 1.5, fc=shade,
                                   ec="white", lw=0.3))
    ax.add_patch(Rectangle((35.2, 26.1), 2.0, 1.6, fc=GRAY_D, lw=0))
    ax.add_patch(Arc((36.2, 27.7), 1.4, 1.8, theta1=0, theta2=180, color=GRAY_D, lw=0.7))
    number(25.5, 28.5, 1)
    ax.text(30.5, 18.6, "ProtT5", ha="center", va="top", fontsize=7)
    arrow(38.5, 43, 24)


    nodes = {(48, 28): [(45.5, 24.5), (50.5, 24.5)],
             (45.5, 24.5): [(44, 21), (47, 21)], (50.5, 24.5): [(49, 21), (52, 21)]}
    for (x0, y0), kids in nodes.items():
        for x1, y1 in kids:
            ax.plot([x0, x1], [y0, y1], color=ACCENT, lw=0.7, zorder=1)
    for x, y in [(48, 28), (45.5, 24.5), (50.5, 24.5), (44, 21), (47, 21), (49, 21), (52, 21)]:
        ax.add_patch(Circle((x, y), 0.75, fc=ACCENT, ec="none", zorder=2))
    number(43.5, 28.5, 2)
    ax.text(48, 18.6, "boosted trees", ha="center", va="top", fontsize=7)
    arrow(53.5, 57.5, 24)

    # PI and disulfide weights use DMS residuals; glycosylation uses a label-free constraint and charge burial is inactive.
    icons_x = [60.5, 65, 71.5, 77]
    for x, ec in zip(icons_x, (ACCENT, ACCENT, ACCENT, GRAY_M)):
        ax.add_patch(Circle((x, 24), 1.8, fc="white", ec=ec, lw=0.7,
                            ls="-" if ec == ACCENT else (0, (2, 1.5))))
    ax.text(icons_x[0], 24, r"$\pm$", ha="center", va="center", fontsize=8, color=ACCENT)

    ax.plot([icons_x[1] - 0.8, icons_x[1] + 0.8], [24, 24], color=ACCENT, lw=0.8)
    for dx in (-0.8, 0.8):
        ax.add_patch(Circle((icons_x[1] + dx, 24), 0.45, fc=ACCENT, ec="none"))
    ax.add_patch(RegularPolygon((icons_x[2], 24), 6, radius=0.9, fc="none", ec=ACCENT, lw=0.7))
    ax.add_patch(Rectangle((icons_x[3] - 1.3, 23.3), 2.6, 1.4, fc="white", ec=GRAY_M, lw=0.4))
    ax.text(icons_x[3], 24, "+", ha="center", va="center", fontsize=7, color=GRAY_M)

    for x0, x1, label in ((58.7, 66.8, "residual fit"), (69.7, 73.3, "constraint")):
        ax.plot([x0, x0, x1, x1], [22.0, 21.5, 21.5, 22.0], color=ACCENT, lw=0.5)
        ax.text((x0 + x1) / 2, 21.1, label, ha="center", va="top", fontsize=5.5,
                color=ACCENT)
    ax.text(icons_x[3], 21.1, "off", ha="center", va="top", fontsize=5.5, color=GRAY_M)
    number(57.5, 28.5, 3)
    ax.text(68.75, 18.6, "rule corrections", ha="center", va="top", fontsize=7)
    arrow(79.5, 83.5, 24)


    grad = np.linspace(0, 1, 64)[None, :]
    ax.imshow(grad, extent=(84, 98, 23, 25), cmap="Greys", aspect="auto", zorder=1)
    ax.add_patch(Rectangle((84, 23), 14, 2, fc="none", ec=GRAY_D, lw=0.5))
    ax.plot([93, 93], [22.2, 25.8], color=ACCENT, lw=1.4)
    ax.text(91, 18.6, "score", ha="center", va="top", fontsize=7)


    bracket(1, 76, 13.8, "Audit")

    ax.add_patch(Rectangle((3, 5), 26, 5, fc="#f2f2f2", ec="none"))
    ax.plot([3, 29], [10, 10], color=GRAY_M, lw=0.6)
    ax.plot([3, 29], [5, 5], color=GRAY_M, lw=0.6)
    yy = np.linspace(1.6, 12.8, 200)
    ax.plot(16 + 1.2 * np.sin((yy - 1.6) * 2.4), yy, color=GRAY_D, lw=1.0)
    for x in (10.5, 13, 19.5, 22):
        ax.text(x, 2.6, "+", ha="center", va="center", fontsize=8, color=GRAY_D, weight="bold")
    ax.text(4, 11.6, "out", fontsize=7, va="center")
    ax.text(4, 2.6, "in", fontsize=7, va="center")
    arrow(31, 37, 7.5)


    for x0, kr, side in ((40, 2.0, "in"), (52, 5.5, "out")):
        ax.add_patch(Rectangle((x0, 3.5), 1.6, kr, fc=GRAY_D, lw=0))
        for i, h in enumerate((3.4, 3.0, 3.7)):
            ax.add_patch(Rectangle((x0 + 2.2 + i * 1.8, 3.5), 1.4, h, fc=GRAY_L, lw=0))
        ax.plot([x0 - 0.5, x0 + 8.0], [3.5, 3.5], color=GRAY_M, lw=0.5)
        ax.text(x0 + 3.7, 1.8, side, ha="center", va="center", fontsize=7)
    ax.text(40.8, 6.3, "K/R", ha="center", va="bottom", fontsize=7, color=GRAY_D,
            rotation=90)
    arrow(62, 68, 7.5)
    ax.text(71, 7.5, r"$\Delta$", ha="center", va="center", fontsize=11)

    note(name, "illustration only; no data values plotted")
    save(fig, name)



def fig2_validation():
    name = "fig2_validation"
    est = pd.read_csv(os.path.join(FORM_A, "estimates.csv"))
    crit = load_json(os.path.join(FORM_A, "criteria.json"))
    bench = load_json(os.path.join(METHOD, "benchmark_results.json"))

    planted = ["P_neg05", "P_neg10", "P_pos10", "P_decay"]
    # Reserve the side-blind S2 model for the refusal panel rather than recovery.
    models = {"S0": "noise-free generator", "S0n": "noisy generator",
              "S3": "side-free chemistry model", "S1": "full model"}
    slope = crit["S1"]["slope_through_origin"]

    fig, axd = plt.subplot_mosaic([["a", "b", "c"]], figsize=(FULL_W, 1.9),
                                  gridspec_kw={"wspace": 0.55})


    ax = axd["a"]
    lim = 0.13
    ax.plot([-lim, lim], [-lim, lim], color=GRAY_L, lw=0.8, zorder=0)
    for scen in planted:
        for m in models:
            rows = est[(est.scenario == scen) & (est.student == m)]
            x = rows.oracle_raw.iloc[0]
            mean = rows.estimate_raw.mean()
            half = T_975[19] * rows.estimate_raw.std(ddof=1) / np.sqrt(len(rows))
            col = ACCENT if m == "S1" else GRAY_M
            ax.errorbar(x, mean, yerr=half, fmt="o", color=col, lw=EBAR,
                        zorder=3 if m == "S1" else 2)
            note(name, f"a {scen} {m} ({models[m]}): planted {x:+.6f}, "
                       f"mean recovered {mean:+.6f}, 95% t-interval of mean "
                       f"[{mean - half:+.6f}, {mean + half:+.6f}], n={len(rows)}")
    ax.text(0.04, 0.96, f"slope {slope:.3f}", transform=ax.transAxes, va="top", fontsize=7)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xticks([-0.1, 0, 0.1])
    ax.set_yticks([-0.1, 0, 0.1])
    ax.set_xlabel("planted magnitude")
    ax.set_ylabel(r"recovered $\hat{\Delta}$")

    s1 = est[(est.student == "S1") & est.scenario.isin(planted)]
    recomputed = (s1.oracle_raw * s1.estimate_raw).sum() / (s1.oracle_raw ** 2).sum()
    note(name, f"a slope through origin, criteria.json: {slope:.10f}; "
               f"recomputed from S1 rows (n={len(s1)}): {recomputed:.10f}")
    check(name, "slope_matches_recomputed", slope, recomputed, 1e-9)
    check(name, "slope", slope, 0.819, 0.0005)


    ax = axd["b"]
    s2 = np.abs(est[est.student == "S2"].estimate_raw.values)
    s1_abs = np.abs(s1.estimate_raw.values)
    s2_max = crit["S2_identity"]["max_abs_estimate"]
    rng = np.random.default_rng(0)
    ax.scatter(rng.uniform(-0.15, 0.15, len(s2)), s2, s=5, color=GRAY_D, lw=0, alpha=0.8)
    ax.scatter(1 + rng.uniform(-0.15, 0.15, len(s1_abs)), s1_abs, s=5, color=ACCENT,
               lw=0, alpha=0.8)
    ratio = s1_abs.min() / s2_max
    exponent = int(np.floor(np.log10(ratio)))
    ax.plot([0.5, 0.5], [s2_max, s1_abs.min()], color=GRAY_D, lw=0.6)
    for yv in (s2_max, s1_abs.min()):
        ax.plot([0.45, 0.55], [yv, yv], color=GRAY_D, lw=0.6)
    ax.text(0.58, np.sqrt(ratio) * s2_max, rf"$>10^{{{exponent}}}\times$", va="center",
            fontsize=7)
    ax.set_yscale("log")
    ax.set_ylim(1e-20, 1)
    ax.set_yticks([1e-18, 1e-12, 1e-6, 1])
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["side-blind", "side-aware"])
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylabel(r"$|\hat{\Delta}|$")
    note(name, f"b side-blind (S2) runs: {len(s2)}; max |Δ̂| criteria.json {s2_max:.4e}; "
               f"max from estimates.csv {s2.max():.4e}; exact zeros {int((s2 == 0).sum())}")
    note(name, f"b side-aware = full model (S1) on planted scenarios, {len(s1_abs)} runs; "
               f"|Δ̂| min {s1_abs.min():.6f}, max {s1_abs.max():.6f}")
    note(name, f"b ratio min side-aware / max side-blind = {ratio:.3e} "
               f"(log10 {np.log10(ratio):.3f}; computed); label '>10^{exponent}x'")
    check(name, "side_blind_max", s2_max, 3.62e-18, 0.005e-18)
    check(name, "side_blind_max_matches_estimates", s2.max(), s2_max, 1e-30)
    check(name, "orders_of_magnitude_floor_matches_paper_fifteen", exponent, 15, 0)


    ax = axd["c"]
    cells = {(c["scenario"], c["method"]): c for c in bench["cells"]}
    scenarios = [("shared_position_confound", "position effects"),
                 ("reference_composition_confound", "composition")]
    methods = [("packaged_within_site_stratified", "full method", ACCENT, -0.15),
               ("naive_charge_arm_only", "without Stage 2", GRAY_D, 0.0),
               ("naive_unstratified_paired", "without Stage 3", GRAY_M, 0.15)]
    ax.axhline(0.05, color=GRAY_L, lw=0.8, zorder=0)
    for i, (scen, _) in enumerate(scenarios):
        for meth, label, col, dx in methods:
            c = cells[(scen, meth)]
            rate, (lo, hi) = c["zero_rejection"], c["zero_rejection_mc95"]
            ax.errorbar(i + dx, rate, yerr=[[rate - lo], [hi - rate]], fmt="o", color=col,
                        lw=EBAR, label=label if i == 0 else None)
            note(name, f"c {scen} | {meth} ({label}): false-rejection {rate:.3f} "
                       f"[{lo:.4f}, {hi:.4f}], {round(rate * c['replicates'])}/"
                       f"{c['replicates']}, bias {c['bias']:+.6f}")
    ax.set_xticks([0, 1])
    ax.set_xticklabels([s for _, s in scenarios])
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(-0.03, 1.05)
    ax.set_yticks([0, 0.5, 1])
    ax.set_ylabel("false-rejection rate")
    fig.legend(loc="lower center", bbox_to_anchor=(0.5, 0.98), ncol=3, handletextpad=0.2,
               columnspacing=1.2, borderaxespad=0)

    sp, cp = scenarios[0][0], scenarios[1][0]
    check(name, "stage2_shortcut_rate", cells[(sp, methods[1][0])]["zero_rejection"], 1.0, 0)
    check(name, "stage3_shortcut_rate", cells[(cp, methods[2][0])]["zero_rejection"], 1.0, 0)
    check(name, "stage2_shortcut_bias", cells[(sp, methods[1][0])]["bias"], 0.600, 0.0005)
    check(name, "stage3_shortcut_bias", cells[(cp, methods[2][0])]["bias"], 0.543, 0.0005)

    for letter in "abc":
        panel_label(axd[letter], letter)
    save(fig, name)



def fig3_forest():
    """All five rows come from one saved run (same bootstrap), w = 10 residues."""
    name = "fig3_forest"
    sens = load_json(os.path.join(RESULTS, "comparator_sensitivity_20261001",
                                  "comparator_sensitivity.json"))
    manifest = {m["model"]: m for m in load_json(os.path.join(METHOD, "figure_manifest.json"))
                ["main_figure_models"]}
    rows = [("all_uncharged", r"$\mathbf{all\ uncharged}$"), ("polar_neutral_QNST", "polar (NQST)"),
            ("hydrophobic_AVLIMF", "hydrophobic (AFILMV)"), ("small_AGS", "small (AGS)"),
            ("aromatic_FWY", "aromatic (FWY)")]
    models = [("alphamissense", "AlphaMissense", AM, 0.0121, 0.00005, (0.0105, 0.0138), 0.0001),
              ("esm1b", "ESM1b", ESM, -0.147, 0.0005, (-0.172, -0.121), 0.001)]
    note(name, f"source: comparator_sensitivity.json, window {sens['window_residues']} residues")

    fig, axd = plt.subplot_mosaic([["a", "b"]], figsize=(FULL_W, 1.6),
                                  gridspec_kw={"wspace": 0.25})
    for letter, (key, label, col, paper_pt, tol, paper_ci, ci_tol) in zip("ab", models):
        ax = axd[letter]
        res = sens["results"][key]
        ends = []
        for k, (comp, comp_label) in enumerate(rows):
            r = res[comp]
            est, (lo, hi) = r["estimate_raw"], r["ci95_raw"]
            ends += [lo, hi]
            declared = comp == "all_uncharged"
            ax.errorbar(est, k, xerr=[[est - lo], [hi - est]], fmt="o", color=col, lw=EBAR,
                        mfc=col if declared else "white", mec=col, zorder=3)
            note(name, f"{letter} {label} | {comp} ({r['comparator']}): {est:+.6f} "
                       f"[{lo:+.6f}, {hi:+.6f}]{' (declared)' if declared else ''}")
        span = max(abs(v) for v in ends) * 1.1
        ax.set_xlim(-span, span)
        ax.axvspan(-span, 0, color=RULE_SHADE, lw=0, zorder=0)
        ax.axvline(0, color=GRAY_L, lw=0.8, zorder=1)
        ax.set_ylim(len(rows) - 0.5, -0.5)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([r[1] for r in rows] if letter == "a" else [])
        ax.set_title(label, loc="left", color=col)
        ax.set_xlabel(r"$\Delta$ (more negative = agrees with rule)")
        ax.locator_params(axis="x", nbins=4)
        panel_label(ax, letter)

        # Use the production figure manifest's audit for the declared comparator point.
        d = res["all_uncharged"]
        m = manifest[key]
        check(name, f"{label}_declared_matches_manifest", d["estimate_raw"], m["estimate_raw"], 1e-12)
        check(name, f"{label}_declared_point", d["estimate_raw"], paper_pt, tol)
        note(name, f"{label} declared contributing sites {m['contributing_sites']:,}, "
                   f"proteins {m['contributing_proteins']:,} (figure_manifest.json)")
        lo, hi = d["ci95_raw"]
        if abs(lo - paper_ci[0]) > ci_tol / 2 or abs(hi - paper_ci[1]) > ci_tol / 2:
            line = (f"TEXT {label} declared interval in this run is [{lo:+.4f}, {hi:+.4f}]; "
                    f"the paper prints [{paper_ci[0]:+.4f}, {paper_ci[1]:+.4f}] from a separate "
                    "bootstrap run of the same audit. Update the text to this run's interval.")
            note(name, line)
            TEXT_ISSUES.append(f"{name}: {line}")
        alts = [res[c]["estimate_raw"] for c, _ in rows[1:]]
        note(name, f"{label} alternatives range {min(alts):+.6f} to {max(alts):+.6f}; "
                   f"intervals excluding zero on the rule side: "
                   f"{sum(res[c]['ci95_raw'][1] < 0 for c, _ in rows)}/5")

    line = ("TEXT paper.tex §5.6 quotes alternative ranges -0.051 to +0.059 (AlphaMissense) "
            "and -0.175 to -0.023 (ESM1b); no saved output contains them. This figure's "
            "ranges are in VALUES_fig3_forest.txt.")
    note(name, line)
    TEXT_ISSUES.append(f"{name}: {line}")
    save(fig, name)



def fig5_accuracy():
    name = "fig5_accuracy"
    per = pd.read_csv(os.path.join(PVALS, "per_protein.csv"))
    per_assay = pd.read_csv(os.path.join(PVALS, "per_assay.csv"))
    main = pd.read_csv(os.path.join(PVALS, "highlighted_configuration_main_comparisons.csv"))
    vs_raw = pd.read_csv(os.path.join(PVALS, "refinements_vs_raw.csv"))
    scheme, refined, three = "guarded_loocv", "positive_inside+disulfide_loss", "PI_DS_constraint_glyco"
    columns = [("GEMME", "GEMME", 3, 0.428), ("VespaG", "VespaG", 3, 0.451),
               (refined, "Ours", 4, 0.4626), (three, "Ours\n+ glyco", 4, 0.4608)]

    per = per[per.scheme == scheme]
    wide = per.pivot(index="gene", columns="model", values="spearman")

    # Require the repair's raw and PI-plus-disulfide columns to reproduce the comparator tables exactly.
    g2 = pd.read_csv(os.path.join(GLYCO, "accuracy_per_protein.csv"))
    g2 = g2[g2.membrane == 1].pivot(index="accession", columns="model", values="spearman")
    gene_of = per.drop_duplicates("accession").set_index("accession").gene   # The run2 tables use HGNC gene symbols.
    g2 = g2.rename(index=gene_of).reindex(wide.index)
    assert not g2[["existing_PI_DS", "raw", three]].isna().any().any(), "protein sets differ"
    same = max(abs(g2["existing_PI_DS"] - wide[refined]).max(), abs(g2["raw"] - wide["raw"]).max())
    check(name, "run2_PI_DS_and_raw_match_comparator_file", same, 0.0, 1e-12)
    wide[three] = g2[three]
    acc = load_json(os.path.join(GLYCO, "accuracy_results.json"))["comparisons"]
    c3 = next(c for c in acc if c["model"] == three and c["reference"] == "existing_PI_DS")
    note(name, f"three-rule minus PI+DS (file, text only): {c3['delta_spearman']:+.6f} "
               f"[{c3['ci95_delta_spearman'][0]:+.6f}, {c3['ci95_delta_spearman'][1]:+.6f}]")
    check(name, "three_rule_minus_PI_DS", c3["delta_spearman"], -0.001762, 0.0000005)

    # Average assay-level correlations equally within each protein.
    pa = per_assay[per_assay.scheme == scheme]
    re_agg = pa.groupby(["gene", "model"]).spearman.mean().unstack()
    worst = max(abs(re_agg[m] - wide[m]).max() for m in (refined, "raw", "GEMME", "VespaG"))
    note(name, f"aggregation check: max |per_protein - mean(per_assay)| = {worst:.2e}")
    check(name, "per_protein_is_mean_of_assays", worst, 0.0, 1e-12)
    check(name, "mean_Ours_raw_(table_only)", wide["raw"].mean(), 0.4623, 0.00005)

    fig, axd = plt.subplot_mosaic([["a", "b"]], figsize=(FULL_W, 2.2),
                                  gridspec_kw={"width_ratios": [0.52, 0.48], "wspace": 0.35})


    ax = axd["a"]
    for gene, row in wide.iterrows():
        ax.plot(range(len(columns)), [row[m] for m, _, _, _ in columns], color=GRAY_L,
                lw=0.6, zorder=1)
    for i, (m, label, ndp, paper) in enumerate(columns):
        if label.startswith("Ours"):
            filled = label == "Ours"
            ax.scatter(np.full(len(wide), i), wide[m], s=9 if filled else 8, zorder=2,
                       color=ACCENT if filled else "white", edgecolor=ACCENT,
                       lw=0 if filled else 0.6)
        else:
            ax.scatter(np.full(len(wide), i), wide[m], color=GRAY_M, s=9, lw=0, zorder=2)
        vals = wide[m].values
        mean = vals.mean()
        half = T_975[12] * vals.std(ddof=1) / np.sqrt(len(vals))
        ax.plot([i + 0.2, i + 0.2], [mean - half, mean + half], color="black", lw=EBAR)
        ax.plot([i + 0.12, i + 0.28], [mean, mean], color="black", lw=2.0)
        label = label.replace("\n", " ")
        check(name, f"mean_{label}", mean, paper, 0.5 * 10 ** -ndp)
        note(name, f"a {label}: n={len(vals)} proteins, mean {mean:.6f}, 95% t-interval "
                   f"[{mean - half:.6f}, {mean + half:.6f}] (computed)")
    for gene, row in wide.iterrows():
        note(name, "a " + gene + ": " + ", ".join(
            f"{label.replace(chr(10), ' ')} {row[m]:.6f}" for m, label, _, _ in columns))
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels([label for _, label, _, _ in columns])
    ax.set_xlim(-0.35, len(columns) - 0.5)
    ax.set_ylim(0, 0.75)
    ax.set_ylabel(r"Spearman $\rho$")

    r = vs_raw[(vs_raw.scheme == scheme) & (vs_raw.model == refined)].iloc[0]
    lo, hi = json.loads(r.ci95_delta_spearman)
    note(name, f"refined minus raw (file, table only): {r.delta_spearman:+.6f} [{lo:+.6f}, {hi:+.6f}]")
    check(name, "refined_minus_raw", r.delta_spearman, 0.000266, 0.0000005)


    ax = axd["b"]
    diff = (wide[refined] - wide["VespaG"]).sort_values()
    y = np.arange(len(diff))
    cols = [ACCENT if d > 0 else GRAY_M for d in diff.values]
    ax.axvline(0, color=GRAY_L, lw=0.8, zorder=0)
    ax.hlines(y, 0, diff.values, color=cols, lw=1.0)
    ax.scatter(diff.values, y, color=cols, s=12, lw=0, zorder=3)
    for gene, d in diff.items():
        note(name, f"b {gene}: Ours - VespaG = {d:+.6f} (computed)")
    t = main[(main.scheme == scheme) & (main.model == refined) & (main.reference == "VespaG")].iloc[0]
    tlo, thi = json.loads(t.ci95_delta_spearman)
    ax.errorbar(t.delta_spearman, len(diff) + 0.5,
                xerr=[[t.delta_spearman - tlo], [thi - t.delta_spearman]],
                fmt="D", ms=3.5, color="black", lw=EBAR)
    ax.set_yticks(list(y) + [len(diff) + 0.5])
    ax.set_yticklabels(list(diff.index) + ["mean"])
    ax.set_ylim(-0.7, len(diff) + 1.1)
    ax.set_xlabel(r"$\rho$(Ours) $-$ $\rho$(VespaG)")
    note(name, f"b paired test (file): delta {t.delta_spearman:+.6f} [{tlo:+.6f}, {thi:+.6f}], "
               f"p {t.p_two_sided:.6f}, wins/losses {t.wins}/{t.losses}")
    check(name, "mean_of_differences_equals_file_delta", diff.mean(), t.delta_spearman, 1e-9)

    for letter in "ab":
        panel_label(axd[letter], letter)
    save(fig, name)



def refinement_rows():
    """Rule contrasts before and after refinement, raw and standardized.

    Each audit saved its estimate divided by the SD of the predictor's scores at
    that audit's eligible sites ("standardized within this scorer"). The SD is
    recovered as estimate / standardized estimate and applied to the interval.
    "After" is the final constrained three-rule model (PI + disulfide + constrained
    glycosylation); its PI and C->S scores equal PI+DS exactly (run2 REPORT.md).
    Glycosylation comes from the constrained repair (run2). Charge burial has no
    active coefficient, so it is one point (before = None): the PI model's value,
    which the paper reports for both the raw and the current three-rule model.
    """
    audits = pd.read_csv(os.path.join(PI_AUDIT, "all_audits.csv"))
    ss = load_json(os.path.join(MULTI, "disulfide_internal_consistency_20261025",
                                "internal_consistency.json"))["results"]
    ex = load_json(os.path.join(MULTI, "label_free_audit_20261025",
                                "existing_refinement_audits.json"))["results"]

    def pi(model):
        r = audits[(audits.model == model) & (audits.window == 10) &
                   (audits.comparator_name == "all_uncharged")].iloc[0]
        return r.estimate_raw, r.group_ci_low, r.group_ci_high, r.penalty_sd

    def ds(model):   # Use C-to-S as the declared disulfide-loss probe.
        r = ss[model]["S"]["measurement"]
        return r["estimate"], r["ci95"][0], r["ci95"][1], r["sd_of_values"]

    gr = load_json(os.path.join(GLYCO, "audit_results.json"))["guarded_group_out"]

    def gly(model, arm):
        r = gr[model]["original_all_motifs"][f"{arm}_noncytosolic"]
        g = r["group_resampling"]
        return g["estimate"], g["ci95"][0], g["ci95"][1], r["legacy_whole_protein"]["sd_of_z"]

    def cb(model):
        r = ex[model]["charge_burial"]
        g = r["grouped_endpoint"]
        sd = r["legacy_results"]["z"]["sd_of_z_all_sites"]
        std = r["legacy_results"]["z"]["adjusted_curve_vs_interface"]["bin4"][
            "delta_standardized_within_this_scorer"]
        assert abs(g["estimate"] / sd - std) < 1e-9   # Reuse the saved SD so before and after values share a standardization scale.
        return g["estimate"], g["ci95"][0], g["ci95"][1], sd


    return [
        ("positive-inside", pi("guarded_oof_raw"), pi("guarded_oof__positive_inside+disulfide_loss"), -1),
        ("disulfide loss", ds("guarded_oof_PI"), ds("guarded_oof_PI_plus_disulfide"), +1),
        ("glycosylation (Pro at X)", gly("existing_PI_DS", "A"),
         gly("PI_DS_constraint_glyco", "A"), +1),
        (r"glycosylation (S$\leftrightarrow$T)", gly("existing_PI_DS", "B"),
         gly("PI_DS_constraint_glyco", "B"), -1),
        ("charge burial (control)", None, cb("positive_inside"), +1),
    ]


def fig6_refinement():
    name = "fig6_refinement"
    rows = refinement_rows()

    # Check raw paper values before standardizing the plotted contrasts.
    (pb, pa), (db, da) = (rows[0][1], rows[0][2]), (rows[1][1], rows[1][2])
    check(name, "PI_before", pb[0], 0.00854, 0.000005)
    check(name, "PI_after", pa[0], -0.02029, 0.000005)
    check(name, "PI_after_ci_low", pa[1], -0.02955, 0.000005)
    check(name, "PI_after_ci_high", pa[2], -0.00484, 0.000005)
    check(name, "DS_before", db[0], 0.060, 0.0005)
    check(name, "DS_before_ci_low", db[1], -0.001, 0.0005)
    check(name, "DS_before_ci_high", db[2], 0.104, 0.0005)
    check(name, "DS_after", da[0], 0.126, 0.0005)
    check(name, "DS_after_ci_low", da[1], 0.071, 0.0005)
    check(name, "DS_after_ci_high", da[2], 0.170, 0.0005)
    for (k, state, i), paper in {
            (2, "before", 0): -0.01140, (2, "before", 1): -0.01862, (2, "before", 2): -0.00379,
            (2, "after", 0): +0.04696, (2, "after", 1): +0.01473, (2, "after", 2): +0.08041,
            (3, "before", 0): +0.01811, (3, "before", 1): +0.00858, (3, "before", 2): +0.03096,
            (3, "after", 0): -0.03977, (3, "after", 1): -0.06310, (3, "after", 2): -0.01592}.items():
        arm = "A" if k == 2 else "B"
        part = ("", "_ci_low", "_ci_high")[i]
        check(name, f"GLY_{arm}_{state}{part}", rows[k][1 if state == "before" else 2][i],
              paper, 0.000005)
    check(name, "CB_current", rows[4][2][0], 0.039452, 0.0000005)
    check(name, "CB_current_ci_low", rows[4][2][1], 0.028827, 0.0000005)
    check(name, "CB_current_ci_high", rows[4][2][2], 0.051985, 0.0000005)

    fig, ax = plt.subplots(figsize=(HALF_W, 1.9))
    oriented = []
    for k, (label, before, after, sign) in enumerate(rows):
        y = len(rows) - 1 - k
        pts = []
        for vals, state in ((before, "before"), (after, "after")):
            if vals is None:
                pts.append(None)
                continue
            est, lo, hi, sd = vals
            # Orient standardized contrasts so movement to the right agrees with the declared rule.
            o = sorted((sign * lo / sd, sign * hi / sd))
            pts.append((sign * est / sd, o[0], o[1]))
            note(name, f"{label} {state}: raw {est:+.6f} [{lo:+.6f}, {hi:+.6f}], site SD "
                       f"{sd:.6f}, oriented standardized {sign * est / sd:+.4f} "
                       f"[{o[0]:+.4f}, {o[1]:+.4f}] (computed; sign {sign:+d})")
        (a, alo, ahi) = pts[1]
        oriented += [alo, ahi]
        if pts[0] is None:   # Show one point when the current model has no active correction.
            ax.errorbar(a, y, xerr=[[a - alo], [ahi - a]], fmt="o", color=ACCENT,
                        lw=EBAR, zorder=4)
            continue
        (b, blo, bhi) = pts[0]
        oriented += [blo, bhi]
        ax.errorbar(b, y + 0.15, xerr=[[b - blo], [bhi - b]], fmt="o", mfc="white",
                    mec=GRAY_M, color=GRAY_M, lw=EBAR, zorder=3)
        ax.errorbar(a, y - 0.15, xerr=[[a - alo], [ahi - a]], fmt="o", color=ACCENT,
                    lw=EBAR, zorder=4)
        ax.annotate("", xy=(a, y - 0.15), xytext=(b, y + 0.15),
                    arrowprops=dict(arrowstyle="-|>", lw=0.5, color=GRAY_M,
                                    mutation_scale=5, shrinkA=2, shrinkB=2), zorder=2)
    span = max(abs(v) for v in oriented) * 1.08
    ax.set_xlim(-span, span)
    ax.axvspan(0, span, color=RULE_SHADE, lw=0, zorder=0)
    ax.axvline(0, color=GRAY_L, lw=0.8, zorder=1)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows][::-1])
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xlabel("standardized contrast (right = agrees with rule)")
    note(name, "orientation: positive-inside and glycosylation S<->T multiplied by -1 "
               "(their expected direction is negative); others unchanged")
    save(fig, name)



def figS1_robustness():
    name = "figS1_robustness"
    audits = pd.read_csv(os.path.join(PI_AUDIT, "all_audits.csv"))
    sub = audits[audits.model == "guarded_oof__positive_inside+disulfide_loss"]
    order = {"all_uncharged": "all uncharged", "polar_neutral": "polar",
             "hydrophobic": "hydrophobic", "small": "small", "aromatic": "aromatic"}
    grid = sub.pivot(index="comparator_name", columns="window", values="estimate_raw")
    grid = grid.loc[list(order)]
    grid.index = list(order.values())
    labels = grid.map(lambda v: f"${v:#.2g}$")

    fig, ax = plt.subplots(figsize=(HALF_W, 1.7))
    sns.heatmap(grid, ax=ax, cmap="Blues_r", annot=labels, fmt="", linewidths=0.5,
                annot_kws={"fontsize": 7}, cbar_kws={"label": r"$\Delta$ after refinement"})
    ax.set_xticklabels([rf"$\leq${w}" for w in grid.columns], rotation=0)
    ax.set_xlabel("flank window (residues)")
    ax.set_ylabel("")
    ax.tick_params(axis="y", rotation=0, length=0)
    ax.tick_params(axis="x", length=0)
    n_neg = int((grid.values < 0).sum())
    verdicts = sub.primary_group_verdict.value_counts().to_dict()
    for comp, row in grid.iterrows():
        note(name, comp + ": " + ", ".join(f"<={w} {v:+.6f}" for w, v in row.items()))
    note(name, f"cells negative: {n_neg}/{grid.size}; verdicts {verdicts}")
    check(name, "all_20_negative", n_neg, 20, 0)
    check(name, "all_20_expected_direction_supported",
          verdicts.get("expected_direction_supported", 0), 20, 0)
    save(fig, name)


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")   # Write UTF-8 to retain delta and minus symbols in exported values.
    print("building figures ->", OUT)
    for build in (fig1_overview, fig2_validation, fig3_forest, fig5_accuracy,
                  fig6_refinement, figS1_robustness):
        print(build.__name__)
        build()
    print("\nnot built here: fig4_casestudy (run figures_casestudy.py); fig6_distance "
          "optional, no matching binned outputs")
    if MISMATCH:
        print("\nMISMATCHES:")
        for fig, lines in MISMATCH.items():
            for line in lines:
                print(f"  {fig}: {line}")
    for line in TEXT_ISSUES:
        print("\n" + line)

"""Generate the paper's accuracy, validation and correction figures from saved results."""

import itertools
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES = os.path.join(PROJ, "results", "beat_am_20261016")
FORM_A = os.path.join(PROJ, "results", "form_a_20261008")
OUT = os.path.join(PROJ, "figures", "paper_20261020")

# Use the Okabe-Ito palette to distinguish loss, gain, baseline and reference series.
C_LOSS = "#D55E00"
C_GAIN = "#0072B2"
C_BASE = "#666666"
C_REF = "#009E73"
C_FAIL = "#CC79A7"


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=3)
    return ax


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=200,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  wrote {name}.png / .svg")



def panel_reserved():
    per = pd.read_csv(os.path.join(RES, "confirmation", "confirmation_per_protein.csv"))
    res = json.load(open(os.path.join(RES, "confirmation", "confirmation_result.json")))

    mem = per[per.membrane == 1].sort_values("difference").reset_index(drop=True)
    non = per[per.membrane == 0].sort_values("difference").reset_index(drop=True)
    m = res["membrane_full_scan"]
    a = res["allhuman_full_scan"]
    dev_mem, dev_non = 0.031136, 0.055883

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.4),
                             gridspec_kw={"width_ratios": [1, 1.25]})

    for ax, d, stat, dev, label, n_lab in (
        (axes[0], mem, m, dev_mem, "membrane (the pre-registered claim)", "7"),
        (axes[1], non, a, dev_non, "non-membrane (reported separately)", "20"),
    ):
        style(ax)
        y = np.arange(len(d))
        colors = [C_GAIN if v > 0 else C_LOSS for v in d.difference]

        # Draw intervals behind markers to keep protein-level observations visible.
        ax.axvspan(stat["ci95"][0], stat["ci95"][1], color=C_BASE, alpha=0.13, lw=0)
        ax.axvline(stat["mean_paired_difference"], color="black", lw=1.6, zorder=2)
        ax.axvline(0, color="black", lw=0.9, ls="-", alpha=0.55, zorder=1)
        ax.axvline(dev, color=C_REF, lw=1.5, ls="--", zorder=2)

        ax.hlines(y, 0, d.difference, color=colors, lw=1.3, alpha=0.6, zorder=3)
        ax.scatter(d.difference, y, s=46, color=colors, zorder=4,
                   edgecolor="white", linewidth=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{g}  ({v:,})" for g, v in zip(d.gene, d.n_variants)],
                           fontsize=8.5)
        ax.set_ylim(-0.8, len(d) - 0.2)
        ax.set_xlabel("paired Spearman difference vs AlphaMissense", fontsize=9.5)
        ax.set_title(
            f"{label}\nn = {n_lab}   mean {stat['mean_paired_difference']:+.4f}   "
            f"95% CI [{stat['ci95'][0]:+.4f}, {stat['ci95'][1]:+.4f}]   "
            f"exact p = {stat['signflip_p']:.3f}",
            fontsize=9.5, loc="left", pad=10)

        pad = 0.035
        lo = min(d.difference.min(), stat["ci95"][0], 0, dev) - pad
        hi = max(d.difference.max(), stat["ci95"][1], 0, dev) + pad
        ax.set_xlim(lo, hi)


    handles = [
        plt.Line2D([], [], color="black", lw=1.6, label="reserved-set mean"),
        plt.Line2D([], [], color=C_BASE, lw=7, alpha=0.3, label="95% bootstrap CI"),
        plt.Line2D([], [], color=C_REF, lw=1.5, ls="--",
                   label="development estimate"),
        plt.Line2D([], [], marker="o", color=C_LOSS, lw=0, markersize=7,
                   label="worse than AlphaMissense"),
        plt.Line2D([], [], marker="o", color=C_GAIN, lw=0, markersize=7,
                   label="better than AlphaMissense"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
               fontsize=8.5, bbox_to_anchor=(0.5, -0.045))
    fig.suptitle(
        "Development said the model won; the reserved proteins said it lost",
        fontsize=12.5, y=1.015, x=0.055, ha="left")
    fig.text(0.055, 0.965,
             "Each point is one protein the model was never fitted on. "
             "Protein counts in parentheses are measured variants.",
             fontsize=8.8, color="#444444", ha="left")
    save(fig, "fig3_reserved_proteins")

    return mem, m, dev_mem



def panel_validation():
    est = pd.read_csv(os.path.join(FORM_A, "estimates.csv"))
    crit = json.load(open(os.path.join(FORM_A, "criteria.json")))
    slope = crit["S1"]["slope_through_origin"]
    s2_max = crit["S2_identity"]["max_abs_estimate"]

    # Use S1 for recovery and reserve the side-blind S2 model for the refusal check.
    s1 = est[est.student == "S1"]
    s2 = est[est.student == "S2"]

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.0),
                             gridspec_kw={"width_ratios": [1.15, 1]})


    ax = style(axes[0])
    lim = max(abs(s1.oracle_raw).max(), abs(s1.estimate_raw).max()) * 1.18
    ax.axhline(0, color="black", lw=0.8, alpha=0.5)
    ax.axvline(0, color="black", lw=0.8, alpha=0.5)
    xs = np.array([-lim, lim])
    ax.plot(xs, xs, color=C_BASE, lw=1.0, ls=":", label="perfect recovery (slope 1)")
    ax.plot(xs, slope * xs, color=C_GAIN, lw=1.8,
            label=f"fitted slope through origin = {slope:.4f}")

    planted = s1[s1.oracle_raw != 0]
    null = s1[s1.oracle_raw == 0]
    ax.scatter(planted.oracle_raw, planted.estimate_raw, s=34, color=C_GAIN,
               alpha=0.75, edgecolor="white", linewidth=0.6,
               label=f"planted rule, {len(planted)} runs")
    ax.scatter(null.oracle_raw, null.estimate_raw, s=34, color=C_BASE,
               alpha=0.65, marker="s", edgecolor="white", linewidth=0.6,
               label=f"no rule planted, {len(null)} runs")
    ax.scatter(s2.oracle_raw, s2.estimate_raw, s=60, color=C_FAIL, marker="X",
               zorder=5, edgecolor="white", linewidth=0.7,
               label=f"side-blind model, {len(s2)} runs")

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("rule magnitude planted in the training data", fontsize=9.5)
    ax.set_ylabel("magnitude the audit recovered\nfrom the trained model", fontsize=9.5)
    ax.set_title("The audit recovers what was planted", fontsize=10.5,
                 loc="left", pad=8)
    # Place the legend in the empty lower-right quadrant to avoid the recovery line.
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    ax.set_aspect("equal")


    ax = style(axes[1])
    groups = [
        ("side-blind model\n(cannot see side)", np.abs(s2.estimate_raw.values), C_FAIL),
        ("full model,\nno rule planted", np.abs(null.estimate_raw.values), C_BASE),
        ("full model,\nrule planted", np.abs(planted.estimate_raw.values), C_GAIN),
    ]
    for i, (lab, vals, col) in enumerate(groups):
        vals = np.maximum(vals, 1e-20)
        ax.scatter(vals, np.full(len(vals), i) + np.random.default_rng(0).normal(0, 0.055, len(vals)),
                   s=30, color=col, alpha=0.7, edgecolor="white", linewidth=0.5)
    ax.set_xscale("log")
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels([g[0] for g in groups], fontsize=8.5)
    ax.set_ylim(-0.6, len(groups) - 0.4)
    ax.set_xlabel("|recovered magnitude|   (log scale)", fontsize=9.5)
    ax.axvline(s2_max, color=C_FAIL, lw=1.2, ls="--")
    ax.annotate(f"largest side-blind\nreading: {s2_max:.2e}",
                xy=(s2_max, 0.42), xytext=(s2_max * 60, 0.62),
                fontsize=8, color=C_FAIL,
                arrowprops=dict(arrowstyle="-", color=C_FAIL, lw=0.8))
    ax.set_title("and reads zero when the mechanism is absent",
                 fontsize=10.5, loc="left", pad=8)

    fig.suptitle("Validating the instrument against planted answers",
                 fontsize=12.5, y=1.045, x=0.055, ha="left")
    fig.text(0.055, 0.985,
             "Rules of known size were planted in training data, models were trained on it, "
             "and the audit read the models' predictions.",
             fontsize=8.8, color="#444444", ha="left")
    save(fig, "fig1_instrument_validation")

    return slope, s2_max, len(planted), len(null), len(s2)



def panel_enforcement():
    frozen = json.load(open(os.path.join(RES, "rule", "beta_frozen_m0.045.json")))
    nominal = json.load(open(os.path.join(RES, "rule", "beta_both_splits.json")))
    cost = json.load(open(os.path.join(RES, "rule",
                                       "accuracy_cost_beta0.017660.json")))

    def sd_units(entry):
        sd = entry["penalty_sd"]
        lo, hi = entry["ci95_raw"]
        return entry["contrast_standardized"], lo / sd, hi / sd

    rows = []
    for split, label in (("sequence", "sequence-similarity split"),
                         ("symbol", "gene-family split")):
        h = frozen["held_out"][split]
        n = nominal["held_out"][split]
        rows.append((label, [
            ("AlphaMissense", sd_units(h["alphamissense"]), C_BASE),
            ("model, uncorrected", sd_units(h["model_pre"]), C_LOSS),
            ("model, corrected\n(nominal margin)", sd_units(n["model_post_joint"]), C_FAIL),
            ("model, corrected\n(frozen, 1.5x margin)", sd_units(h["model_post_joint"]), C_GAIN),
        ], h["alphamissense"]["n_proteins"], h["alphamissense"]["n_sites"]))

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.6), sharex=True)
    for panel, (ax, (label, entries, npro, nsite)) in enumerate(zip(axes, rows)):
        style(ax)
        ax.axvspan(-0.22, -0.03, color=C_REF, alpha=0.085, lw=0)
        ax.axvline(0, color="black", lw=0.9, alpha=0.6)
        ax.axvline(-0.03, color=C_REF, lw=1.3, ls="--")

        for i, (name, (pt, lo, hi), col) in enumerate(entries):
            y = len(entries) - 1 - i
            ax.plot([lo, hi], [y, y], color=col, lw=2.4, solid_capstyle="round",
                    alpha=0.85)
            ax.scatter([pt], [y], s=62, color=col, zorder=4,
                       edgecolor="white", linewidth=0.9)
            ok = "pass" if hi < 0 and pt <= -0.03 else "fail"
            ax.text(hi + 0.006, y, f"{pt:+.4f}  {ok}", va="center", fontsize=8,
                    color=col if ok == "pass" else "#777777")

        ax.set_yticks(range(len(entries)))

        labels = [n for n, _, _ in entries][::-1] if panel == 0 else [""] * len(entries)
        ax.set_yticklabels(labels, fontsize=8.5)
        ax.set_ylim(-0.6, len(entries) - 0.22)
        ax.set_xlabel("audited contrast, in the scorer's own SD unit", fontsize=9.5)
        ax.set_title(f"{label}\n{npro:,} held-out proteins, {nsite:,} sites",
                     fontsize=9.8, loc="left", pad=8)
        ax.set_xlim(-0.17, 0.29)
        ax.annotate("rule satisfied", xy=(-0.10, 3.48), fontsize=8,
                    color=C_REF, ha="center")
        ax.annotate("rule not satisfied", xy=(0.12, 3.48), fontsize=8,
                    color="#777777", ha="center")

    overall = cost["scopes"]["all"]
    fig.suptitle("Enforcement moves the contrast across zero on both splits, "
                 "for no measurable accuracy cost",
                 fontsize=12.0, y=1.08, x=0.055, ha="left")
    fig.text(0.055, 1.012,
             f"Accuracy change from enforcement, {overall['n']} development proteins: "
             f"{overall['post_mean_vs_am'] - overall['pre_mean_vs_am']:+.5f} Spearman "
             f"(rank fidelity {frozen['fidelity_full_panel']['spearman']:.5f}). "
             "A negative contrast agrees with the rule.",
             fontsize=8.8, color="#444444", ha="left")
    save(fig, "fig2_enforcement")

    return rows, overall



def signflip_p(diffs, signs):
    """Exact two-sided sign-flip p for one vector of paired differences."""
    obs = abs(diffs.mean())
    means = np.abs(signs @ np.abs(diffs)) / len(diffs)
    return (means >= obs - 1e-15).mean()


def power_at(delta, sd, n=7, reps=4000, nboot=400, seed=11):
    """Power of the declared three-part win condition.

    All three conditions are invariant to a positive rescaling, so power
    depends only on d = delta / sd.
    """
    rng = np.random.default_rng(seed)
    signs = np.array(list(itertools.product([-1, 1], repeat=n)), dtype=float)
    d = delta / sd
    wins = 0
    for _ in range(reps):
        x = rng.normal(d, 1.0, size=n)
        if x.mean() <= 0:
            continue
        if signflip_p(x, signs) >= 0.05:
            continue
        idx = rng.integers(0, n, size=(nboot, n))
        if np.quantile(x[idx].mean(axis=1), 0.025) > 0:
            wins += 1
    return wins / reps


def panel_power(mem, mstat, dev_mem):
    sd_assumed, sd_observed = 0.025, 0.0760
    deltas = np.linspace(0.0, 0.14, 29)

    curves = {}
    for lab, sd in (("pre-registered SD = 0.025", sd_assumed),
                    ("observed SD = 0.0760", sd_observed)):
        curves[lab] = [power_at(dl, sd) for dl in deltas]

    # Use more repetitions for annotated values to stabilize their printed precision.
    p_assumed = power_at(dev_mem, sd_assumed, reps=20000)
    p_observed = power_at(dev_mem, sd_observed, reps=20000)

    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.6),
                            gridspec_kw={"width_ratios": [1.3, 1]})


    ax = style(axes[0])
    for (lab, ys), col, ls in zip(curves.items(), (C_REF, C_LOSS), ("-", "-")):
        ax.plot(deltas, ys, color=col, lw=2.0, ls=ls, label=lab)
    ax.axvline(dev_mem, color=C_BASE, lw=1.2, ls="--")
    ax.scatter([dev_mem, dev_mem], [p_assumed, p_observed],
               s=70, color=[C_REF, C_LOSS], zorder=5,
               edgecolor="white", linewidth=0.9)
    ax.annotate(f"{p_assumed:.2f}", (dev_mem, p_assumed), xytext=(11, -3),
                textcoords="offset points", fontsize=9.5, color=C_REF, weight="bold")
    ax.annotate(f"{p_observed:.2f}", (dev_mem, p_observed), xytext=(11, -4),
                textcoords="offset points", fontsize=9.5, color=C_LOSS, weight="bold")
    ax.annotate("development\neffect +0.0313", (dev_mem, 0.99),
                xytext=(-7, -2), textcoords="offset points", ha="right",
                va="top", fontsize=8.2, color="#444444")
    ax.axhline(0.8, color=C_BASE, lw=0.8, ls=":", alpha=0.8)
    ax.text(0.138, 0.815, "80%", fontsize=7.8, color=C_BASE, ha="right")
    ax.set_xlabel("true effect size (paired Spearman difference)", fontsize=9.5)
    ax.set_ylabel("power of the declared win condition\nat n = 7", fontsize=9.5)
    ax.set_ylim(0, 1.02)
    ax.set_xlim(0, 0.14)
    ax.legend(fontsize=8.5, frameon=False, loc="lower right")
    ax.set_title("Power of the declared win condition at n = 7",
                 fontsize=10.2, loc="left", pad=8)


    ax = style(axes[1])
    se = sd_observed / np.sqrt(7)
    xs = np.linspace(dev_mem - 4.3 * se, dev_mem + 3.0 * se, 400)
    ys = np.exp(-0.5 * ((xs - dev_mem) / se) ** 2)
    ax.plot(xs, ys, color=C_REF, lw=1.8)
    ax.fill_between(xs, 0, ys, color=C_REF, alpha=0.12)
    obs = mstat["mean_paired_difference"]
    z = (dev_mem - obs) / se
    ax.axvline(obs, color=C_LOSS, lw=2.0)
    ax.axvline(0, color="black", lw=0.9, alpha=0.5)
    ax.annotate(f"observed\n{obs:+.4f}\n{z:.1f} SE below", (obs, 0.72),
                xytext=(10, 0), textcoords="offset points", fontsize=8.6,
                color=C_LOSS, weight="bold")
    ax.annotate("sampling distribution if the\ndevelopment effect were real",
                (dev_mem + 1.15 * se, 0.62), fontsize=8.2, color="#2a6e5c")
    ax.set_yticks([])
    ax.set_ylim(0, 1.25)
    ax.set_xlabel("mean paired difference at n = 7", fontsize=9.5)
    ax.set_title("Low power does not explain the loss", fontsize=10.2,
                 loc="left", pad=8)

    fig.suptitle("Development understated the variance, so the study was far "
                 "weaker than it claimed",
                 fontsize=12.5, y=1.10, x=0.055, ha="left")
    fig.text(0.055, 1.035,
             "Power is simulated under the declared three-part win condition and depends only on effect / SD.\n"
             "The pre-registration recorded 0.74 at the assumed SD; this independent recomputation gives 0.749 +- 0.003.",
             fontsize=8.8, color="#444444", ha="left", va="top")
    save(fig, "fig4_variance_and_power")

    return p_assumed, p_observed, z


if __name__ == "__main__":
    print("building figures ->", OUT)
    slope, s2_max, n_pl, n_nu, n_s2 = panel_validation()
    rows, overall = panel_enforcement()
    mem, mstat, dev_mem = panel_reserved()
    p_a, p_o, z = panel_power(mem, mstat, dev_mem)

    print("\nvalues drawn, for cross-checking against the draft:")
    print(f"  fig3  membrane mean {mstat['mean_paired_difference']:+.6f} "
          f"CI [{mstat['ci95'][0]:+.6f}, {mstat['ci95'][1]:+.6f}] "
          f"p {mstat['signflip_p']} improved {mstat['fraction_improved']}")
    print(f"  fig1  slope {slope:.7f}  side-blind max |est| {s2_max:.4e}  "
          f"runs planted/null/side-blind {n_pl}/{n_nu}/{n_s2}")
    for label, entries, npro, nsite in rows:
        pts = "  ".join(f"{n.splitlines()[0]}={v[0]:+.4f}" for n, v, _ in entries)
        print(f"  fig2  {label}: {pts}")
    print(f"  fig2  enforcement cost (all, n={overall['n']}): "
          f"{overall['post_mean_vs_am'] - overall['pre_mean_vs_am']:+.6f}")
    print(f"  fig4  power at +0.0313: SD 0.025 -> {p_a:.2f}; "
          f"SD 0.0760 -> {p_o:.2f};  observed is {z:.2f} SE below")

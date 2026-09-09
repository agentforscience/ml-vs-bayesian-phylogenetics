"""Figures for the ML-vs-Bayesian outbreak-phylogenetics report.

Design rules followed (validated categorical palette, `--pairs all`, light mode):
  * at most three categorical hues anywhere; the scientific contrast that carries
    colour is BETWEEN-method vs WITHIN-method null.  The four different nulls are
    separated by POSITION, not by four more hues.
  * one axis per panel, no dual scales; recessive grid; legend whenever >= 2
    series are drawn; direct labels where a number matters.
  * aqua (#1baf7a) is below 3:1 on the light surface, so wherever it is used the
    mark is accompanied by a visible label (the palette's "relief rule").
"""
import collections, json, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

FIGDIR = "figures"
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8985"
GRID = "#e3e2dd"
BETWEEN = "#2a78d6"     # slot 1 blue   -- between-method (ML vs Bayesian)
NULL = "#eb6834"        # slot 2 orange -- within-method null
THIRD = "#1baf7a"       # slot 3 aqua   -- third series, always directly labelled

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "font.size": 9,
    "axes.edgecolor": GRID, "axes.labelcolor": INK, "axes.titlesize": 10,
    "axes.titleweight": "bold", "axes.titlecolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.axisbelow": True, "legend.frameon": False, "figure.dpi": 130,
    "axes.spines.top": False, "axes.spines.right": False,
})

# information gradient from datasets/_dataset_stats.json (parsimony-informative
# sites per taxon) -- the covariate the whole study is ordered by
def pis_per_taxon():
    stats = json.load(open("datasets/_dataset_stats.json"))
    return {s["dataset"]: s["alignment"]["informative_sites_per_taxon"] for s in stats}


SHORT = {"sarscov2_usa_early2020": "SARS-CoV-2\nUSA 2020",
         "fmdv_uk_2007": "FMDV\nUK 2007",
         "mers_korea_2015": "MERS\nKorea 2015",
         "ebov_sierraleone_2014": "EBOV\nSL 2014",
         "h5n1_dairy_cattle_2024": "H5N1\ncattle 2024",
         "ebov_westafrica_makona": "EBOV\nW.Africa",
         "fmdv_uk_2001": "FMDV\nUK 2001",
         "zika_americas_2015_2016": "Zika\nAmericas"}


def load():
    topo = {r["dataset"]: r for r in json.load(open("results/topology_summary.json"))
            if "error" not in r}
    epi = {r["dataset"]: r for r in json.load(open("results/epi_summary.json"))
           if "error" not in r}
    pis = pis_per_taxon()
    order = sorted(topo, key=lambda d: pis[d])   # low information -> high
    return topo, epi, pis, order


def _finish(fig, name, caption):
    os.makedirs(FIGDIR, exist_ok=True)
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/{name}.png", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {FIGDIR}/{name}.png  -- {caption}")


# ---------------------------------------------------------------- figure 1
def fig1_scale(topo, pis, order):
    """The headline: between-method distance against four within-method nulls."""
    fig, ax = plt.subplots(figsize=(10, 5.2))
    xs = np.arange(len(order))
    w = 0.17
    # null bands drawn as vertical ranges at offset positions
    lanes = [("null_ml_replicate", "re-run IQ-TREE\n(new seed)", -1.5),
             ("null_mcmc_replicate", "re-run MrBayes\n(chain 1 vs 2)", -0.5),
             ("null_bootstrap_spread", "ML tree vs its own\nUFBoot replicates", 0.5),
             ("null_posterior_spread", "Bayes consensus vs its\nown posterior trees", 1.5)]
    for key, lab, off in lanes:
        lo, hi, md = [], [], []
        for d in order:
            v = topo[d].get(key)
            if not v:
                lo.append(np.nan); hi.append(np.nan); md.append(np.nan); continue
            if "nrf" in v and "median" not in v:            # scalar null (N1)
                lo.append(v["nrf"]); hi.append(v["nrf"]); md.append(v["nrf"])
            else:
                lo.append(v["q025"]); hi.append(v["q975"]); md.append(v["median"])
        x = xs + off * w
        ax.vlines(x, lo, hi, color=NULL, lw=3.2, alpha=0.35,
                  solid_capstyle="round", zorder=2)
        ax.plot(x, md, "o", ms=5, color=NULL, mec=SURFACE, mew=1.2, zorder=3)
    bet = [topo[d]["between"]["nrf"] for d in order]
    ax.plot(xs, bet, "D", ms=9, color=BETWEEN, mec=SURFACE, mew=1.4, zorder=5,
            label="BETWEEN methods: ML tree vs Bayesian consensus")
    for x, y in zip(xs, bet):
        ax.annotate(f"{y:.2f}", (x, y), textcoords="offset points", xytext=(0, 11),
                    ha="center", fontsize=8, color=BETWEEN, fontweight="bold")
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{SHORT[d]}\n{pis[d]:.2f} PIS/taxon" for d in order], fontsize=7.5)
    ax.set_ylabel("normalised Robinson–Foulds distance")
    ax.set_ylim(0, max(1.02, max(bet) * 1.25))
    ax.set_title("Between-method topological distance vs four within-method nulls\n"
                 "(datasets ordered by phylogenetic information per taxon, left = least)")
    ax.legend(handles=[
        Line2D([], [], marker="D", ls="", color=BETWEEN, ms=8, label="BETWEEN methods (ML vs Bayesian)"),
        Line2D([], [], marker="o", ls="-", color=NULL, lw=3, alpha=0.6, ms=5,
               label="WITHIN-method nulls: median and 95% range\n(4 lanes per dataset, left to right: "
                     "IQ-TREE re-seed · MrBayes chain 1 vs 2 · UFBoot spread · posterior spread)")],
        loc="upper center", bbox_to_anchor=(0.5, -0.22), fontsize=8)
    _finish(fig, "fig1_between_vs_within_null",
            "between-method RF against within-method nulls")


# ---------------------------------------------------------------- figure 2
def fig2_two_scales(topo, pis, order):
    """Major clade structure vs fine-scale branching order."""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.4))
    xs = np.arange(len(order))
    w = 0.36
    ml = [topo[d]["major_fine"]["strong_ml_recovered_by_mb"] or 0 for d in order]
    mb = [topo[d]["major_fine"]["strong_mb_recovered_by_ml"] or 0 for d in order]
    a1.bar(xs - w / 2, ml, w * 0.94, color=BETWEEN, label="UFBoot ≥ 95 splits found in the Bayesian tree")
    a1.bar(xs + w / 2, mb, w * 0.94, color=NULL, label="PP ≥ 0.95 splits found in the ML tree")
    a1.axhline(0.85, color=INK, lw=1.2, ls=(0, (4, 3)))
    a1.annotate("0.85 — the hypothesised agreement floor", (len(order) - 0.4, 0.855),
                ha="right", va="bottom", fontsize=8, color=INK)
    a1.set_ylim(0, 1.12); a1.set_ylabel("fraction recovered by the other method")
    a1.set_title("MAJOR clade structure\n(confidently supported splits)")
    a1.legend(loc="lower center", bbox_to_anchor=(0.5, -0.55), fontsize=8)

    n = [topo[d]["n_taxa"] for d in order]
    fine = [topo[d]["major_fine"]["fine_scale_rf"] / (2 * (t - 3))
            for d, t in zip(order, n)]
    a2.bar(xs, fine, 0.6, color=BETWEEN)
    for x, y in zip(xs, fine):
        a2.annotate(f"{y:.2f}", (x, y), textcoords="offset points", xytext=(0, 3),
                    ha="center", fontsize=8, color=INK2)
    a2.set_ylim(0, max(fine) * 1.25 if fine else 1)
    a2.set_ylabel("weakly-supported splits that differ\n(as a fraction of max RF)")
    a2.set_title("FINE-SCALE branching order\n(everything not confidently supported)")
    for ax in (a1, a2):
        ax.set_xticks(xs)
        ax.set_xticklabels([SHORT[d] for d in order], fontsize=7.5)
    fig.suptitle("The two scales the hypothesis distinguishes behave completely differently",
                 fontsize=11, fontweight="bold", y=1.03)
    _finish(fig, "fig2_two_scales", "major vs fine-scale agreement")


# ---------------------------------------------------------------- figure 3
def fig3_calibration(topo, order):
    """Bootstrap proportion vs posterior probability, one panel per dataset."""
    ncol = 4
    nrow = math.ceil(len(order) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(11, 2.8 * nrow), squeeze=False)
    for i, d in enumerate(order):
        ax = axes[i // ncol][i % ncol]
        pts = topo[d].get("split_support_pairs", [])
        if pts:
            u = np.array([p[0] for p in pts]); p = np.array([p[1] for p in pts])
            ax.plot([0, 1], [0, 1], color=MUTED, lw=1, ls=(0, (4, 3)), zorder=1)
            ax.scatter(u, p, s=9, color=BETWEEN, alpha=0.45, lw=0, zorder=2)
            above = float(np.mean(p > u))
            ax.annotate(f"PP > UFBoot for {above:.0%} of splits\nn = {len(pts)}",
                        (0.04, 0.96), xycoords="axes fraction", va="top",
                        fontsize=7.5, color=INK2)
        ax.set_xlim(-0.03, 1.03); ax.set_ylim(-0.03, 1.03)
        ax.set_title(SHORT[d].replace("\n", " "), fontsize=9)
        if i % ncol == 0: ax.set_ylabel("posterior probability")
        if i // ncol == nrow - 1: ax.set_xlabel("UFBoot proportion")
    for j in range(len(order), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle("Split support: the two methods do not measure the same quantity\n"
                 "(dashed line = equality; points above it are splits the posterior is more confident about)",
                 fontsize=10.5, fontweight="bold", y=1.0)
    _finish(fig, "fig3_support_calibration", "PP vs UFBoot per split")


# ---------------------------------------------------------------- figure 4
def fig4_conclusions(epi, order):
    """Conclusion-change matrix: does the between-method flip beat the null?"""
    comps = [("BETWEEN_ml_vs_bayes", "ML vs Bayesian\n(the tested comparison)"),
             ("PRIOR_gammadir_vs_exponential", "Bayesian: compound-Dirichlet\nvs exponential brlen prior"),
             ("NULL_mcmc_replicate", "NULL: MrBayes chain 1\nvs chain 2"),
             ("__mlnull__", "NULL: IQ-TREE re-seed\n(mean over 10 pairs)")]
    rules = [("tier1_who", "WHO infected whom\n(index unit or edge reversed)"),
             ("tier2_how_many", "HOW MANY introductions"),
             ("tier3_scope", "cluster scope"),
             ("ANY", "ANY conclusion change")]
    grid = np.full((len(rules), len(comps), len(order)), np.nan)
    for k, d in enumerate(order):
        e = epi[d]
        no_trait = e.get("primary_trait") is None
        for j, (ck, _) in enumerate(comps):
            for i, (rk, _) in enumerate(rules):
                if no_trait and rk in ("tier1_who", "tier2_how_many"):
                    continue          # not applicable: no epidemiological grouping variable
                if ck == "__mlnull__":
                    keys = [x for x in e["comparisons"] if x.startswith("NULL_ml_reseed")]
                    grid[i, j, k] = np.mean([e["comparisons"][x][rk] for x in keys])
                elif ck in e["comparisons"]:
                    grid[i, j, k] = float(e["comparisons"][ck][rk])
    fig, axes = plt.subplots(len(rules), 1, figsize=(10, 1.55 * len(rules) + 1.2),
                             sharex=True)
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("seq", ["#f4f4f1", BETWEEN])  # single hue
    for i, (rk, rlab) in enumerate(rules):
        ax = axes[i]
        im = ax.imshow(grid[i], cmap=cmap, vmin=0, vmax=1, aspect="auto")
        ax.set_yticks(range(len(comps)))
        ax.set_yticklabels([c[1] for c in comps], fontsize=7.5)
        ax.set_ylabel(rlab, fontsize=8, rotation=0, ha="right", va="center",
                      labelpad=8, color=INK, fontweight="bold")
        ax.grid(False)
        for j in range(len(comps)):
            for k in range(len(order)):
                v = grid[i, j, k]
                if np.isnan(v):
                    ax.text(k, j, "n/a", ha="center", va="center", fontsize=7, color=MUTED)
                else:
                    ax.text(k, j, ("yes" if v == 1 else "no" if v == 0 else f"{v:.1f}"),
                            ha="center", va="center", fontsize=7.5,
                            color=(SURFACE if v > 0.55 else INK2),
                            fontweight=("bold" if v == 1 else "normal"))
    axes[-1].set_xticks(range(len(order)))
    axes[-1].set_xticklabels([SHORT[d] for d in order], fontsize=7.5)
    fig.suptitle("Does the epidemiological conclusion change?  Between-method vs within-method null\n"
                 "(cells are labelled, so colour is never the only encoding)",
                 fontsize=10.5, fontweight="bold", y=1.0)
    _finish(fig, "fig4_conclusion_changes", "conclusion-change matrix")


# ---------------------------------------------------------------- figure 5
def fig5_distributions(epi, order):
    """Introduction-count and index-unit distributions: ML bootstrap vs posterior."""
    have = [d for d in order if epi[d].get("distributions", {}).get("intro_total")]
    if not have:
        print("  (fig5 skipped: no trait-based distributions)"); return
    ncol = min(4, len(have)); nrow = math.ceil(len(have) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 2.9 * nrow), squeeze=False)
    for i, d in enumerate(have):
        ax = axes[i // ncol][i % ncol]
        dd = epi[d]["distributions"]["intro_total"]
        pa, pb = dd["ufboot"]["pmf"], dd["posterior"]["pmf"]
        keys = sorted({int(k) for k in pa} | {int(k) for k in pb})
        x = np.arange(len(keys)); w = 0.4
        ax.bar(x - w / 2, [pa.get(str(k), 0) for k in keys], w * 0.94,
               color=BETWEEN, label="ML (1000 UFBoot)")
        ax.bar(x + w / 2, [pb.get(str(k), 0) for k in keys], w * 0.94,
               color=NULL, label="Bayesian (posterior)")
        ax.set_xticks(x); ax.set_xticklabels(keys, fontsize=7)
        ax.set_title(f'{SHORT[d].replace(chr(10)," ")}\noverlap = {dd["overlap"]:.2f}, '
                     f'p = {dd["test"]["p"]:.1e}', fontsize=8.5)
        if i % ncol == 0: ax.set_ylabel("probability")
        ax.set_xlabel(f'independent introductions ({epi[d]["primary_trait"]})', fontsize=7.5)
        if i == 0: ax.legend(fontsize=7.5, loc="upper right")
    for j in range(len(have), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle("Propagating uncertainty collapses the point-estimate disagreement:\n"
                 "the two methods' DISTRIBUTIONS of the introduction count",
                 fontsize=10.5, fontweight="bold", y=1.02)
    _finish(fig, "fig5_introduction_distributions", "introduction-count distributions")


# ---------------------------------------------------------------- figure 6
def fig6_gradient(topo, pis, order):
    """Does phylogenetic information predict between-method divergence?"""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.1))
    x = np.array([pis[d] for d in order])
    b = np.array([topo[d]["between"]["nrf"] for d in order])
    nml = np.array([topo[d]["null_ml_replicate"]["mean"] for d in order])
    a1.scatter(x, b, s=70, marker="D", color=BETWEEN, zorder=3,
               label="BETWEEN methods (ML vs Bayesian)")
    a1.scatter(x, nml, s=55, marker="o", color=NULL, zorder=3,
               label="WITHIN ML (mean over re-seed pairs)")
    for xi, yi, d in zip(x, b, order):
        a1.annotate(SHORT[d].replace("\n", " "), (xi, yi), fontsize=6.5,
                    color=INK2, textcoords="offset points", xytext=(5, 5))
    from scipy import stats as st
    r1 = st.spearmanr(x, b); r2 = st.spearmanr(x, nml)
    a1.set_xscale("log"); a1.set_xlabel("parsimony-informative sites per taxon (log)")
    a1.set_ylabel("normalised RF")
    a1.set_title(f"Divergence vs information\nSpearman ρ = {r1.statistic:.2f} "
                 f"(p = {r1.pvalue:.2f}) between; {r2.statistic:.2f} within")
    a1.legend(fontsize=8, loc="best")

    ratio = b / np.maximum(nml, 1e-9)
    a2.scatter(x, ratio, s=70, color=BETWEEN, zorder=3)
    a2.axhline(1.0, color=INK, lw=1.2, ls=(0, (4, 3)))
    a2.annotate("ratio = 1: model choice moves the tree\nno more than re-running one method",
                (x.min(), 1.03), fontsize=7.5, color=INK, va="bottom")
    for xi, yi, d in zip(x, ratio, order):
        a2.annotate(SHORT[d].replace("\n", " "), (xi, yi), fontsize=6.5,
                    color=INK2, textcoords="offset points", xytext=(5, 4))
    a2.set_xscale("log"); a2.set_xlabel("parsimony-informative sites per taxon (log)")
    a2.set_ylabel("between-method RF ÷ within-ML null RF")
    a2.set_title("Effect size relative to the null")
    _finish(fig, "fig6_information_gradient", "divergence vs information content")


# ---------------------------------------------------------------- figure 7
def fig7_edge_support(epi, order):
    """Directed transmission-edge support: bootstrap vs posterior."""
    have = [d for d in order if epi[d].get("distributions", {}).get("edge_support")]
    if not have:
        print("  (fig7 skipped)"); return
    ncol = min(4, len(have)); nrow = math.ceil(len(have) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 2.9 * nrow), squeeze=False)
    for i, d in enumerate(have):
        ax = axes[i // ncol][i % ncol]
        tab = epi[d]["distributions"]["edge_support"]["table"]
        u = np.array([v[0] for v in tab.values()]); p = np.array([v[1] for v in tab.values()])
        ax.plot([0, 1], [0, 1], color=MUTED, lw=1, ls=(0, (4, 3)))
        ax.scatter(u, p, s=16, color=BETWEEN, alpha=0.6, lw=0)
        ax.set_xlim(-0.03, 1.03); ax.set_ylim(-0.03, 1.03)
        ax.set_title(f'{SHORT[d].replace(chr(10)," ")}\n{len(tab)} directed edges', fontsize=8.5)
        if i % ncol == 0: ax.set_ylabel("posterior probability")
        if i // ncol == nrow - 1: ax.set_xlabel("UFBoot proportion")
    for j in range(len(have), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle("Support for each directed who-infected-whom edge, ML bootstrap vs Bayesian posterior",
                 fontsize=10.5, fontweight="bold", y=1.01)
    _finish(fig, "fig7_edge_support", "transmission edge support")


if __name__ == "__main__":
    topo, epi, pis, order = load()
    print(f"figures for {len(order)} datasets")
    fig1_scale(topo, pis, order)
    fig2_two_scales(topo, pis, order)
    fig3_calibration(topo, order)
    fig4_conclusions(epi, order)
    fig5_distributions(epi, order)
    fig6_gradient(topo, pis, order)
    fig7_edge_support(epi, order)
    print("done")

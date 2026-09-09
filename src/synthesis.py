"""Cross-dataset synthesis: the three clauses of the hypothesis, tested.

The hypothesis has three quantitative clauses.  Each is stated here as a null,
tested, and reported with an interval -- plus, where the clause is not
operationally defined in the original wording, under every reasonable reading
rather than the one that happens to agree.

  C1  "agree on major clade structure in over 85% of cases"
  C2  "disagree on fine-scale branching order in 30-50% of outbreak datasets"
  C3  "disagreements change the inferred transmission chain in at least 2 of 8
       datasets in ways that would alter public health response recommendations"

C3 additionally gets the null the amended plan requires: how many datasets are
conclusion-changing when the SAME method is merely re-run.

Writes results/synthesis.json and results/synthesis_tables.md.
"""
import collections, itertools, json, math, os, statistics, sys
import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def wilson(k, n, z=1.959963985):
    """Wilson score interval -- correct for proportions near 0 or 1, which is
    exactly where these recovery fractions sit."""
    if n == 0:
        return (None, None)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (round(max(0.0, c - h), 4), round(min(1.0, c + h), 4))


def load():
    topo = {r["dataset"]: r for r in json.load(open("results/topology_summary.json"))
            if "error" not in r}
    epi = {r["dataset"]: r for r in json.load(open("results/epi_summary.json"))
           if "error" not in r}
    pis = {s["dataset"]: s["alignment"]["informative_sites_per_taxon"]
           for s in json.load(open("datasets/_dataset_stats.json"))}
    return topo, epi, pis


# --------------------------------------------------------------------- C1
def clause1(topo):
    """Major clade structure: >85% agreement?"""
    per = {}
    kml = kmb = nml = nmb = 0
    for d, r in topo.items():
        mf = r["major_fine"]
        a, b = mf["strong_ml_recovered_by_mb"], mf["strong_mb_recovered_by_ml"]
        per[d] = dict(n_strong_ml=mf["n_strong_ml"], n_strong_mb=mf["n_strong_mb"],
                      ml_recovered_by_mb=a, mb_recovered_by_ml=b,
                      strong_conflicts=mf["strong_conflict"],
                      both_ge_085=(a is not None and b is not None
                                   and a >= 0.85 and b >= 0.85))
        if a is not None:
            kml += round(a * mf["n_strong_ml"]); nml += mf["n_strong_ml"]
        if b is not None:
            kmb += round(b * mf["n_strong_mb"]); nmb += mf["n_strong_mb"]
    k, n = kml + kmb, nml + nmb
    # one-sided exact binomial test of H0: p <= 0.85
    p_bin = stats.binomtest(int(k), int(n), 0.85, alternative="greater").pvalue if n else None
    ds = sum(v["both_ge_085"] for v in per.values())
    return dict(
        reading_A_pooled_splits=dict(
            k=int(k), n=int(n), proportion=round(k / n, 4) if n else None,
            ci95=wilson(k, n), p_one_sided_vs_0_85=float(p_bin) if p_bin is not None else None,
            ml_direction=dict(k=int(kml), n=int(nml),
                              proportion=round(kml / nml, 4) if nml else None,
                              ci95=wilson(kml, nml)),
            mb_direction=dict(k=int(kmb), n=int(nmb),
                              proportion=round(kmb / nmb, 4) if nmb else None,
                              ci95=wilson(kmb, nmb))),
        reading_B_datasets=dict(k=ds, n=len(per), proportion=round(ds / len(per), 4),
                                ci95=wilson(ds, len(per))),
        total_strong_split_conflicts=sum(v["strong_conflicts"] for v in per.values()),
        per_dataset=per)


# --------------------------------------------------------------------- C2
def clause2(topo):
    """Fine-scale branching order: disagreement in 30-50% of datasets?

    The original wording gives no threshold for "disagree", so every reading is
    reported: any difference at all, and normalised fine-scale RF above a series
    of cut-offs."""
    per = {}
    for d, r in topo.items():
        n = r["n_taxa"]; mx = 2 * (n - 3)
        mf = r["major_fine"]
        per[d] = dict(fine_scale_rf=mf["fine_scale_rf"],
                      fine_scale_rf_normalised=round(mf["fine_scale_rf"] / mx, 4),
                      n_weak_ml=mf["n_weak_ml"], n_weak_mb=mf["n_weak_mb"],
                      weak_ml_recovered_by_mb=mf["weak_ml_recovered_by_mb"],
                      weak_mb_recovered_by_ml=mf["weak_mb_recovered_by_ml"],
                      overall_nrf=r["between"]["nrf"])
    N = len(per)
    thresholds = {}
    for t in (0.0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50):
        k = sum(1 for v in per.values() if v["fine_scale_rf_normalised"] > t)
        thresholds[f">{t:.2f}"] = dict(k=k, n=N, proportion=round(k / N, 4),
                                       ci95=wilson(k, N),
                                       inside_30_50_window=(0.30 <= k / N <= 0.50))
    return dict(per_dataset=per, by_threshold=thresholds,
                n_with_any_fine_scale_difference=sum(
                    1 for v in per.values() if v["fine_scale_rf"] > 0))


# --------------------------------------------------------------------- C3
def clause3(epi):
    """Conclusion change: >=2 of 8 datasets -- and does it beat the null?"""
    rows, null_rows = {}, {}
    for d, r in epi.items():
        c = r["comparisons"]
        between = c["BETWEEN_ml_vs_bayes"]
        nullkeys = [k for k in c if k.startswith("NULL_ml_reseed")]
        n_null_any = sum(c[k]["ANY"] for k in nullkeys)
        n_null_t1 = sum(c[k]["tier1_who"] for k in nullkeys)
        mcmc = c.get("NULL_mcmc_replicate")
        prior = c.get("PRIOR_gammadir_vs_exponential")
        # Monte-Carlo p-value: how often does merely re-seeding ONE method
        # produce a change at least as strong as the between-method change?
        p_any = (n_null_any + 1) / (len(nullkeys) + 1) if between["ANY"] else 1.0
        p_t1 = (n_null_t1 + 1) / (len(nullkeys) + 1) if between["tier1_who"] else 1.0
        rows[d] = dict(
            primary_trait=r["primary_trait"],
            between_ANY=between["ANY"], between_tier1=between["tier1_who"],
            between_tier2=between["tier2_how_many"], between_tier3=between["tier3_scope"],
            index_ml=r["point_readouts"]["ml"].get("index"),
            index_bayes=r["point_readouts"]["mb"].get("index"),
            index_differs=between.get("index_group_differs"),
            reversed_edges=between.get("reversed_edges", []),
            intro_ml=r["point_readouts"]["ml"].get("intro_total"),
            intro_bayes=r["point_readouts"]["mb"].get("intro_total"),
            null_ml_reseed_ANY=f"{n_null_any}/{len(nullkeys)}",
            null_ml_reseed_tier1=f"{n_null_t1}/{len(nullkeys)}",
            null_mcmc_ANY=(mcmc or {}).get("ANY"),
            null_mcmc_tier1=(mcmc or {}).get("tier1_who"),
            prior_arm_ANY=(prior or {}).get("ANY"),
            prior_arm_tier1=(prior or {}).get("tier1_who"),
            mc_p_ANY=round(p_any, 4), mc_p_tier1=round(p_t1, 4))
        null_rows[d] = dict(rate_any=n_null_any / len(nullkeys),
                            rate_tier1=n_null_t1 / len(nullkeys))
    N = len(rows)
    k_any = sum(v["between_ANY"] for v in rows.values())
    k_t1 = sum(v["between_tier1"] for v in rows.values())
    k_null_any = sum(1 for v in null_rows.values() if v["rate_any"] > 0)
    k_null_t1 = sum(1 for v in null_rows.values() if v["rate_tier1"] > 0)
    k_mcmc = sum(1 for v in rows.values() if v["null_mcmc_ANY"])
    k_prior = sum(1 for v in rows.values() if v["prior_arm_ANY"])
    # McNemar on the paired indicator (between vs "re-seeding ML ever changes it")
    b = sum(1 for d in rows if rows[d]["between_ANY"] and null_rows[d]["rate_any"] == 0)
    c_ = sum(1 for d in rows if not rows[d]["between_ANY"] and null_rows[d]["rate_any"] > 0)
    mcn = (stats.binomtest(b, b + c_, 0.5).pvalue if (b + c_) > 0 else 1.0)
    return dict(
        n_datasets=N,
        between_ANY=dict(k=k_any, n=N, ci95=wilson(k_any, N),
                         meets_at_least_2_of_8=(k_any >= 2)),
        between_tier1_who=dict(k=k_t1, n=N, ci95=wilson(k_t1, N)),
        null_ml_reseed_datasets_ever_changing=dict(k=k_null_any, n=N, ci95=wilson(k_null_any, N)),
        null_ml_reseed_tier1_datasets=dict(k=k_null_t1, n=N),
        null_mcmc_replicate_datasets=dict(k=k_mcmc, n=N),
        prior_arm_datasets=dict(k=k_prior, n=N),
        mcnemar=dict(between_only=b, null_only=c_, p=float(mcn)),
        per_dataset=rows)


# ------------------------------------------------------- distance vs null
def distance_tests(topo, pis):
    ds = sorted(topo, key=lambda d: pis[d])
    b = np.array([topo[d]["between"]["nrf"] for d in ds])
    nml = np.array([topo[d]["null_ml_replicate"]["mean"] for d in ds])
    nmc = np.array([topo[d].get("null_mcmc_replicate", {}).get("nrf", np.nan) for d in ds])
    nps = np.array([topo[d]["null_posterior_spread"]["mean"] for d in ds])
    nbs = np.array([topo[d]["null_bootstrap_spread"]["mean"] for d in ds])
    out = dict(datasets=ds, between=b.round(4).tolist(),
               null_ml_reseed=nml.round(4).tolist(),
               null_mcmc_replicate=np.round(nmc, 4).tolist(),
               null_posterior_spread=nps.round(4).tolist(),
               null_bootstrap_spread=nbs.round(4).tolist())
    for name, null in (("vs_ml_reseed", nml), ("vs_posterior_spread", nps),
                       ("vs_bootstrap_spread", nbs)):
        w = stats.wilcoxon(b, null, alternative="greater")
        diff = b - null
        out[name] = dict(
            median_difference=round(float(np.median(diff)), 4),
            mean_ratio=round(float(np.mean(b / np.maximum(null, 1e-9))), 4),
            wilcoxon_stat=float(w.statistic), p_one_sided_greater=float(w.pvalue),
            n_datasets_between_exceeds_null=int((b > null).sum()))
    # per-dataset percentile of the between distance inside each null distribution
    out["percentile_in_null"] = {
        d: {k: topo[d]["calibration"].get(k) for k in
            ("pct_of_ml_replicate_null", "pct_of_posterior_spread",
             "pct_of_bootstrap_spread", "pct_of_posterior_pairs")} for d in ds}
    # does information content predict divergence?
    x = np.array([pis[d] for d in ds])
    out["information_gradient"] = dict(
        spearman_between=dict(rho=round(float(stats.spearmanr(x, b).statistic), 4),
                              p=float(stats.spearmanr(x, b).pvalue)),
        spearman_ratio=dict(
            rho=round(float(stats.spearmanr(x, b / np.maximum(nml, 1e-9)).statistic), 4),
            p=float(stats.spearmanr(x, b / np.maximum(nml, 1e-9)).pvalue)))
    return out


# ------------------------------------------------------ support calibration
def support_calibration(topo):
    out = {}
    pooled_u, pooled_p = [], []
    for d, r in topo.items():
        pts = r.get("split_support_pairs", [])
        if not pts:
            continue
        u = np.array([p[0] for p in pts]); p = np.array([p[1] for p in pts])
        pooled_u += list(u); pooled_p += list(p)
        w = stats.wilcoxon(p, u, alternative="greater") if len(u) > 5 else None
        out[d] = dict(n_splits=len(pts),
                      mean_ufboot=round(float(u.mean()), 4),
                      mean_pp=round(float(p.mean()), 4),
                      mean_pp_minus_ufboot=round(float((p - u).mean()), 4),
                      frac_pp_greater=round(float((p > u).mean()), 4),
                      spearman=round(float(stats.spearmanr(u, p).statistic), 4),
                      wilcoxon_p_pp_greater=(float(w.pvalue) if w else None),
                      resolution_ufboot=r["split_support"]["resolution_ufboot"],
                      resolution_pp=r["split_support"]["resolution_pp"])
    u = np.array(pooled_u); p = np.array(pooled_p)
    out["POOLED"] = dict(n_splits=len(u), mean_ufboot=round(float(u.mean()), 4),
                         mean_pp=round(float(p.mean()), 4),
                         mean_pp_minus_ufboot=round(float((p - u).mean()), 4),
                         frac_pp_greater=round(float((p > u).mean()), 4),
                         spearman=round(float(stats.spearmanr(u, p).statistic), 4))
    return out


# ---------------------------------------------------------- distributions
def distribution_tests(epi):
    out = {}
    for d, r in epi.items():
        dd = r.get("distributions", {})
        row = {}
        for key in ("intro_total", "n_transitions", "n_clusters_standard"):
            if key in dd and dd[key].get("ufboot"):
                row[key] = dict(
                    ufboot_median=dd[key]["ufboot"]["median"],
                    ufboot_95=[dd[key]["ufboot"]["lo95"], dd[key]["ufboot"]["hi95"]],
                    posterior_median=dd[key]["posterior"]["median"],
                    posterior_95=[dd[key]["posterior"]["lo95"], dd[key]["posterior"]["hi95"]],
                    overlap=dd[key].get("overlap"),
                    p=dd[key]["test"]["p"], effect=dd[key]["test"]["effect"])
        if "index_group" in dd:
            row["index_group"] = dict(
                ufboot_top=dd["index_group"]["ufboot_top"],
                ufboot_top_prob=dd["index_group"]["ufboot_top_prob"],
                posterior_top=dd["index_group"]["posterior_top"],
                posterior_top_prob=dd["index_group"]["posterior_top_prob"],
                total_variation=dd["index_group"]["tv"],
                overlap=dd["index_group"]["overlap"],
                point_estimates_agree=(dd["index_group"]["ufboot_top"] ==
                                       dd["index_group"]["posterior_top"]))
        if "edge_support" in dd:
            row["edge_support"] = {k: dd["edge_support"][k] for k in
                                   ("n_edges", "max_abs_diff", "n_edges_diff_gt_0_5",
                                    "n_confident_both", "n_confident_conflict")}
        out[d] = row
    return out


def convergence(topo):
    return {d: {t: r.get(f"mcmc_{t}") for t in ("mb", "mbexp")} for d, r in topo.items()}


# ------------------------------------------------------------------ tables
def markdown_tables(S, topo, epi, pis):
    ds = S["distance_tests"]["datasets"]
    L = []
    A = L.append
    A("# Synthesis tables\n")
    A("## T1 — Convergence gate (all downstream results are conditional on this)\n")
    A("| dataset | ngen | ASDSF | max PSRF | min ESS | converged | runtime (s) |")
    A("|---|---|---|---|---|---|---|")
    for d in ds:
        m = topo[d].get("mcmc_mb") or {}
        A(f'| `{d}` | {m.get("ngen")} | {m.get("asdsf_final")} | {m.get("max_psrf")} | '
          f'{m.get("min_ess")} | {"**yes**" if m.get("converged") else "NO"} | {m.get("seconds")} |')
    A("\n## T2 — Two scales of agreement (clause C1 and C2)\n")
    A("| dataset | PIS/taxon | strong ML splits | recovered by Bayes | strong Bayes splits | "
      "recovered by ML | strong conflicts | fine-scale nRF | overall nRF |")
    A("|---|---|---|---|---|---|---|---|---|")
    for d in ds:
        c1 = S["clause1"]["per_dataset"][d]; c2 = S["clause2"]["per_dataset"][d]
        A(f'| `{d}` | {pis[d]:.2f} | {c1["n_strong_ml"]} | {c1["ml_recovered_by_mb"]} | '
          f'{c1["n_strong_mb"]} | {c1["mb_recovered_by_ml"]} | {c1["strong_conflicts"]} | '
          f'{c2["fine_scale_rf_normalised"]} | {c2["overall_nrf"]} |')
    A("\n## T3 — Between-method distance vs the within-method nulls\n")
    A("| dataset | BETWEEN nRF | null: IQ-TREE re-seed | null: MrBayes chain1 vs 2 | "
      "null: UFBoot spread | null: posterior spread | between ÷ re-seed | percentile in posterior spread |")
    A("|---|---|---|---|---|---|---|---|")
    dt = S["distance_tests"]
    for i, d in enumerate(ds):
        pc = dt["percentile_in_null"][d]
        ratio = dt["between"][i] / max(dt["null_ml_reseed"][i], 1e-9)
        A(f'| `{d}` | **{dt["between"][i]}** | {dt["null_ml_reseed"][i]} | '
          f'{dt["null_mcmc_replicate"][i]} | {dt["null_bootstrap_spread"][i]} | '
          f'{dt["null_posterior_spread"][i]} | {ratio:.2f} | {pc["pct_of_posterior_spread"]} |')
    A("\n## T4 — Epidemiological conclusion change (clause C3), with the null\n")
    A("| dataset | trait | BETWEEN any | tier1 who | index ML → Bayes | introductions ML/Bayes | "
      "NULL re-seed (any) | NULL re-seed (tier1) | NULL MrBayes chains | prior arm | MC p (any) |")
    A("|---|---|---|---|---|---|---|---|---|---|---|")
    for d in ds:
        r = S["clause3"]["per_dataset"][d]
        A(f'| `{d}` | {r["primary_trait"]} | {"**YES**" if r["between_ANY"] else "no"} | '
          f'{"**YES**" if r["between_tier1"] else "no"} | {r["index_ml"]} → {r["index_bayes"]} | '
          f'{r["intro_ml"]}/{r["intro_bayes"]} | {r["null_ml_reseed_ANY"]} | '
          f'{r["null_ml_reseed_tier1"]} | {r["null_mcmc_ANY"]} | {r["prior_arm_ANY"]} | '
          f'{r["mc_p_ANY"]} |')
    A("\n## T5 — Support calibration (UFBoot vs posterior probability)\n")
    A("| dataset | n splits | mean UFBoot | mean PP | mean PP − UFBoot | % splits PP > UFBoot | "
      "Spearman ρ | resolution UFBoot≥.95 | resolution PP≥.95 |")
    A("|---|---|---|---|---|---|---|---|---|")
    for d in ds + ["POOLED"]:
        s = S["support_calibration"].get(d)
        if not s: continue
        A(f'| `{d}` | {s["n_splits"]} | {s["mean_ufboot"]} | {s["mean_pp"]} | '
          f'{s["mean_pp_minus_ufboot"]} | {s["frac_pp_greater"]:.1%} | {s["spearman"]} | '
          f'{s.get("resolution_ufboot","–")} | {s.get("resolution_pp","–")} |')
    A("\n## T6 — Uncertainty propagation: distributions, not point estimates\n")
    A("| dataset | introductions ML median [95%] | introductions Bayes median [95%] | overlap | p | "
      "index unit ML (prob) | index unit Bayes (prob) | TV distance |")
    A("|---|---|---|---|---|---|---|---|")
    for d in ds:
        r = S["distribution_tests"].get(d, {})
        it = r.get("intro_total"); ig = r.get("index_group")
        if not it and not ig: continue
        A(f'| `{d}` | '
          f'{it["ufboot_median"] if it else "–"} {it["ufboot_95"] if it else ""} | '
          f'{it["posterior_median"] if it else "–"} {it["posterior_95"] if it else ""} | '
          f'{it["overlap"] if it else "–"} | {f"{it['p']:.2g}" if it else "–"} | '
          f'{ig["ufboot_top"]+" ("+str(ig["ufboot_top_prob"])+")" if ig else "–"} | '
          f'{ig["posterior_top"]+" ("+str(ig["posterior_top_prob"])+")" if ig else "–"} | '
          f'{ig["total_variation"] if ig else "–"} |')
    open("results/synthesis_tables.md", "w").write("\n".join(L) + "\n")
    return "\n".join(L)


if __name__ == "__main__":
    topo, epi, pis = load()
    S = dict(n_datasets=len(topo),
             clause1=clause1(topo), clause2=clause2(topo), clause3=clause3(epi),
             distance_tests=distance_tests(topo, pis),
             support_calibration=support_calibration(topo),
             distribution_tests=distribution_tests(epi),
             convergence=convergence(topo))
    json.dump(S, open("results/synthesis.json", "w"), indent=1, default=str)
    markdown_tables(S, topo, epi, pis)
    print(json.dumps({k: v for k, v in S.items()
                      if k in ("clause1", "clause2", "clause3")}, indent=1,
                     default=str)[:4000])
    print("\nwrote results/synthesis.json, results/synthesis_tables.md")

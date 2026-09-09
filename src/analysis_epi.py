"""EXPERIMENT 3-6 -- do ML and Bayesian trees imply different public-health answers?

E3 (direction D2, pre-registered).  Four readouts a response team acts on, taken
from the ML point estimate and the Bayesian consensus:
      introductions   how many independent introductions per region/unit
      index_group     which unit sits at the root (the putative index case)
      transitions     the directed unit-to-unit who-infected-whom edge set
      clusters        genetic-distance transmission clusters at 3 thresholds
    A dataset is "conclusion-changing" under the rule pre-declared in
    planning.md: a change in the COUNT of introductions or of clusters, a change
    in the INDEX group, or a REVERSAL of a directed transmission edge.

E4 (ADDED).  The same rule applied to WITHIN-method pairs -- two IQ-TREE runs
    that differ only in random seed, and the two independent MrBayes chains.
    This is the null the ">=2 of 8 datasets" claim has to beat.

E5 (direction D3).  The readouts computed over whole tree DISTRIBUTIONS (1000
    UFBoot replicates vs the posterior sample) instead of point estimates, so the
    comparison is distribution-vs-distribution rather than tip-of-the-iceberg.

E6.  Branch-length-prior arm (compound Dirichlet vs exponential(10)) and a
    rooting-sensitivity arm (independent root-to-tip rooting vs a common fixed
    root), to attribute any difference to its actual cause.

Usage:  python src/analysis_epi.py [dataset ...]
Writes: results/epi/<dataset>.json, results/epi_summary.json
"""
import collections, itertools, json, math, os, random, statistics, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "code"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import dendropy
import treesets as T
import traitclean
import fastroot
from mlbayes import epi
from mlbayes.config import DATASETS

SEED = 42
POST_LIMIT = 1000       # posterior trees per dataset for the distribution arm
UFBOOT_LIMIT = 1000
CLUSTER_SUBSAMPLE = 300  # clustering is O(n^2) per tree -> subsample the sets
THRESHOLDS = {"tight": 1e-4, "standard": 5e-4, "loose": 1e-3}

# The grouping variable a public-health decision would actually turn on.
PRIMARY_TRAIT = {
    "fmdv_uk_2001":            "host",         # species jumps -> species-specific control
    "fmdv_uk_2007":            "epi_group",    # infected premises -> which farm to cull
    "ebov_sierraleone_2014":   "division",     # district -> where to send teams
    "ebov_westafrica_makona":  "geo_country",  # cross-border -> border control
    "zika_americas_2015_2016": "geo_country",  # country introductions
    "sarscov2_usa_early2020":  "division",     # state introductions -> travel measures
    "h5n1_dairy_cattle_2024":  "division",     # interstate cattle movement restrictions
    "mers_korea_2015":         None,           # single country: clusters only
}


# ------------------------------------------------------------------ readouts
def readouts(rooted, trait, with_clusters=True):
    """All epidemiological quantities for ONE rooted tree."""
    out = {}
    if trait:
        intro = epi.introductions(rooted, trait)
        trans = epi.group_transitions(rooted, trait)
        idx = epi.index_group(rooted, trait)
        out.update(intro_total=intro["total_introductions"],
                   intro_per_region=intro["introductions"],
                   parsimony_score=intro["parsimony_score"],
                   index=idx["root_state"],
                   index_ambiguous=idx["ambiguous"],
                   index_states=idx["root_states"],
                   n_transitions=trans["n_transitions"],
                   edge_set=trans["edge_set"])
    if with_clusters:
        for name, thr in THRESHOLDS.items():
            c = epi.transmission_clusters(rooted, max_patristic=thr)
            out[f"n_clusters_{name}"] = c["n_clusters"]
            out[f"n_clustered_{name}"] = c["n_clustered"]
            out[f"partition_{name}"] = sorted(tuple(x) for x in c["clusters"])
    return out


def _labels_from_partition(partition, all_taxa):
    """Cluster partition -> per-taxon label, singletons given unique labels, so
    two partitions can be compared with the adjusted Rand index."""
    lab, i = {}, 0
    for cl in partition:
        for t in cl:
            lab[t] = i
        i += 1
    for t in all_taxa:
        if t not in lab:
            lab[t] = i; i += 1
    return [lab[t] for t in sorted(all_taxa)]


def adjusted_rand(x, y):
    """Adjusted Rand index between two labellings of the same items."""
    from math import comb
    n = len(x)
    if n < 2:
        return 1.0
    ct = collections.Counter(zip(x, y))
    a = collections.Counter(x); b = collections.Counter(y)
    s_ij = sum(comb(v, 2) for v in ct.values())
    s_a = sum(comb(v, 2) for v in a.values())
    s_b = sum(comb(v, 2) for v in b.values())
    tot = comb(n, 2)
    exp = s_a * s_b / tot
    mx = 0.5 * (s_a + s_b)
    return round((s_ij - exp) / (mx - exp), 6) if mx != exp else 1.0


def conclusion_change(a, b, trait_present=True):
    """The rule pre-declared in planning.md, component by component.

    Binary flags alone saturate: with a 17-level trait over 120 taxa almost any
    two trees give a different introduction count, so `ANY` is True everywhere
    and carries no information.  Each flag is therefore paired with a CONTINUOUS
    effect size (`d_*`, `jaccard_edges`, `ari_*`), which is what the
    between-method value is actually compared against in its null.
    """
    f = {}
    if trait_present:
        f["introduction_count_differs"] = (a.get("intro_total") != b.get("intro_total")
                                           or a.get("intro_per_region") != b.get("intro_per_region"))
        f["index_group_differs"] = a.get("index") != b.get("index")
        ea, eb = set(a.get("edge_set", [])), set(b.get("edge_set", []))
        only_a, only_b = ea - eb, eb - ea
        rev = sorted(e for e in only_a
                     if "->".join(reversed(e.split("->"))) in only_b)
        f["reversed_edge"] = len(rev) > 0
        f["reversed_edges"] = rev
        f["edges_only_a"] = sorted(only_a)
        f["edges_only_b"] = sorted(only_b)
        # --- continuous effect sizes
        f["d_intro_total"] = abs((a.get("intro_total") or 0) - (b.get("intro_total") or 0))
        ra, rb = a.get("intro_per_region") or {}, b.get("intro_per_region") or {}
        f["d_intro_L1"] = sum(abs(ra.get(k, 0) - rb.get(k, 0))
                              for k in set(ra) | set(rb))
        f["n_regions_differing"] = sum(1 for k in set(ra) | set(rb)
                                       if ra.get(k, 0) != rb.get(k, 0))
        uni = ea | eb
        f["jaccard_edges"] = round(1 - len(ea & eb) / len(uni), 6) if uni else 0.0
        f["d_n_transitions"] = abs((a.get("n_transitions") or 0) - (b.get("n_transitions") or 0))
    f["cluster_count_differs"] = any(
        a.get(f"n_clusters_{n}") != b.get(f"n_clusters_{n}") for n in THRESHOLDS)
    f["cluster_partition_differs"] = any(
        a.get(f"partition_{n}") != b.get(f"partition_{n}") for n in THRESHOLDS)
    f["ANY"] = bool(f.get("introduction_count_differs") or f.get("index_group_differs")
                    or f.get("reversed_edge") or f["cluster_count_differs"])
    # severity tiers: what kind of action would actually change
    f["tier1_who"] = bool(f.get("index_group_differs") or f.get("reversed_edge"))
    f["tier2_how_many"] = bool(f.get("introduction_count_differs"))
    f["tier3_scope"] = bool(f["cluster_count_differs"] or f["cluster_partition_differs"])
    for nm in THRESHOLDS:
        ka, kb = a.get(f"n_clusters_{nm}"), b.get(f"n_clusters_{nm}")
        if ka is not None and kb is not None:
            f[f"d_n_clusters_{nm}"] = abs(ka - kb)
        pa, pb = a.get(f"partition_{nm}"), b.get(f"partition_{nm}")
        if pa is not None and pb is not None:
            taxa = {t for cl in pa for t in cl} | {t for cl in pb for t in cl}
            if taxa:
                f[f"ari_{nm}"] = adjusted_rand(_labels_from_partition(pa, taxa),
                                               _labels_from_partition(pb, taxa))
    return f


# ------------------------------------------------- distribution comparison
def _dist(vals):
    if not vals: return None
    c = collections.Counter(vals)
    v = sorted(vals); n = len(v)
    return dict(n=n, mean=round(statistics.fmean(v), 4), median=v[n // 2],
                sd=round(statistics.pstdev(v), 4) if n > 1 else 0.0,
                lo95=v[int(0.025 * n)], hi95=v[min(int(0.975 * n), n - 1)],
                mode=c.most_common(1)[0][0],
                pmf={str(k): round(cnt / n, 5) for k, cnt in sorted(c.items())})


def tv_distance(pmf_a, pmf_b):
    """Total-variation distance between two categorical distributions."""
    keys = set(pmf_a) | set(pmf_b)
    return round(0.5 * sum(abs(pmf_a.get(k, 0.0) - pmf_b.get(k, 0.0)) for k in keys), 6)


def overlap(pmf_a, pmf_b):
    keys = set(pmf_a) | set(pmf_b)
    return round(sum(min(pmf_a.get(k, 0.0), pmf_b.get(k, 0.0)) for k in keys), 6)


def mann_whitney(x, y):
    """Two-sided Mann-Whitney U with normal approximation + tie correction, and
    the rank-biserial effect size.  Used for integer-valued readouts
    (introduction counts, cluster counts) where normality cannot be assumed."""
    from scipy import stats
    if not x or not y or (len(set(x)) == 1 and set(x) == set(y)):
        return dict(U=None, p=1.0, effect=0.0, note="degenerate")
    try:
        u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
    except ValueError:
        return dict(U=None, p=1.0, effect=0.0, note="degenerate")
    n1, n2 = len(x), len(y)
    return dict(U=float(u), p=float(p),
                effect=round(2 * float(u) / (n1 * n2) - 1, 4),   # rank-biserial
                median_diff=float(statistics.median(x) - statistics.median(y)))


def edge_support(readout_list):
    """Bootstrap proportion / posterior probability of each directed edge."""
    n = len(readout_list)
    c = collections.Counter()
    for r in readout_list:
        c.update(set(r.get("edge_set", [])))
    return {e: v / n for e, v in c.items()}, n


# ------------------------------------------------------------------ driver
def analyse(dsid, outdir="results/epi"):
    t0 = time.time()
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(SEED)
    tns = T.new_namespace()
    dts = T.dates(dsid)
    all_traits = T.traits(dsid)
    ptrait_name = PRIMARY_TRAIT.get(dsid)
    ptrait, clean_report = None, None
    if ptrait_name and ptrait_name in all_traits:
        # normalise synonyms and mixed granularity BEFORE any inference readout,
        # otherwise annotation style is scored as epidemiology (see traitclean.py)
        ptrait, clean_report = traitclean.normalise(dsid, ptrait_name,
                                                   all_traits[ptrait_name])
    res = dict(dataset=dsid, primary_trait=ptrait_name,
               trait_cleaning=clean_report,
               n_trait_levels=(len(set(ptrait.values())) if ptrait else 0),
               traits_available=sorted(all_traits))

    # ---------------- point estimates ------------------------------------
    kinds = ["ml", "mlrep2", "mlrep3", "mlrep4", "mlrep5", "mb"]
    if os.path.exists(f"results/{dsid}/{dsid}_mbexp.con.tre"):
        kinds.append("mbexp")
    if os.path.exists(f"results/{dsid}/{dsid}_mlgtr.treefile"):
        kinds.append("mlgtr")   # matched-model ML arm: paradigm without model
    trees = {k: T.point_tree(dsid, k, tns) for k in kinds}

    # The two MrBayes chains, summarised separately -> the MCMC-replicate null.
    # Built with greedy (allcompat) consensus and MEAN split branch lengths, so
    # they are constructed exactly like MrBayes' own .con.tre; a 50%-majority
    # tree would be only partly resolved and would not be comparable, and the
    # epidemiological readouts need branch lengths anyway (rooting, clustering).
    # full post-burn-in sample for the consensus trees (sub-sampling before
    # building a consensus was measured to move it as much as the effect under
    # test); a fixed-size sub-sample is used for the distribution readouts below
    post_full, post_runs_full = T.posterior(dsid, tns, tag="mb", limit=None)
    for r, tl in post_runs_full.items():
        trees[f"mbrun{r}"] = T.greedy_consensus(tl, tns)
    rng_p = random.Random(SEED)
    post_all = (rng_p.sample(post_full, POST_LIMIT)
                if len(post_full) > POST_LIMIT else list(post_full))

    rooted, rttinfo = {}, {}
    for k, t in trees.items():
        rt = fastroot.rtt_root(t, dts)
        rooted[k] = rt
        rttinfo[k] = dict(r2=(round(rt.rtt_r2, 4) if rt.rtt_r2 is not None else None),
                          rate=(f"{rt.rtt_rate:.3e}" if rt.rtt_rate else None),
                          n_dated=rt.rtt_n)
    res["rooting"] = rttinfo

    # ---- external anchor: the earliest-sampled group ---------------------
    # There is no case-level ground truth for these outbreaks, but the index
    # unit must be detectable no later than any unit it seeded.  The
    # earliest-sampled group is therefore a weak but INDEPENDENT check on the
    # inferred index unit: an index call on a group first sampled well after
    # another group requires unsampled intermediates to be plausible.
    # (Note: `papers/cottam2008_fmdv_transmission_pathways.pdf`, catalogued in
    # papers/README.md as ground truth for fmdv_uk_2007, is in fact an analysis
    # of the 2001 UK outbreak, so it is NOT used as ground truth here.)
    if ptrait:
        first = {}
        for lab, g in ptrait.items():
            if lab in dts:
                first[g] = min(first.get(g, 1e9), dts[lab])
        order_g = sorted(first, key=first.get)
        res["earliest_group"] = dict(
            group=order_g[0] if order_g else None,
            first_sampled={g: round(first[g], 4) for g in order_g},
            gap_to_second=(round(first[order_g[1]] - first[order_g[0]], 4)
                           if len(order_g) > 1 else None))
    # the root clade of the ML tree, used to impose a COMMON rooting later
    ml_root_clade = sorted(lf.taxon.label for lf in
                           min(rooted["ml"].seed_node.child_nodes(),
                               key=lambda c: len(list(c.leaf_iter()))).leaf_iter())
    res["ml_root_clade_size"] = len(ml_root_clade)

    pr = {k: readouts(v, ptrait) for k, v in rooted.items()}
    if ptrait and res.get("earliest_group", {}).get("group"):
        eg = res["earliest_group"]["group"]
        res["index_matches_earliest_group"] = {
            k: (v.get("index") == eg) for k, v in pr.items()}
    res["point_readouts"] = {k: {kk: vv for kk, vv in v.items()
                                 if not kk.startswith("partition_")} for k, v in pr.items()}

    # ---------------- E3 + E4: the rule, between and within ---------------
    tp = ptrait is not None
    comps = {"BETWEEN_ml_vs_bayes": ("ml", "mb")}
    if "mbexp" in pr:
        comps["PRIOR_gammadir_vs_exponential"] = ("mb", "mbexp")
        comps["BETWEEN_ml_vs_bayes_expprior"] = ("ml", "mbexp")
    if "mbrun1" in pr and "mbrun2" in pr:
        comps["NULL_mcmc_replicate"] = ("mbrun1", "mbrun2")
    if "mlgtr" in pr:
        # paradigm alone (same substitution model on both sides), and the pure
        # substitution-model effect inside ML
        comps["BETWEEN_mlGTR_vs_bayes_PARADIGM_ONLY"] = ("mlgtr", "mb")
        comps["MODEL_ml_modelfinder_vs_mlGTR"] = ("ml", "mlgtr")
    mlpairs = list(itertools.combinations(["ml", "mlrep2", "mlrep3", "mlrep4", "mlrep5"], 2))
    for a, b in mlpairs:
        comps[f"NULL_ml_reseed_{a}_vs_{b}"] = (a, b)

    res["comparisons"] = {name: conclusion_change(pr[a], pr[b], tp)
                          for name, (a, b) in comps.items()}
    # headline summary: between-method vs the within-ML null rate
    nullkeys = [k for k in res["comparisons"] if k.startswith("NULL_ml_reseed")]
    res["null_ml_reseed_change_rate"] = round(
        sum(res["comparisons"][k]["ANY"] for k in nullkeys) / len(nullkeys), 4)
    res["null_ml_reseed_tier1_rate"] = round(
        sum(res["comparisons"][k]["tier1_who"] for k in nullkeys) / len(nullkeys), 4)

    # ---------------- E6: rooting sensitivity -----------------------------
    # impose the ML tree's root clade on the Bayesian consensus (and vice versa)
    fixed = {}
    for k in ("ml", "mb"):
        rt, ok = fastroot.root_at_leafset(trees[k], ml_root_clade, dts)
        fixed[k] = (rt, ok)
    res["rooting_sensitivity"] = dict(
        common_root_applied={k: v[1] for k, v in fixed.items()},
        under_common_root=conclusion_change(
            readouts(fixed["ml"][0], ptrait), readouts(fixed["mb"][0], ptrait), tp),
    )
    mid = {k: epi.root_midpoint(trees[k]) for k in ("ml", "mb")}
    res["rooting_sensitivity"]["under_midpoint_root"] = conclusion_change(
        readouts(mid["ml"], ptrait), readouts(mid["mb"], ptrait), tp)

    # ---------------- E5: distributions ------------------------------------
    ub = T.ufboot(dsid, tns, limit=UFBOOT_LIMIT, seed=SEED)
    sets = {"ufboot": ub, "posterior": post_all}
    dist_read, dist_read_clu = {}, {}
    for name, tl in sets.items():
        rt = [fastroot.rtt_root(t, dts) for t in tl]
        dist_read[name] = [readouts(t, ptrait, with_clusters=False) for t in rt]
        sub = rt if len(rt) <= CLUSTER_SUBSAMPLE else rng.sample(rt, CLUSTER_SUBSAMPLE)
        dist_read_clu[name] = [readouts(t, None, with_clusters=True) for t in sub]

    dsum = {}
    if ptrait:
        for key in ("intro_total", "n_transitions", "parsimony_score"):
            a = [r[key] for r in dist_read["ufboot"]]
            b = [r[key] for r in dist_read["posterior"]]
            dsum[key] = dict(ufboot=_dist(a), posterior=_dist(b),
                             test=mann_whitney(a, b),
                             overlap=overlap(_dist(a)["pmf"], _dist(b)["pmf"]),
                             tv=tv_distance(_dist(a)["pmf"], _dist(b)["pmf"]))
        ia = [r["index"] for r in dist_read["ufboot"]]
        ib = [r["index"] for r in dist_read["posterior"]]
        pa = {k: v / len(ia) for k, v in collections.Counter(ia).items()}
        pb = {k: v / len(ib) for k, v in collections.Counter(ib).items()}
        eg = (res.get("earliest_group") or {}).get("group")
        dsum["index_group"] = dict(
            earliest_sampled_group=eg,
            ufboot_prob_of_earliest=round(pa.get(eg, 0.0), 4) if eg else None,
            posterior_prob_of_earliest=round(pb.get(eg, 0.0), 4) if eg else None,
            ufboot_pmf={k: round(v, 4) for k, v in sorted(pa.items(), key=lambda x: -x[1])},
            posterior_pmf={k: round(v, 4) for k, v in sorted(pb.items(), key=lambda x: -x[1])},
            tv=tv_distance(pa, pb), overlap=overlap(pa, pb),
            ufboot_top=max(pa, key=pa.get), posterior_top=max(pb, key=pb.get),
            ufboot_top_prob=round(max(pa.values()), 4),
            posterior_top_prob=round(max(pb.values()), 4))
        ea, na = edge_support(dist_read["ufboot"])
        eb, nb = edge_support(dist_read["posterior"])
        keys = sorted(set(ea) | set(eb))
        dsum["edge_support"] = dict(
            n_edges=len(keys), n_ufboot=na, n_posterior=nb,
            table={k: [round(ea.get(k, 0.0), 4), round(eb.get(k, 0.0), 4)] for k in keys},
            max_abs_diff=(round(max(abs(ea.get(k, 0.0) - eb.get(k, 0.0)) for k in keys), 4)
                          if keys else None),
            n_edges_diff_gt_0_5=sum(1 for k in keys
                                    if abs(ea.get(k, 0.0) - eb.get(k, 0.0)) > 0.5),
            n_confident_both=sum(1 for k in keys
                                 if ea.get(k, 0) >= 0.95 and eb.get(k, 0) >= 0.95),
            n_confident_conflict=sum(1 for k in keys
                                     if (ea.get(k, 0) >= 0.95) != (eb.get(k, 0) >= 0.95)))
    for name in THRESHOLDS:
        a = [r[f"n_clusters_{name}"] for r in dist_read_clu["ufboot"]]
        b = [r[f"n_clusters_{name}"] for r in dist_read_clu["posterior"]]
        dsum[f"n_clusters_{name}"] = dict(ufboot=_dist(a), posterior=_dist(b),
                                          test=mann_whitney(a, b),
                                          overlap=overlap(_dist(a)["pmf"], _dist(b)["pmf"]))
    res["distributions"] = dsum

    res["analysis_seconds"] = round(time.time() - t0, 1)
    json.dump(res, open(os.path.join(outdir, f"{dsid}.json"), "w"), indent=1)
    return res


if __name__ == "__main__":
    ids = sys.argv[1:] or DATASETS
    summary = []
    for d in ids:
        print(f"--- {d}", flush=True)
        try:
            r = analyse(d)
            c = r["comparisons"]["BETWEEN_ml_vs_bayes"]
            print(f'    BETWEEN change={c["ANY"]} tier1={c["tier1_who"]} '
                  f'intro={c["introduction_count_differs"] if r["primary_trait"] else "n/a"} '
                  f'| NULL ml-reseed rate={r["null_ml_reseed_change_rate"]} '
                  f'tier1={r["null_ml_reseed_tier1_rate"]} [{r["analysis_seconds"]}s]', flush=True)
        except Exception as e:
            import traceback; traceback.print_exc()
            r = dict(dataset=d, error=str(e))
        summary.append(r)
        json.dump(summary, open("results/epi_summary.json", "w"), indent=1)
    print("wrote results/epi_summary.json")

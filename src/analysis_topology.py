"""EXPERIMENT 1 + 2 -- topological concordance, with a within-method null.

E1 (direction D1, pre-registered).  For each dataset, quantify agreement at the
two scales the hypothesis distinguishes:
    MAJOR   fraction of confidently supported splits (UFBoot >= 95 / PP >= 0.95)
            of one method that are present in the other method's point estimate
    FINE    normalised Robinson-Foulds over the remaining, weakly supported
            splits; quartet distance; the Bayesian topological credible set

E2 (ADDED after the pilot -- see STATE.md).  A between-method difference is only
evidence about MODEL CHOICE if it exceeds what the SAME method produces when you
merely re-run it.  The pilot showed the posterior never samples the same topology
twice (max topology probability = 1/n_samples), so a raw ML-vs-Bayesian RF has no
interpretable scale.  E2 supplies that scale with four null distributions:

    N1  MCMC replicate    consensus(MrBayes run 1) vs consensus(run 2)
    N2  ML search replicate  all pairs among 5 IQ-TREE runs differing only in seed
    N3  posterior spread   Bayesian consensus vs each posterior sample tree
    N4  bootstrap spread   ML tree vs each UFBoot replicate
    P   prior sensitivity  compound-Dirichlet vs exponential(10) branch-length prior

Usage:  python src/analysis_topology.py [dataset ...]
Writes: results/topology/<dataset>.json  and results/topology_summary.json
"""
import collections, itertools, json, math, os, random, statistics, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "code"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import dendropy
import treesets as T
from mlbayes import compare
from mlbayes.config import DATASETS

SEED = 42
POST_LIMIT = 2000      # sub-sample size for per-tree SPREAD statistics only
N_CONSENSUS_BOOT = 20  # resampled-consensus pairs for the point-estimate null
UFBOOT_LIMIT = 1000    # all of them
N_QUARTETS = 20000
STRONG_UFB, STRONG_PP = 95.0, 0.95
ML_REPS = ["ml", "mlrep2", "mlrep3", "mlrep4", "mlrep5"]


def _stats(vals):
    if not vals:
        return None
    v = sorted(vals)
    n = len(v)
    return dict(n=n, mean=round(statistics.fmean(v), 6),
                sd=round(statistics.pstdev(v), 6) if n > 1 else 0.0,
                min=round(v[0], 6), q025=round(v[int(0.025 * n)], 6),
                median=round(v[n // 2], 6),
                q975=round(v[min(int(0.975 * n), n - 1)], 6), max=round(v[-1], 6))


def _percentile_of(x, vals):
    """Fraction of the null distribution that is <= x."""
    if not vals:
        return None
    return round(sum(1 for v in vals if v <= x) / len(vals), 6)


def analyse(dsid, outdir="results/topology"):
    t0 = time.time()
    os.makedirs(outdir, exist_ok=True)
    tns = T.new_namespace()
    res = dict(dataset=dsid)

    # ---------------- point estimates -----------------------------------
    ml = T.point_tree(dsid, "ml", tns)
    mb = T.point_tree(dsid, "mb", tns)
    reps = {k: T.point_tree(dsid, k, tns) for k in ML_REPS[1:]}
    mbexp = None
    if os.path.exists(f"results/{dsid}/{dsid}_mbexp.con.tre"):
        mbexp = T.point_tree(dsid, "mbexp", tns)
    # matched-model ML arm (GTR+F+I+G4 == the MrBayes block), so the
    # substitution model can be separated from the inference paradigm
    mlgtr = None
    if os.path.exists(f"results/{dsid}/{dsid}_mlgtr.treefile"):
        mlgtr = T.point_tree(dsid, "mlgtr", tns)

    n = len(tns)
    res["n_taxa"] = n
    maxrf = 2 * (n - 3)
    res["max_rf"] = maxrf

    s_ml, s_mb = T.split_set(ml), T.split_set(mb)
    s_reps = {k: T.split_set(v) for k, v in reps.items()}
    s_exp = T.split_set(mbexp) if mbexp is not None else None
    s_gtr = T.split_set(mlgtr) if mlgtr is not None else None

    # ---------------- E1: two-scale agreement ---------------------------
    # reuse the reference support parsers (IQ-TREE 'SH/UFB' labels, MrBayes PP)
    ml2 = T.point_tree(dsid, "ml", tns); ml2.encode_bipartitions(suppress_unifurcations=False)
    mb2 = T.point_tree(dsid, "mb", tns); mb2.encode_bipartitions(suppress_unifurcations=False)
    agree = compare.split_agreement(ml2, mb2, STRONG_UFB, STRONG_PP)
    paired = agree.pop("paired_support")
    res["major_fine"] = agree
    res["paired_support"] = paired
    res["between"] = dict(rf=T.rf(s_ml, s_mb),
                          nrf=round(T.rf(s_ml, s_mb) / maxrf, 6))

    try:
        res["quartet"] = compare.quartet_distance(ml2, mb2, tns, max_quartets=N_QUARTETS)
    except Exception as e:
        res["quartet"] = {"error": str(e)}

    # ---------------- N2: ML search replicates ---------------------------
    allml = {"ml": s_ml, **s_reps}
    pair_nrf = [T.rf(allml[a], allml[b]) / maxrf
                for a, b in itertools.combinations(ML_REPS, 2)]
    res["null_ml_replicate"] = _stats(pair_nrf)
    res["null_ml_replicate"]["pairs"] = {f"{a}|{b}": round(T.rf(allml[a], allml[b]) / maxrf, 6)
                                         for a, b in itertools.combinations(ML_REPS, 2)}
    # how many of the MAIN run's strongly-supported splits survive a re-seed?
    sup = compare.ml_support(ml2)
    strong_ml_splits = {b.split_bitmask for b, (sh, u) in sup.items()
                        if u is not None and u >= STRONG_UFB and b is not None
                        and not b.is_trivial()}
    # ml2 and ml are separate Tree objects; their bitmasks are only comparable if
    # both were encoded unrooted on the same namespace.  Fail loudly if not.
    assert strong_ml_splits <= s_ml, (
        "support-labelled splits are not a subset of the ML split set -- "
        "bipartition encodings are inconsistent")
    res["n_strong_ufb_splits"] = len(strong_ml_splits)
    res["strong_ml_recovered_by_ml_reseed"] = _stats(
        [len(strong_ml_splits & s) / len(strong_ml_splits) for s in s_reps.values()]
    ) if strong_ml_splits else None

    # ---------------- distributions --------------------------------------
    random.seed(SEED)
    ub = T.ufboot(dsid, tns, limit=UFBOOT_LIMIT, seed=SEED)
    s_ub = T.split_sets(ub)
    # Load the FULL post-burn-in posterior.  Sub-sampling it before building a
    # consensus was measured to move the consensus by nRF 0.43 on fmdv_uk_2007 --
    # i.e. sub-sampling noise alone was as large as the effect under test.  With
    # the full sample our greedy consensus reproduces MrBayes' own .con.tre
    # exactly (nRF = 0); that identity is recorded below as a validation field.
    post_full, post_runs_full = T.posterior(dsid, tns, tag="mb", limit=None)
    s_post_full = T.split_sets(post_full)
    s_post_runs = {r: T.split_sets(tr) for r, tr in post_runs_full.items()}
    # a fixed-size sub-sample is still used for the SPREAD statistics, which are
    # per-tree distances and do not need every tree
    rng0 = random.Random(SEED)
    s_post = (rng0.sample(s_post_full, POST_LIMIT)
              if len(s_post_full) > POST_LIMIT else list(s_post_full))

    # N4 bootstrap spread, N3 posterior spread
    res["null_bootstrap_spread"] = _stats([T.rf(s_ml, s) / maxrf for s in s_ub])
    res["null_posterior_spread"] = _stats([T.rf(s_mb, s) / maxrf for s in s_post])
    # cross terms: is the ML tree an atypical member of the posterior cloud?
    ml_to_post = [T.rf(s_ml, s) / maxrf for s in s_post]
    mb_to_ub = [T.rf(s_mb, s) / maxrf for s in s_ub]
    res["ml_to_posterior"] = _stats(ml_to_post)
    res["bayes_to_ufboot"] = _stats(mb_to_ub)
    # within-posterior pair distance = the Bayesian method's intrinsic scale
    rng = random.Random(SEED)
    pp = [T.rf(*rng.sample(s_post, 2)) / maxrf for _ in range(2000)] if len(s_post) > 2 else []
    res["null_posterior_pairs"] = _stats(pp)
    uu = [T.rf(*rng.sample(s_ub, 2)) / maxrf for _ in range(2000)] if len(s_ub) > 2 else []
    res["null_ufboot_pairs"] = _stats(uu)

    # N1a: summarise each MrBayes chain separately.
    # The consensus MUST be built the way MrBayes builds its own (contype=allcompat,
    # i.e. greedy / extended majority rule, fully resolved).  A 50%-majority-rule
    # tree is only partly resolved, and RF between two trees of different
    # resolution is not commensurable with RF between two resolved trees.
    if len(s_post_runs) == 2:
        cons = {r: T.split_set(T.greedy_consensus(ss, tns))
                for r, ss in s_post_runs.items()}
        r1, r2 = sorted(cons)
        s_pooled = T.split_set(T.greedy_consensus(s_post_full, tns))
        res["null_mcmc_replicate"] = dict(
            rf=T.rf(cons[r1], cons[r2]),
            nrf=round(T.rf(cons[r1], cons[r2]) / maxrf, 6),
            n_trees_run1=len(s_post_runs[r1]), n_trees_run2=len(s_post_runs[r2]),
            nrf_run1_vs_pooled=round(T.rf(cons[r1], s_pooled) / maxrf, 6),
            nrf_run2_vs_pooled=round(T.rf(cons[r2], s_pooled) / maxrf, 6),
            # VALIDATION: our pooled greedy consensus should reproduce MrBayes'
            # own contype=allcompat .con.tre exactly, i.e. nRF = 0
            nrf_pooled_vs_mrbayes_contre=round(T.rf(s_pooled, s_mb) / maxrf, 6),
            n_splits_run1=len(cons[r1]), n_splits_run2=len(cons[r2]),
            n_splits_pooled=len(s_pooled), n_splits_mrbayes=len(s_mb),
            n_internal_max=n - 3)

    # N1b: Monte-Carlo instability of the Bayesian POINT ESTIMATE at the actual
    # sample size.  Resample the pooled posterior WITH REPLACEMENT to the same
    # size and rebuild the consensus; the distance between two such consensus
    # trees is how much the published tree would move if the chain had merely
    # drawn a different sample.  This is the fairest null for the between-method
    # distance, because both sides then rest on the same number of trees.
    rngb = random.Random(SEED + 7)
    m = len(s_post_full)
    boot = []
    for _ in range(N_CONSENSUS_BOOT):
        a = T.split_set(T.greedy_consensus(
            [s_post_full[rngb.randrange(m)] for _ in range(m)], tns))
        b = T.split_set(T.greedy_consensus(
            [s_post_full[rngb.randrange(m)] for _ in range(m)], tns))
        boot.append(T.rf(a, b) / maxrf)
    res["null_consensus_resample"] = _stats(boot)

    # the same instability on the ML side: resample the UFBoot replicates and
    # rebuild the bootstrap consensus
    rngc = random.Random(SEED + 8)
    mu = len(s_ub)
    bootml = []
    for _ in range(N_CONSENSUS_BOOT):
        a = T.split_set(T.greedy_consensus(
            [s_ub[rngc.randrange(mu)] for _ in range(mu)], tns))
        b = T.split_set(T.greedy_consensus(
            [s_ub[rngc.randrange(mu)] for _ in range(mu)], tns))
        bootml.append(T.rf(a, b) / maxrf)
    res["null_ufboot_consensus_resample"] = _stats(bootml)

    # ---------------- model-vs-paradigm decomposition ---------------------
    if s_gtr is not None:
        best = None
        lg = f"results/{dsid}/{dsid}_ml.log"
        if os.path.exists(lg):
            for line in open(lg, errors="ignore"):
                if line.startswith("Best-fit model"):
                    best = line.split(":", 1)[1].split(" chosen")[0].strip(); break
        res["model_control"] = dict(
            modelfinder_model=best, matched_model="GTR+F+I+G4",
            nrf_paradigm_plus_model=round(T.rf(s_ml, s_mb) / maxrf, 6),
            nrf_paradigm_only=round(T.rf(s_gtr, s_mb) / maxrf, 6),
            nrf_model_only_within_ml=round(T.rf(s_ml, s_gtr) / maxrf, 6))

    # ---------------- P: branch-length prior arm --------------------------
    if s_exp is not None:
        res["prior_arm"] = dict(
            nrf_gammadir_vs_exponential=round(T.rf(s_mb, s_exp) / maxrf, 6),
            nrf_ml_vs_exponential=round(T.rf(s_ml, s_exp) / maxrf, 6),
            nrf_ml_vs_gammadir=round(T.rf(s_ml, s_mb) / maxrf, 6))

    # ---------------- where does the between-method distance sit? ---------
    b = res["between"]["nrf"]
    res["calibration"] = dict(
        between_nrf=b,
        pct_of_ml_replicate_null=_percentile_of(b, pair_nrf),
        pct_of_posterior_spread=_percentile_of(b, [T.rf(s_mb, s) / maxrf for s in s_post]),
        pct_of_bootstrap_spread=_percentile_of(b, [T.rf(s_ml, s) / maxrf for s in s_ub]),
        pct_of_posterior_pairs=_percentile_of(b, pp),
        pct_of_consensus_resample=_percentile_of(b, boot),
        ratio_between_over_ml_replicate=(round(b / res["null_ml_replicate"]["mean"], 4)
                                         if res["null_ml_replicate"]["mean"] else None),
        ratio_between_over_posterior_spread=(
            round(b / res["null_posterior_spread"]["mean"], 4)
            if res["null_posterior_spread"] and res["null_posterior_spread"]["mean"] else None),
    )

    # ---------------- split-level support comparison ----------------------
    fb, nb = T.split_frequencies(s_ub)
    fp, npp = T.split_frequencies(s_post_full)
    allsplits = set(fb) | set(fp)
    pairs = [(fb.get(s, 0.0), fp.get(s, 0.0)) for s in allsplits]
    res["split_support"] = dict(
        n_splits_seen=len(allsplits),
        n_ufboot_ge95=sum(1 for s in allsplits if fb.get(s, 0) >= 0.95),
        n_pp_ge95=sum(1 for s in allsplits if fp.get(s, 0) >= 0.95),
        n_both_ge95=sum(1 for s in allsplits if fb.get(s, 0) >= 0.95 and fp.get(s, 0) >= 0.95),
        resolution_ufboot=round(sum(1 for s in allsplits if fb.get(s, 0) >= 0.95) / (n - 3), 4),
        resolution_pp=round(sum(1 for s in allsplits if fp.get(s, 0) >= 0.95) / (n - 3), 4),
        mean_pp_minus_ufboot=round(statistics.fmean([p - u for u, p in pairs]), 6),
        n_bootstrap_trees=nb, n_posterior_trees=npp,
    )
    res["split_support_pairs"] = [[round(u, 4), round(p, 4)] for u, p in pairs]

    # ---------------- credible set (reference implementation) -------------
    tp = f"results/{dsid}/{dsid}_mb.trprobs"
    if os.path.exists(tp):
        try:
            res["credible_set"] = compare.credible_set(
                tp, f"results/{dsid}/{dsid}_ml.treefile")
        except Exception as e:
            res["credible_set"] = {"error": str(e)}
    # sample-size-free alternative: how concentrated is the posterior really?
    res["credible_set"] = res.get("credible_set", {})
    res["credible_set"]["distinct_topologies_in_sample"] = len(set(s_post_full))
    res["credible_set"]["n_sampled"] = len(s_post_full)
    res["credible_set"]["max_topology_freq_in_sample"] = (
        round(max(collections.Counter(s_post_full).values()) / len(s_post_full), 6)
        if s_post_full else None)

    # ---------------- MCMC convergence gate -------------------------------
    for tag in ("mb", "mbexp"):
        f = f"results/{dsid}/{dsid}_{tag}.mbstats.json"
        if os.path.exists(f):
            m = json.load(open(f))
            res[f"mcmc_{tag}"] = {k: m.get(k) for k in
                                  ("asdsf_final", "max_psrf", "min_ess", "ngen", "seconds")}
            res[f"mcmc_{tag}"]["converged"] = bool(
                m.get("asdsf_final") is not None and m["asdsf_final"] < 0.01 and
                (m.get("max_psrf") is None or m["max_psrf"] < 1.02) and
                (m.get("min_ess") is None or m["min_ess"] > 200))

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
            print(f"    between nRF={r['between']['nrf']:.3f}  "
                  f"ML-reseed null mean={r['null_ml_replicate']['mean']:.3f}  "
                  f"MCMC-replicate null={r.get('null_mcmc_replicate',{}).get('nrf')}  "
                  f"post spread={r['null_posterior_spread']['mean']:.3f}  "
                  f"[{r['analysis_seconds']}s]", flush=True)
        except Exception as e:
            import traceback; traceback.print_exc()
            r = dict(dataset=d, error=str(e))
        summary.append(r)
        json.dump(summary, open("results/topology_summary.json", "w"), indent=1)
    print("wrote results/topology_summary.json")

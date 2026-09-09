"""Loading tree point-estimates and tree DISTRIBUTIONS, and fast split-set
distances between them.

Every tree for a given dataset is loaded onto ONE shared TaxonNamespace, which
is what makes split bitmasks comparable across files (IQ-TREE newick, MrBayes
NEXUS consensus, MrBayes .t posterior samples, IQ-TREE .ufboot replicates).

Distances are computed on normalised split-bitmask SETS rather than through
dendropy's tree-comparison API: the analysis needs O(10^5)-O(10^6) pairwise
distances (posterior spreads, bootstrap spreads, all-pairs replicate matrices)
and the set-based form is ~2 orders of magnitude faster.  `test_treesets.py`
checks it against dendropy.treecompare on real trees.
"""
import collections, csv, gzip, os, random, re, statistics
import dendropy

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# --------------------------------------------------------------- metadata
def metadata(dsid):
    """taxon_map.tsv keyed by GenBank accession (= the tip label in the trees)."""
    path = os.path.join(ROOT, "datasets", dsid, "taxon_map.tsv")
    return {r["id"]: r for r in csv.DictReader(open(path), delimiter="\t")}


def dates(dsid):
    m = metadata(dsid)
    return {k: float(v["decimal_date"]) for k, v in m.items() if v.get("decimal_date")}


def traits(dsid, min_levels=2):
    """Candidate epidemiological grouping variables with >1 level."""
    m = metadata(dsid)
    out = {}
    for col in ("geo_country", "division", "epi_group", "host"):
        vals = {k: v[col] for k, v in m.items() if v.get(col)}
        if len(set(vals.values())) >= min_levels and len(vals) >= 0.5 * len(m):
            out[col] = vals
    return out


# ------------------------------------------------------------------ loading
def new_namespace():
    return dendropy.TaxonNamespace()


def _get(path, schema, tns, rooted=False):
    return dendropy.Tree.get(path=path, schema=schema, taxon_namespace=tns,
                             preserve_underscores=True,
                             rooting="force-unrooted" if not rooted else "default-rooted")


def point_tree(dsid, kind, tns):
    """kind: 'ml' | 'mlrep2'..'mlrep5' | 'mb' | 'mbexp'  (Bayesian -> consensus)."""
    d = os.path.join(ROOT, "results", dsid)
    if kind.startswith("ml"):
        return _get(os.path.join(d, f"{dsid}_{kind}.treefile"), "newick", tns)
    return _get(os.path.join(d, f"{dsid}_{kind}.con.tre"), "nexus", tns)


def ufboot(dsid, tns, limit=None, seed=0):
    """The 1000 UFBoot replicate trees = the ML method's own topology distribution."""
    path = os.path.join(ROOT, "results", dsid, f"{dsid}_ml.ufboot")
    tl = dendropy.TreeList.get(path=path, schema="newick", taxon_namespace=tns,
                               preserve_underscores=True, rooting="force-unrooted")
    trees = list(tl)
    if limit and len(trees) > limit:
        trees = random.Random(seed).sample(trees, limit)
    return trees


def posterior(dsid, tns, tag="mb", burnin=0.25, limit=None, seed=0, runs=(1, 2)):
    """Post-burn-in posterior tree samples from MrBayes .run<N>.t files.

    Returns (all_trees, {run_index: [trees]}) so that the two independent runs
    can also be used separately as within-Bayesian replicates.
    """
    per_run, allt = {}, []
    for r in runs:
        path = os.path.join(ROOT, "results", dsid, f"{dsid}_{tag}.run{r}.t")
        if not os.path.exists(path):
            continue
        tl = dendropy.TreeList.get(path=path, schema="nexus", taxon_namespace=tns,
                                   preserve_underscores=True, rooting="force-unrooted")
        keep = list(tl)[int(len(tl) * burnin):]
        per_run[r] = keep
        allt.extend(keep)
    if limit:
        rng = random.Random(seed)
        if len(allt) > limit:
            allt = rng.sample(allt, limit)
        for r in per_run:
            if len(per_run[r]) > limit:
                per_run[r] = rng.sample(per_run[r], limit)
    return allt, per_run


def consensus(trees, tns, min_freq=0.5):
    """Majority-rule consensus of a tree list (used to build a Bayesian point
    estimate from an arbitrary subset of the posterior)."""
    tl = dendropy.TreeList(taxon_namespace=tns)
    tl.extend(trees)
    return tl.consensus(min_freq=min_freq)


# ------------------------------------------------------------ split sets
def split_set(tree):
    """Normalised, non-trivial split bitmasks of an UNROOTED tree."""
    tree.is_rooted = False
    tree.encode_bipartitions(suppress_unifurcations=True)
    return frozenset(b.split_bitmask for b in tree.bipartition_encoding
                     if not b.is_trivial())


def split_sets(trees):
    return [split_set(t) for t in trees]


def rf(a, b):
    """Unrooted Robinson-Foulds distance between two split sets."""
    return len(a ^ b)


def nrf(a, b, n_taxa):
    """RF normalised by its maximum, 2(n-3), for a fully resolved pair."""
    m = 2 * (n_taxa - 3)
    return rf(a, b) / m if m > 0 else 0.0


def split_frequencies(split_list):
    """{split -> frequency} over a list of split sets (posterior probability /
    bootstrap proportion, depending on which distribution was passed in)."""
    from collections import Counter
    c = Counter()
    for s in split_list:
        c.update(s)
    n = len(split_list)
    return {k: v / n for k, v in c.items()}, n


def mean_split_lengths(trees):
    """{split bitmask -> mean edge length over the trees that contain it}.

    This is how MrBayes labels an `allcompat` consensus, and the epidemiological
    readouts need it: root-to-tip rooting and patristic-distance clustering are
    both functions of branch length, so a consensus without lengths is unusable
    downstream.
    """
    acc = collections.defaultdict(lambda: [0.0, 0])
    for t in trees:
        t.is_rooted = False
        t.encode_bipartitions(suppress_unifurcations=True)
        for e in t.postorder_edge_iter():
            if e.tail_node is None or e.length is None or e.bipartition is None:
                continue
            a = acc[e.bipartition.split_bitmask]
            a[0] += e.length; a[1] += 1
    return {k: v[0] / v[1] for k, v in acc.items() if v[1]}


def greedy_consensus(trees_or_splitsets, tns, n_taxa=None, with_lengths=True):
    """Greedy (extended majority-rule / 'allcompat') consensus from a tree sample.

    WHY NOT `TreeList.consensus(min_freq=0.5)`: MrBayes summarises the posterior
    with `contype=allcompat`, which is FULLY RESOLVED.  A 50%-majority-rule tree
    built from the same sample is not, and Robinson-Foulds between two trees of
    different resolution is not comparable to RF between two resolved trees --
    it would make the within-Bayesian null (chain 1 vs chain 2) incommensurable
    with the between-method distance it is supposed to calibrate.

    This builds every Bayesian point estimate the same way MrBayes does: sort
    splits by sample frequency, add each one that is compatible with everything
    already accepted, stop at n-3 internal splits.

    Two unrooted splits A, B over the full taxon set F are compatible iff any of
    A&B, A&~B, ~A&B, ~A&~B is empty.
    """
    if trees_or_splitsets and isinstance(trees_or_splitsets[0], frozenset):
        sets = trees_or_splitsets
    else:
        sets = split_sets(trees_or_splitsets)
    n = n_taxa or len(tns)
    full = (1 << n) - 1
    freq, _ = split_frequencies(sets)
    accepted = []
    for s, f in sorted(freq.items(), key=lambda kv: (-kv[1], kv[0])):
        ok = True
        for a in accepted:
            if (s & a) and (s & ~a & full) and (~s & a & full) and (~s & ~a & full):
                ok = False
                break
        if ok:
            accepted.append(s)
        if len(accepted) >= n - 3:
            break
    t = dendropy.Tree.from_split_bitmasks(
        split_bitmasks=list(accepted) + [1 << i for i in range(n)],
        taxon_namespace=tns, is_rooted=False)
    t.is_rooted = False
    # frequency of each accepted split == its posterior probability / bootstrap
    t.split_support = {s: freq[s] for s in accepted}
    if with_lengths and not (trees_or_splitsets and
                             isinstance(trees_or_splitsets[0], frozenset)):
        lens = mean_split_lengths(trees_or_splitsets)
        t.encode_bipartitions(suppress_unifurcations=True)
        default = statistics.fmean(lens.values()) if lens else 1e-8
        for e in t.postorder_edge_iter():
            if e.tail_node is None:
                continue
            e.length = lens.get(e.bipartition.split_bitmask, default)
    return t

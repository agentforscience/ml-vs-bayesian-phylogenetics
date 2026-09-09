"""O(n) root-to-tip-regression rooting.

Why this exists
---------------
`mlbayes.epi.root_by_tip_dates` finds the TempEst rooting by cloning the tree
and re-rooting it once per candidate edge: O(n^2) tree copies.  That is fine for
a handful of point estimates but hopeless for the D3 analysis, which needs the
same rooting applied to 1000 UFBoot replicates and thousands of posterior trees
per dataset.

This module computes the regression sufficient statistics for EVERY candidate
root position in two tree traversals, using the standard up/down dynamic
program, and only materialises the single winning rooted tree.

The regression
--------------
For a root placed at distance `x` from node `v` along edge (u,v) of length L,
the root-to-tip distance of tip i is

    y_i = p_i + x            if i is below v      (p_i = dist(v, i))
    y_i = q_i + (L - x)      otherwise            (q_i = dist(u, i))

so the regression of y on sampling date d needs only these per-side sums:

    n, sum(p), sum(p^2), sum(d), sum(d*p)

DOWN(v) holds them for the leaves below v; UP(v) for the leaves outside.  Both
are obtained in one postorder and one preorder pass.  Shifting a side's
distances by a constant c updates the sums in closed form (`_shift`).

R^2(x) is then a ratio of low-order polynomials in x, maximised on a grid over
[0, L] refined by golden-section search.  Rootings implying a NEGATIVE clock
rate are rejected, exactly as `epi.root_by_tip_dates` does.
"""
import math
import dendropy

# a "side" is the 5-tuple (n, sum_p, sum_p2, sum_d, sum_dp)
_ZERO = (0, 0.0, 0.0, 0.0, 0.0)


def _merge(a, b):
    return tuple(x + y for x, y in zip(a, b))


def _shift(s, c):
    """Distances of every leaf on this side increase by c."""
    n, sp, sp2, sd, sdp = s
    return (n, sp + n * c, sp2 + 2 * c * sp + n * c * c, sd, sdp + c * sd)


def _sub(total, part):
    return tuple(x - y for x, y in zip(total, part))


def _r2_slope(sideA, sideB, x, L, sdd_total, d_mean, n_tot):
    """R^2 and slope of the regression of root-to-tip distance on date, for a
    root at distance x from the head of the edge (sideA below it)."""
    a = _shift(sideA, x)
    b = _shift(sideB, L - x)
    n, sp, sp2, sd, sdp = _merge(a, b)
    if n != n_tot or n < 3:
        return None
    syy = sp2 - sp * sp / n
    sdy = sdp - sd * sp / n
    if sdd_total <= 0 or syy <= 0:
        return None
    return (sdy * sdy) / (sdd_total * syy), sdy / sdd_total


def _best_x(sideA, sideB, L, sdd, d_mean, n_tot, grid=17, refine=25):
    """Maximise R^2 over the root position on one edge."""
    if L <= 0:
        r = _r2_slope(sideA, sideB, 0.0, 0.0, sdd, d_mean, n_tot)
        return (r[0], 0.0, r[1]) if r else None
    best = None
    for i in range(grid + 1):
        x = L * i / grid
        r = _r2_slope(sideA, sideB, x, L, sdd, d_mean, n_tot)
        if r and (best is None or r[0] > best[0]):
            best = (r[0], x, r[1])
    if best is None:
        return None
    # golden-section refinement around the grid winner
    step = L / grid
    lo, hi = max(0.0, best[1] - step), min(L, best[1] + step)
    phi = (math.sqrt(5) - 1) / 2
    c, d = hi - phi * (hi - lo), lo + phi * (hi - lo)
    fc = _r2_slope(sideA, sideB, c, L, sdd, d_mean, n_tot)
    fd = _r2_slope(sideA, sideB, d, L, sdd, d_mean, n_tot)
    for _ in range(refine):
        if fc is None or fd is None:
            break
        if fc[0] > fd[0]:
            hi, d, fd = d, c, fc
            c = hi - phi * (hi - lo)
            fc = _r2_slope(sideA, sideB, c, L, sdd, d_mean, n_tot)
        else:
            lo, c, fc = c, d, fd
            d = lo + phi * (hi - lo)
            fd = _r2_slope(sideA, sideB, d, L, sdd, d_mean, n_tot)
        cand = max([f for f in (fc, fd) if f], key=lambda f: f[0], default=None)
        if cand and cand[0] > best[0]:
            best = (cand[0], c if cand is fc else d, cand[1])
    return best


def rtt_root(tree, dates, require_positive_slope=True, inplace=False):
    """Root `tree` at the position maximising root-to-tip vs date R^2.

    Parameters
    ----------
    tree   : dendropy.Tree (unrooted or arbitrarily rooted); NOT modified
             unless inplace=True.
    dates  : {taxon_label: decimal_date}.  Tips missing from `dates` are
             excluded from the regression but kept in the tree.

    Returns a NEW rooted dendropy.Tree carrying attributes
    `rtt_r2`, `rtt_rate`, `rtt_n`.  Falls back to midpoint rooting when fewer
    than 3 tips are dated or no rooting yields a positive rate.
    """
    t = tree if inplace else tree.clone(depth=1)
    t.is_rooted = True

    dated = [lf for lf in t.leaf_node_iter()
             if lf.taxon is not None and lf.taxon.label in dates]
    n_tot = len(dated)
    if n_tot < 3:
        t.reroot_at_midpoint(update_bipartitions=True, suppress_unifurcations=False)
        t.rtt_r2, t.rtt_rate, t.rtt_n = None, None, n_tot
        return t

    ds = [float(dates[lf.taxon.label]) for lf in dated]
    d_mean = sum(ds) / n_tot
    sdd = sum((d - d_mean) ** 2 for d in ds)

    # ---- postorder: DOWN(v) = stats over dated leaves below v -------------
    down = {}
    for nd in t.postorder_node_iter():
        if nd.is_leaf():
            if nd.taxon is not None and nd.taxon.label in dates:
                d = float(dates[nd.taxon.label])
                down[nd] = (1, 0.0, 0.0, d, 0.0)
            else:
                down[nd] = _ZERO
            continue
        acc = _ZERO
        for ch in nd.child_nodes():
            acc = _merge(acc, _shift(down[ch], ch.edge.length or 0.0))
        down[nd] = acc
    total = down[t.seed_node]
    assert total[0] == n_tot, (total[0], n_tot)

    # ---- preorder: UP(v) = stats over dated leaves OUTSIDE v's subtree ----
    up = {t.seed_node: _ZERO}
    for nd in t.preorder_node_iter():
        if nd is t.seed_node:
            continue
        par = nd.parent_node
        sibs = _ZERO
        for ch in par.child_nodes():
            if ch is nd:
                continue
            sibs = _merge(sibs, _shift(down[ch], ch.edge.length or 0.0))
        # everything outside `par`, plus the siblings' subtrees, all measured
        # from `par`, then lifted across the edge par->nd
        up[nd] = _shift(_merge(up[par], sibs), nd.edge.length or 0.0)

    # ---- scan every edge --------------------------------------------------
    best = None   # (r2, node, x, slope)
    for nd in t.preorder_node_iter():
        if nd is t.seed_node or nd.edge.length is None:
            continue
        L = nd.edge.length
        # the root sits on the edge above nd; sideA = below nd (already at nd),
        # sideB = the rest measured from nd's parent
        sideA = down[nd]
        # up[nd] holds the outside leaves measured FROM nd, i.e. already
        # including this edge; shift back by -L to measure from nd's parent
        sideB = _shift(up[nd], -L)
        if sideA[0] == 0 or sideB[0] == 0:
            continue
        r = _best_x(sideA, sideB, L, sdd, d_mean, n_tot)
        if r is None:
            continue
        r2, x, slope = r
        if require_positive_slope and slope <= 0:
            continue
        if best is None or r2 > best[0]:
            best = (r2, nd, x, slope)

    if best is None:
        t.reroot_at_midpoint(update_bipartitions=True, suppress_unifurcations=False)
        t.rtt_r2, t.rtt_rate, t.rtt_n = None, None, n_tot
        return t

    r2, nd, x, slope = best
    L = nd.edge.length
    # dendropy: length1 = new edge from root to nd, length2 = root to parent
    t.reroot_at_edge(nd.edge, length1=max(x, 0.0), length2=max(L - x, 0.0),
                     update_bipartitions=True)
    t.rtt_r2, t.rtt_rate, t.rtt_n = r2, slope, n_tot
    return t


def root_at_leafset(tree, leafset, dates=None):
    """Root at the branch subtending `leafset` if it is monophyletic.

    Used to hold rooting CONSTANT across a tree distribution: the root clade is
    fixed once from a point estimate, then imposed on every replicate.  Returns
    (rooted_tree, True) on success, or (rtt/midpoint-rooted tree, False) when
    `leafset` is not a clade in this tree.
    """
    t = tree.clone(depth=1)
    t.is_rooted = True
    want = frozenset(leafset)
    labels = {lf.taxon.label for lf in t.leaf_node_iter()}
    other = labels - want
    for target in (want, other):
        if not target or target == labels:
            continue
        taxa = [tx for tx in t.taxon_namespace if tx.label in target]
        try:
            mrca = t.mrca(taxa=taxa)
        except Exception:
            continue
        if mrca is None or mrca.parent_node is None:
            continue
        if frozenset(lf.taxon.label for lf in mrca.leaf_iter()) != target:
            continue
        t.reroot_at_edge(mrca.edge, update_bipartitions=True)
        return t, True
    if dates:
        return rtt_root(tree, dates), False
    t2 = tree.clone(depth=1)
    t2.reroot_at_midpoint(update_bipartitions=True, suppress_unifurcations=False)
    return t2, False

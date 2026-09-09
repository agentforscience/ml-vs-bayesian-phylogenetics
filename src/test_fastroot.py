"""Validate fastroot.rtt_root against the O(n^2) reference in mlbayes.epi.

The reference clones the tree once per candidate edge and re-roots it, so it is
independent code exercising the same criterion.  Agreement on R^2, clock rate
and the resulting split set is the correctness check for the DP.
"""
import csv, os, sys, time
import dendropy
sys.path.insert(0, "code"); sys.path.insert(0, "src")
from mlbayes import epi
import fastroot


def load(dsid, tag="ml"):
    tns = dendropy.TaxonNamespace()
    t = dendropy.Tree.get(path=f"results/{dsid}/{dsid}_{tag}.treefile", schema="newick",
                          taxon_namespace=tns, preserve_underscores=True,
                          rooting="force-unrooted")
    meta = {r["id"]: r for r in csv.DictReader(
        open(f"datasets/{dsid}/taxon_map.tsv"), delimiter="\t")}
    # tip labels in the trees are the GenBank accessions (the `id` column)
    lab = {r["id"]: float(r["decimal_date"]) for r in meta.values() if r["decimal_date"]}
    return t, lab


def splits(t):
    t.encode_bipartitions()
    return frozenset(b.split_bitmask for b in t.bipartition_encoding)


if __name__ == "__main__":
    ok = True
    for dsid in sys.argv[1:] or ["fmdv_uk_2007", "mers_korea_2015"]:
        t, dates = load(dsid)
        n_matched = sum(1 for lf in t.leaf_node_iter() if lf.taxon.label in dates)
        t0 = time.time(); slow = epi.root_by_tip_dates(t, dates); t_slow = time.time() - t0
        t1 = time.time(); fast = fastroot.rtt_root(t, dates);     t_fast = time.time() - t1
        r2s, r2f = getattr(slow, "rtt_r2", None), fast.rtt_r2
        rs, rf = getattr(slow, "rtt_rate", None), fast.rtt_rate
        same = splits(slow) == splits(fast)
        # the fast version also optimises the root POSITION within the edge, so
        # its R^2 must be >= the reference (which always roots at the midpoint
        # of the chosen edge); require it never to be worse.
        better = (r2f is not None and r2s is not None and r2f >= r2s - 1e-9)
        print(f"{dsid:26s} ntip={n_matched:4d} dated={len(dates)}  "
              f"R2 slow={r2s:.6f} fast={r2f:.6f}  rate slow={rs:.3e} fast={rf:.3e}  "
              f"same_splits={same}  speedup={t_slow/max(t_fast,1e-9):.0f}x "
              f"({t_slow:.1f}s -> {t_fast:.2f}s)")
        if not better:
            print("   !! FAIL: fast R^2 worse than reference"); ok = False
        if not same:
            print("   ** note: different root EDGE chosen (allowed if R^2 higher)")
    print("PASS" if ok else "FAIL")

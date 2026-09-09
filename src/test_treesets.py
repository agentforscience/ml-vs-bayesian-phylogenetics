"""Check the fast split-set RF against dendropy.treecompare.symmetric_difference."""
import sys, random
sys.path.insert(0, "src")
import dendropy
from dendropy.calculate import treecompare
import treesets as T

if __name__ == "__main__":
    dsid = sys.argv[1] if len(sys.argv) > 1 else "fmdv_uk_2007"
    tns = T.new_namespace()
    trees = [T.point_tree(dsid, k, tns) for k in ("ml", "mlrep2", "mlrep3", "mlrep4", "mlrep5")]
    ss = T.split_sets(trees)
    bad = 0
    for i in range(len(trees)):
        for j in range(i + 1, len(trees)):
            a = trees[i].clone(depth=1); b = trees[j].clone(depth=1)
            a.encode_bipartitions(); b.encode_bipartitions()
            ref = treecompare.symmetric_difference(a, b)
            mine = T.rf(ss[i], ss[j])
            flag = "OK " if ref == mine else "BAD"
            if ref != mine: bad += 1
            print(f"{flag} pair({i},{j}) dendropy={ref:4d} fast={mine:4d}")
    print("PASS" if bad == 0 else f"FAIL ({bad} mismatches)")

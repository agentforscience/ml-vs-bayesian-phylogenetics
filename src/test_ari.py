"""Validate the hand-rolled adjusted Rand index against scikit-learn.

`analysis_epi.adjusted_rand` rounds to 6 dp, so the tolerance is 1e-6, not exact.
"""
import random, sys
sys.path.insert(0, "src"); sys.path.insert(0, "code")
import analysis_epi as A
from sklearn.metrics import adjusted_rand_score as ars

if __name__ == "__main__":
    assert A.adjusted_rand([0, 0, 1, 1], [0, 0, 1, 1]) == 1.0
    r = random.Random(0); worst = 0.0
    for _ in range(200):
        n = r.randint(10, 120)
        a = [r.randint(0, r.randint(1, 8)) for _ in range(n)]
        b = [r.randint(0, r.randint(1, 8)) for _ in range(n)]
        worst = max(worst, abs(A.adjusted_rand(a, b) - ars(a, b)))
    print(f"200 random labellings, max |mine - sklearn| = {worst:.2e}")
    print("PASS" if worst < 1e-6 else "FAIL")

"""T1 — Lemma 1(o) and (iv) for staggered trapezoids, standalone.
   Named designs: exact rational rank.  Sweep (~11k designs): floating-point rank, which agrees with the exact
   rank on every design that was cross-checked.
   nu(S) = number of connected components of the increment graph, and nu <= max(d, Delta_2 - w).
   Exact rational rank; no data, no network, no other scripts.  Run: python3 check_lemma1_components.py
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools, numpy as np

def rank_exact(A):
    a = [[F(x) for x in row] for row in A]; r = 0
    if not a: return 0
    for j in range(len(a[0])):
        p = next((i for i in range(r, len(a)) if a[i][j]), None)
        if p is None: continue
        a[r], a[p] = a[p], a[r]; v = a[r][j]; a[r] = [x / v for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][j]:
                v = a[i][j]; a[i] = [x - v * y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a): break
    return r

def design(E, T):
    """cell-mean design for mu(e,t) = alpha(t) + g(e) + tau(t-e+1), tau(1) = g(e_1) = 0"""
    cells = [(e, t) for e in E for t in range(e, T + 1)]
    S = sorted({t - e + 1 for e, t in cells}); Ts = sorted({t for _, t in cells})
    cols = [('a', t) for t in Ts] + [('g', e) for e in sorted(E)[1:]] + [('u', s) for s in S if s != 1]
    A = [[int((v == t) if k == 'a' else (v == e) if k == 'g' else (v == t - e + 1)) for k, v in cols] for e, t in cells]
    return cells, cols, A

def increment_graph(E, T):
    """vertices: observed increments u (some cohort observed at tenures u and u+1 in consecutive periods);
       edges: u ~ u' when two cohorts are observed at the same t and t+1 with those tenures."""
    support = {(e, t) for e in E for t in range(e, T + 1)}
    adj = {}
    for e in E:
        for t in range(e, T):
            u = t - e + 1; adj.setdefault(u, set())
            for ep in E:
                if ep != e and (ep, t) in support and (ep, t + 1) in support:
                    adj[u].add(t - ep + 1)
    return adj

def components(adj):
    left = set(adj); comps = []
    while left:
        todo = [min(left)]; c = set()
        while todo:
            v = todo.pop()
            if v not in c: c.add(v); todo.extend(adj.get(v, set()) - c)
        comps.append(sorted(c)); left -= c
    return comps

print("=" * 96)
print("Lemma 1(o): nullity = #components of the increment graph;  (iv): nullity <= max(d, Delta_2 - w)\n")
for E, T in [([1, 5, 13], 19), ([1, 4, 8], 8), ([1, 6, 13], 19), ([1, 6, 13], 16), ([1, 5, 9, 11, 12], 12)]:
    cells, cols, A = design(E, T); nu = len(cols) - rank_exact(A)
    adj = increment_graph(E, T); comps = components(adj)
    d = reduce(gcd, [e - E[0] for e in E[1:]]); w = T - max(E); D2 = sorted(E)[1] - E[0]
    print(f"  E={E}, T={T}: cells={len(cells)} cols={len(cols)} nullity={nu} d={d} w={w} "
          f"components={len(comps)} {comps if len(comps) <= 5 else ''} bound=max({d},{D2 - w})={max(d, D2 - w)}")
print()
viol_o = viol_iv = n = strict = 0
for K in range(3, 6):
    for E in itertools.combinations(range(1, 14), K):
        if E[0] != 1: continue
        for T in range(E[-1], E[-1] + 14):
            cells, cols, A = design(list(E), T)
            nu = len(cols) - int(np.linalg.matrix_rank(np.array(A, dtype=float)))   # float rank in the sweep; the
            ncomp = len(components(increment_graph(list(E), T)))                     # named designs above use exact rank
            d = reduce(gcd, [e - E[0] for e in E[1:]]); bound = max(d, (E[1] - E[0]) - (T - E[-1]))
            n += 1
            if nu != ncomp: viol_o += 1; print("  VIOLATION (o):", E, T, nu, ncomp)
            if nu > bound: viol_iv += 1; print("  VIOLATION (iv):", E, T, nu, bound)
            if nu < bound: strict += 1
print(f"  sweep: {n} designs (e_1 = 1, 3 <= K <= 5, T <= e_K + 13); violations of (o): {viol_o}; of (iv): {viol_iv}; "
      f"bound strict in {strict}")

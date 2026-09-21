#!/usr/bin/env python3
"""Counterexample to Theorem 1 of T1 as currently stated, and the exact repair.

Found 2026-09-20 after a second reader reported "a counterexample to the central
identification theorem".  This script locates it, characterises it exactly, and
delimits the damage.

MODEL      mu(e,t) = alpha(t) + g(e) + tau(s),  s = t - e + 1,
           normalisations tau(1) = 0, g(e0) = 0.

THEOREM 1 (current) claims the identified set for tau is the ONE-dimensional
family {tau(s) + c(s-1)}.

THE GAP    The converse proof argues: Delta delta_tau(s) = Delta delta_tau(s')
           for tenures s, s' co-observed at a common calendar time, "hence
           Delta delta_tau is a constant".  That last step needs the
           co-observation graph on tenure values to be CONNECTED.  Tenures
           co-observed at time t differ by differences of entry times, so the
           graph has d = gcd{e - e' : e, e' in E} components and
           Delta delta_tau is constant only on residue classes mod d.

RESULT     dim(identified set) = d, not 1.  Verified on 14/14 designs below.
           JLPS (entries at waves 1, 5, 13) has d = 4.

DAMAGE     Corollary 2 ("all second differences of tau are identified") is FALSE
           whenever d > 1 -- the paper's positive result.
           Section 6 (Theorems 3-5) SURVIVES, for a reason the current proof does
           not give: d divides every entry difference, so all entries are
           congruent mod d, so s mod d is a function of t alone and is absorbed
           by the time effects.  The identified set is
               {affine in s}  +  {functions of s mod d}
           and two-way fixed effects annihilate BOTH.  That is a strictly
           stronger absorption theorem than the one in the paper.
"""
import numpy as np
from math import gcd
from functools import reduce

def identified_set(entries, Tmax):
    """Null space of the cell-mean map, returned as tau-parts indexed by tenure."""
    cells = sorted({(e, t) for e in entries for t in range(e, Tmax + 1)})
    E = sorted({e for e, _ in cells}); T = sorted({t for _, t in cells})
    S = sorted({t - e + 1 for e, t in cells}); e0 = E[0]
    cols = [('a', t) for t in T] + [('g', e) for e in E if e != e0] \
         + [('u', s) for s in S if s != 1]
    idx = {c: k for k, c in enumerate(cols)}
    M = np.zeros((len(cells), len(cols)))
    for r, (e, t) in enumerate(cells):
        s = t - e + 1
        M[r, idx[('a', t)]] = 1
        if e != e0: M[r, idx[('g', e)]] = 1
        if s != 1:  M[r, idx[('u', s)]] = 1
    rank = np.linalg.matrix_rank(M, tol=1e-10)
    _, _, Vt = np.linalg.svd(M)
    B = Vt[rank:].T
    out = np.zeros((B.shape[1], max(S) + 1))
    for k in range(B.shape[1]):
        for s in S:
            if s != 1: out[k, s] = B[idx[('u', s)], k]
    return out, S

def twoway_annihilator(entries, Tmax, n_per=40):
    rows = []
    u = 0
    for e in entries:
        for _ in range(n_per):
            rows += [(u, t, t - e + 1) for t in range(e, Tmax + 1)]
            u += 1
    a = np.array(rows); i, t, s = a[:, 0], a[:, 1], a[:, 2]
    un, tu = np.unique(i), np.unique(t)
    W = np.hstack([np.eye(len(un))[np.searchsorted(un, i)],
                   np.eye(len(tu))[np.searchsorted(tu, t)][:, 1:]])
    return W, s

if __name__ == '__main__':
    designs = [([1,2],10), ([1,2,3],10), ([1,3],10), ([1,3,5],12), ([1,4,7],14),
               ([1,5,13],19), ([1,5,9],16), ([1,2,5],12), ([1,6],14), ([1,7,13],20),
               ([1,5,13],25), ([1,2,5,13],19), ([1,3,7],14), ([1,9],20)]
    print(f"{'entries':<16}{'Tmax':>5}{'d=gcd':>7}{'dim':>5}{'dim==d':>8}"
          f"{'max|D2 tau|':>13}{'Cor.2':>8}{'max|M.dir|':>12}{'Thm 3':>8}")
    for entries, Tmax in designs:
        d = reduce(gcd, [abs(a - b) for a in entries for b in entries if a != b])
        B, S = identified_set(entries, Tmax)
        W, s = twoway_annihilator(entries, Tmax)
        d2, ann = 0.0, 0.0
        for k in range(B.shape[0]):
            v = B[k] / (np.abs(B[k]).max() or 1)
            d2  = max(d2, np.abs(np.diff(np.diff(v[1:]))).max())
            col = v[s].reshape(-1, 1)
            ann = max(ann, np.abs(col - W @ np.linalg.lstsq(W, col, rcond=None)[0]).max())
        print(f"{str(entries):<16}{Tmax:>5}{d:>7}{B.shape[0]:>5}{str(B.shape[0]==d):>8}"
              f"{d2:>13.3f}{('holds' if d2<1e-9 else 'FAILS'):>8}"
              f"{ann:>12.1e}{('holds' if ann<1e-8 else 'FAILS'):>8}")

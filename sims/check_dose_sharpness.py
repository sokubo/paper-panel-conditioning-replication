"""T1 v0.8 — sharpness of the interrupted-participation kernel (Theorem 2(b)), exact arithmetic.

A counterexample refuted the v0.7 statement, whose condition (P) asked only that every
spacing in dZ ∩ [1, j_K] be REALIZED by some cohort pair.  The corrected sufficient condition (P') is
structural:

  (P')  entry cohorts form an arithmetic progression e_0, e_0 + d, ..., e_0 + Ld  (continuous recruitment
        at stride d); every scheduled cell is observed (field period through e_0 + Ld + j_K); the stride is a
        spacing of the pattern, d ∈ J − J; and L·d ≥ Γ(J, d), the largest gap between consecutive members of J
        that lie in the same residue class modulo d.

Under (P') the tau-projection of the kernel is exactly {m j(k) + rho(j(k))}, of dimension = number of residue
classes mod d represented in J.  This script (i) reproduces the counterexample (tau-kernel dimension 2, not 1),
(ii) sweeps patterns, strides and window lengths and confirms that whenever (P') holds the dimension equals the
prediction, (iii) shows (P') is not necessary, and (iv) confirms the CPS threshold of ten monthly cohorts.
Run: python3 check_dose_sharpness.py      (about a minute)
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools

def rank_and_kernel(A, ncol):
    a = [[F(x) for x in row] for row in A]; piv = []; r = 0
    for j in range(ncol):
        p = next((i for i in range(r, len(a)) if a[i][j]), None)
        if p is None: continue
        a[r], a[p] = a[p], a[r]; v = a[r][j]; a[r] = [x / v for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][j]:
                v = a[i][j]; a[i] = [x - v * y for x, y in zip(a[i], a[r])]
        piv.append(j); r += 1
        if r == len(a): break
    free = [j for j in range(ncol) if j not in piv]; basis = []
    for f in free:
        v = [F(0)] * ncol; v[f] = F(1)
        for i, j in enumerate(piv): v[j] = -a[i][f]
        basis.append(v)
    return r, basis

def rank_of(vectors):
    if not vectors: return 0
    return rank_and_kernel([list(v) for v in vectors], len(vectors[0]))[0]

def dose_design(E, J, T):
    cells = [(e, e + j, k + 1) for e in E for k, j in enumerate(J) if e + j <= T]
    Ts = sorted({t for _, t, _ in cells}); Ks = sorted({k for _, _, k in cells}); e0 = min(E)
    cols = [('a', t) for t in Ts] + [('g', e) for e in sorted(E) if e != e0] + [('u', k) for k in Ks if k != 1]
    idx = {c: i for i, c in enumerate(cols)}; A = []
    for e, t, k in cells:
        row = [0] * len(cols); row[idx[('a', t)]] = 1
        if e != e0: row[idx[('g', e)]] = 1
        if k != 1: row[idx[('u', k)]] = 1
        A.append(row)
    return cells, cols, A, Ks

def tau_dim(E, J, T):
    """full nullity and the dimension of the tau-PROJECTION of the kernel (rank of projected basis)."""
    cells, cols, A, Ks = dose_design(E, J, T)
    r, ker = rank_and_kernel(A, len(cols))
    tau_idx = [i for i, (kind, _) in enumerate(cols) if kind == 'u']
    proj = [[v[i] for i in tau_idx] for v in ker]
    return len(cols) - r, rank_of(proj), len(cells), len(cols), r

def classes_represented(J, d):
    return len({j % d for j in J}) if d > 0 else 1

def Gamma(J, d):
    """largest gap between consecutive members of J within a residue class mod d (0 if every class is a singleton)"""
    g = 0
    for rcls in {j % d for j in J}:
        js = sorted(j for j in J if j % d == rcls)
        for a, b in zip(js, js[1:]): g = max(g, b - a)
    return g

def P_prime(J, d, L):
    return (d in {b - a for a in J for b in J if b > a}) and (L * d >= Gamma(J, d))

if __name__ == "__main__":
    print("(i) the counterexample to the v0.7 statement: E = {1,2,4}, J = {0,1,3}, T = 7")
    nul, td, nc, ncol, r = tau_dim([1, 2, 4], [0, 1, 3], 7)
    print(f"    cells={nc} cols={ncol} rank={r} full nullity={nul}  tau-projection dimension={td}  "
          f"(v0.7 asserted 1; the cohort set has a gap, (P') fails)")
    assert td == 2

    print("\n(ii) sweep: arithmetic-progression cohorts, full field period; (P') predicts dim = #classes mod d")
    patterns = [[0,1,3], [0,1,2,3,12,13,14,15], [0,1,4,5], [0,2,4,6,8,10], [0,1,2,3], [0,3], [0,1,5], [0,2,3,7],
                [0,1,2,6,7,8], [0,4,8,9], [0,1,3,4,10,11]]
    ok_pred = fail_pred = 0; not_nec = 0; strict_fail = 0; worst = None
    for J in patterns:
        jK = max(J)
        for d in (1, 2, 3):
            for L in range(1, 17):
                E = [1 + d * i for i in range(L + 1)]
                T = max(E) + jK
                nul, td, *_ = tau_dim(E, J, T)
                pred = classes_represented(J, d)
                if P_prime(J, d, L):
                    if td == pred: ok_pred += 1
                    else:
                        fail_pred += 1; worst = (J, d, L, td, pred)
                else:
                    if td == pred: not_nec += 1
                    elif td > pred: strict_fail += 1
    print(f"    (P') holds and dim == #classes : {ok_pred}")
    print(f"    (P') holds and dim != #classes : {fail_pred}   <- must be 0 for the theorem")
    print(f"    (P') fails but dim == #classes : {not_nec}   (the condition is sufficient, not necessary)")
    print(f"    (P') fails and dim  > #classes : {strict_fail}   (the condition is not vacuous)")
    assert fail_pred == 0, worst

    print("\n(iii) gapped cohort sets with the same pattern can exceed the prediction even when spacings are realized")
    for E in ([1, 2, 4], [1, 2, 3, 4], [1, 3, 4], [1, 2, 4, 5]):
        nul, td, *_ = tau_dim(E, [0, 1, 3], 8)
        print(f"    E={E}: tau-dim={td} (prediction under (P') would be 1; (P') {'holds' if E == [1,2,3,4] else 'fails'})")

    print("\n(iv) CPS pattern, monthly entry: Gamma = 9, so (P') needs L >= 9, i.e. ten cohorts")
    Jc = [0,1,2,3,12,13,14,15]
    print(f"    Gamma(J,1) = {Gamma(Jc,1)};  1 in J-J: {1 in {b-a for a in Jc for b in Jc if b>a}}")
    for L in (7, 8, 9, 10, 12):
        E = list(range(1, L + 2)); nul, td, *_ = tau_dim(E, Jc, max(E) + 15)
        print(f"    {L+1:2d} cohorts: tau-dim={td}  (P') {'holds' if P_prime(Jc,1,L) else 'fails'}")
        assert (td == 1) == P_prime(Jc, 1, L)
    print("\n(v) the pattern's difference set omits 4..8, which (P) as written required and the CPS cannot realize")
    print("    J-J =", sorted({b-a for a in Jc for b in Jc if b>a}))

    print("\n(vi) two facts cited in the text: (c) for d > 1 is not spanned by within-class divided differences, and")
    print("     (d) for d > 1 an unequally spaced pattern can still contain the dose-linear path")
    from fractions import Fraction as F
    E = [1, 3, 5]; J = [0, 1, 2, 3]; nul, td, *_ = tau_dim(E, J, max(E) + 3)
    print(f"    J={J}, d=2: tau-dim={td}, identified functionals={(len(J)-1)-td} (one cross-class contrast; no class has three offsets)")
    assert td == 2
    E = [1, 4, 7]; J = [0, 1, 3, 4]; cells, cols, A, Ks = dose_design(E, J, max(E) + 4); r, ker = rank_and_kernel(A, len(cols))
    ti = [i for i, (k, _) in enumerate(cols) if k == 'u']; proj = [[v[i] for i in ti] for v in ker]
    lin = [F(k - 1) for k in Ks if k != 1]; in_ker = rank_of([list(v) for v in proj] + [lin]) == rank_of(proj)
    print(f"    J={J}, d=3: tau-dim={rank_of(proj)}, dose-linear path in kernel: {in_ker}  (k-1 = (2/3) j + rho, rho = 1/3 on the class of 1)")
    assert in_ker
    print("\nall assertions passed")

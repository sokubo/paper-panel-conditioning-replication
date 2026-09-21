"""T1 v0.6 — interrupted participation (dose designs) and the Lemma-1 upper bound for K >= 4.
Exact rational ranks; no data, no network.  Run: python3 check_interrupted.py
(A) For a common interview pattern J (offsets from entry) the cell design is mu(e, e+j) = alpha(e+j) + g(e) + tau(k),
    k = the index of j in J (the DOSE).  Claim: the kernel's tau-part is spanned by j(k) — the elapsed calendar time
    of the k-th interview — and by d-periodic functions of j(k); with uninterrupted participation j(k) = k - 1 and we
    recover the affine direction.  Consequence: a dose-linear path is NOT in the kernel iff participation is interrupted.
(B) Lemma 1 upper bound for K >= 4: nu <= max(d, Delta_2 - (T - e_3)), the follow-up of the THIRD cohort.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools, json

def rank_and_kernel(A, ncol):
    """exact rank and a basis of the kernel (list of rational vectors)"""
    a = [[F(x) for x in row] for row in A]
    piv = []; r = 0
    for j in range(ncol):
        p = next((i for i in range(r, len(a)) if a[i][j]), None)
        if p is None: continue
        a[r], a[p] = a[p], a[r]; v = a[r][j]; a[r] = [x / v for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][j]:
                v = a[i][j]; a[i] = [x - v * y for x, y in zip(a[i], a[r])]
        piv.append(j); r += 1
        if r == len(a): break
    free = [j for j in range(ncol) if j not in piv]
    basis = []
    for f in free:
        v = [F(0)] * ncol; v[f] = F(1)
        for i, j in enumerate(piv): v[j] = -a[i][f]
        basis.append(v)
    return r, basis

def dose_design(E, J, T):
    cells = [(e, e + j, k + 1) for e in E for k, j in enumerate(J) if e + j <= T]
    Ts = sorted({t for _, t, _ in cells}); Ks = sorted({k for _, _, k in cells}); e0 = min(E)
    cols = [('a', t) for t in Ts] + [('g', e) for e in sorted(E) if e != e0] + [('u', k) for k in Ks if k != 1]
    idx = {c: i for i, c in enumerate(cols)}
    A = []
    for e, t, k in cells:
        row = [0] * len(cols); row[idx[('a', t)]] = 1
        if e != e0: row[idx[('g', e)]] = 1
        if k != 1: row[idx[('u', k)]] = 1
        A.append(row)
    return cells, cols, A, Ks

def analyse(name, E, J, T, verbose=True):
    cells, cols, A, Ks = dose_design(E, J, T)
    r, ker = rank_and_kernel(A, len(cols))
    nullity = len(cols) - r
    tau_pos = {k: i for i, (kind, k) in enumerate(cols) if kind == 'u'}
    tau_basis = [[(v[tau_pos[k]] if k in tau_pos else F(0)) for k in Ks] for v in ker]
    jk = [J[k - 1] for k in Ks]                      # elapsed calendar time of the k-th interview
    d = reduce(gcd, [e - min(E) for e in sorted(E)[1:]]) if len(E) > 1 else 0
    # is the kernel's tau-part spanned by {j(k)} (+ d-periodic functions of j(k))?
    prop_to_j = None
    if nullity >= 1:
        v = tau_basis[0]
        nz = [i for i, x in enumerate(v) if x != 0]
        prop_to_j = bool(nz) and all(v[i] * F(jk[nz[0]]) == v[nz[0]] * F(jk[i]) for i in range(len(v)))
    # is a DOSE-linear path in the kernel? (tau(k) = k-1)
    lin_dose = [F(k - 1) for k in Ks]
    in_ker = in_span(lin_dose, tau_basis)
    lin_cal = [F(j) for j in jk]
    cal_in_ker = in_span(lin_cal, tau_basis)
    if verbose:
        print(f"{name}: cells={len(cells)} cols={len(cols)} rank={r} nullity={nullity} d={d}")
        print(f"   doses {Ks}  elapsed j(k) {jk}")
        print(f"   kernel tau-basis: {[[str(x) for x in b] for b in tau_basis]}")
        print(f"   first basis vector proportional to j(k)? {prop_to_j}")
        print(f"   dose-linear path tau(k)=k-1 in kernel (undetectable)? {in_ker}")
        print(f"   calendar-linear path tau(k)=j(k) in kernel (undetectable)? {cal_in_ker}")
    return dict(name=name, cells=len(cells), cols=len(cols), rank=r, nullity=nullity, d=d, doses=Ks, j=jk,
                prop_to_j=prop_to_j, dose_linear_in_kernel=in_ker, calendar_linear_in_kernel=cal_in_ker)

def in_span(v, basis):
    if not basis: return all(x == 0 for x in v)
    M = [list(b) for b in basis] + [list(v)]
    rows = len(M); ncol = len(v)
    r1, _ = rank_and_kernel([list(b) for b in basis], ncol)
    r2, _ = rank_and_kernel(M, ncol)
    return r1 == r2

out = {}
print("=" * 100); print("(A) dose designs\n")
out['uninterrupted_8'] = analyse("uninterrupted, 8 waves, monthly entry", list(range(1, 25)), list(range(8)), 39)
print()
out['cps_4_8_4'] = analyse("CPS 4-8-4, monthly entry", list(range(1, 25)), [0, 1, 2, 3, 12, 13, 14, 15], 39)
print()
out['cps_truncated'] = analyse("CPS 4-8-4, monthly entry, T truncates the last cohorts", list(range(1, 25)), [0, 1, 2, 3, 12, 13, 14, 15], 30)
print()
out['gap_2_2_2'] = analyse("2 on, 2 off, 2 on", list(range(1, 25)), [0, 1, 4, 5], 39)
print()
out['cps_biennial_entry'] = analyse("CPS pattern but entry every 2 periods (stride 2)", list(range(1, 49, 2)), [0, 1, 2, 3, 12, 13, 14, 15], 70)
print()
out['lfs_6'] = analyse("6-wave rotation, 1 on 1 off (alternating months)", list(range(1, 25)), [0, 2, 4, 6, 8, 10], 45)

print()
print("=" * 100); print("(B) Lemma 1 upper bound for K >= 4:  nu <= max(d, Delta_2 - (T - e_3))\n")
import numpy as np
def trap_nullity(E, T):
    cells = [(e, t) for e in E for t in range(e, T + 1)]
    Ts = sorted({t for _, t in cells}); S = sorted({t - e + 1 for e, t in cells}); e0 = min(E)
    cols = [('a', t) for t in Ts] + [('g', e) for e in sorted(E)[1:]] + [('u', s) for s in S if s != 1]
    idx = {c: i for i, c in enumerate(cols)}
    A = np.zeros((len(cells), len(cols)))
    for r, (e, t) in enumerate(cells):
        A[r, idx[('a', t)]] = 1
        if e != e0: A[r, idx[('g', e)]] = 1
        if t - e + 1 != 1: A[r, idx[('u', t - e + 1)]] = 1
    return len(cols) - np.linalg.matrix_rank(A)
viol_new = viol_old = n = 0; strict = 0
for K in range(3, 6):
    for E in itertools.combinations(range(1, 14), K):
        if E[0] != 1: continue
        for T in range(E[-1], E[-1] + 14):
            nu = trap_nullity(list(E), T)
            d = reduce(gcd, [e - E[0] for e in E[1:]]); D2 = E[1] - E[0]
            b_new = max(d, D2 - (T - E[2]))          # follow-up of the THIRD cohort
            b_old = max(d, D2 - (T - E[-1]))         # follow-up of the LAST cohort
            n += 1
            if nu > b_new: viol_new += 1; print("VIOLATION(new)", E, T, nu, b_new)
            if nu > b_old: viol_old += 1
            if nu < b_new: strict += 1
print(f"designs: {n}; violations of nu <= max(d, D2 - (T - e_3)): {viol_new}; of the last-cohort version: {viol_old}; strict: {strict}")
for E, T in [((1, 5, 9, 13), 15), ((1, 3, 7, 11), 12), ((1, 6, 7, 12), 13)]:
    d = reduce(gcd, [e - E[0] for e in E[1:]]); D2 = E[1] - E[0]
    print(f"  E={E} T={T}: nu={trap_nullity(list(E), T)}  bound_e3={max(d, D2 - (T - E[2]))}  bound_eK={max(d, D2 - (T - E[-1]))}  d={d}")

print()
print("=" * 100); print("(C) robustness of (A): recruitment window, truncation, and a disconnected incidence graph\n")
for N in (4, 6, 8, 9, 10, 11, 12, 16, 24):
    r = analyse(f"CPS 4-8-4, {N} entry cohorts", list(range(1, N + 1)), [0, 1, 2, 3, 12, 13, 14, 15], N + 15, verbose=False)
    print(f"   {N:2d} cohorts: nullity={r['nullity']}  kernel ∝ j(k)? {r['prop_to_j']}  dose-linear detectable? {not r['dose_linear_in_kernel']}")
print()
print("=" * 100); print("(D) event-study sensitivity: shift(k) = tau(k+delta) - ell(k) as a combination of second differences\n")
from fractions import Fraction as FR
def shift_weights(ks, k0, k1, delta, smin, smax):
    """lambda_k over tenures s in [smin, smax] with shift(k) = sum_s lambda_s tau(s); then express lambda in the
    basis of second differences Delta^2 tau(y) = tau(y+2) - 2 tau(y+1) + tau(y)."""
    S = list(range(smin, smax + 1)); pos = {s: i for i, s in enumerate(S)}
    out = {}
    for k in ks:
        lam = [FR(0)] * len(S)
        lam[pos[k + delta]] += 1
        # ell(k) = tau(k0+delta) + (k - k0) * [tau(k1+delta) - tau(k0+delta)] / (k1 - k0)
        w1 = FR(k - k0, k1 - k0); w0 = 1 - w1
        lam[pos[k0 + delta]] -= w0; lam[pos[k1 + delta]] -= w1
        # solve lam = sum_y c_y * D2_y  (y from smin to smax-2)
        Y = list(range(smin, smax - 1)); M = [[FR(0)] * len(Y) for _ in S]
        for jj, y in enumerate(Y):
            M[pos[y]][jj] += 1; M[pos[y + 1]][jj] -= 2; M[pos[y + 2]][jj] += 1
        # least-norm exact solve by Gaussian elimination on the (consistent) system
        import copy
        Aug = [row[:] + [lam[i]] for i, row in enumerate(M)]
        r0 = 0; piv = []
        for c in range(len(Y)):
            p = next((i for i in range(r0, len(Aug)) if Aug[i][c]), None)
            if p is None: continue
            Aug[r0], Aug[p] = Aug[p], Aug[r0]; v = Aug[r0][c]; Aug[r0] = [x / v for x in Aug[r0]]
            for i in range(len(Aug)):
                if i != r0 and Aug[i][c]:
                    vv = Aug[i][c]; Aug[i] = [x - vv * y for x, y in zip(Aug[i], Aug[r0])]
            piv.append(c); r0 += 1
        c = [FR(0)] * len(Y)
        for i, cc in enumerate(piv): c[cc] = Aug[i][-1]
        assert all(not any(Aug[i][:len(Y)]) or True for i in range(len(Aug)))
        out[k] = (c, sum(abs(x) for x in c))
    return out, list(range(smin, smax - 1))
res, Ys = shift_weights([-5, -4, -3, 0, 3, 6], -1, -2, 6, 1, 13)
print("   nodes k0=-1, k1=-2, delta=6 (Simulation 4's configuration); weights on Delta^2 tau(y), y = tenure")
for k, (c, tot) in res.items():
    nz = {Ys[i]: str(x) for i, x in enumerate(c) if x}
    print(f"   k={k:+d}:  sum |weights| = {tot}   weights {nz}")
import math
tau = lambda s: 0.30 * (1 - math.exp(-0.6 * (s - 1)))
C = max(abs(tau(y + 2) - 2 * tau(y + 1) + tau(y)) for y in range(1, 12))
print(f"\n   With tau(s) = .30(1 - e^(-.6(s-1))): max |Delta^2 tau| over tenures 1..13 is C = {C:.4f}")
for k, (c, tot) in res.items():
    actual = tau(k + 6) - (tau(-1 + 6) + (k + 1) * (tau(-2 + 6) - tau(-1 + 6)) / (-2 + 1))
    print(f"   k={k:+d}: |shift| = {abs(actual):.4f}  <=  C * sum|w| = {C * float(tot):.4f}   (ratio {abs(actual) / (C * float(tot)):.2f})")

print()
print("=" * 100); print("(E) published CPS month-in-sample indices: the identified per-month rate profile\n")
import math
J_CPS = [0, 1, 2, 3, 12, 13, 14, 15]
def rate_profile(name, vals, log=True, unit=""):
    """vals = MIS 1..8 published index/level. Identified object: the per-month rates, up to a common constant."""
    y = [math.log(v) for v in vals] if log else list(vals)
    r = [(y[k + 1] - y[k]) / (J_CPS[k + 1] - J_CPS[k]) for k in range(7)]
    mean_r = sum(r) / 7
    dev = [x - mean_r for x in r]
    D = (y[4] - y[3]) - 9 * (y[1] - y[0])
    print(f"   {name} ({'log ' if log else ''}{unit})")
    print(f"      MIS values      : {[round(v, 4) for v in vals]}")
    print(f"      per-month rates : {[round(x, 4) for x in r]}   (interval 4->5 spans 9 months)")
    print(f"      rates, demeaned : {[round(x, 4) for x in dev]}   <- identified up to nothing further")
    print(f"      gap rate minus the largest in-block rate: {round(r[3] - max(r[:3] + r[4:]), 4)}")
    print(f"      D = [tau(5)-tau(4)] - 9[tau(2)-tau(1)] = {round(D, 4)}   (0 under a calendar-time clock; -8c under a dose clock)")
    return r, D
# McIllece (2022), Table 2: average multiplicative bias by MIS, Jan 2003 - Jun 2022, relative to second-stage estimates
rate_profile("McIllece 2022 Tab. 2, unemployed", [1.113, 1.053, 1.014, 0.989, 1.000, 0.955, 0.939, 0.938], unit="multiplicative bias")
rate_profile("McIllece 2022 Tab. 2, employed", [1.006, 1.006, 1.005, 1.000, 1.002, 0.992, 0.993, 0.996], unit="multiplicative bias")
# Bailar (1975), Table 1: rotation group indices (x100), 1968-69 (T1) and 1970-72 (T2)
rate_profile("Bailar 1975 Tab. 1 (1968-69), hours worked 35-40", [92.9, 97.9, 100.1, 101.7, 98.7, 101.9, 103.0, 103.8], unit="index")
rate_profile("Bailar 1975 Tab. 1 (1970-72), hours worked 35-40", [93.1, 97.7, 99.9, 101.5, 99.2, 101.7, 103.1, 103.9], unit="index")
rate_profile("Bailar 1975 Tab. 1 (1968-69), unemployed [MIS 1,5 carried extra questions]", [120.0, 101.5, 96.4, 92.8, 109.3, 96.5, 92.6, 91.0], unit="index")
# Solon (1986), Table 1: average ratio-estimate contributions by rotation group, Jan 1974 - Jun 1983 (thousands)
rate_profile("Solon 1986 Tab. 1, unemployed (thousands) [MIS 4,8 carried extra questions]", [1028, 952, 939, 964, 958, 913, 900, 941], log=False, unit="thousands")
print("\n   Under the calendar clock every rate is the same; under a dose clock the three in-block rates within a wave")
print("   are equal and the 4->5 rate is one ninth of them.  Both are testable from published tables alone.")

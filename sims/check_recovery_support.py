#!/usr/bin/env python3
"""T1 v1.0 — support conditions behind the recovery results (Section 5), the reference-category premise of
Theorem 5, the aliasing caveat of Theorem 4, and the population widths of Simulation 3's design G4.
Exact rational arithmetic (fractions.Fraction) throughout this file; standard library only; no data, no
network, no other scripts.  Run:  python3 check_recovery_support.py     (exit status 1 on any failure)

(A) Full kernel versus tenure projection.  Condition C_d characterises the TENURE PROJECTION K_tau of the
    kernel of the cell design.  When the cohort–period incidence graph is disconnected the full kernel is
    larger by one dimension per extra component:  dim K = dim K_tau + (c - 1).  The six-cell support
    E = {1, 4, 6}, cells (e, e) and (e, e+1), has d = 1, a single observed increment (C_1 holds), incidence
    components c = 3, full nullity 3 and projected nullity 1.  On it, exact knowledge of the cohort effects g
    does NOT pin tau(2): the vector h_tau(2) = a, h_alpha(e+1) = -a, h_g = 0 preserves every cell mean and both
    normalisations.  Propositions 2 and 3 therefore assume a connected incidence graph (condition (CG) in
    v0.9), which every staggered trapezoid satisfies; Proposition 1 does not, because it uses only K_tau.  (CG) is
    sufficient, not necessary: on cohorts {1, 2, 5} with the same two-cell schedule the graph has two components,
    yet every kernel vector moves g and knowledge of g pins tau(2) -- block (A) checks this too.  The same cohorts on the common-end trapezoid
    (T = 7) are connected and knowledge of g pins the whole path.  A sweep over trapezoids with up to two cells
    removed verifies dim K = dim K_tau + (c - 1) on every design.
(B) Theorem 5's parenthetical premise.  With homogeneous step effects 1{k >= 0}, zero untreated outcomes, no
    anticipation and full rank, omitting the categories k0 = -1 and k1 = +1 (a POST-event reference whose
    effect is nonzero) gives a pre-period coefficient of +1/2 at k = -2; omitting -1 and -2 gives zero pre-period
    coefficients.  The measurement-shift identity is unaffected; the sufficient condition for zero unconditioned
    pre-period coefficients must require the effect to vanish at both reference categories.
(C) Theorem 4's condition (R).  When D = 1{s >= 4} on a staggered panel, D lies in the span of the fixed
    effects and the retained tenure indicators; the joint regression has no unique coefficient on D, and two
    exact fits with different coefficients on D are exhibited.  Which column a software package reports as
    aliased depends on column order, so (R) has to be checked by a rank computation.
(D) Corollary 1's residue classes for d = 2, and the centred plateau boundary of Proposition 1.
(E) Simulation 3, design G4: population width of the identified set for m at M = M0 with both true drifts
    equal to M0 is 2 M0 - (max d - min d) = 2 M0 = .050, and the widths for tau(13) and tau(17) are .600 and .800.
(G) Proposition 3 without (CG) (round-5 review, N1).  The identified set for tau(2) under the drift bound M is
    computed EXACTLY by Fourier-Motzkin elimination over the kernel coefficients, and compared with the interval
    Proposition 3 displays for a given representative.  On cohorts {1, 2, 5}, two-cell schedule, zero cell means,
    zero representative: identified set [-M, M] = displayed interval, although the incidence graph has two
    components -- (CG) is sufficient for sharpness, not necessary.  Same support, representative with
    g_hat(5) = 2.7 M: displayed [-M, 0.1 M], identified set still [-M, M] -- without (CG) the displayed interval
    depends on the representative and can be strictly narrower than the identified set.  Six-cell support
    E = {1, 4, 6}: displayed [-M, M], identified set unbounded for every M, including M = 0.  A connected
    trapezoid is included as a control: displayed = identified.
"""
from fractions import Fraction as F
import itertools, sys

# ----------------------------------------------------------------------------------------------------------
def rank_exact(rows):
    """rank over the rationals by Gaussian elimination; rows = list of lists"""
    A = [[F(x) for x in r] for r in rows]
    if not A: return 0
    m, n = len(A), len(A[0]); r = 0
    for c in range(n):
        p = next((i for i in range(r, m) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][c]; A[r] = [x / pv for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        r += 1
        if r == m: break
    return r

def solve_exact(rows, rhs):
    """one exact solution of A x = b (A full column rank assumed for uniqueness; consistent system assumed)"""
    A = [[F(x) for x in r] + [F(b)] for r, b in zip(rows, rhs)]
    m, n = len(A), len(A[0]) - 1; r = 0; piv = []
    for c in range(n):
        p = next((i for i in range(r, m) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][c]; A[r] = [x / pv for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
        if r == m: break
    for i in range(r, m):
        assert A[i][n] == 0, "inconsistent system"
    x = [F(0)] * n
    for i, c in enumerate(piv): x[c] = A[i][n]
    return x, len(piv) == n

def cell_design(cells):
    """cells = set of (e, t); columns alpha(t), g(e) e != e0, tau(s) s != 1; returns (rows, col labels)"""
    E = sorted({e for e, _ in cells}); Ts = sorted({t for _, t in cells}); S = sorted({t - e + 1 for e, t in cells})
    e0 = E[0]
    cols = [('a', t) for t in Ts] + [('g', e) for e in E[1:]] + [('u', s) for s in S if s != 1]
    idx = {c: i for i, c in enumerate(cols)}
    rows = []
    for e, t in sorted(cells):
        r = [0] * len(cols); r[idx[('a', t)]] = 1
        if e != e0: r[idx[('g', e)]] = 1
        if t - e + 1 != 1: r[idx[('u', t - e + 1)]] = 1
        rows.append(r)
    return rows, cols

def nullities(cells):
    rows, cols = cell_design(cells)
    n = len(cols); rk = rank_exact(rows); full = n - rk
    nuis = [i for i, c in enumerate(cols) if c[0] != 'u']
    rk_nuis = rank_exact([[r[i] for i in nuis] for r in rows])
    proj = (n - len(nuis)) - rk + rk_nuis          # dim ker[W Z] - dim ker W
    return full, proj, rk, n

def incidence_components(cells):
    """connected components of the bipartite cohort–period graph"""
    nodes = {('e', e) for e, _ in cells} | {('t', t) for _, t in cells}
    adj = {v: set() for v in nodes}
    for e, t in cells: adj[('e', e)].add(('t', t)); adj[('t', t)].add(('e', e))
    seen, c = set(), 0
    for v in nodes:
        if v in seen: continue
        c += 1; stack = [v]
        while stack:
            u = stack.pop()
            if u in seen: continue
            seen.add(u); stack.extend(adj[u] - seen)
    return c

def increment_graph(cells):
    """(observed increments, components of the increment graph) as in Section 3"""
    U = sorted({t - e + 1 for e, t in cells if (e, t + 1) in cells})
    parent = {u: u for u in U}
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for (e, t) in cells:
        for (e2, t2) in cells:
            if t2 == t and e2 != e and (e, t + 1) in cells and (e2, t + 1) in cells:
                a, b = find(t - e + 1), find(t - e2 + 1)
                if a != b: parent[a] = b
    return U, len({find(u) for u in U})

def trapezoid(E, T): return {(e, t) for e in E for t in range(e, T + 1)}

fails = 0
def check(cond, msg):
    global fails
    print(("   ok   " if cond else "   FAIL ") + msg)
    if not cond: fails += 1

# ----------------------------------------------------------------------------------------------------------
print("=" * 100); print("(A) full kernel vs tenure projection; connectedness in Propositions 2-3\n")
six = {(e, e) for e in (1, 4, 6)} | {(e, e + 1) for e in (1, 4, 6)}
full, proj, rk, n = nullities(six)
U, ncomp = increment_graph(six); c = incidence_components(six)
print(f"   six-cell support E={{1,4,6}}, cells (e,e),(e,e+1): {len(six)} cells, {n} columns, rank {rk}, full nullity {full}, projected nullity {proj}")
print(f"   observed increments {U}, increment-graph components {ncomp}, incidence-graph components {c}")
check((full, proj, rk, n) == (3, 1, 6, 9), "rank 6 of 9, full nullity 3, projected nullity 1 (the six-cell example)")
check(U == [1] and ncomp == 1, "C_1 holds: single observed increment, one component")
check(c == 3, "incidence graph has three components")
check(full == proj + c - 1, "dim K = dim K_tau + (c - 1)")
# g known: drop the g columns (fix them at their true values) and ask whether tau(2) is still free
rows, cols = cell_design(six)
keep = [i for i, cc in enumerate(cols) if cc[0] != 'g']
rk_g = rank_exact([[r[i] for i in keep] for r in rows])
check(len(keep) - rk_g == 1, "with g known exactly, the (alpha, tau) design still has nullity 1: tau(2) is not pinned")
# the explicit kernel vector
a = F(7, 3)
h = {('a', t): (-a if t in (2, 5, 7) else F(0)) for _, t in six}; h.update({('g', 4): F(0), ('g', 6): F(0), ('u', 2): a})
check(all(sum(h.get(cc, F(0)) * r[i] for i, cc in enumerate(cols)) == 0 for r in rows), "h_tau(2)=a, h_alpha(e+1)=-a, h_g=0 preserves every cell mean")
# (CG) is sufficient, not necessary: cohorts {1,2,5} on the same schedule -- two components, but g known pins tau(2)
five = {(e, e) for e in (1, 2, 5)} | {(e, e + 1) for e in (1, 2, 5)}
full, proj, rk, n = nullities(five); c = incidence_components(five)
rows, cols = cell_design(five); keep = [i for i, cc in enumerate(cols) if cc[0] != 'g']
print(f"   cohorts {{1,2,5}}, same schedule: full nullity {full}, projected {proj}, incidence components {c}; nullity with g known {len(keep) - rank_exact([[r[i] for i in keep] for r in rows])}")
check(c == 2 and full == 2 and proj == 1 and len(keep) - rank_exact([[r[i] for i in keep] for r in rows]) == 0,
      "disconnected (c = 2) yet g known pins tau(2): (CG) is sufficient for Proposition 2, not necessary")
# the same cohorts on the common-end trapezoid are connected and g pins the path
tz = trapezoid([1, 4, 6], 7)
full, proj, rk, n = nullities(tz); c = incidence_components(tz)
rows, cols = cell_design(tz); keep = [i for i, cc in enumerate(cols) if cc[0] != 'g']
check(c == 1 and full == proj, "trapezoid E={1,4,6}, T=7: connected, full nullity = projected nullity")
check(len(keep) - rank_exact([[r[i] for i in keep] for r in rows]) == 0, "on the trapezoid, g known pins tau entirely (d = 1)")
# sweep: trapezoids with <= 2 cells removed, E subset of {1..6}, K in {2,3}, T <= 8
ndes = nviol = 0; ndisc = 0
for K in (2, 3):
    for E in itertools.combinations(range(1, 7), K):
        for T in range(E[-1], 9):
            base = trapezoid(list(E), T); cl = sorted(base)
            for k in (0, 1, 2):
                for rem in itertools.combinations(cl, k):
                    S = base - set(rem)
                    if any(e not in {e2 for e2, _ in S} for e in E): continue      # keep every cohort present
                    if not any(t - e + 1 == 1 for e, t in S): continue           # tenure 1 must be observed for tau(1) = 0 to bind
                    full, proj, rk, n = nullities(S); c = incidence_components(S)
                    ndes += 1; ndisc += (c > 1)
                    if full != proj + c - 1: nviol += 1
print(f"   sweep: {ndes} supports (trapezoids with up to two cells removed; {ndisc} disconnected): {nviol} violations of dim K = dim K_tau + (c - 1)")
check(nviol == 0 and ndisc > 0, "identity holds on every support in the sweep, including disconnected ones")

# ----------------------------------------------------------------------------------------------------------
print(); print("=" * 100); print("(B) Theorem 5: the reference categories must carry zero effect\n")
def event_study(refs, effect):
    events = list(range(3, 8)); window = list(range(-2, 3))
    obs = [(i, e + k, k) for i, e in enumerate(events) for k in window]
    periods = sorted({t for _, t, _ in obs}); kept = [k for k in window if k not in refs]
    X = [[int(i == u) for u in range(len(events))] + [int(t == p) for p in periods[1:]] + [int(k == kk) for kk in kept] for i, t, k in obs]
    y = [effect(k) for _, _, k in obs]
    rk = rank_exact(X); x, uniq = solve_exact(X, y)
    coef = {kk: x[len(events) + len(periods) - 1 + j] for j, kk in enumerate(kept)}
    return rk, len(X[0]), uniq, coef
step = lambda k: 1 if k >= 0 else 0
rk, ncol, uniq, coef = event_study([-1, 1], step)
print(f"   events 3..7, window -2..2, effect 1{{k>=0}}, references {{-1,+1}}: rank {rk} of {ncol}; coefficients {dict((k, str(v)) for k, v in coef.items())}")
check(rk == ncol == 16 and uniq, "full column rank 16 of 16")
check(coef[-2] == F(1, 2), "pre-period coefficient at k=-2 equals +1/2 although effects are homogeneous and untreated outcomes are zero")
rk, ncol, uniq, coef = event_study([-1, -2], step)
print(f"   same design, references {{-1,-2}}: rank {rk} of {ncol}; coefficients {dict((k, str(v)) for k, v in coef.items())}")
check(rk == ncol and all(coef[k] == 0 for k in coef if k < 0) and all(coef[k] == 1 for k in coef if k >= 0), "with both references pre-event the pre-period coefficients are zero and the post-period ones equal the effect")

# ----------------------------------------------------------------------------------------------------------
print(); print("=" * 100); print("(C) Theorem 4, condition (R): tenure-determined exposure has no unique coefficient\n")
E = [1, 5, 13]; T = 19
units = [(e, u) for e in E for u in range(2)]                       # two units per cohort suffice for the algebra
obs = [(i, t) for i, (e, _) in enumerate(units) for t in range(e, T + 1)]
periods = sorted({t for _, t in obs}); S = sorted({t - units[i][0] + 1 for i, t in obs})
def design(order_D_first):
    W = [[int(i == u) for u in range(len(units))] + [int(t == p) for p in periods[1:]] for i, t in obs]
    Z = [[int(t - units[i][0] + 1 == s) for s in S if s != 1] for i, t in obs]
    D = [[int(t - units[i][0] + 1 >= 4)] for i, t in obs]
    return [w + (d + z if order_D_first else z + d) for w, z, d in zip(W, Z, D)]
Xd = design(True); rkW = rank_exact([r[:len(units) + len(periods) - 1] for r in Xd])
rk_all = rank_exact(Xd); rk_noD = rank_exact([r[:len(units) + len(periods) - 1] + r[len(units) + len(periods):] for r in Xd])
print(f"   staggered panel E={E}, T={T}: rank[W, D, Z_tau] = {rk_all}, rank[W, Z_tau] = {rk_noD} -> D adds nothing")
check(rk_all == rk_noD, "D = 1{s >= 4} lies in span[W, Z_tau]: condition (R) fails")
# two exact fits of the same outcome with different coefficients on D
y = [F(t - units[i][0] + 1 >= 4) for i, t in obs]                   # outcome = D itself
beta_candidates = []
for beta in (F(0), F(1)):
    resid = [yy - beta * r[len(units) + len(periods) - 1] for yy, r in zip(y, Xd)]
    x, _ = solve_exact([r[:len(units) + len(periods) - 1] + r[len(units) + len(periods):] for r in Xd], resid)
    beta_candidates.append(beta)
check(len(beta_candidates) == 2, "exact fits exist with beta = 0 and with beta = 1: the coefficient on D is not determined by the data")

# ----------------------------------------------------------------------------------------------------------
print(); print("=" * 100); print("(D) Corollary 1 at d = 2; plateau boundary\n")
tz = trapezoid([1, 3, 5], 9)
rows, cols = cell_design(tz); S2 = sorted({t - e + 1 for e, t in tz})
# tau-kernel basis: exact null space restricted to tau coordinates via solving; use the characterisation and test functionals
def identified(lam):
    """lambda' tau identified iff lambda lies in the row space of the design restricted to tau columns after
    projecting out nuisance columns: equivalently rank([rows; lam_row]) == rank(rows) with lam placed in tau cols"""
    row = [0] * len(cols)
    for s, w in lam.items():
        if s != 1: row[cols.index(('u', s))] = w
    return rank_exact(rows + [row]) == rank_exact(rows)
check(not identified({3: 1, 4: -2, 5: 1}), "d = 2: ordinary second difference Delta^2 tau(4) is NOT identified (weights in two residue classes)")
check(identified({2: 1, 4: -2, 6: 1}), "d = 2: centred lag-2 second difference Delta_2^2 tau(4) = tau(6) - 2 tau(4) + tau(2) IS identified")
check(not identified({4: 1, 2: -1}), "d = 2: the increment tau(4) - tau(2) is not identified")
tau = {1: F(0), 2: F(1), 3: F(1), 4: F(1)}
check(tau[3] - 2 * tau[2] + tau[1] == -1 and tau[4] - 2 * tau[3] + tau[2] == 0,
      "plateau from S* = 2 with d = 1: Delta^2 tau(2) = -1 (centre at S*), Delta^2 tau(3) = 0 (all three points on the plateau)")

# ----------------------------------------------------------------------------------------------------------
print(); print("=" * 100); print("(E) Simulation 3, design G4: population widths at M = M0\n")
M0 = F(1, 40); drifts = [M0, M0]
w_m = 2 * M0 - (max(drifts) - min(drifts))
print(f"   drifts (+M0, +M0): width of the set for m = 2 M0 - (max - min) = {w_m} = {float(w_m):.3f}; tau(13): {float(12 * w_m):.3f}; tau(17): {float(16 * w_m):.3f}")
check(w_m == F(1, 20) and 12 * w_m == F(3, 5) and 16 * w_m == F(4, 5), "widths .050, .600, .800")
drifts3 = [M0, -M0]; w3 = 2 * M0 - (max(drifts3) - min(drifts3))
check(w3 == 0, "design G3 (+M0, -M0): the set for m is a point at M = M0")

# ----------------------------------------------------------------------------------------------------------
print(); print("=" * 100); print("(F) Table 1, uninterrupted rows, by exact rank (the interrupted rows are in check_interrupted.py / check_dose_sharpness.py)\n")
for name, E, T, want in [("adjacent entry, long follow-up: E={1,2,3}, T=8", [1, 2, 3], 8, 1),
                         ("entry every d=3 waves, long follow-up: E={1,4,7}, T=12", [1, 4, 7], 12, 3),
                         ("JLPS: E={1,5,13}, T=19", [1, 5, 13], 19, 4),
                         ("E={1,4,8}, T=8 (no follow-up of the last cohort)", [1, 4, 8], 8, 3),
                         ("single refreshment at wave 5: E={1,5}, T=12", [1, 5], 12, 4)]:
    tz = trapezoid(E, T); full, proj, rk, n = nullities(tz); c = incidence_components(tz)
    print(f"   {name}: {len(tz)} cells, {n} parameters, rank {rk}, nullity {full} (projected {proj}, components {c})")
    check(full == proj == want and c == 1, f"dim Theta_tau = {want}")
tz = trapezoid([1, 5, 13], 19); full, proj, rk, n = nullities(tz)
check((len(tz), n, rk) == (41, 39, 35), "JLPS design: 41 cells, 39 free parameters, rank 35 (as quoted in Sections 3 and 8 and Appendix C)")
tz = trapezoid([1, 4, 8], 8); full, proj, rk, n = nullities(tz)
check((len(tz), n, rk) == (14, 17, 14), "E={1,4,8}, T=8: 14 cells, 17 parameters, rank 14 (as quoted after Theorem 1)")

# ----------------------------------------------------------------------------------------------------------
print(); print("=" * 100); print("(G) Proposition 3 without (CG): identified set for tau(2) under |drift| <= M, by exact Fourier-Motzkin\n")
def kernel_basis(cells):
    """rational basis of ker(cell design); returns (basis vectors, cols)"""
    rows, cols = cell_design(cells); n = len(cols)
    A = [[F(x) for x in r] for r in rows]; piv = []; rr = 0
    for c in range(n):
        pr = next((i for i in range(rr, len(A)) if A[i][c] != 0), None)
        if pr is None: continue
        A[rr], A[pr] = A[pr], A[rr]; pv = A[rr][c]; A[rr] = [x / pv for x in A[rr]]
        for i in range(len(A)):
            if i != rr and A[i][c] != 0:
                f = A[i][c]; A[i] = [a - f * b for a, b in zip(A[i], A[rr])]
        piv.append(c); rr += 1
    basis = []
    for fc in [c for c in range(n) if c not in piv]:
        v = [F(0)] * n; v[fc] = F(1)
        for i, c in enumerate(piv): v[c] = -A[i][fc]
        basis.append(v)
    return basis, cols

def fm_bounds(ineqs, nvar, keep):
    """ineqs: list of (coef list, rhs) meaning coef.x <= rhs; eliminate every variable except `keep`;
    returns (lo, hi) for x[keep], None meaning unbounded on that side"""
    cur = [(list(c), r) for c, r in ineqs]
    for v in range(nvar):
        if v == keep: continue
        pos = [(c, r) for c, r in cur if c[v] > 0]; neg = [(c, r) for c, r in cur if c[v] < 0]
        new = [(c, r) for c, r in cur if c[v] == 0]
        for cp, rp in pos:
            for cn, rn in neg:
                a, b = cp[v], -cn[v]
                new.append(([b * x + a * y for x, y in zip(cp, cn)], b * rp + a * rn))
        cur = new
    lo = hi = None
    for c, r in cur:
        a = c[keep]
        if a > 0: hi = r / a if hi is None else min(hi, r / a)
        elif a < 0: lo = r / a if lo is None else max(lo, r / a)
        elif r < 0: raise ValueError("infeasible")
    return lo, hi

def identified_set_tau2(cells, M, rep_g):
    """identified set for tau(2) = rep_tau(2) + h_tau(2) over all kernel vectors h with |drift(rep + h)| <= M;
    rep_g maps e -> representative's g_hat(e) (default 0); rep_tau(2) = 0 throughout"""
    basis, cols = kernel_basis(cells); k = len(basis); idx = {c: i for i, c in enumerate(cols)}
    E = sorted({e for e, _ in cells}); e0 = E[0]
    def gvec(e):
        if e == e0: return [F(0)] * k, F(0)
        return [b[idx[('g', e)]] for b in basis], F(rep_g.get(e, 0))
    ineqs = []
    for ea, eb in zip(E, E[1:]):
        va, ca = gvec(ea); vb, cb = gvec(eb); sp = eb - ea
        d = [(y - x) / sp for x, y in zip(va, vb)]; cd = (cb - ca) / sp
        ineqs.append((d + [F(0)], M - cd)); ineqs.append(([-x for x in d] + [F(0)], M + cd))
    tv = [b[idx[('u', 2)]] for b in basis]
    ineqs.append((tv + [F(-1)], F(0))); ineqs.append(([-x for x in tv] + [F(1)], F(0)))   # z = tau(2)
    return fm_bounds(ineqs, k + 1, k)

def displayed_interval(cells, M, rep_g):
    """the interval Proposition 3 displays for tau(2) - tau(1) at the representative: [m_lo, m_hi] x (s - 1)"""
    E = sorted({e for e, _ in cells})
    dh = [(F(rep_g.get(eb, 0)) - F(rep_g.get(ea, 0))) / (eb - ea) for ea, eb in zip(E, E[1:])]
    return -M - min(dh), M - max(dh)

M = F(1)
five = {(e, e) for e in (1, 2, 5)} | {(e, e + 1) for e in (1, 2, 5)}
lo, hi = identified_set_tau2(five, M, {}); dlo, dhi = displayed_interval(five, M, {})
print(f"   cohorts {{1,2,5}}, two-cell schedule, zero means, zero representative: identified set [{lo}, {hi}] M, displayed [{dlo}, {dhi}] M")
check((lo, hi) == (-M, M) and (dlo, dhi) == (-M, M) and incidence_components(five) == 2,
      "two components, yet the displayed interval [-M, M] IS the identified set: (CG) is sufficient for Proposition 3, not necessary")
lo, hi = identified_set_tau2(five, M, {5: F(27, 10)}); dlo, dhi = displayed_interval(five, M, {5: F(27, 10)})
print(f"   same support, representative g_hat(5) = 2.7 M: identified set [{lo}, {hi}] M, displayed [{dlo}, {dhi}] M")
check((lo, hi) == (-M, M) and (dlo, dhi) == (-M, F(1, 10)),
      "without (CG) the displayed interval depends on the representative and here is strictly narrower than the identified set")
for Mv in (F(1), F(0)):
    lo, hi = identified_set_tau2(six, Mv, {}); dlo, dhi = displayed_interval(six, Mv, {})
    print(f"   six-cell E={{1,4,6}}, M = {Mv}: identified set [{lo}, {hi}] (None = unbounded), displayed [{dlo}, {dhi}]")
    check(lo is None and hi is None, f"tau(2) unbounded at M = {Mv} on the six-cell support while the displayed interval is [{dlo}, {dhi}]")
tz = trapezoid([1, 2, 3], 5)
lo, hi = identified_set_tau2(tz, M, {}); dlo, dhi = displayed_interval(tz, M, {})
check((lo, hi) == (dlo, dhi) == (-M, M) and incidence_components(tz) == 1, "control: connected trapezoid E={1,2,3}, T=5 -- displayed interval = identified set = [-M, M]")

print(); print("=" * 100)
print(f"check_recovery_support.py: {fails} failure(s)")
sys.exit(1 if fails else 0)

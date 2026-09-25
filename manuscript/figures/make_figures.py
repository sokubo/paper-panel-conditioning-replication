#!/usr/bin/env python3
"""Figures 1 and 2 of "Panel Conditioning in Fixed-Effects Models: Identification and Bias Propagation".
Data-free; the PNG files regenerate byte-identically (the PDF files differ only in their creation date): every plotted quantity is design geometry stated in the text (the Japanese panel's
entries at waves 1, 5 and 13, observed through wave 19; the CPS 4-8-4 offsets 0, 1, 2, 3, 12, 13, 14, 15).
Run:  python3 make_figures.py   (Python 3 with NumPy and Matplotlib) -> fig1_jlps_support.{pdf,png},
fig2_cps_rates.{pdf,png} next to this script. Before drawing, the script checks the quantities the text states
(41 cells; 39 free parameters, rank 35, kernel dimension 4; tenure class = period class in every cell;
D = 0 for a constant per-month path and D = -8c for a constant per-interview path) and exits non-zero if one fails."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"svg.hashsalt": "t1", "font.family": "serif", "font.serif": ["Latin Modern Roman", "DejaVu Serif"],
                     "mathtext.fontset": "cm", "font.size": 10, "axes.linewidth": .6,
                     "xtick.major.width": .6, "ytick.major.width": .6, "pdf.fonttype": 42})
fails = []
def check(what, ok):
    print(("ok    " if ok else "FAIL  ") + what)
    if not ok: fails.append(what)

# ---------------- Figure 1: JLPS support, tenure residue classes ----------------
E, T, d = [1, 5, 13], 19, 4
cells = [(e, t) for e in E for t in range(e, T + 1)]
check("JLPS support has 41 cells", len(cells) == 41)
# tenure residue equals period residue in every cell (the reason periodic directions are absorbed by alpha(t))
check("s mod 4 == (t - e1 + 1) mod 4 in every cell", all(((t - e + 1) % d) == ((t - E[0] + 1) % d) for e, t in cells))
# rank of the cell design: 41 cells; free params alpha(19) + g(2) + tau(18) = 39; rank 35; kernel dim 4
S = sorted(set(t - e + 1 for e, t in cells))
cols = [("a", t) for t in range(1, T + 1)] + [("g", e) for e in E[1:]] + [("s", s) for s in S if s != 1]
X = np.zeros((len(cells), len(cols)))
idx = {c: i for i, c in enumerate(cols)}
for r, (e, t) in enumerate(cells):
    X[r, idx[("a", t)]] = 1
    if e != E[0]: X[r, idx[("g", e)]] = 1
    s = t - e + 1
    if s != 1: X[r, idx[("s", s)]] = 1
rk = np.linalg.matrix_rank(X)
check("39 free parameters, rank 35, kernel dimension 4", len(cols) == 39 and rk == 35 and len(cols) - rk == 4)

shades = {1: "#ffffff", 2: "#d9d9d9", 3: "#969696", 0: "#525252"}
txtcol = {1: "black", 2: "black", 3: "white", 0: "white"}
fig, ax = plt.subplots(figsize=(6.5, 2.3))
for i, e in enumerate(E):
    y = len(E) - 1 - i
    for t in range(1, T + 1):
        if t < e:
            continue
        s = t - e + 1
        ax.add_patch(Rectangle((t - .5, y - .4), 1, .8, facecolor=shades[s % d], edgecolor="black", lw=.4))
        ax.text(t, y, str(s), ha="center", va="center", fontsize=9, color=txtcol[s % d])
ax.set_xlim(.4, T + .6); ax.set_ylim(-.6, len(E) - .4)
ax.set_xticks(range(1, T + 1)); ax.set_xlabel("Wave (calendar period $t$)")
ax.set_yticks(range(len(E))); ax.set_yticklabels(["entered wave %d" % e for e in reversed(E)])
for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
ax.tick_params(axis="y", length=0)
handles = [Rectangle((0, 0), 1, 1, facecolor=shades[r], edgecolor="black", lw=.4) for r in (1, 2, 3, 0)]
ax.legend(handles, [r"$s \equiv %d$" % r for r in (1, 2, 3, 0)], title=r"tenure $s$ mod 4", ncol=4, fontsize=9,
          title_fontsize=9, frameon=False, loc="lower center", bbox_to_anchor=(.5, 1.0), handlelength=1.2, columnspacing=1.2)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(HERE, "fig1_jlps_support." + ext), dpi=300, bbox_inches="tight")
plt.close(fig)

# ---------------- Figure 2: CPS 4-8-4, two constant-rate paths ----------------
j = np.array([0, 1, 2, 3, 12, 13, 14, 15]); k = np.arange(1, 9)
c = 1.0
per_interview = c * (k - 1)
per_month = (7.0 / 15.0) * j          # same total change over the eight interviews, for display
D = lambda tau: (tau[4] - tau[3]) - 9 * (tau[1] - tau[0])
check("D = 0 for a constant per-month path", abs(D(per_month)) < 1e-12)
check("D = -8c for a constant per-interview path", abs(D(per_interview) + 8 * c) < 1e-12)
fig, ax = plt.subplots(figsize=(4.6, 2.7))
ax.axvspan(3, 12, color="#efefef", lw=0)
ax.text(7.5, 7.2, "eight months\nwithout interview", ha="center", va="top", fontsize=9)
ax.plot(j, per_month, "o--", color="#737373", ms=4, lw=1, label="constant rate per month")
ax.plot(j, per_interview, "s-", color="black", ms=4, lw=1, label="constant rate per interview")
for kk, jj, yy in zip(k, j, per_interview):
    ax.annotate(str(kk), (jj, yy), textcoords="offset points", xytext=(-9, 3) if kk <= 4 else (4, -12), fontsize=9)
ax.set_xlabel("Months since first interview, $j(k)$"); ax.set_ylabel(r"Conditioning $\tau(k)$ (arbitrary units)")
ax.set_xticks([0, 1, 2, 3, 12, 13, 14, 15]); ax.set_ylim(-.4, 7.6)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
ax.legend(frameon=False, fontsize=9, loc="lower right")
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(HERE, "fig2_cps_rates." + ext), dpi=300, bbox_inches="tight")
plt.close(fig)

print("%d check(s) failed" % len(fails) if fails else "all figure checks passed")
sys.exit(1 if fails else 0)

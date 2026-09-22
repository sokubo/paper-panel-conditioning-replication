# Replication archive: Panel Conditioning in Fixed-Effects Models: Identification and Bias Propagation

Code and outputs reproducing every table, appendix number and design computation of the paper.
Preprint: arXiv:XXXX.XXXXX (to be filled at posting). Author: Shoki Okubo (Toyo University).

## Contents
- manuscript
- manuscript/latex
- sims

`FILE_MANIFEST.txt` lists every file in this snapshot except itself, the `.git` history and the
release-check records (`RELEASE_CHECK*`), which are written after the manifest.

## Checked commit and release record
Computational commit checked against the manuscript: `a3b36521410184b0ce27815b6f8edef9ba09aba1` — see `RELEASE_CHECK.md` (with `RELEASE_CHECK_run.log` and `RELEASE_CHECK_sessionInfo.txt` when a clean-copy run was made). Later commits change documentation and the release record only — `git diff --stat a3b36521410184b0ce27815b6f8edef9ba09aba1 HEAD` lists them — so the scripts and outputs are those of the checked commit; after any change to code or outputs the release check is rerun and this line is regenerated.
Tag matching this version of the manuscript: `paper-v0.9`. Tag matching the posted preprint version: to be added at posting (`arxiv-<id>v<n>`).
Third-party reproduction: none. The release record is the author's own re-execution of the published snapshot in a clean copy.

## Data
No data are distributed with this archive and none are needed. Every number in the paper comes from a seeded simulation that generates its own data, from a deterministic algebraic computation on a stated design, or from month-in-sample indices published in Bailar (1975, Tab. 1), Solon (1986, Tab. 1) and McIllece (2022, Tab. 2), which are reproduced inside `sims/check_interrupted.py`. The Japanese Life Course Panel Surveys (JLPS) microdata are **not** used anywhere in this paper: Section 8 is a design illustration whose only inputs are the entry waves 1, 5, 13 and the final wave 19.

## Environment

- R 4.3.3 (also re-run independently under R 4.6.0, arm64, with byte-identical results up to 1e-13); packages `data.table` 1.14.10, `sandwich` 3.1.0 (Simulation 4 only).
- Python 3 with `numpy` (support and design checks only; developed under 3.11 / numpy 2.4, and the release check of 21 September 2026 ran Python 3.14 / numpy 2.5.3 on macOS arm64 with the same integer ranks and the same sweep counts).
- Quarto 1.6.42 + XeLaTeX (TeX Live 2023) for the manuscript; `manuscript/build_latex.py` reproduces `main.pdf` and `latex/main.tex`.

## Execution order and runtimes (single core)

| Step | Script | Output | Runtime | Paper location |
|---|---|---|---|---|
| 1 | `sims/sim1_normalizations.R` | `sim1_results.csv` | ~1 min | Appendix A.1 (table; noiseless identities) |
| 2 | `sims/sim2_direct_contrast.R` | `sim2_results.csv` | ~1 min | Appendix A.2 |
| 3 | `sims/sim3_drift_bounds.R` (v4) | `sim3_results.csv` | ~18 min | Appendix A.3 (three DGPs G2, G3, G4; adaptive and fixed-z coverage) |
| 4 | `sims/sim4_downstream.R` (v3) | `sim4_results.csv`, `sim4_pretrend_coefs.csv`, `sim4_mc.csv`, `sim4_manuscript_values.csv`, `sim4_results.txt` | ~30 min | Appendix C (identities, alignment checks, Monte Carlo table) |
| 5 | `sims/check_manuscript_values.R` | console; exit status | seconds | **asserts** every number quoted in Appendices A and C against its named output column (71 assertions, including key uniqueness and row counts); each selector must return exactly one row, each vector must have its expected length and only finite values, and floating-point identities are asserted at the scale-aware tolerance 1e-10; exits 1 on any failure |
| 5b | `sims/check_manuscript_values.R --selftest` | console; exit status | seconds | negative controls: perturbed copies of the outputs (deleted row, duplicated row, altered value, non-finite value, deleted quantity, identity residual 1e-6) must each make the checker fail, and an unmodified copy and a residual of 5e-11 must pass |
| 6 | `sims/design_rank.R` | console | seconds | Theorem 1(d), Lemma 1(o) examples; the {1,4,8} design |
| 7 | `sims/check_thm1_support.py` | console | seconds | nullity = gcd on 14 C_d designs |
| 8 | `sims/check_lemma1_trapezoid.py` | console | ~10 s | Lemma 1 sweep (1,547 exact three-cohort designs; 9,401 designs for the bound) |
| 9 | `sims/check_lemma1_components.py` | console | ~10 s | Lemma 1(o) and (iv): exact rank for five named designs, then a sweep of 10,934 trapezoids |
| 10 | `sims/check_interrupted.py` | console | ~5 s | Theorem 2(a),(d), Table 1 (dose designs); Lemma 1(iv) sweep (10,934 designs; the 72 failures reported there are of the withdrawn third-cohort conjecture, not of the lemma); Proposition 5's weights and the constant-curvature equality test; §8.1's descriptive Table 2 (Appendix D) |
| 11 | `sims/check_dose_sharpness.py` | console | ~10 s | Theorem 2(b): the earlier counterexample (cohorts {1,2,4}, offsets {0,1,3}) and condition (P') verified by exact arithmetic over 496 designs; the CPS ten-cohort threshold; the two d > 1 examples cited in §3.4 (Appendix D) |
| 12 | `sims/check_recovery_support.py` | console; exit status | ~15 s | exact arithmetic: the six-cell disconnected support (C_1 holds, full nullity 3, projected nullity 1, tau(2) free with g known) and the identity dim K = dim K_tau + (c-1) on 9,183 supports (Theorem 1(d), condition (CG) of §5); the post-event-reference example for Theorem 5's premise (rank 16, pre-coefficient +1/2); the failure of (R) under D = 1{s>=4} with two exact fits; Corollary 1 at d = 2; the plateau boundary; the population widths of G4; the uninterrupted rows of Table 1 (Appendix D) |

Seeds are fixed inside each script (sim1–2: as in file; sim3: 20260830; sim4: 20260916). Scripts are run from `sims/` and write into it. The whole sequence takes about 50 minutes, of which Simulation 4's Monte Carlo (five designs, B = 100) is 30 and Simulation 3 is 18; steps 5–12 together take about a minute. Steps 6–12 are deterministic — no randomness and no data — so they either reproduce exactly or fail. Steps 11 and 12, the named designs of step 9, and blocks (A) and (D') of step 10 use exact rational arithmetic; steps 6, 7, 8 and the two 10,934-design sweeps (steps 9 and 10) use floating-point rank (numpy SVD at tolerance 1e-10, R's `qr`), which for these small 0–1 matrices returns the same integers as the exact code on every named design.

## How to read the release check

The release script (`release_check.sh` in the author's workflow, whose record is `RELEASE_CHECK.md`) exits with a
nonzero status — and says so in its last line — if any step of the sequence fails, if any file listed in
`FILE_MANIFEST.txt` is missing from the anonymous download or any unlisted file is present, if any restricted-data
or review-material pattern matches, if any regenerated numeric value differs from the shipped one by more than
`TOL` (default 1e-8), or if any output cannot be compared at all. A tag is placed only on a snapshot whose check
passed.


`RELEASE_CHECK.md` records a re-execution of the documented sequence in a clean copy of the published
snapshot, downloaded anonymously, with the shipped outputs set aside and the regenerated ones compared
token by token. Two kinds of difference are expected and are not errors.

*Floating point.* The shipped outputs were produced under R 4.3.3 on x86-64 Linux; the release check of
21 September 2026 ran R 4.6.0 on macOS arm64, with a different BLAS and LAPACK. Differences of order
1e-13 in regression coefficients and simulation summaries follow from that and from nothing else. The
deterministic identities the paper relies on are all checked *inside* the scripts against their own
analytic formulas, to tolerances stated in Appendix C, so a platform difference of this size cannot
disturb any claim; `sims/check_manuscript_values.R` asserts every number quoted in Appendices A and C
from its named output column, at a scale-aware tolerance of 1e-10 for the floating-point identities (the archived residuals of order 1e-14 to 1e-16 are platform-specific and are not themselves required to recur).

*Print format.* `sim4_results.txt` is a transcript of console output produced with `sink()`. It is now
written with base R's data-frame printer at a fixed, very wide console, so its layout does not depend on
the installed `data.table` version: 1.15 and later add a type-annotation row under each header, and the
column widths — hence where a wide table wraps into a second block — have changed between versions, either
of which makes a token-by-token comparison fail on another machine for no substantive reason. (The check of
22 September 2026 hit exactly that: under `data.table` 1.18.4 the ten-row event-study table wrapped into two
blocks, adding a header and ten row labels, while the shipped transcript, written under 1.14.10, printed it
in one. The transcript was regenerated in the stable format; it carries the same 124 numbers, and the only
tokens that disappeared are `data.table`'s row labels `1:` to `10:`.)

The comparator keeps two safety nets for this file. It drops type-annotation rows, and if a transcript's
token stream still differs it compares the **real numbers of the two streams in order** — every token
containing a decimal point or an exponent — requiring the same count and agreement within the release
tolerance, and prints the first divergence so that a genuine change is never hidden behind the word
"format". Such a file is reported as `format differs; numbers match`, with the largest numerical
difference; one whose real numbers differ in count, or by more than the tolerance, is a failure and stops
the release. Integers in the transcript (ranks, counts) are not compared that way, but each is also
produced by a deterministic script that fails loudly, and every number the paper quotes is asserted by
`check_manuscript_values.R` against a named column of a CSV, which is compared token by token.

The deterministic checks — steps 6 to 12 — carry no such caveat. They use no randomness, so they either
reproduce exactly or fail. In the 21 September 2026 check
they reproduced exactly, including Lemma 1's sweep over 10,934 designs (no violation of either part; the
bound strict in 1,188 of them) and every entry of Table 2, on a different operating system, a different
R, a different Python and a different `numpy` from those used to write the paper.

## Map from manuscript objects to code

- Theorem 1(a)–(d), Corollary 1: `design_rank.R` (`design_rank()`, `increment_graph()`, `is_identified()`), `check_thm1_support.py`; Theorem 1(d)'s full-kernel/projection identity, condition (CG) of §5 and the six-cell example, Corollary 1 at d = 2, Proposition 1's plateau boundary: `check_recovery_support.py`.
- Theorem 2(a),(d), Corollary 2, Table 1 (interrupted participation / dose designs) and Proposition 5's sensitivity weights and equality test: `check_interrupted.py`.
- Theorem 2(b), condition (P'), the counterexample to the earlier statement, and the d > 1 examples of §3.4: `check_dose_sharpness.py`.
- Lemma 1 (o)–(iv): `check_lemma1_components.py`, `check_lemma1_trapezoid.py`, `check_interrupted.py`, `design_rank.R`.
- Proposition 1 (plateau test), normalizations: `sim1_normalizations.R`.
- Corollary 1 readouts: `sim2_direct_contrast.R`.
- Proposition 3 (drift bounds, IM intervals): `sim3_drift_bounds.R`.
- Theorem 3 (absorption, kernel invariance), Theorem 4 (joint regression, (R)), Theorem 5 (pre-trend identity, rank-deficient design, exact-alignment rank), Proposition 4 (partial-alignment counterexample), Corollary 3(a): `sim4_downstream.R`; Theorem 4's aliasing caveat and Theorem 5's reference-category premise (exact examples): `check_recovery_support.py`.
- Simulation 3's population widths (G2, G3, G4): `check_manuscript_values.R` and `check_recovery_support.py` (E).

## Manuscript

`manuscript/main.qmd` (source), `references.bib`, `latex/preamble_extra_T1.tex`, `build_latex.py`. Render: `quarto render main.qmd --to html` and `python3 build_latex.py`.

## Citation
Okubo, S. (2026). Panel Conditioning in Fixed-Effects Models: Identification and Bias Propagation. Working paper. arXiv:XXXX.XXXXX.

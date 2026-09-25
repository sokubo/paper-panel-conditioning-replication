# Replication archive: Panel Conditioning in Fixed-Effects Models: Identification and Bias Propagation

Code and outputs reproducing every table, appendix number and design computation of the paper.
Preprint: arXiv:2609.28871. Author: Shoki Okubo (Toyo University).

## Contents
- manuscript
- manuscript/figures
- manuscript/latex
- sims

`FILE_MANIFEST.txt` lists every file in this snapshot except itself, the `.git` history and the
release-check records (`RELEASE_CHECK*`), which are written after the manifest.

## Checked commit and release record
Computational commit checked against the numerical results of the manuscript: not yet checked for this version. The `RELEASE_CHECK*` files in this snapshot are the record of the previously checked commit `f06283e3039a4693bbe94e71060715191b1d4d48`, whose scripts or outputs differ from these; the release check is rerun on this version before it is tagged, and this line is then regenerated.
Tag of the checked release: `paper-v1.1`. Tag matching the posted preprint version: `arxiv-2609.28871v1`.
Third-party reproduction: none. The release record is the author's own re-execution of the published snapshot in a clean copy.

## Data
No data are distributed with this archive and none are needed. Every number in the paper comes from a seeded simulation that generates its own data, from a deterministic algebraic computation on a stated design, or from month-in-sample indices published in Bailar (1975, Tab. 1), Solon (1986, Tab. 1) and McIllece (2022, Tab. 2), which are reproduced inside `sims/check_interrupted.py`. The Japanese Life Course Panel Surveys (JLPS) microdata are **not** used anywhere in this paper: Section 8 is a design illustration whose only inputs are the entry waves 1, 5, 13 and the final wave 19.

## Environment

- R 4.3.3 (also re-run under R 4.6.0, arm64, with results identical to 1e-13); packages `data.table` 1.14.10, `sandwich` 3.1.0 (Simulation 4 only).
- Python 3 with `numpy` (support and design checks only; developed under 3.11 / numpy 2.4, and the release check of 21 September 2026 ran Python 3.14 / numpy 2.5.3 on macOS arm64 with the same integer ranks and the same sweep counts).
- Quarto 1.6.42 (Pandoc 3.4) + XeLaTeX (TeX Live 2023) for the manuscript; `manuscript/build_latex.py` reproduces `main.pdf` and
  `latex/main.tex`. The PDF engine (XeLaTeX) and the citation style (`manuscript/latex/chicago-author-date.csl`, Pandoc's
  built-in Chicago author-date style) are pinned in `main.qmd`. The numerical scripts do not depend on this toolchain;
  exact pagination and citation formatting do, so a build with another Quarto or TeX Live version can differ in
  page breaks and reference layout without any change to a statement, number or formula.

## Execution order and runtimes (single core)

| Step | Script | Output | Runtime | Paper location |
|---|---|---|---|---|
| 1 | `sims/sim1_normalizations.R` | `sim1_results.csv` | ~1 min | Appendix A.1 (table; noiseless identities) |
| 2 | `sims/sim2_direct_contrast.R` | `sim2_results.csv` | ~1 min | Appendix A.2 |
| 3 | `sims/sim3_drift_bounds.R` (v4) | `sim3_results.csv` | ~18 min | Appendix A.3 (three DGPs G2, G3, G4; adaptive and fixed-z coverage) |
| 4 | `sims/sim4_downstream.R` (v3) | `sim4_results.csv`, `sim4_pretrend_coefs.csv`, `sim4_mc.csv`, `sim4_manuscript_values.csv`, `sim4_results.txt` | ~30 min | Appendix C (identities, alignment checks, Monte Carlo table) |
| 5 | `sims/check_manuscript_values.R` | console; exit status | seconds | **asserts** every number quoted in Appendices A and C against its named output column (73 assertions, including key uniqueness and row counts); each selector must return exactly one row, each vector must have its expected length and only finite values, and floating-point identities are asserted at the scale-aware tolerance 1e-10; exits 1 on any failure |
| 5b | `sims/check_manuscript_values.R --selftest` | console; exit status | seconds | negative controls: perturbed copies of the outputs (deleted row, duplicated row, altered value, non-finite value, deleted quantity, identity residual 1e-6) must each make the checker fail, and an unmodified copy and a residual of 5e-11 must pass |
| 6 | `sims/design_rank.R` | console | seconds | Theorem 1(d), Lemma 1(o) examples; the {1,4,8} design |
| 7 | `sims/check_thm1_support.py` | console | seconds | nullity = gcd on 14 C_d designs |
| 8 | `sims/check_lemma1_trapezoid.py` | console | ~10 s | Lemma 1 sweep (1,547 exact three-cohort designs; 9,401 designs for the bound) |
| 9 | `sims/check_lemma1_components.py` | console | ~10 s | Lemma 1(o) and (iv): exact rank for five named designs, then a sweep of 10,934 trapezoids |
| 10 | `sims/check_interrupted.py` | console | ~5 s | Theorem 2(a),(d), Table 1 (dose designs); Lemma 1(iv) sweep (10,934 designs; the 72 failures reported there are of the withdrawn third-cohort conjecture, not of the lemma); Proposition 5's weights and the constant-curvature equality test; §8.1's descriptive Table 2 (Appendix D) |
| 11 | `sims/check_dose_sharpness.py` | console | ~10 s | Theorem 2(b): the earlier counterexample (cohorts {1,2,4}, offsets {0,1,3}) and condition (P') verified by exact arithmetic over 496 designs; the CPS ten-cohort threshold; the two d > 1 examples cited in §3.4 (Appendix D) |
| 12 | `sims/check_recovery_support.py` | console; exit status | ~15 s | exact arithmetic: the six-cell disconnected support (C_1 holds, full nullity 3, projected nullity 1, tau(2) free with g known) and the identity dim K = dim K_tau + (c-1) on 9,183 supports (Theorem 1(d), condition (CG) of §5); the identified sets for tau(2) under a drift bound on the {1,2,5} and six-cell supports by exact Fourier–Motzkin elimination, against the interval Proposition 3 displays (block G: (CG) is sufficient for that interval to be the identified set, not necessary; without it the displayed interval depends on the representative and can be strictly narrower); the post-event-reference example for Theorem 5's premise (rank 16, pre-coefficient +1/2); the failure of (R) under D = 1{s>=4} with two exact fits; Corollary 1 at d = 2; the plateau boundary; the population widths of G4; the uninterrupted rows of Table 1 (Appendix D) |
| 13 | `sims/compare_outputs.py --selftest .` | console; exit status | seconds | self-test of the output comparator used by the release check, on corrupted copies of the shipped outputs in `sims/`: 21 corruptions (changed label, changed integer count with and without a label change, perturbed value, value replaced by NA, NA replaced by a value, appended record of empty fields, duplicated / deleted / swapped row, dropped / renamed column, missing file, changed transcript number, changed transcript verdict, deleted transcript line, and four controls added on 22 September 2026) must each fail, and 6 benign variants (unmodified copy, perturbation of 1e-12, NA written as NaN in a CSV and in the transcript, appended blank lines, trailing whitespace with CRLF endings) must pass |

Seeds are fixed inside each script (sim1–2: as in file; sim3: 20260830; sim4: 20260916). Scripts are run from `sims/` and write into it. The whole sequence (14 commands, counting 5b) takes about 50 minutes, of which Simulation 4's Monte Carlo (five designs, B = 100) is 30 and Simulation 3 is 18; steps 5–13 together take about a minute. Steps 6–13 are deterministic — no randomness and no data — so they either reproduce exactly or fail. Steps 11 and 12, the named designs of step 9, and blocks (A) and (D') of step 10 use exact rational arithmetic; steps 6, 7, 8 and the two 10,934-design sweeps (steps 9 and 10) use floating-point rank (numpy SVD at tolerance 1e-10, R's `qr`), which for these small 0–1 matrices returns the same integers as the exact code on every named design.

## How to read the release check

The release script (`release_check.sh` in the author's workflow, whose record is `RELEASE_CHECK.md`) exits with a
nonzero status — and says so in its last line — if any step of the sequence fails, if any file listed in
`FILE_MANIFEST.txt` is missing from the anonymous download or any unlisted file is present, or if
`sims/compare_outputs.py` reports a failure on any output file. A tag is
placed only on a snapshot whose check passed.

`RELEASE_CHECK.md` records a re-execution of the documented sequence in a clean copy of the published
snapshot, downloaded anonymously, with the shipped outputs set aside and the regenerated ones compared by
`sims/compare_outputs.py`, which is part of this archive so that the rule applied is the rule a reader can
inspect. The rule is file by file and cell by cell, with no fallback:

- every `.csv` and `.txt` output present in either set is compared, and a file present in only one set is a failure;
- for a CSV, the column names must be identical and in the same order, the numbers of rows and columns must be
  identical, a missing cell (empty, `NA` or `NaN`) must be missing in both, a cell that is a number in both must agree
  to within `TOL` (default 1e-8) in absolute value — integer counts included, so `B = 100` against `101` fails — and
  any other cell (a scenario label, a design name) must be byte-identical. (In the comparator shipped with tags
  `paper-v1.0` and `arxiv-2609.28871v1`, a record that contains only delimiters and empty fields was dropped before the row counts were
  compared, so such an appended record was not caught by this function alone; the named-value checker of step 5
  (`sims/check_manuscript_values.R`) rejects the resulting row count. From tag `paper-v1.1` on, such a record is counted as a row and the case is part of the self-test; only physical blank lines are ignored.)
- for the console transcript `sim4_results.txt`, the number of non-blank lines must be identical after trailing
  whitespace and line endings are normalised, each line must split into the same number of whitespace-separated
  tokens, and each token is compared by the same cell rule, so a changed rank (`rank=35` against `rank=34`) or a
  changed verdict (`not identified` against `identified`) fails.

Beyond whitespace (line endings, trailing and inter-token spacing, blank lines) and numeric spellings that agree to
within `TOL`, the only print-format difference tolerated is therefore an empty cell, `NA` and `NaN` for a missing
entry, which R's printers emit inconsistently across versions for a mean over an empty set (the re-executions under R 4.6.0 on 22 September
2026 printed `NaN` for exactly three such entries — the tenure-determined design's corrected bias, RMSE and coverage,
whose CSV entries are missing and whose unidentifiable-fit share is one — where the shipped transcript, written under
R 4.3.3, prints `NA`). Step 13 of the sequence runs the comparator's self-test on corrupted copies of
the shipped outputs and fails if any of its 20 corruptions is accepted or any of its 5 benign variants is rejected.

*History.* Until 22 September 2026 the release check used a comparator embedded in `release_check.sh` that fell
back, after any token mismatch, to comparing the real-number tokens of the two files in order. That fallback was
introduced to survive a `data.table` print-format change in the transcript, but it applied to CSVs as well and
discarded text and integer tokens, so — as controlled tests showed — a changed scenario label
or a changed integer count could be accepted as `format differs; numbers match`. The transcript is now written with
base R's data-frame printer at a fixed, very wide console, so that its layout does not depend on the installed
`data.table` version, and the comparator no longer has a fallback. The release records of the earlier tags
(`paper-v0.7`, `paper-v0.9`) were produced with the old comparator and do not carry the guarantee described above;
the records from `paper-v1.0` on re-execute the computational scripts of their snapshot and compare the outputs with
the comparator of that snapshot.

*Floating point.* The shipped outputs were produced under R 4.3.3 on x86-64 Linux; the release checks of
21–22 September 2026 ran R 4.6.0 on macOS arm64, with a different BLAS and LAPACK. Differences of order
1e-13 in regression coefficients and simulation summaries follow from that and from nothing else. The
deterministic identities the paper relies on are all checked *inside* the scripts against their own
analytic formulas, to tolerances stated in Appendix C, so a platform difference of this size cannot
disturb any claim; `sims/check_manuscript_values.R` asserts every number quoted in Appendices A and C
from its named output column, at a scale-aware tolerance of 1e-10 for the floating-point identities (the archived residuals of order 1e-14 to 1e-16 are platform-specific and are not themselves required to recur). What `check_manuscript_values.R` checks is the set of quantities the paper quotes, together with the row counts and key uniqueness of every output file; it is not a schema validator for the outputs as files — that is what `compare_outputs.py` is for.

The deterministic checks — steps 6 to 13 — carry no such caveat. They use no randomness, so they either
reproduce exactly or fail. In the 21 and 22 September 2026 checks
they reproduced exactly, including Lemma 1's sweep over 10,934 designs (no violation of either part; the
bound strict in 1,188 of them) and every entry of Table 2, on a different operating system, a different
R, a different Python and a different `numpy` from those used to write the paper.

## Map from manuscript objects to code

- Theorem 1(a)–(d), Corollary 1: `design_rank.R` (`design_rank()`, `increment_graph()`, `is_identified()`), `check_thm1_support.py`; Theorem 1(d)'s full-kernel/projection identity, condition (CG) of §5, the six-cell example and the identified sets of Proposition 3 with and without (CG) (block G), Corollary 1 at d = 2, Proposition 1's plateau boundary: `check_recovery_support.py`.
- Theorem 2(a),(d), Corollary 2, Table 1 (interrupted participation / dose designs) and Proposition 5's sensitivity weights and equality test: `check_interrupted.py`.
- Theorem 2(b), condition (P'), the counterexample to the earlier statement, and the d > 1 examples of §3.4: `check_dose_sharpness.py`.
- Lemma 1 (o)–(iv): `check_lemma1_components.py`, `check_lemma1_trapezoid.py`, `check_interrupted.py`, `design_rank.R`.
- Proposition 1 (plateau test), normalizations: `sim1_normalizations.R`.
- Corollary 1 readouts: `sim2_direct_contrast.R`.
- Proposition 3 (drift bounds, IM intervals): `sim3_drift_bounds.R`.
- Theorem 3 (absorption, kernel invariance), Theorem 4 (joint regression, (R)), Theorem 5 (pre-trend identity, rank-deficient design, exact-alignment rank), Proposition 4 (partial-alignment counterexample), Corollary 3(a): `sim4_downstream.R`; Theorem 4's aliasing caveat and Theorem 5's reference-category premise (exact examples): `check_recovery_support.py`.
- Simulation 3's population widths (G2, G3, G4): `check_manuscript_values.R` and `check_recovery_support.py` (E).
- Release verification (regenerated versus shipped outputs): `compare_outputs.py`; its `--selftest` is step 13.
- Figures 1 and 2 (the Japanese panel's support shaded by tenure class; the two constant-rate paths of Corollary 2(iii) under the CPS 4–8–4 pattern): `manuscript/figures/make_figures.py`, which draws them from the stated designs without data and first checks the quantities the text states (41 cells, rank 35 and kernel dimension 4; D = 0 and D = −8c). Added with the arXiv version (dated 24 September 2026); no numerical result changes.

## Manuscript

`manuscript/main.qmd` (source), `references.bib`, `latex/preamble_extra_T1.tex`, `latex/chicago-author-date.csl`, `build_latex.py`, `figures/` (Figures 1 and 2 and `make_figures.py`). Render: `python3 figures/make_figures.py`, then `quarto render main.qmd --to html` and `python3 build_latex.py`.

*Versions.* The numerical results are those of tag `paper-v1.0` (computational commit `f06283e`). The manuscript in that tag is the 73-page version of 23 September 2026. The arXiv version (dated 24 September 2026) adds presentation and reference updates and Figures 1 and 2, and changes no numerical result; its source, PDF and figure script are in `manuscript/` from the tag `arxiv-2609.28871v1` on, and the scripts and outputs outside `manuscript/` in that tag are those of `paper-v1.0`. Tag `paper-v1.1` changes, outside `manuscript/`, only the output comparator (a record of empty fields is counted as a row, with the two self-test cases above) and comments or a printed label in four check scripts; no script that produces a number of the paper changed, its release check re-executes the documented sequence with the new comparator, and the numerical results are those of `paper-v1.0`.

## Citation
Okubo, S. (2026). Panel Conditioning in Fixed-Effects Models: Identification and Bias Propagation. arXiv preprint arXiv:2609.28871.

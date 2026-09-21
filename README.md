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
Computational commit checked against the manuscript: not yet checked (a release check is run after the first publication; see `RELEASE_CHECK.md` once present)
Tag matching this version of the manuscript: `paper-v0.7`. Tag matching the posted preprint version: to be added at posting (`arxiv-<id>v<n>`).
Third-party reproduction: none. The release record is the author's own re-execution of the published snapshot in a clean copy.

## Data
No data are distributed with this archive and none are needed. Every number in the paper comes from a seeded simulation that generates its own data, from a deterministic algebraic computation on a stated design, or from month-in-sample indices published in Bailar (1975, Tab. 1), Solon (1986, Tab. 1) and McIllece (2022, Tab. 2), which are reproduced inside `sims/check_interrupted.py`. The Japanese Life Course Panel Surveys (JLPS) microdata are **not** used anywhere in this paper: Section 8 is a design illustration whose only inputs are the entry waves 1, 5, 13 and the final wave 19.

## Environment

- R 4.3.3 (also re-run independently under R 4.6.0, arm64, with byte-identical results up to 1e-13); packages `data.table` 1.14.10, `sandwich` 3.1.0 (Simulation 4 only).
- Python 3 with `numpy` 2.4 (support checks only).
- Quarto 1.6.42 + XeLaTeX (TeX Live 2023) for the manuscript; `manuscript/build_latex.py` reproduces `main.pdf` and `latex/main.tex`.

## Execution order and runtimes (single core)

| Step | Script | Output | Runtime | Paper location |
|---|---|---|---|---|
| 1 | `sims/sim1_normalizations.R` | `sim1_results.csv` | ~1 min | Appendix A.1 (table; noiseless identities) |
| 2 | `sims/sim2_direct_contrast.R` | `sim2_results.csv` | ~1 min | Appendix A.2 |
| 3 | `sims/sim3_drift_bounds.R` (v3) | `sim3_results.csv` | ~12 min | Appendix A.3 (both DGPs; IM and fixed-z coverage) |
| 4 | `sims/sim4_downstream.R` (v3) | `sim4_results.csv`, `sim4_pretrend_coefs.csv`, `sim4_mc.csv`, `sim4_manuscript_values.csv`, `sim4_results.txt` | ~30 min | Appendix C (identities, alignment checks, Monte Carlo table) |
| 5 | `sims/check_manuscript_values.R` | console | seconds | prints every number quoted in Appendices A and C from named columns |
| 6 | `sims/design_rank.R` | console | seconds | Theorem 1(d), Lemma 1(o) examples; the {1,4,8} design |
| 7 | `sims/check_thm1_support.py` | console | seconds | nullity = gcd on 14 C_d designs |
| 8 | `sims/check_lemma1_trapezoid.py` | console | ~5 min | Lemma 1 sweep (1,547 exact three-cohort designs; 9,401 designs for the bound) |
| 9 | `sims/check_lemma1_components.py` | console | ~10 s | Lemma 1(o) and (iv): exact rank for five named designs, then a sweep of 10,934 trapezoids |
| 10 | `sims/check_interrupted.py` | console | ~5 min | Theorem 2, Corollary 2, Table 1 (dose designs); Lemma 1(iv) sweep (10,934 designs); Proposition 5's weights; §8.1's Table 2 from published CPS indices (Appendix D) |

Seeds are fixed inside each script (sim1–2: as in file; sim3: 20260830; sim4: 20260916). Scripts are run from `sims/` and write into it.

## Map from manuscript objects to code

- Theorem 1(a)–(d), Corollary 1: `design_rank.R` (`design_rank()`, `increment_graph()`, `is_identified()`), `check_thm1_support.py`.
- Theorem 2, Corollary 2, Table 1 (interrupted participation / dose designs) and Proposition 5's sensitivity weights: `check_interrupted.py`.
- Lemma 1 (o)–(iv): `check_lemma1_components.py`, `check_lemma1_trapezoid.py`, `check_interrupted.py`, `design_rank.R`.
- Proposition 1 (plateau test), normalizations: `sim1_normalizations.R`.
- Corollary 1 readouts: `sim2_direct_contrast.R`.
- Proposition 3 (drift bounds, IM intervals): `sim3_drift_bounds.R`.
- Theorem 3 (absorption, kernel invariance), Theorem 4 (joint regression, (R)), Theorem 5 (pre-trend identity, rank-deficient design, exact-alignment rank), Proposition 4 (partial-alignment counterexample), Corollary 3(a): `sim4_downstream.R`.

## Manuscript

`manuscript/main.qmd` (source), `references.bib`, `latex/preamble_extra_T1.tex`, `build_latex.py`. Render: `quarto render main.qmd --to html` and `python3 build_latex.py`.

## Citation
Okubo, S. (2026). Panel Conditioning in Fixed-Effects Models: Identification and Bias Propagation. Working paper. arXiv:XXXX.XXXXX.

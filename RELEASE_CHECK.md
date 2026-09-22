# Release check — paper-panel-conditioning-replication

Date: 2026-09-22T09:56:34Z. Snapshot downloaded anonymously (no credentials, no gh CLI) from `https://codeload.github.com/sokubo/paper-panel-conditioning-replication/tar.gz/ea044e655ced00e05f2a2442ab0fc19c607c5d77`.

- ref: `ea044e655ced00e05f2a2442ab0fc19c607c5d77`; commit: `ea044e655ced00e05f2a2442ab0fc19c607c5d77`
- archive SHA-256: `7b06d28b5ac4b957d55aca80065bd2e90c2e10fd46f19684da5b25cf0ae63bf0`
- files in snapshot (excluding FILE_MANIFEST.txt and RELEASE_CHECK*): 30; listed in FILE_MANIFEST.txt: 30; missing from snapshot: 0; not listed in manifest: 0
- restricted-data / review-material scan of the published snapshot: 0 file(s) matched
- clean run: documented sequence executed in a clean copy with shipped outputs set aside (1475s); log and sessionInfo kept; comparison below
- third-party reproduction: none; this record is the author's own re-execution.

## Environment of the clean run

```
R version 4.6.0 (2026-04-24)
Platform: aarch64-apple-darwin23
Running under: macOS Sequoia 15.7.3

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] ja_JP.UTF-8/ja_JP.UTF-8/ja_JP.UTF-8/C/ja_JP.UTF-8/ja_JP.UTF-8

time zone: Asia/Tokyo
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] sandwich_3.1-1    data.table_1.18.4

loaded via a namespace (and not attached):
[1] zoo_1.8-15     compiler_4.6.0 grid_4.6.0     lattice_0.22-9

Python side:
Python 3.14.7
numpy 2.5.3
```

## Regenerated vs shipped outputs

```
                            file                        status max_abs_diff
           sims/sim1_results.csv                     identical 0.000000e+00
           sims/sim2_results.csv                     identical 0.000000e+00
           sims/sim3_results.csv            numeric difference 9.992007e-15
 sims/sim4_manuscript_values.csv            numeric difference 1.540434e-13
                sims/sim4_mc.csv            numeric difference 9.992007e-16
    sims/sim4_pretrend_coefs.csv            numeric difference 3.589351e-13
           sims/sim4_results.csv            numeric difference 3.497203e-15
           sims/sim4_results.txt format differs; numbers match 4.330000e-15

files compared: 8; identical: 2; numeric difference: 5; format differs, numbers match: 1; max abs difference: 3.59e-13; not comparable: 0

format diagnostics (printed whenever a token stream differs, so that a genuine change is never hidden):
  sims/sim4_results.txt: regenerated 501 tokens vs shipped 501; first divergence at token 147
      regenerated: tau_sat 0.5682 0.5937 2.556e-02 2.556e-02 5.516e-15 4.530e-14 6.217e-15 NA NA NA
      shipped:     tau_sat 0.5682 0.5937 2.556e-02 2.556e-02 7.140e-15 4.441e-14 5.440e-15 NA NA NA 
```

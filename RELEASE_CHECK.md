# Release check — paper-panel-conditioning-replication

Date: 2026-09-22T08:47:50Z. Snapshot downloaded anonymously (no credentials, no gh CLI) from `https://codeload.github.com/sokubo/paper-panel-conditioning-replication/tar.gz/a3b36521410184b0ce27815b6f8edef9ba09aba1`.

- ref: `a3b36521410184b0ce27815b6f8edef9ba09aba1`; commit: `a3b36521410184b0ce27815b6f8edef9ba09aba1`
- archive SHA-256: `77059eaa8956c66dbfb341f7d3a6249b4e5030a2c7ea41aef08863919d40a4b1`
- files in snapshot (excluding FILE_MANIFEST.txt and RELEASE_CHECK*): 30; listed in FILE_MANIFEST.txt: 30; missing from snapshot: 0; not listed in manifest: 0
- restricted-data / review-material scan of the published snapshot: 0 file(s) matched
- clean run: documented sequence executed in a clean copy with shipped outputs set aside (1433s); log and sessionInfo kept; comparison below
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
           sims/sim4_results.txt format differs; numbers match 8.488000e-01

files compared: 8; identical: 2; numeric difference: 5; format differs, numbers match: 1; max abs difference: 8.49e-01; not comparable: 0

format diagnostics (printed whenever a token stream differs, so that a genuine change is never hidden):
  sims/sim4_results.txt: lengths 581 vs 571; first divergence at token 132
      shipped:     2: A_absorption tau_lin 0.5682 0.5682 7.772e-15 3: A_absorption tau_mix 0.5682 0.5937
      regenerated: 2: A_absorption tau_lin 0.5682 0.5682 8.771e-15 3: A_absorption tau_mix 0.5682 0.5937 
```

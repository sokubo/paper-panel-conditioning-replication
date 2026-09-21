# Release check — paper-panel-conditioning-replication

Date: 2026-09-21T07:42:15Z. Snapshot downloaded anonymously (no credentials, no gh CLI) from `https://codeload.github.com/sokubo/paper-panel-conditioning-replication/tar.gz/a48bd805f765ac4fa8fd282437affedf2af57e31`.

- ref: `a48bd805f765ac4fa8fd282437affedf2af57e31`; commit: `a48bd805f765ac4fa8fd282437affedf2af57e31`
- archive SHA-256: `382e205c7c615451d0d90edaaefdbf899b80ad6808e2dbba0462e8902d56d199`
- files in snapshot (excluding FILE_MANIFEST.txt and RELEASE_CHECK*): 28; listed in FILE_MANIFEST.txt: 27; missing from snapshot: 0; not listed in manifest: 1
- restricted-data / review-material scan of the published snapshot: 0 file(s) matched

Not listed:
- .gitignore
- clean run: documented sequence executed in a clean copy with shipped outputs set aside (1393s); log and sessionInfo kept; comparison below
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
```

## Regenerated vs shipped outputs

```
                            file                status max_abs_diff
           sims/sim1_results.csv             identical 0.000000e+00
           sims/sim2_results.csv             identical 0.000000e+00
           sims/sim3_results.csv    numeric difference 9.992007e-15
 sims/sim4_manuscript_values.csv    numeric difference 1.540434e-13
                sims/sim4_mc.csv    numeric difference 9.992007e-16
    sims/sim4_pretrend_coefs.csv    numeric difference 3.589351e-13
           sims/sim4_results.csv    numeric difference 3.497203e-15
           sims/sim4_results.txt different token count           NA

files compared: 8; identical: 2; numeric difference: 5 (max 3.59e-13); other: 1
```

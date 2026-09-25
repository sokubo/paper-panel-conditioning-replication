# Release check — paper-panel-conditioning-replication

Date: 2026-09-25T03:59:38Z. Snapshot downloaded anonymously (no credentials, no gh CLI) from `https://codeload.github.com/sokubo/paper-panel-conditioning-replication/tar.gz/d8561e30d141f0dd7fe171459fa95f28678e7e0d`.

- ref: `d8561e30d141f0dd7fe171459fa95f28678e7e0d`; commit: `d8561e30d141f0dd7fe171459fa95f28678e7e0d`
- archive SHA-256: `0eb6fce8384bef65f66fc2edc646406819ebe06ddfbe8b7b4a566c6f9701b10f`
- files in snapshot (excluding FILE_MANIFEST.txt and RELEASE_CHECK*): 37; listed in FILE_MANIFEST.txt: 37; missing from snapshot: 0; not listed in manifest: 0
- content scan of the published snapshot: 0 file(s) matched
- clean run: documented sequence executed in a clean copy with shipped outputs set aside (1434s); log and sessionInfo kept; comparison below
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
file                        max |diff|          cells/tokens  status
sim1_results.csv            0.000e+00           -             identical
sim2_results.csv            0.000e+00           -             identical
sim3_results.csv            9.992e-15           250           within tolerance
sim4_manuscript_values.csv  1.540e-13           52            within tolerance
sim4_mc.csv                 9.992e-16           60            within tolerance
sim4_pretrend_coefs.csv     3.589e-13           66            within tolerance
sim4_results.csv            3.497e-15           250           within tolerance
sim4_results.txt            4.330e-15           491           within tolerance

files compared: 8; identical: 2; within tolerance: 6; FAILED: 0; largest numeric difference: 3.589e-13; tol: 1e-08
comparison passed
```

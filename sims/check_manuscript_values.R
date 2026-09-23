# ============================================================
# check_manuscript_values.R (v4, 2026-09-22) — every number quoted in Appendices A and C of T1 (v1.0) is ASSERTED
# against the named column of the simulation output it comes from. Exit status 1 on any failure.
#
# v3 hardening (after the fourth review round):
#   * every row selector must return EXACTLY ONE row (a deleted or duplicated row is a failure, not a vacuous pass:
#     v2 let all(numeric(0) == x) succeed on an empty selection);
#   * every compared vector must have the expected length and only finite values;
#   * floating-point IDENTITIES are asserted at a scale-aware tolerance (IDENT_TOL = 1e-10, five orders of magnitude
#     below any displayed precision) instead of the archived machine-error residuals, which vary across platforms;
#   * `--selftest` perturbs copies of the outputs (deleted row, duplicated row, altered value, non-finite value) and
#     requires the checker to FAIL on each; it is run by the release check.
# Run from sims/ after sim1–sim4:   Rscript check_manuscript_values.R            (checks the outputs in ./)
#                                   Rscript check_manuscript_values.R --selftest (negative controls on temp copies)
#                                   Rscript check_manuscript_values.R <dir>      (checks the outputs in <dir>)
# ============================================================
options(width = 120)
IDENT_TOL <- 1e-10

run_checks <- function(dir = ".", quiet = FALSE) {
  rd <- function(f) { p <- file.path(dir, f); if (!file.exists(p)) stop("missing output file: ", p); read.csv(p, stringsAsFactors = FALSE) }
  s1 <- rd("sim1_results.csv"); s2 <- rd("sim2_results.csv"); s3 <- rd("sim3_results.csv")
  s4 <- rd("sim4_mc.csv"); mv <- rd("sim4_manuscript_values.csv")
  ok <- 0L; bad <- 0L
  say <- function(...) if (!quiet) cat(sprintf(...))
  chk <- function(what, got, want, digits = NULL, tol = NULL) {
    g <- suppressWarnings(as.numeric(got)); w <- as.numeric(want)
    reason <- NULL
    if (length(g) != length(w) || length(w) == 0L) reason <- sprintf("length %d, expected %d", length(g), length(w))
    else if (any(!is.finite(g))) reason <- "non-finite value"
    else if (!is.null(digits)) { if (!all(round(g, digits) == round(w, digits))) reason <- "value differs at the displayed rounding" }
    else if (!all(abs(g - w) <= tol)) reason <- sprintf("differs by more than tol=%g", tol)
    hit <- is.null(reason)
    if (hit) ok <<- ok + 1L else bad <<- bad + 1L
    say("%-66s %-18s %-18s %s\n", what, paste(signif(g, 6), collapse = " "), paste(w, collapse = " "),
        if (hit) "ok" else paste("MISMATCH:", reason))
  }
  one <- function(df, cond, what) {                          # exactly one row, else a recorded failure and an NA row
    r <- df[cond & !is.na(cond), , drop = FALSE]
    if (nrow(r) != 1L) { bad <<- bad + 1L; say("%-66s %s\n", what, sprintf("MISMATCH: %d rows selected, expected 1", nrow(r)))
      r <- df[0, , drop = FALSE]; r[1, ] <- NA }
    r
  }
  val <- function(q) { v <- mv$value[mv$quantity == q]; if (length(v) != 1L) NA_real_ else as.numeric(v) }
  say("%-66s %-18s %-18s\n", "quantity", "from output", "quoted"); say("%s\n", strrep("-", 112))

  ## ---- structural checks: key uniqueness and completeness of every output file --------------------------------
  chk("keys unique in sim1 (scenario x cohort_g)", anyDuplicated(s1[, c("scenario", "cohort_g")]) == 0, 1, tol = 0)
  chk("sim1 has 8 rows (S1-S4 x G0,G1)", nrow(s1), 8, tol = 0)
  chk("keys unique in sim3 (dgp x s_eval x M_mult)", anyDuplicated(s3[, c("dgp", "s_eval", "M_mult")]) == 0, 1, tol = 0)
  chk("sim3 has 24 rows (G2,G3,G4 x s=13,17 x 4 M)", nrow(s3), 24, tol = 0)
  chk("keys unique in sim4_mc (design)", anyDuplicated(s4$design) == 0, 1, tol = 0)
  chk("keys unique in sim4_manuscript_values (quantity)", anyDuplicated(mv$quantity) == 0, 1, tol = 0)
  chk("all values in sim4_manuscript_values finite", all(is.finite(suppressWarnings(as.numeric(mv$value)))), 1, tol = 0)

  ## ---- Appendix A.1 (sim1, G0 rows) ------------------------------------------------------------------------
  tau2 <- function(s) 0.2 * (1 - exp(-0.6 * (s - 1)))
  chk("A.1 closed form Delta_4^2 tau(5) = tau(9)-2tau(5)+tau(1), S2", tau2(9) - 2 * tau2(5) + tau2(1), -0.1654, digits = 4)
  chk("A.1 closed form Delta_4^2 tau(9) = tau(13)-2tau(9)+tau(5), S2", tau2(13) - 2 * tau2(9) + tau2(5), -0.0150, digits = 4)
  tab <- list(S1 = c(0.000, 0.000, 0.000, 0.001, .064, 0.001, .058),
              S2 = c(-0.051, -0.102, -0.165, 0.004, .061, -0.006, .048),
              S3 = c(-0.089, -0.140, -0.165, 0.001, .065, -0.097, .054),
              S4 = c(-0.173, -0.124, -0.129, 0.005, .060, -0.408, 1.000))
  g0 <- s1[s1$cohort_g == "G0", ]
  for (sc in names(tab)) {
    r <- one(s1, s1$cohort_g == "G0" & s1$scenario == sc, paste("A.1 row", sc, "(G0)"))
    chk(paste("A.1 table row", sc), c(r$bias_lev_N1, r$bias_lev_N2, r$true_d4, r$bias_d4_contrast, r$sd_d4_contrast, r$bias_taubar_N3, r$rej_overid), tab[[sc]], digits = 3)
  }
  chk("A.1 contrast MC standard error <= .003 (sd/sqrt(500))", length(g0$sd_d4_contrast) == 4 && max(g0$sd_d4_contrast) / sqrt(500) <= 0.003, 1, tol = 0)
  chk("A.1 plateau-constrained bias under S3", one(s1, s1$cohort_g == "G0" & s1$scenario == "S3", "A.1 S3")$bias_taubar_N3, -0.097, digits = 3)
  chk("A.1 overid rejection rates S1-S4 (G0), quoted in the text (.058 .048 .054 1)",
      sapply(c("S1", "S2", "S3", "S4"), function(sc) one(s1, s1$cohort_g == "G0" & s1$scenario == sc, paste("A.1", sc))$rej_overid), c(.058, .048, .054, 1), digits = 3)
  chk("A.1 overid rejection rate S4 (G1), quoted as 499 of 500", one(s1, s1$cohort_g == "G1" & s1$scenario == "S4", "A.1 S4 G1")$rej_overid, 499 / 500, digits = 3)

  ## ---- Appendix A.2 (sim2) ---------------------------------------------------------------------------------
  chk("sim2 has 2 rows (S2, S3)", nrow(s2), 2, tol = 0)
  r2 <- one(s2, s2$scenario == "S2", "A.2 S2"); r3 <- one(s2, s2$scenario == "S3", "A.2 S3")
  chk("A.2 theta true", r2$theta_true, -0.037, digits = 3)
  chk("A.2 table S2", c(r2$bias_contrast, r2$sd_contrast, r2$bias_readout, r2$sd_readout, r2$d2_3_true_unidentified, r2$mean_d2_3_readout),
      c(0.006, .110, 0.000, .080, -0.022, -0.075), digits = 3)
  chk("A.2 table S3", c(r3$bias_contrast, r3$sd_contrast, r3$bias_readout, r3$sd_readout, r3$d2_3_true_unidentified, r3$mean_d2_3_readout),
      c(0.003, .107, 0.003, .081, -0.022, -0.072), digits = 3)
  chk("A.2 MC standard error <= .005", length(s2$sd_contrast) == 2 && max(s2$sd_contrast) / sqrt(500) <= 0.005, 1, tol = 0)

  ## ---- Appendix A.3 (sim3 v4: G2, G3, G4) ------------------------------------------------------------------
  row3 <- function(d, s, m) one(s3, s3$dgp == d & s3$s_eval == s & s3$M_mult == m, sprintf("A.3 row %s s=%d M/M0=%g", d, s, m))
  a3 <- list(
    list("G2", 1,   c(.957, .957, 0.69, 0.48, 1.64), c(.957, .957, 0.91)),
    list("G2", 2,   c(1, 1, 1.29, 1.08, 1.64),       c(1, 1, 1.71)),
    list("G2", 4,   c(1, 1, 2.49, 2.28, 1.64),       c(1, 1, 3.31)),
    list("G3", 1,   c(.953, .900, 0.25, 0.00, 1.93), c(.940, .910, 0.33)),
    list("G3", 2,   c(1, 1, 0.82, 0.60, 1.64),       c(1, 1, 1.08)),
    list("G3", 4,   c(1, 1, 2.02, 1.80, 1.64),       c(1, 1, 2.68)),
    list("G4", 1,   c(.880, .880, 0.74, 0.55, 1.64), c(.890, .890, 0.99)),
    list("G4", 2,   c(1, 1, 1.34, 1.15, 1.64),       c(1, 1, 1.79)),
    list("G4", 4,   c(1, 1, 2.54, 2.35, 1.64),       c(1, 1, 3.39)))
  for (e in a3) {
    r13 <- row3(e[[1]], 13, e[[2]]); r17 <- row3(e[[1]], 17, e[[2]])
    chk(sprintf("A.3 %s M/M0=%g, tau(13): cov IM/fixed, width, set width, c", e[[1]], e[[2]]),
        c(r13$coverage, r13$coverage_fixed_z, r13$med_width, r13$med_set_width, r13$med_c), e[[3]], digits = 2)
    chk(sprintf("A.3 %s M/M0=%g, tau(17): cov IM/fixed, width", e[[1]], e[[2]]),
        c(r17$coverage, r17$coverage_fixed_z, r17$med_width), e[[4]], digits = 2)
  }
  # violated-bound rows and the empty rates quoted as ranges
  chk("A.3 G2 0.4: coverage tau(13) IM/fixed", c(row3("G2",13,.4)$coverage, row3("G2",13,.4)$coverage_fixed_z), c(.127, .123), digits = 3)
  chk("A.3 G2 0.4: coverage tau(17) IM/fixed", c(row3("G2",17,.4)$coverage, row3("G2",17,.4)$coverage_fixed_z), c(.107, .107), digits = 3)
  chk("A.3 G3 0.4: coverage tau(13) IM/fixed", c(row3("G3",13,.4)$coverage, row3("G3",13,.4)$coverage_fixed_z), c(.020, .007), digits = 3)
  chk("A.3 G3 0.4: coverage tau(17) IM/fixed", c(row3("G3",17,.4)$coverage, row3("G3",17,.4)$coverage_fixed_z), c(.023, .007), digits = 3)
  chk("A.3 G4 0.4: coverage tau(13) IM/fixed", c(row3("G4",13,.4)$coverage, row3("G4",13,.4)$coverage_fixed_z), c(.017, .017), digits = 3)
  chk("A.3 G4 0.4: coverage tau(17) IM/fixed", c(row3("G4",17,.4)$coverage, row3("G4",17,.4)$coverage_fixed_z), c(.013, .013), digits = 3)
  er <- function(d, m) c(row3(d, 13, m)$empty_rate, row3(d, 17, m)$empty_rate)
  chk("A.3 G2 0.4 empty rate within .05-.07", all(is.finite(er("G2", .4))) && all(er("G2", .4) >= .045 & er("G2", .4) <= .075), 1, tol = 0)
  chk("A.3 G3 1.0 empty rate within .47-.50", all(is.finite(er("G3", 1))) && all(er("G3", 1) >= .465 & er("G3", 1) <= .505), 1, tol = 0)
  chk("A.3 G3 0.4 empty rate = 1", er("G3", .4), c(1, 1), tol = 0)
  g4e <- s3$empty_rate[s3$dgp == "G4"]
  chk("A.3 G4 empty rate = 0 at every M (8 rows)", length(g4e) == 8 && all(g4e == 0), 1, tol = 0)
  chk("A.3 MC standard error of a coverage of .95 with B=300 ~ .013", sqrt(.95 * .05 / 300), .013, digits = 3)
  chk("A.3 MC standard error at .880 and .890 with B=300 ~ .019, .018", sqrt(c(.88 * .12, .89 * .11) / 300), c(.019, .018), digits = 3)
  # population widths of the drift-bound sets (Proposition 3): 2M - (max d - min d), scaled by (s - 1)
  M0 <- 0.025; drift <- list(G2 = c(.025, .015), G3 = c(.025, -.025), G4 = c(.025, .025))
  wm <- sapply(drift, function(d) 2 * M0 - (max(d) - min(d)))
  chk("A.3 population width of the set for m at M=M0: G2, G3, G4", wm, c(.040, 0, .050), tol = 1e-12)
  chk("A.3 G4 population widths for tau(13), tau(17) at M=M0: .600, .800", c(12, 16) * wm[["G4"]], c(.600, .800), tol = 1e-12)
  chk("A.3 G4 median estimated set width for tau(13) below its population width .600 (biased-down maximum)", row3("G4",13,1)$med_set_width < .600, 1, tol = 0)

  ## ---- Appendix C (sim4) -----------------------------------------------------------------------------------
  chk("C: TWFE bias, saturating path", val("bias_TWFE_tau_sat"), 0.0256, digits = 4)
  chk("C: |bias - formula| within identity tolerance", val("max_abs_diff_formula"), 0, tol = IDENT_TOL)
  chk("C: affine-direction invariance gap within identity tolerance", val("invariance_gap_affine_max"), 0, tol = IDENT_TOL)
  chk("C: periodic-direction invariance gap within identity tolerance", val("invariance_gap_periodic_max"), 0, tol = IDENT_TOL)
  chk("C: single-cohort bias within identity tolerance", val("single_cohort_bias"), 0, tol = IDENT_TOL)
  chk("C: pre-trend identity error within identity tolerance", val("pretrend_max_abs_err"), 0, tol = IDENT_TOL)
  chk("C: pre-period SHIFTS at k=-5,-4,-3", c(val("pretrend_SHIFT_k-5"), val("pretrend_SHIFT_k-4"), val("pretrend_SHIFT_k-3")), c(-0.183, -0.070, -0.018), digits = 3)
  chk("C: pre-period COEFS on reports", c(val("pretrend_COEF_reported_k-5"), val("pretrend_COEF_reported_k-4"), val("pretrend_COEF_reported_k-3")), c(0.068, 0.061, 0.014), digits = 3)
  chk("C: pre-period COEFS on the unconditioned outcome", c(val("pretrend_COEF_unconditioned_k-5"), val("pretrend_COEF_unconditioned_k-4"), val("pretrend_COEF_unconditioned_k-3")), c(0.252, 0.131, 0.032), digits = 3)
  chk("C: partial-alignment shifts at k=-5,-4,-3", c(val("partial_alignment_shift_k-5"), val("partial_alignment_shift_k-4"), val("partial_alignment_shift_k-3")), c(-0.1419, -0.0544, -0.0142), digits = 4)
  chk("6.4: Simulation-4 shifts quoted after Proposition 5 (.18, .07, .02)", abs(c(val("pretrend_SHIFT_k-5"), val("pretrend_SHIFT_k-4"), val("pretrend_SHIFT_k-3"))), c(.18, .07, .02), digits = 2)
  chk("6.4: max centered curvature of the saturating path = .061", 0.30 * (1 - exp(-0.6))^2, .061, digits = 3)
  mc <- list(baseline = c(0.025, .067, .95, 0.002, .063, .98), zero_conditioning = c(0.010, .063, .99, 0.012, .065, .98),
             heterogeneous_effect = c(0.018, .059, .99, -0.006, .058, .99), cohort_adoption = c(0.042, .074, .92, -0.001, .064, .98))
  for (dn in names(mc)) { r <- one(s4, s4$design == dn, paste("C MC row", dn))
    chk(paste("C MC table", dn), c(r$bias_reported, r$rmse_reported, r$cover_reported, r$bias_corrected, r$rmse_corrected, r$cover_corrected), mc[[dn]], digits = 3) }
  r <- one(s4, s4$design == "tenure_adoption", "C MC row tenure_adoption")
  chk("C MC table tenure_adoption (reported)", c(r$bias_reported, r$rmse_reported, r$cover_reported), c(0.149, .164, .46), digits = 3)
  chk("C MC table tenure_adoption: corrected estimator not identified", r$share_corrected_NA, 1, tol = 0)

  say("%s\n", strrep("-", 112))
  say("checked %d assertions (quoted quantities and structural checks): %d reproduced, %d MISMATCHED\n", ok + bad, ok, bad)
  list(ok = ok, bad = bad)
}

## ---- negative controls ----------------------------------------------------------------------------------------
selftest <- function(src = ".") {
  files <- c("sim1_results.csv", "sim2_results.csv", "sim3_results.csv", "sim4_mc.csv", "sim4_manuscript_values.csv")
  mk <- function(mut) { d <- tempfile("chk_"); dir.create(d); for (f in files) file.copy(file.path(src, f), d); mut(d); d }
  cases <- list(
    "unmodified copy (must PASS)" = list(mut = function(d) NULL, expect_fail = FALSE),
    "delete the G4/s=13/M=M0 row of sim3 (must FAIL)" = list(mut = function(d) { p <- file.path(d, "sim3_results.csv"); x <- read.csv(p)
      write.csv(x[!(x$dgp == "G4" & x$s_eval == 13 & x$M_mult == 1), ], p, row.names = FALSE) }, expect_fail = TRUE),
    "duplicate the G4/s=13/M=M0 row of sim3 (must FAIL)" = list(mut = function(d) { p <- file.path(d, "sim3_results.csv"); x <- read.csv(p)
      write.csv(rbind(x, x[x$dgp == "G4" & x$s_eval == 13 & x$M_mult == 1, ]), p, row.names = FALSE) }, expect_fail = TRUE),
    "set that row's coverage to .50 (must FAIL)" = list(mut = function(d) { p <- file.path(d, "sim3_results.csv"); x <- read.csv(p)
      x$coverage[x$dgp == "G4" & x$s_eval == 13 & x$M_mult == 1] <- .50; write.csv(x, p, row.names = FALSE) }, expect_fail = TRUE),
    "set that row's coverage to NA (must FAIL)" = list(mut = function(d) { p <- file.path(d, "sim3_results.csv"); x <- read.csv(p)
      x$coverage[x$dgp == "G4" & x$s_eval == 13 & x$M_mult == 1] <- NA; write.csv(x, p, row.names = FALSE) }, expect_fail = TRUE),
    "delete the baseline row of sim4_mc (must FAIL)" = list(mut = function(d) { p <- file.path(d, "sim4_mc.csv"); x <- read.csv(p)
      write.csv(x[x$design != "baseline", ], p, row.names = FALSE) }, expect_fail = TRUE),
    "delete the pretrend_max_abs_err quantity (must FAIL)" = list(mut = function(d) { p <- file.path(d, "sim4_manuscript_values.csv"); x <- read.csv(p)
      write.csv(x[x$quantity != "pretrend_max_abs_err", ], p, row.names = FALSE) }, expect_fail = TRUE),
    "set an identity residual to 1e-6 (must FAIL)" = list(mut = function(d) { p <- file.path(d, "sim4_manuscript_values.csv"); x <- read.csv(p)
      x$value[x$quantity == "invariance_gap_affine_max"] <- 1e-6; write.csv(x, p, row.names = FALSE) }, expect_fail = TRUE),
    "set an identity residual to 5e-11 (must PASS: within the scale-aware tolerance)" = list(mut = function(d) { p <- file.path(d, "sim4_manuscript_values.csv"); x <- read.csv(p)
      x$value[x$quantity == "invariance_gap_affine_max"] <- 5e-11; write.csv(x, p, row.names = FALSE) }, expect_fail = FALSE),
    "delete the S3/G0 row of sim1 (must FAIL)" = list(mut = function(d) { p <- file.path(d, "sim1_results.csv"); x <- read.csv(p)
      write.csv(x[!(x$scenario == "S3" & x$cohort_g == "G0"), ], p, row.names = FALSE) }, expect_fail = TRUE))
  allok <- TRUE
  cat("self-test: negative controls on perturbed copies of the outputs\n")
  for (nm in names(cases)) {
    d <- mk(cases[[nm]]$mut); res <- run_checks(d, quiet = TRUE); unlink(d, recursive = TRUE)
    failed <- res$bad > 0; good <- failed == cases[[nm]]$expect_fail; allok <- allok && good
    cat(sprintf("   %-82s -> %d mismatch(es): %s\n", nm, res$bad, if (good) "as expected" else "UNEXPECTED"))
  }
  cat(if (allok) "self-test passed: every perturbation is detected and the unmodified copy passes\n" else "SELF-TEST FAILED\n")
  allok
}

args <- commandArgs(TRUE)
if (length(args) && args[1] == "--selftest") {
  src <- if (length(args) > 1) args[2] else "."
  quit(status = if (selftest(src)) 0 else 1)
} else {
  res <- run_checks(if (length(args)) args[1] else ".")
  if (res$bad > 0) quit(status = 1)
}

# ============================================================
# T1 — Simulation 4 (v3, 2026-09-20 evening): downstream propagation of tenure drift
#   v3 (revision): attrition is MONOTONE (a unit that leaves does not return; hazard increasing in
#   tenure), consistent with D2; block C adds the exact-alignment rank check [MX, MZ_tau] and the
#   partial-alignment design (delta_i in {6,7}, composition balanced, nonzero shift) with the corrected condition
#   X'M Z_tau; the pre-trend table reports coefficients AND shifts in named columns; a manuscript-value
#   check writes every quoted number to sim4_manuscript_values.csv.
#   (A) Absorption (Theorem 3, v0.3): TWFE bias on reported outcomes equals (D'MD)^{-1} D'M tau(s) exactly and is
#       invariant along the ENTIRE identified set of Theorem 1 — the affine direction AND the d-periodic directions
#       (here d = 4): both lie in col(W) because s = t + (1 - e) and, all entries being congruent mod d, s mod d is a
#       function of t alone.
#   (B) Identified correction (Theorem 4, v0.3): the corrected estimator is the coefficient on D in the JOINT
#       regression Y ~ D + unit FE + time FE + tenure dummies (aliased tenure columns dropped). No separate
#       calibration of tau from the outcome (v1 calibrated tau from cell means that contained
#       0.5 D, letting the genuine effect leak into the tenure component). Estimability condition: [MD, MZ_tau]
#       has full column rank after dropping aliased columns; fails when adoption is a deterministic function of
#       tenure (design "tenure-adoption" below).
#   (C) Spurious pre-trends (Theorem 5): exact shift tau(k + delta) - ell(k) under exact alignment; plus the
#       a rank-deficient example showing that the withdrawn parenthetical condition is not sufficient.
#   (D) Single-cohort protection.
#   (E) Monte Carlo over five designs: baseline; zero conditioning; heterogeneous treatment effects;
#       cohort-dependent adoption; tenure-determined adoption. Bias, RMSE, coverage (cluster-robust by unit).
# Output: sim4_results.csv, sim4_pretrend_coefs.csv, sim4_mc.csv, sim4_results.txt. Runs in ~3 min.
# ============================================================
suppressMessages({library(data.table); library(sandwich)})
set.seed(20260916)
out_dir <- if (dir.exists("sims")) "sims" else "."
tau_sat <- function(s) 0.30 * (1 - exp(-0.6 * (s - 1)))
tau_lin <- function(s) 0.05 * (s - 1)
tau_mix <- function(s) tau_sat(s) + 0.03 * (s - 1)
dgcd <- 4; rho4 <- function(s) c(0, 0.30, -0.10, 0.20)[((s - 1) %% dgcd) + 1]   # a 4-periodic null direction, rho(1)=0

build_panel <- function(E = c(1, 5, 13), Tmax = 19, n = 200, attrit = 0, adopt = "random", beta_sd = 0) {
  d <- rbindlist(lapply(seq_along(E), function(c) { e <- E[c]; CJ(id = paste0("c", c, "_", seq_len(n)), t = e:Tmax)[, e := e] }))
  d[, s := t - e + 1]; d[, u := rnorm(1), by = id]; lam <- rnorm(Tmax); d[, lam := lam[t]]
  if (adopt == "random") {            # adoption at a random calendar time after entry; 30% never
    d[, A := { if (runif(1) < .3) Inf else as.numeric(sample((e[1] + 1):Tmax, 1)) }, by = id]
  } else if (adopt == "cohort") {     # adoption probability and timing depend on entry cohort
    d[, A := { pnever <- c(.5, .3, .1)[match(e[1], E)]; if (runif(1) < pnever) Inf else as.numeric(sample((e[1] + 1):Tmax, 1)) }, by = id]
  } else if (adopt == "tenure") {     # adoption is a deterministic function of tenure: D = 1{s >= 4}
    d[, A := e + 3]
  }
  d[, D := as.integer(t >= A)]
  d[, beta_i := 0.5 + beta_sd * rnorm(1), by = id]
  d[, y := u + lam + beta_i * D + rnorm(.N, 0, 1)]
  if (attrit > 0) {                    # monotone attrition: per-wave exit hazard attrit*(s-1)/Tmax; once out, out
    setorder(d, id, t)
    d[, keep := { h <- attrit * (s - 1) / Tmax; u <- runif(.N); first_exit <- which(u < h)[1]
                  if (is.na(first_exit)) rep(TRUE, .N) else seq_len(.N) < first_exit }, by = id]
    d <- d[keep == TRUE]; d[, keep := NULL]
  }
  d[]
}
twfe <- function(y, D, id, t) coef(lm(y ~ D + factor(id) + factor(t)))["D"]
bias_formula <- function(v, D, id, t) { Dd <- resid(lm(D ~ factor(id) + factor(t))); sum(Dd * v) / sum(Dd^2) }
## corrected estimator: joint regression with tenure dummies (aliased ones dropped by lm); returns coef, cluster-robust se, rank info
corrected <- function(d, yvar = "ytil") {
  ## D enters LAST so that, when it is collinear with [W, Z_tau] (condition (R) fails), lm aliases D rather than a tenure dummy
  d <- copy(d); d[, fs := factor(s)]
  f <- lm(as.formula(paste(yvar, "~ factor(id) + factor(t) + fs + D")), data = d)
  cf <- coef(f); nd <- sum(is.na(cf[grep("^fs", names(cf))]))
  if (is.na(cf["D"])) return(c(est = NA_real_, se = NA_real_, n_tenure_dropped = nd))
  ## explicit residual-rank check of (R): residual of MD after projecting on the retained tenure indicators
  Dd <- resid(lm(D ~ factor(id) + factor(t), data = d)); Zr <- model.matrix(~ fs, d)[, -1, drop = FALSE]; Zr <- Zr[, !is.na(cf[colnames(Zr)]), drop = FALSE]
  Zd <- resid(lm(Zr ~ factor(id) + factor(t), data = d)); rr <- sum(resid(lm(Dd ~ Zd - 1))^2) / sum(Dd^2)
  if (rr < 1e-8) return(c(est = NA_real_, se = NA_real_, n_tenure_dropped = nd))
  V <- vcovCL(f, cluster = d$id); c(est = unname(cf["D"]), se = sqrt(V["D", "D"]), n_tenure_dropped = nd)
}
naive_se <- function(d, yvar = "ytil") { f <- lm(as.formula(paste(yvar, "~ D + factor(id) + factor(t)")), data = d); sqrt(vcovCL(f, cluster = d$id)["D", "D"]) }

res <- list(); txt <- c()
## ---- (A) absorption over the whole identified set ------------------------------------------------------
d <- build_panel(n = 200, attrit = 0.15)
for (nm in c("tau_sat", "tau_lin", "tau_mix")) {
  tf <- get(nm); d[, ytil := y + tf(s)]
  b_true <- twfe(d$y, d$D, d$id, d$t); b_rep <- twfe(d$ytil, d$D, d$id, d$t); bf <- bias_formula(tf(d$s), d$D, d$id, d$t)
  d[, ytil2 := y + tf(s) + 0.7 + 0.2 * (s - 1)]; b_rep2 <- twfe(d$ytil2, d$D, d$id, d$t)            # affine direction
  d[, ytil3 := y + tf(s) + rho4(s)];              b_rep3 <- twfe(d$ytil3, d$D, d$id, d$t)            # periodic direction
  res[[length(res) + 1]] <- data.table(block = "A_absorption", tau = nm, twfe_true = b_true, twfe_reported = b_rep,
    bias_observed = b_rep - b_true, bias_formula = bf, abs_diff_formula = abs((b_rep - b_true) - bf),
    invariance_gap_affine = abs(b_rep2 - b_rep), invariance_gap_periodic = abs(b_rep3 - b_rep))
}
## cell-mean design of this support: columns / rank / nullity (Appendix D count)
cm <- unique(d[, .(e, t, s)]); X <- cbind(model.matrix(~ factor(t), cm), model.matrix(~ factor(e), cm)[, -1, drop = FALSE], model.matrix(~ factor(s), cm)[, -1, drop = FALSE])
txt <- c(txt, sprintf("Cell-mean design (entries 1,5,13; T=19; tau(1)=g(1)=0): cells=%d columns=%d rank=%d nullity=%d", nrow(cm), ncol(X), qr(X)$rank, ncol(X) - qr(X)$rank))

## ---- (B) identified correction by joint regression ------------------------------------------------------
d[, ytil := y + tau_mix(s)]
b_true <- twfe(d$y, d$D, d$id, d$t); b_rep <- twfe(d$ytil, d$D, d$id, d$t); cc <- corrected(d)
res[[length(res) + 1]] <- data.table(block = "B_correction", tau = "tau_mix", twfe_true = b_true, twfe_reported = b_rep,
  bias_observed = b_rep - b_true, corrected_twfe = cc["est"], corrected_minus_true = cc["est"] - b_true, corrected_se = cc["se"],
  n_tenure_dropped = cc["n_tenure_dropped"])
## estimability failure: adoption a deterministic function of tenure
dt <- build_panel(n = 200, attrit = 0, adopt = "tenure"); dt[, ytil := y + tau_mix(s)]; cct <- corrected(dt)
res[[length(res) + 1]] <- data.table(block = "B_tenure_adoption", tau = "tau_mix", twfe_reported = twfe(dt$ytil, dt$D, dt$id, dt$t),
  corrected_twfe = cct["est"], corrected_se = cct["se"], n_tenure_dropped = cct["n_tenure_dropped"])
txt <- c(txt, sprintf("Tenure-determined adoption (D = 1{s>=4}): corrected coefficient on D = %s (NA = D aliased with tenure dummies: not identified)", format(cct["est"])))

## ---- (C) spurious pre-trends -------------------------------------------------------------------------------
dd <- build_panel(E = 1:8, Tmax = 20, n = 60, attrit = 0); delta <- 6; dd[, Ev := e + delta - 1]; dd[, k := t - Ev]; dd <- dd[k >= -5 & k <= 6]
dd[, D := NULL]; dd[, y := u + lam + rnorm(.N)]; dd[, ytil := y + tau_sat(s)]
ks <- sort(unique(dd$k)); k0 <- -1; k1 <- -2; kk <- setdiff(ks, c(k0, k1))
for (k in kk) dd[[paste0("ev_", k)]] <- as.integer(dd$k == k)
rhs <- paste(sprintf("`ev_%d`", kk), collapse = "+")
cf_rep <- coef(lm(as.formula(paste("ytil ~", rhs, "+ factor(id) + factor(t)")), dd)); cf_true <- coef(lm(as.formula(paste("y ~", rhs, "+ factor(id) + factor(t)")), dd))
names(cf_rep) <- gsub("`", "", names(cf_rep)); names(cf_true) <- gsub("`", "", names(cf_true)); stopifnot(!anyNA(cf_rep[paste0("ev_", kk)]))
get <- function(cf, k) unname(cf[paste0("ev_", k)]); ell <- function(k) tau_sat(k0 + delta) + (tau_sat(k1 + delta) - tau_sat(k0 + delta)) / (k1 - k0) * (k - k0)
pre <- rbindlist(lapply(kk, function(k) data.table(k = k, beta_reported = get(cf_rep, k), beta_true_outcome = get(cf_true, k), predicted_shift = tau_sat(k + delta) - ell(k))))
pre[, observed_shift := beta_reported - beta_true_outcome]; pre[, abs_err := abs(observed_shift - predicted_shift)]
res[[length(res) + 1]] <- data.table(block = "C_pretrend", tau = "tau_sat", n_coefs = nrow(pre), max_abs_err = max(pre$abs_err), max_pre_shift = max(abs(pre[k < 0, predicted_shift])), n_pre_coefs = sum(pre$k < 0))
## rank-deficient event-study design: E in {0,2,4}, t in 4..7, window k in 2..5, omit k = 2, 3 -> 8 columns, rank 7
rv <- CJ(E = c(0, 2, 4), t = 4:7)[, k := t - E][k >= 2 & k <= 5]; Xr <- cbind(model.matrix(~ factor(E), rv), model.matrix(~ factor(t), rv)[, -1], rv[, .(ev4 = as.integer(k == 4), ev5 = as.integer(k == 5))])
txt <- c(txt, sprintf("Rank-deficient event-study design (E in {0,2,4}, t 4..7, window 2..5): columns=%d rank=%d (all E even: k mod 2 = t mod 2)", ncol(as.matrix(Xr)), qr(as.matrix(Xr))$rank))
## (C2) exact alignment: X_k = Z_{k+delta}; the within-sample joint regression cannot separate them ((R) fails)
Mres <- function(V, id, t) { V <- as.matrix(V); apply(V, 2, function(v) resid(lm(v ~ factor(id) + factor(t)))) }
Xev <- as.matrix(dd[, paste0("ev_", kk), with = FALSE]); Zten <- model.matrix(~ factor(s), dd)[, -1, drop = FALSE]
MX <- Mres(Xev, dd$id, dd$t); MZ <- Mres(Zten, dd$id, dd$t)
r_MX <- qr(MX)$rank; r_MZ <- qr(MZ)$rank; r_MXZ <- qr(cbind(MX, MZ))$rank
txt <- c(txt, sprintf("Exact alignment: rank(MX)=%d, rank(MZ_tau)=%d, rank([MX, MZ_tau])=%d of %d columns -> (R) fails; no within-sample correction", r_MX, r_MZ, r_MXZ, ncol(MX) + ncol(MZ)))
res[[length(res) + 1]] <- data.table(block = "C2_exact_alignment_rank", tau = "tau_sat", rank_MX = r_MX, rank_MZ = r_MZ, rank_MX_MZ = r_MXZ, cols_MX_MZ = ncol(MX) + ncol(MZ))
## (C3) partial alignment design: E = 6..13, two units per date with delta in {6, 7}, k in [-5, 6],
##      references -1, -2; tenure composition of every event-time category is balanced across periods, yet the
##      coefficient shift (X'MX)^{-1} X'M tau(s) is nonzero. Deterministic (no noise): shift = coefficients of tau(s).
pa <- CJ(Ev = 6:13, delta = c(6, 7), k = -5:6); pa[, id := paste0("u", Ev, "_", delta)]; pa[, t := Ev + k]; pa[, s := k + delta]
for (k in kk) pa[[paste0("ev_", k)]] <- as.integer(pa$k == k)
cf_pa <- coef(lm(as.formula(paste("tau_sat(s) ~", rhs, "+ factor(id) + factor(t)")), pa)); names(cf_pa) <- gsub("`", "", names(cf_pa))
Xpa <- as.matrix(pa[, paste0("ev_", kk), with = FALSE]); Zpa <- model.matrix(~ factor(s), pa)[, -1, drop = FALSE]
MXpa <- Mres(Xpa, pa$id, pa$t); XMZ <- crossprod(MXpa, Zpa)
comp <- pa[, .(share_lo = mean(delta == 6)), by = .(k, t)]
txt <- c(txt, sprintf("Partial alignment (delta_i in {6,7}): composition balanced (share at delta=6 is %s in every (k,t)); rank([W0,X]) = %d of %d; shifts at k=-5,-4,-3: %.10f %.10f %.10f; max|X'M Z_tau| = %.3f",
  paste(unique(round(comp$share_lo, 3)), collapse = "/"), qr(cbind(model.matrix(~ factor(id) + factor(t), pa), Xpa))$rank, ncol(model.matrix(~ factor(id) + factor(t), pa)) + ncol(Xpa),
  cf_pa["ev_-5"], cf_pa["ev_-4"], cf_pa["ev_-3"], max(abs(XMZ))))
pa_out <- data.table(k = kk, shift = unname(cf_pa[paste0("ev_", kk)]))
res[[length(res) + 1]] <- data.table(block = "C3_partial_alignment", tau = "tau_sat", shift_k_minus5 = cf_pa["ev_-5"], shift_k_minus4 = cf_pa["ev_-4"], shift_k_minus3 = cf_pa["ev_-3"], max_abs_XMZ = max(abs(XMZ)))

## ---- (D) single-cohort protection ---------------------------------------------------------------------------
d1 <- build_panel(E = 1, Tmax = 19, n = 600, attrit = 0.15); d1[, ytil := y + tau_sat(s)]
res[[length(res) + 1]] <- data.table(block = "D_single_cohort", tau = "tau_sat", twfe_true = twfe(d1$y, d1$D, d1$id, d1$t), twfe_reported = twfe(d1$ytil, d1$D, d1$id, d1$t),
  bias_observed = twfe(d1$ytil, d1$D, d1$id, d1$t) - twfe(d1$y, d1$D, d1$id, d1$t))

## ---- (E) Monte Carlo over designs ------------------------------------------------------------------------------
B <- 100
designs <- list(baseline = list(tf = tau_mix, adopt = "random", beta_sd = 0), zero_conditioning = list(tf = function(s) 0 * s, adopt = "random", beta_sd = 0),
                heterogeneous_effect = list(tf = tau_mix, adopt = "random", beta_sd = 0.2), cohort_adoption = list(tf = tau_mix, adopt = "cohort", beta_sd = 0),
                tenure_adoption = list(tf = tau_mix, adopt = "tenure", beta_sd = 0))
mc <- rbindlist(lapply(names(designs), function(nm) { g <- designs[[nm]]
  rbindlist(lapply(seq_len(B), function(b) {
    d <- build_panel(n = 120, attrit = 0.15, adopt = g$adopt, beta_sd = g$beta_sd); d[, ytil := y + g$tf(s)]
    bt <- twfe(d$y, d$D, d$id, d$t); br <- twfe(d$ytil, d$D, d$id, d$t); ser <- naive_se(d); cc <- corrected(d)
    data.table(design = nm, b = b, twfe_true = bt, twfe_reported = br, se_reported = ser, twfe_corrected = cc["est"], se_corrected = cc["se"])
  })) }))
mcs <- mc[, .(B = .N, bias_reported = mean(twfe_reported - 0.5), rmse_reported = sqrt(mean((twfe_reported - 0.5)^2)),
              cover_reported = mean(abs(twfe_reported - 0.5) <= 1.96 * se_reported),
              bias_corrected = mean(twfe_corrected - 0.5, na.rm = TRUE), rmse_corrected = sqrt(mean((twfe_corrected - 0.5)^2, na.rm = TRUE)),
              cover_corrected = mean(abs(twfe_corrected - 0.5) <= 1.96 * se_corrected, na.rm = TRUE), share_corrected_NA = mean(is.na(twfe_corrected)),
              mcse_reported = sd(twfe_reported) / sqrt(.N)), by = design]
R <- rbindlist(res, fill = TRUE)
fwrite(R, file.path(out_dir, "sim4_results.csv")); fwrite(pre, file.path(out_dir, "sim4_pretrend_coefs.csv")); fwrite(mcs, file.path(out_dir, "sim4_mc.csv"))
## The transcript is written with base R's data.frame printer at a fixed, very wide console, so that its
## layout does not depend on the installed data.table version. data.table 1.15 and later add a type row
## under each header, and the column widths -- hence where a wide table wraps into a second block -- have
## changed between versions; either makes a token-by-token comparison of this file fail on another machine
## for no substantive reason. Every number in it is also in the CSVs written above, which are compared
## value by value by the release check and asserted by check_manuscript_values.R.
print_tab <- function(x) print(as.data.frame(x), digits = 4, row.names = FALSE)
old_width <- getOption("width"); options(width = 10000)
sink(file.path(out_dir, "sim4_results.txt"))
cat("Simulation 4 v2 (T1 downstream) — seed 20260916\n\n"); cat(paste(txt, collapse = "\n"), "\n\n"); print_tab(R)
cat("\nEvent-study coefficients (block C):\n"); print_tab(pre); cat("\nMonte Carlo by design (block E):\n"); print_tab(mcs)
sink()
options(width = old_width)
cat(paste(txt, collapse = "\n"), "\n\n"); print(R, digits = 4); print(pre, digits = 4); print(mcs, digits = 4)
## manuscript values (Appendix C), each from a named column
mv <- rbind(
  data.table(quantity = "bias_TWFE_tau_sat", value = R[block == "A_absorption" & tau == "tau_sat", bias_observed]),
  data.table(quantity = "max_abs_diff_formula", value = max(R[block == "A_absorption", abs_diff_formula])),
  data.table(quantity = "invariance_gap_affine_max", value = max(R[block == "A_absorption", invariance_gap_affine])),
  data.table(quantity = "invariance_gap_periodic_max", value = max(R[block == "A_absorption", invariance_gap_periodic])),
  data.table(quantity = "single_cohort_bias", value = R[block == "D_single_cohort", bias_observed]),
  data.table(quantity = "pretrend_max_abs_err", value = max(pre$abs_err)),
  data.table(quantity = paste0("pretrend_SHIFT_k", pre[k < 0, k]), value = pre[k < 0, predicted_shift]),
  data.table(quantity = paste0("pretrend_COEF_reported_k", pre[k < 0, k]), value = pre[k < 0, beta_reported]),
  data.table(quantity = paste0("pretrend_COEF_unconditioned_k", pre[k < 0, k]), value = pre[k < 0, beta_true_outcome]),
  data.table(quantity = paste0("partial_alignment_shift_k", pa_out$k), value = pa_out$shift))
fwrite(mv, file.path(out_dir, "sim4_manuscript_values.csv")); print(mv, digits = 6)

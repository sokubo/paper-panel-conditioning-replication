# ============================================================
# T1 Simulation 1 (v2, 2026-09-20): normalizations, identified functionals, plateau test
#   Design: JLPS-calibrated staggered support, entries E = (1, 5, 13), T = 19 -> d = gcd(4, 8, 12) = 4.
#   Under Theorem 1 (v0.3) the identified set of tau has dimension d = 4: tau + m(s-1) + rho(s),
#   rho 4-periodic with rho(1) = 0. Consequences tested here:
#   (A) NOISELESS IDENTITY TESTS: adding the periodic null direction leaves every cell mean unchanged,
#       moves the ordinary second difference D2 tau(3) (NOT identified when d > 1), and leaves the
#       lag-d second difference D4^2 tau(1) = tau(9) - 2 tau(5) + tau(1) unchanged (identified).
#   (B) MONTE CARLO: level estimates depend on the normalization; the lag-d curvature estimated by an
#       explicit contrast is unbiased; the plateau *observable implications* test (cohort-gap constancy
#       within the plateau window) has correct size, no power against pure linear drift (Theorem 1),
#       and power against nonlinear violations.
#   Spec (definitive; text follows code): cell noise sd 1/sqrt(1000) (n = 1000 respondents per cell,
#   outcome sd 1); S2: tau(s) = 0.2(1 - exp(-0.6(s-1))); S3 = S2 + 0.01(s-1); S4 = S2 + 0.004 max(s-6,0)^2;
#   plateau window s >= 6; B = 500; seed 20260825.
# ============================================================
set.seed(20260825)
B <- 500; TT <- 19; ENTRY <- c(1, 5, 13); SIG <- 1 / sqrt(1000); dgcd <- 4
tau_fn <- list(S1 = function(s) 0 * s,
               S2 = function(s) 0.2 * (1 - exp(-0.6 * (s - 1))),
               S3 = function(s) 0.2 * (1 - exp(-0.6 * (s - 1))) + 0.01 * (s - 1),
               S4 = function(s) 0.2 * (1 - exp(-0.6 * (s - 1))) + 0.004 * pmax(s - 6, 0)^2)
g_fn <- list(G0 = c(0, 0, 0), G1 = c(0, 0.15, 0))
cells <- do.call(rbind, lapply(seq_along(ENTRY), function(k)
  data.frame(e = k, t = ENTRY[k]:TT, s = (ENTRY[k]:TT) - ENTRY[k] + 1)))
cells$fs <- factor(cells$s, levels = 1:TT); cells$fe <- factor(cells$e); cells$ft <- factor(cells$t, levels = 1:TT)
get_mu <- function(df, e, t) df$mu[df$e == e & df$t == t]

## ---- (A) noiseless identity tests ------------------------------------------------------------------
cat("== (A) noiseless identity tests (d = 4) ==\n")
alpha0 <- 0.05 * sin(1:TT); g0 <- c(0, 0.10, 0.22); tau0 <- tau_fn$S2(1:TT)
mu0 <- alpha0[cells$t] + g0[cells$e] + tau0[cells$s]
## periodic null direction: rho(s) with period 4, rho(1) = 0; compensated by alpha(t) <- alpha(t) - rho(t - e0 + 1), e0 = 1
rho <- function(s) c(0, 0.30, -0.10, 0.20)[((s - 1) %% dgcd) + 1]
mu1 <- (alpha0[cells$t] - rho(cells$t)) + g0[cells$e] + (tau0[cells$s] + rho(cells$s))
cat(sprintf("  max |cell mean change| under tau -> tau + rho(s), alpha -> alpha - rho(t):  %.2e  (must be 0)\n", max(abs(mu1 - mu0))))
D2 <- function(tau, s) tau[s + 1] - 2 * tau[s] + tau[s - 1]                       # ordinary second difference at s
D2d <- function(tau, s) tau[s + 2 * dgcd] - 2 * tau[s + dgcd] + tau[s]            # lag-d second difference at s
cat(sprintf("  ordinary D2 tau(3):     before %+.4f  after %+.4f  -> moves by %+.4f (NOT identified)\n",
            D2(tau0, 3), D2(tau0 + rho(1:TT), 3), D2(tau0 + rho(1:TT), 3) - D2(tau0, 3)))
cat(sprintf("  lag-4 D4^2 tau(1):      before %+.4f  after %+.4f  -> moves by %+.1e (identified)\n",
            D2d(tau0, 1), D2d(tau0 + rho(1:TT), 1), D2d(tau0 + rho(1:TT), 1) - D2d(tau0, 1)))
## explicit-contrast estimator of the lag-d second difference from cell means, cancelling alpha and g:
## D(t) := mu(1,t) - mu(2,t) = tau(t) - tau(t-4) + g(1) - g(2)  for t >= 5;  D(t+4) - D(t) = tau(t+4) - 2 tau(t) + tau(t-4) = D4^2 tau(t-4)
contrast_D4 <- function(df, s) { t <- s + dgcd; (get_mu(df, 1, t + dgcd) - get_mu(df, 2, t + dgcd)) - (get_mu(df, 1, t) - get_mu(df, 2, t)) }
cat(sprintf("  explicit contrast for D4^2 tau(1) from noiseless cells: %+.4f  (truth %+.4f, diff %.1e)\n",
            contrast_D4(cbind(cells, mu = mu0), 1), D2d(tau0, 1), contrast_D4(cbind(cells, mu = mu0), 1) - D2d(tau0, 1)))
cat(sprintf("  same contrast after adding the periodic direction:       %+.4f  (diff %.1e)\n\n",
            contrast_D4(cbind(cells, mu = mu1), 1), contrast_D4(cbind(cells, mu = mu1), 1) - D2d(tau0, 1)))

## ---- (B) Monte Carlo ------------------------------------------------------------------------------
run_one <- function(tau_f, g_v) {
  alpha <- cumsum(rnorm(TT, 0, 0.05))
  df <- cells; df$mu <- alpha[df$t] + g_v[df$e] + tau_f(df$s) + rnorm(nrow(df), 0, SIG)
  ## N1: reference-level normalization (lm drops aliased columns: 4 of them here, one per null direction)
  cf <- coef(lm(mu ~ ft + fe + fs, data = df)); tn <- paste0("fs", 2:TT)
  tau1 <- unname(c(0, cf[tn])); tau1[is.na(tau1)] <- NA
  ## N2: zero-mean normalization along the affine direction only (another representative of the same fit)
  ok <- !is.na(tau1); cvec <- -sum(tau1[ok]) / sum((which(ok) - 1)); tau2 <- tau1 + cvec * (seq_len(TT) - 1)
  ## identified functional by explicit contrast
  d4 <- contrast_D4(df, 1)
  ## plateau (s >= 6): constrained fit gives tau_bar - tau(1) under the restriction; observable-implication test =
  ## cohort-gap constancy over time within the window (an identified functional: differences of lag-4 differences)
  df$fsc <- factor(pmin(df$s, 6)); cf3 <- coef(lm(mu ~ ft + fe + fsc, data = df)); tau_bar <- unname(cf3["fsc6"])
  dfp <- df[df$s >= 6, ]; dfp$ft <- droplevels(dfp$ft); dfp$fe <- droplevels(dfp$fe)
  av <- tryCatch(anova(lm(mu ~ ft + fe, dfp), lm(mu ~ ft + fe + fe:t, dfp)), error = function(e) NULL)
  p_over <- if (!is.null(av) && !is.na(av$`Pr(>F)`[2])) av$`Pr(>F)`[2] else NA_real_
  c(lev_N1 = tau1[5], lev_N2 = tau2[5], d4_contrast = d4, tau_bar = tau_bar, p_over = p_over)
}
res <- list()
for (ts in names(tau_fn)) for (gs in names(g_fn)) {
  out <- t(replicate(B, run_one(tau_fn[[ts]], g_fn[[gs]])))
  tf <- tau_fn[[ts]]; true5 <- tf(5) - tf(1); true_d4 <- tf(9) - 2 * tf(5) + tf(1); true_bar <- tf(10) - tf(1)
  res[[paste(ts, gs)]] <- data.frame(scenario = ts, cohort_g = gs,
    bias_lev_N1 = mean(out[, "lev_N1"], na.rm = TRUE) - true5, bias_lev_N2 = mean(out[, "lev_N2"], na.rm = TRUE) - true5,
    true_d4 = true_d4, bias_d4_contrast = mean(out[, "d4_contrast"]) - true_d4, sd_d4_contrast = sd(out[, "d4_contrast"]),
    bias_taubar_N3 = mean(out[, "tau_bar"]) - true_bar, rej_overid = mean(out[, "p_over"] < 0.05, na.rm = TRUE))
}
tab <- do.call(rbind, res); tab[, -(1:2)] <- round(tab[, -(1:2)], 4)
print(tab, row.names = FALSE)
write.csv(tab, "sim1_results.csv", row.names = FALSE)

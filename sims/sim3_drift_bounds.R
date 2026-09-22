# ============================================================
# T1 sim3 (v4, 2026-09-21): Proposition 3 (bounded cohort drift) — bounds validity & coverage
# v4: a third DGP (G4) with EQUAL true drifts (+M0, +M0). In G3 the two endpoints of the identified set for m
#     coincide at M = M0 but each is a min or max over DISTINCT values (a smooth function of the cell means);
#     G4 makes the max over cohort drifts a TIE, the non-smooth case, while the set for m stays an interval of
#     width 2*M0 - (max d - min d) = 2*M0 = .050 at M = M0 (truth at its upper endpoint; v0.8's header
#     wrote M0 in error). Neither DGP satisfies Imbens–Manski (2004) Assumption 1(iii) at a point-identified
#     parameter: the raw set width (s-1)(2M - |d_hat_1 - d_hat_2|) is continuously distributed, so the published
#     uniform-coverage lemma does not apply to the adaptive critical value as used here; the coverage numbers are
#     a finite-sample experiment, not a validity proof. Interval convention: the critical value c is computed
#     from the width TRUNCATED at zero; the interval itself is [lo - c*se_lo, hi + c*se_hi] with the RAW endpoints,
#     which may cross (an empty estimated set) — 'empty_rate' records how often.
# v3: the interval uses the Imbens–Manski (2004) critical value c solving
#     Phi(c + Delta_hat/sigma_hat) - Phi(-c) = 1 - alpha (c -> 1.645 for a wide set, c -> 1.96 at a point), instead of a
#     fixed z = 1.645; coverage with the fixed z is reported alongside for comparison. A second DGP (G3) with true
#     drifts (+M0, -M0) makes both constraints bind at M = M0, so the identified set for m is a point there.
# v2: under Theorem 1 (v0.3) the identified set on this support (d = 4) is tau + m(s-1) + rho(s) with rho
#     4-periodic, rho(1) = 0. The periodic part vanishes at s = 1 (mod 4), so the drift bounds of Prop. 3
#     apply to tau(s) - tau(1) exactly for s in {5, 9, 13, 17}; the linear part of a representative is
#     read off that subgrid. Evaluation tenures: 13 and 17 (v1 used 19, which is 3 mod 4: invalid).
# 目的: 修正版Prop.3(一階差分=ドリフト制約)の
#   (a) シャープ区間が真値を被覆するか(母集団レベルでの妥当性)
#   (b) 推定版+Imbens-Manski拡張の有限標本被覆
#   (c) 仕様検定(区間空)のサイズ/検出力
#   (d) 対照実験: 二階差分制約はcを制約しない(理論どおり全実数線)の確認
# 実行: Rscript sims/sim3_drift_bounds.R   (TP1ルートから)
# ============================================================
set.seed(20260830)

## --- design (JLPS calibration) ---------------------------------------------
E <- c(1, 5, 13)          # entry times (3 cohorts)
Tmax <- 19
n_cell <- 500             # respondents per cohort-wave cell
sdy <- 1
tau_true <- function(s) 0.25 * (1 - exp(-(s - 1) / 2))    # saturating, tau(1)=0
G_SETS <- list(G2 = c(0, 0.10, 0.22), G3 = c(0, 0.10, -0.10), G4 = c(0, 0.10, 0.30))   # drifts: G2 (.025,.015); G3 (+.025,-.025) coincident endpoints at M0; G4 (+.025,+.025) tied max
g_true <- G_SETS$G2                                        # g(e), g(1)=0 (population check uses G2)
alpha_true <- function(t) 0.05 * sin(t)
d_true <- diff(g_true) / diff(E)                          # per-period drifts .025, .015
M_true <- max(abs(d_true))

cells <- do.call(rbind, lapply(seq_along(E), function(k) {
  t <- E[k]:Tmax
  data.frame(e = E[k], t = t, s = t - E[k] + 1, ge = g_true[k])
}))

## --- estimator: representative triple + drift bounds ------------------------
fit_representative <- function(mu_hat) {
  d <- cells; d$mu <- mu_hat
  d$fs <- factor(d$s, levels = 1:Tmax)
  d$fe <- factor(d$e, levels = E)
  d$ft <- factor(d$t, levels = 1:Tmax)
  # saturated cell-mean model; lm drops one aliased column (the linear dependency)
  m <- lm(mu ~ ft + fe + fs, data = d)
  cf <- coef(m); cf[is.na(cf)] <- 0
  tau_hat <- c(0, cf[paste0("fs", 2:Tmax)]); tau_hat[is.na(tau_hat)] <- 0
  g_hat <- c(0, cf[paste0("fe", E[-1])]); g_hat[is.na(g_hat)] <- 0
  # re-normalize to zero-linear-part tau representative:
  s <- 1:Tmax
  sub <- which((s - 1) %% 4 == 0)                    # residue class of 1 mod d: periodic part is zero there
  chat <- coef(lm(tau_hat[sub] ~ I(s[sub] - 1)))[2]  # linear part of this solution, read on the valid subgrid
  tau_perp <- tau_hat - chat * (s - 1)
  g_perp <- g_hat + (-chat) * (E - E[1]) * (-1)     # g shifts with SAME c as tau
  # (tau,g) shift together: tau+c(s-1), g+c(e-e0). Removing chat from tau means
  # applying c = -chat, so g_perp = g_hat + (-chat)*(E - E[1]).
  g_perp <- g_hat + (-chat) * (E - E[1])
  list(tau_perp = as.numeric(tau_perp),
       d_perp = diff(as.numeric(g_perp)) / diff(E))
}

bounds_at <- function(rep_fit, M, s_eval) {
  c_lo <- -M - min(rep_fit$d_perp)
  c_hi <-  M - max(rep_fit$d_perp)
  list(lo = rep_fit$tau_perp[s_eval] + c_lo * (s_eval - 1),
       hi = rep_fit$tau_perp[s_eval] + c_hi * (s_eval - 1),
       empty = c_lo > c_hi)
}

## --- (a)+(d) population-level check ----------------------------------------
mu_pop <- alpha_true(cells$t) + cells$ge + tau_true(cells$s)
rep_pop <- fit_representative(mu_pop)
cat("== population check ==\n")
cat("d_perp (implied drifts of representative):", round(rep_pop$d_perp, 4), "\n")
for (s_eval in c(13, 17)) for (M in M_true * c(0.4, 1, 2, 4)) {
  b <- bounds_at(rep_pop, M, s_eval)
  cat(sprintf("s=%d M=%.4f: tau(s) in [%.4f, %.4f] empty=%s | truth %.4f inside=%s\n",
              s_eval, M, b$lo, b$hi, b$empty, tau_true(s_eval),
              !b$empty && tau_true(s_eval) >= b$lo - 1e-9 && tau_true(s_eval) <= b$hi + 1e-9))
}
cat("(identity) representative's tau_perp at s = 1 mod 4 vs truth minus its own linear part:\n")
sub <- which((1:Tmax - 1) %% 4 == 0); tt <- tau_true(1:Tmax); ct <- coef(lm(tt[sub] ~ I(sub - 1)))[2]
cat("   max abs diff on the subgrid:", format(max(abs(rep_pop$tau_perp[sub] - (tt[sub] - ct * (sub - 1)))), digits = 3), "\n")
cat("(d) second-difference bound leaves c unconstrained: the transformation\n",
    "   g+c(e-e0) is linear, so Delta^2 g is invariant — no finite interval exists.\n",
    "   [analytic fact; no computation required]\n\n")

## --- (b)+(c) finite-sample coverage ----------------------------------------
B_outer <- 300; B_boot <- 150; z_fixed <- 1.645
c_im <- function(width, sig, alpha = 0.05) {          # Imbens–Manski (2004) critical value
  if (!is.finite(sig) || sig <= 0) return(qnorm(1 - alpha / 2))
  r <- max(width, 0) / sig
  uniroot(function(c) pnorm(c + r) - pnorm(-c) - (1 - alpha), c(0.5, 4))$root
}
Ms <- M_true * c(0.4, 1, 2, 4)
all_summ <- list()
for (gname in names(G_SETS)) {
g_true <- G_SETS[[gname]]; cells$ge <- g_true[match(cells$e, E)]
d_true <- diff(g_true) / diff(E); stopifnot(abs(max(abs(d_true)) - M_true) < 1e-12)
mu_pop <- alpha_true(cells$t) + cells$ge + tau_true(cells$s)
for (s_eval in c(13, 17)) {
res <- array(NA_real_, dim = c(B_outer, length(Ms), 6),
             dimnames = list(NULL, paste0("M", c(0.4, 1, 2, 4)),
                             c("cover", "width", "empty", "cover_fixed", "c_used", "set_width")))
for (b_ in 1:B_outer) {
  mu_hat <- mu_pop + rnorm(nrow(cells), 0, sdy / sqrt(n_cell))
  rf <- fit_representative(mu_hat)
  # parametric bootstrap for endpoint SEs
  lo_b <- hi_b <- matrix(NA_real_, B_boot, length(Ms))
  for (r in 1:B_boot) {
    mu_r <- mu_hat + rnorm(nrow(cells), 0, sdy / sqrt(n_cell))
    rf_r <- fit_representative(mu_r)
    for (m_ in seq_along(Ms)) {
      br <- bounds_at(rf_r, Ms[m_], s_eval)
      lo_b[r, m_] <- br$lo; hi_b[r, m_] <- br$hi
    }
  }
  for (m_ in seq_along(Ms)) {
    b0 <- bounds_at(rf, Ms[m_], s_eval)
    res[b_, m_, "empty"] <- as.numeric(b0$empty)
    se_lo <- sd(lo_b[, m_]); se_hi <- sd(hi_b[, m_]); cc <- c_im(b0$hi - b0$lo, max(se_lo, se_hi))
    lo_im <- b0$lo - cc * se_lo; hi_im <- b0$hi + cc * se_hi
    res[b_, m_, "cover"] <- as.numeric(tau_true(s_eval) >= lo_im & tau_true(s_eval) <= hi_im)
    res[b_, m_, "width"] <- hi_im - lo_im; res[b_, m_, "c_used"] <- cc; res[b_, m_, "set_width"] <- b0$hi - b0$lo
    lo_f <- b0$lo - z_fixed * se_lo; hi_f <- b0$hi + z_fixed * se_hi
    res[b_, m_, "cover_fixed"] <- as.numeric(tau_true(s_eval) >= lo_f & tau_true(s_eval) <= hi_f)
  }
}
summ <- data.frame(dgp = gname, s_eval = s_eval,
  M_mult = c(0.4, 1, 2, 4), M = Ms,
  coverage = apply(res[, , "cover"], 2, mean, na.rm = TRUE),
  coverage_fixed_z = apply(res[, , "cover_fixed"], 2, mean, na.rm = TRUE),
  med_width = apply(res[, , "width"], 2, median, na.rm = TRUE),
  med_set_width = apply(res[, , "set_width"], 2, median, na.rm = TRUE),
  med_c = apply(res[, , "c_used"], 2, median, na.rm = TRUE),
  empty_rate = apply(res[, , "empty"], 2, mean, na.rm = TRUE))
cat("== finite-sample (B=300, n_cell=500, IM 95%) —", gname, "tau(", s_eval, "), truth", round(tau_true(s_eval), 3), "==\n")
print(summ, row.names = FALSE); all_summ[[length(all_summ) + 1]] <- summ
}
}
out <- do.call(rbind, all_summ)
write.csv(out, if (dir.exists("sims")) "sims/sim3_results.csv" else "sim3_results.csv", row.names = FALSE)
cat("\nsaved sim3_results.csv\n")

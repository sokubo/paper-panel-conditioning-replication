# ============================================================
# T1 Simulation 2 (v2, 2026-09-20): identified vs unidentified curvature readouts
#   v1 targeted "D2 tau(3)" and mis-centred its contrast (corrected in revision): the implemented estimator was
#   -D2D(6) = D2 tau(2) - D2 tau(6), not D2 tau(3). Under Theorem 1 (v0.3) with d = 4, D2 tau(3) is NOT
#   identified on this support, while theta := D2 tau(2) - D2 tau(6) IS (its weights sum to zero within
#   every residue class mod 4 and are orthogonal to s - 1). The v1 "boundary-dummy contamination" story is
#   withdrawn: a saturated regression returns an arbitrary representative of the identified set, so reading
#   an unidentified functional off it is arbitrary, and reading an identified functional off it is exact.
#   (A) noiseless identity tests; (B) Monte Carlo: bias of theta by explicit contrast and by regression
#   readout (both unbiased), and the spread of the "D2 tau(3)" readout across representatives.
#   Spec as sim1 (definitive): SIG = 1/sqrt(1000), S2 = 0.2(1 - exp(-0.6(s-1))), S3 = S2 + 0.01(s-1), B = 500.
# ============================================================
set.seed(20260825)
B <- 500; TT <- 19; ENTRY <- c(1, 5, 13); SIG <- 1 / sqrt(1000); dgcd <- 4
tau_fn <- list(S2 = function(s) 0.2 * (1 - exp(-0.6 * (s - 1))),
               S3 = function(s) 0.2 * (1 - exp(-0.6 * (s - 1))) + 0.01 * (s - 1))
cells <- do.call(rbind, lapply(seq_along(ENTRY), function(k)
  data.frame(e = k, t = ENTRY[k]:TT, s = (ENTRY[k]:TT) - ENTRY[k] + 1)))
cells$fs <- factor(cells$s, levels = 1:TT); cells$fe <- factor(cells$e); cells$ft <- factor(cells$t, levels = 1:TT)
get_mu <- function(df, e, t) df$mu[df$e == e & df$t == t]
D2 <- function(v, s) v[s + 1] - 2 * v[s] + v[s - 1]
## explicit contrast: D(t) = mu(1,t) - mu(2,t) = tau(t) - tau(t-4) + const;  -D2D(6) = D2 tau(2) - D2 tau(6)
theta_contrast <- function(df) { D <- sapply(5:7, function(t) get_mu(df, 1, t) - get_mu(df, 2, t)); -(D[3] - 2 * D[2] + D[1]) }
theta_true <- function(tf) D2(tf(1:TT), 2) - D2(tf(1:TT), 6)
fit_tau <- function(df) { cf <- coef(lm(mu ~ ft + fe + fs, data = df)); v <- unname(c(0, cf[paste0("fs", 2:TT)])); v[is.na(v)] <- 0; v }

cat("== (A) noiseless identity tests ==\n")
alpha0 <- 0.05 * sin(1:TT); g0 <- c(0, 0.15, 0); tf <- tau_fn$S2
df0 <- cells; df0$mu <- alpha0[df0$t] + g0[df0$e] + tf(df0$s)
cat(sprintf("  theta = D2 tau(2) - D2 tau(6): truth %+.5f | explicit contrast %+.5f | diff %.1e\n",
            theta_true(tf), theta_contrast(df0), theta_contrast(df0) - theta_true(tf)))
th <- fit_tau(df0)
cat(sprintf("  regression readout of theta from an arbitrary representative: %+.5f (diff %.1e)  <- identified: exact\n",
            D2(th, 2) - D2(th, 6), D2(th, 2) - D2(th, 6) - theta_true(tf)))
cat(sprintf("  regression readout of D2 tau(3): %+.5f vs truth %+.5f  <- NOT identified: representative-dependent\n",
            D2(th, 3), D2(tf(1:TT), 3)))
rho <- function(s) c(0, 0.30, -0.10, 0.20)[((s - 1) %% dgcd) + 1]
df1 <- cells; df1$mu <- (alpha0[df1$t] - rho(df1$t)) + g0[df1$e] + tf(df1$s) + rho(df1$s)
cat(sprintf("  after adding the periodic null direction: cells change by %.1e; contrast theta %+.5f (unchanged)\n\n",
            max(abs(df1$mu - df0$mu)), theta_contrast(df1)))

cat("== (B) Monte Carlo ==\n")
run_one <- function(tf, g_v = c(0, 0.15, 0)) {
  alpha <- cumsum(rnorm(TT, 0, 0.05)); df <- cells
  df$mu <- alpha[df$t] + g_v[df$e] + tf(df$s) + rnorm(nrow(df), 0, SIG)
  th <- fit_tau(df)
  c(theta_contrast = theta_contrast(df), theta_readout = D2(th, 2) - D2(th, 6), d2_3_readout = D2(th, 3))
}
res <- list()
for (ts in names(tau_fn)) {
  tf <- tau_fn[[ts]]; out <- t(replicate(B, run_one(tf)))
  res[[ts]] <- data.frame(scenario = ts, theta_true = theta_true(tf),
    bias_contrast = mean(out[, "theta_contrast"]) - theta_true(tf), sd_contrast = sd(out[, "theta_contrast"]),
    bias_readout = mean(out[, "theta_readout"]) - theta_true(tf), sd_readout = sd(out[, "theta_readout"]),
    d2_3_true_unidentified = D2(tf(1:TT), 3), mean_d2_3_readout = mean(out[, "d2_3_readout"]))
}
tab <- do.call(rbind, res); tab[, -1] <- round(tab[, -1], 4); print(tab, row.names = FALSE)
write.csv(tab, "sim2_results.csv", row.names = FALSE)

# ============================================================
# check_manuscript_values.R — reproduces every number quoted in Appendices A and C of T1 (v0.5) from the NAMED
# columns of the simulation output files, so that prose and outputs cannot drift apart (a safeguard added in revision).
# Run from sims/ after sim1–sim4:  Rscript check_manuscript_values.R
# ============================================================
options(width = 120)
s1 <- read.csv("sim1_results.csv"); s2 <- read.csv("sim2_results.csv"); s3 <- read.csv("sim3_results.csv")
s4 <- read.csv("sim4_mc.csv"); mv <- read.csv("sim4_manuscript_values.csv")
tau2 <- function(s) 0.2 * (1 - exp(-0.6 * (s - 1)))
cat("== Appendix A.1 (Simulation 1) ==\n")
cat(sprintf("  Delta_4^2 tau(1) under S2 (closed form): %.4f ; [mu(1,9)-mu(5,9)]-[mu(1,5)-mu(5,5)] = %.4f ; Delta_4^2 tau(5) = %.4f\n",
            tau2(9) - 2 * tau2(5) + tau2(1), (tau2(9) - tau2(5)) - (tau2(5) - tau2(1)), tau2(13) - 2 * tau2(9) + tau2(5)))
tab1 <- s1[, c("scenario", "cohort_g", "bias_lev_N1", "bias_lev_N2", "true_d4", "bias_d4_contrast", "sd_d4_contrast", "bias_taubar_N3", "rej_overid")]
print(tab1, row.names = FALSE, digits = 3)
cat("\n== Appendix A.2 (Simulation 2) ==\n"); print(s2, row.names = FALSE, digits = 3)
cat("\n== Appendix A.3 (Simulation 3) ==\n"); print(s3, row.names = FALSE, digits = 3)
cat("\n== Appendix C (Simulation 4) — deterministic identities and named values ==\n"); print(mv, row.names = FALSE, digits = 6)
cat("\n== Appendix C — Monte Carlo table ==\n")
print(s4[, c("design", "bias_reported", "rmse_reported", "cover_reported", "bias_corrected", "rmse_corrected", "cover_corrected", "share_corrected_NA", "mcse_reported")], row.names = FALSE, digits = 3)

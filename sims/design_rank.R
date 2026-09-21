# ============================================================
# design_rank(): what does a staggered panel design identify about the conditioning path?
#   Builds the cell-mean design X_S for mu(e,t) = alpha(t) + g(e) + tau(s), s = t - e + 1,
#   with normalisations tau(1) = 0, g(e_1) = 0, on the support S = {(e,t): e in entries, e <= t <= Tmax}
#   (or an arbitrary support given as a data.frame of (e,t) cells), and reports:
#     - cells, free parameters, rank, nullity  (Theorem 1(d) of T1)
#     - the stride d = gcd of entry spacings and whether nullity == d
#     - Lemma 1 checks (trapezoid): w = Tmax - e_K, the sufficient condition w >= (e_2 - e_1) - d,
#       and, for three cohorts, the exact value max(d, e_2 - e_1 - w)
#     - the increment graph (vertices = observed increments u with (u, u+1) observed for some cohort; edges between
#       increments of two cohorts observed at the same t and t+1), its number of connected components, and whether
#       condition C_d holds (components = residue classes mod d; increment-complete).  For a staggered trapezoid the
#       number of components equals the nullity (Lemma 1(o) of T1 v0.5).
#     - a basis of the tau-part of the null space (the directions of the identified set)
#     - a test of whether a given linear functional lambda'tau is identified (orthogonal to the null space)
#   Usage:  source("sims/design_rank.R"); r <- design_rank(c(1, 5, 13), 19); r$summary
#           is_identified(r, lambda = c(1, -2, 1, rep(0, 16)))          # ordinary D2 tau(2): FALSE when d = 4
#           is_identified(r, lambda = replace(numeric(19), c(1, 5, 9), c(1, -2, 1)))  # lag-4 D4^2 tau(1): TRUE
# ============================================================
design_rank <- function(entries, Tmax = NULL, cells = NULL, tol = 1e-9) {
  if (is.null(cells)) {
    stopifnot(!is.null(Tmax), all(diff(entries) > 0), Tmax >= max(entries))
    cells <- do.call(rbind, lapply(entries, function(e) data.frame(e = e, t = e:Tmax)))
  }
  cells <- unique(cells[order(cells$e, cells$t), c("e", "t")]); cells$s <- cells$t - cells$e + 1
  E <- sort(unique(cells$e)); Tt <- sort(unique(cells$t)); S <- sort(unique(cells$s)); e0 <- E[1]
  cols <- c(paste0("a", Tt), if (length(E) > 1) paste0("g", E[-1]), paste0("u", S[S != 1]))
  X <- matrix(0, nrow(cells), length(cols), dimnames = list(NULL, cols))
  for (r in seq_len(nrow(cells))) {
    X[r, paste0("a", cells$t[r])] <- 1
    if (cells$e[r] != e0) X[r, paste0("g", cells$e[r])] <- 1
    if (cells$s[r] != 1)  X[r, paste0("u", cells$s[r])] <- 1
  }
  qrX <- qr(X, tol = tol); rk <- qrX$rank; nullity <- ncol(X) - rk
  d <- if (length(E) > 1) Reduce(function(a, b) { while (b) { t <- b; b <- a %% b; a <- t }; a }, abs(diff(E))) else NA_integer_
  ## null-space basis via SVD (right singular vectors for zero singular values)
  sv <- svd(X, nu = 0, nv = ncol(X)); Nb <- if (nullity > 0) sv$v[, (rk + 1):ncol(X), drop = FALSE] else matrix(0, ncol(X), 0)
  tau_rows <- grep("^u", cols); tau_basis <- matrix(0, length(S), ncol(Nb), dimnames = list(paste0("tau", S), NULL))
  if (ncol(Nb) > 0) tau_basis[S != 1, ] <- Nb[tau_rows, , drop = FALSE]
  w <- max(Tt) - max(E); D2 <- if (length(E) > 1) E[2] - E[1] else NA
  ig <- increment_graph(cells)
  summary <- data.frame(cells = nrow(cells), parameters = ncol(X), rank = rk, nullity = nullity, stride_d = d,
                        nullity_equals_d = nullity == d, increment_components = ig$n_components, increment_complete = ig$complete,
                        C_d = ig$complete && !is.na(d) && ig$n_components == d && all(sapply(ig$components, function(cc) length(unique(cc %% d)) == 1)),
                        w_last_cohort = w, lemma1_condition = if (!is.na(D2)) w >= D2 - d else NA,
                        three_cohort_formula = if (length(E) == 3) max(d, D2 - w) else NA)
  structure(list(summary = summary, tenures = S, tau_null_basis = tau_basis, X = X, cells = cells, cols = cols, increment_graph = ig), class = "design_rank")
}
## increment graph of a support (any set of (e,t) cells)
increment_graph <- function(cells) {
  key <- paste(cells$e, cells$t); has <- function(e, t) paste(e, t) %in% key
  inc <- unique(cells$s[has(cells$e, cells$t + 1)])                      # observed increments u: (u, u+1) both observed
  edges <- list()
  for (r in seq_len(nrow(cells))) { e <- cells$e[r]; t <- cells$t[r]
    if (!has(e, t + 1)) next
    others <- unique(cells$e[cells$t == t & cells$e != e]); others <- others[has(others, t + 1)]
    for (ep in others) edges[[length(edges) + 1]] <- c(t - e + 1, t - ep + 1) }
  comp <- setNames(seq_along(inc), inc)                                    # union-find by relabelling
  for (ed in edges) { a <- comp[as.character(ed[1])]; b <- comp[as.character(ed[2])]; if (a != b) comp[comp == b] <- a }
  comps <- split(as.integer(names(comp)), comp)
  S <- sort(unique(cells$s)); complete <- all(sapply(S, function(s) all(seq_len(s - 1) %in% inc)))
  list(increments = sort(inc), components = unname(comps), n_components = length(comps), complete = complete)
}
print.design_rank <- function(x, ...) { print(x$summary, row.names = FALSE); invisible(x) }
## is a linear functional lambda'tau (lambda indexed by the observed tenures, in order) identified?
is_identified <- function(dr, lambda, tol = 1e-8) {
  stopifnot(length(lambda) == length(dr$tenures))
  if (ncol(dr$tau_null_basis) == 0) return(TRUE)
  max(abs(crossprod(dr$tau_null_basis, lambda))) < tol
}
if (sys.nframe() == 0) {
  r <- design_rank(c(1, 5, 13), 19); print(r)
  lam_d2 <- numeric(19); lam_d2[c(2, 3, 4)] <- c(1, -2, 1)          # ordinary second difference at s = 3
  lam_d4 <- numeric(19); lam_d4[c(1, 5, 9)] <- c(1, -2, 1)          # lag-4 second difference at s = 1
  lam_th <- numeric(19); lam_th[c(1, 2, 3)] <- c(1, -2, 1); lam_th[c(5, 6, 7)] <- lam_th[c(5, 6, 7)] - c(1, -2, 1)  # D2 tau(2) - D2 tau(6)
  cat("identified?  D2 tau(3):", is_identified(r, lam_d2), " | D4^2 tau(1):", is_identified(r, lam_d4), " | D2tau(2)-D2tau(6):", is_identified(r, lam_th), "\n")
  for (E in list(c(1, 2, 3), c(1, 3, 5), c(1, 7, 10), c(1, 6, 9), c(1, 4, 8), c(1, 6, 13))) for (Tm in c(max(E), max(E) + 3, max(E) + 8)) {
    rr <- design_rank(E, Tm)$summary
    cat(sprintf("entries %-12s T=%2d: nullity=%d d=%d components=%d C_d=%s lemma1=%s three-cohort formula=%s\n", paste(E, collapse = ","), Tm,
                rr$nullity, rr$stride_d, rr$increment_components, rr$C_d, rr$lemma1_condition, rr$three_cohort_formula)) }
  cat("\nDesign {1,4,8}, T=8 (d = 1 but C_1 fails):\n"); print(design_rank(c(1, 4, 8), 8)); print(design_rank(c(1, 4, 8), 8)$increment_graph$components)
}

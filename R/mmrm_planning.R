# Sample size for the MMRM comparator (hypothetical estimand under MAR).
#
# With monotone dropout and a correctly specified MMRM fitted by likelihood,
# the variance of the estimated final-visit mean in an arm is
#   sigma^2 / n * sum_j tau_j / r_j,
# where r_j is the share of the arm still observed at visit j and tau_j is the
# share of the final-visit variance first explained by visit j:
#   tau_j = V_j - V_{j-1},  V_j = Var(E[Y_K | Y_1..Y_j]) / Var(Y_K),  V_0 = 0.
# Each increment of information is estimated from the patients still observed
# at that visit (Lu, Luo and Chen 2008). With r_j = 1 for all j this is the
# complete-data variance; with no correlation between visits it is the
# completers-only variance.

# Correlation of PANSS total change between weeks 1 to 6, pooled within arm and
# trial over the 16 six-week placebo-controlled NEWMEDS trials (5,157 patients,
# pairwise complete observations).
.newmeds_corr <- matrix(c(
  1.00, 0.72, 0.59, 0.54, 0.49, 0.51,
  0.72, 1.00, 0.79, 0.69, 0.66, 0.64,
  0.59, 0.79, 1.00, 0.83, 0.76, 0.73,
  0.54, 0.69, 0.83, 1.00, 0.86, 0.83,
  0.49, 0.66, 0.76, 0.86, 1.00, 0.90,
  0.51, 0.64, 0.73, 0.83, 0.90, 1.00), 6, 6)

.mmrm_corr <- function(visits, structure, rho) {
  if (structure == "newmeds") {
    if (visits != 6)
      stop("The NEWMEDS correlation matrix is for 6 weekly visits; use ",
           "structure = \"ar1\" or \"cs\" for another number of visits.")
    return(.newmeds_corr)
  }
  if (is.null(rho) || rho < 0 || rho >= 1)
    stop("'rho' must be in [0, 1) for structure = \"", structure, "\".")
  i <- seq_len(visits)
  if (structure == "ar1") return(rho^abs(outer(i, i, "-")))
  R <- matrix(rho, visits, visits); diag(R) <- 1; R
}

.mmrm_retention <- function(p, visits, pattern) {
  f <- seq_len(visits) / visits
  f <- switch(pattern, even = f, early = 1 - (1 - f)^2, late = f^2)
  1 - (1 - p) * f
}

.mmrm_phi <- function(R, ret) {
  K <- nrow(R)
  V <- vapply(seq_len(K), function(j) {
    s <- R[K, seq_len(j)]
    drop(s %*% solve(R[seq_len(j), seq_len(j), drop = FALSE], s))
  }, numeric(1))
  tau <- diff(c(0, V))
  sum(tau / ret)
}

cs_mmrm_design <- function(d, pA, pB, visits = 6,
                           structure = c("newmeds", "ar1", "cs"), rho = NULL,
                           pattern = c("even", "early", "late"), r = 1,
                           alpha = 0.025) {
  structure <- match.arg(structure); pattern <- match.arg(pattern)
  .check_prob(pA, "pA"); .check_prob(pB, "pB")
  if (!is.numeric(d) || length(d) != 1 || d == 0) stop("'d' must be a non-zero number.")
  if (visits < 1 || visits != round(visits)) stop("'visits' must be a positive integer.")
  if (r <= 0) stop("'r' must be positive.")
  if (alpha <= 0 || alpha >= 0.5) stop("'alpha' must be in (0, 0.5).")
  R <- .mmrm_corr(visits, structure, rho)
  phiA <- .mmrm_phi(R, .mmrm_retention(pA, visits, pattern))
  phiB <- .mmrm_phi(R, .mmrm_retention(pB, visits, pattern))
  structure(list(d = d, pA = pA, pB = pB, visits = visits, structure = structure,
                 rho = rho, pattern = pattern, r = r, alpha = alpha,
                 phiA = phiA, phiB = phiB),
            class = "cs_mmrm_design")
}

cs_mmrm_power <- function(design, nA) {
  with(design, pnorm(abs(d) * sqrt(nA) / sqrt(phiA + phiB / r) - qnorm(1 - alpha)))
}

cs_mmrm_sample_size <- function(design, power = 0.90) {
  with(design, {
    nA <- ceiling((qnorm(1 - alpha) + qnorm(power))^2 * (phiA + phiB / r) / d^2)
    nB <- ceiling(r * nA)
    list(nA = nA, nB = nB, total = nA + nB, power = cs_mmrm_power(design, nA),
         phiA = phiA, phiB = phiB)
  })
}

print.cs_mmrm_design <- function(x, ...) {
  cat("MMRM comparator (hypothetical estimand, assumes missing at random)\n")
  cat(sprintf("  Effect size at the final visit: d = %g; completion: pA = %.3f, pB = %.3f\n",
              x$d, x$pA, x$pB))
  cat(sprintf("  %d visits, correlation %s%s, dropout %s over visits\n",
              x$visits, x$structure,
              if (x$structure == "newmeds") "" else sprintf(" (rho = %g)", x$rho),
              x$pattern))
  cat(sprintf("  Variance inflation over complete data: A %.3f, B %.3f\n",
              x$phiA, x$phiB))
  invisible(x)
}

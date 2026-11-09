# Simulation check of cs_mmrm_sample_size().
#
# With monotone dropout, the likelihood of an MMRM with unstructured covariance
# factors into a sequence of regressions (Y_1, then Y_j on Y_1..Y_{j-1} among
# patients still observed at visit j), and the ML estimate of each arm's
# final-visit mean is obtained from them in closed form. That makes thousands of
# simulated trials fast. A small number of trials is also fitted with the full
# MMRM (cs_mmrm(), nlme REML) to show that the two estimates agree.
#
# Checks: (1) the empirical SD of the estimated difference against the formula
# SD; (2) power, using the formula SE as the denominator; (3) agreement between
# the sequential-regression and full MMRM estimates and SEs.
# Run from the package root: Rscript inst/validation/mmrm_sample_size_check.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(20261107)
K <- 6; d <- 0.40; pA <- 0.65; pB <- 0.55
des <- cs_mmrm_design(d, pA, pB, visits = K, structure = "newmeds", pattern = "even")
ss <- cs_mmrm_sample_size(des, power = 0.90); n <- ss$nA
se_formula <- sqrt((des$phiA + des$phiB) / n)
cat(sprintf("Computed size: %d per arm (total %d); formula SE of the difference %.4f; analytic power %.3f\n",
            n, ss$total, se_formula, ss$power))
L <- chol(compositeSS:::.newmeds_corr)

one_arm <- function(n, mean_k, p, mar) {
  Y <- matrix(rnorm(n * K), n, K) %*% L
  Y <- sweep(Y, 2, mean_k * seq_len(K) / K, "+")
  ret <- compositeSS:::.mmrm_retention(p, K, "even")
  prev <- c(1, ret[-K])
  last <- rep(K, n); alive <- rep(TRUE, n)
  for (j in seq_len(K)) {
    q <- 1 - ret[j] / prev[j]
    lp <- qlogis(q) + if (mar && j > 1) -0.8 * (Y[, j - 1] - mean_k * (j - 1) / K) else 0
    drop <- alive & (runif(n) < plogis(lp))
    last[drop] <- j - 1; alive <- alive & !drop
  }
  list(Y = Y, last = last)
}

# ML estimate of the final-visit mean from the factored likelihood
seq_mean <- function(a) {
  mu <- numeric(K)
  o <- a$last >= 1
  mu[1] <- mean(a$Y[o, 1])
  for (j in 2:K) {
    o <- a$last >= j
    b <- lm.fit(cbind(1, a$Y[o, 1:(j - 1), drop = FALSE]), a$Y[o, j])$coefficients
    mu[j] <- b[1] + sum(b[-1] * mu[1:(j - 1)])
  }
  mu[K]
}

to_long <- function(a, b) {
  f <- function(x, arm, off) do.call(rbind, lapply(which(x$last >= 1), function(i)
    data.frame(id = off + i, visit = seq_len(x$last[i]), arm = arm,
               y = x$Y[i, seq_len(x$last[i])])))
  rbind(f(a, "A", 0), f(b, "B", n))
}

run <- function(reps, mar, full = 0) {
  est <- numeric(reps); cA <- cB <- numeric(reps); cmp <- NULL
  for (r in seq_len(reps)) {
    a <- one_arm(n, d, pA, mar); b <- one_arm(n, 0, pB, mar)
    est[r] <- seq_mean(a) - seq_mean(b)
    cA[r] <- mean(a$last == K); cB[r] <- mean(b$last == K)
    if (r <= full) {
      f <- cs_mmrm(to_long(a, b), k = K, active = "A")
      cmp <- rbind(cmp, c(seq = est[r], mmrm = f$estimate, se = f$se))
    }
  }
  pw <- mean(est / se_formula > qnorm(0.975))
  cat(sprintf("\n%s, %d trials: completion A %.3f, B %.3f\n",
              if (mar) "Dropout depends on the last observed value (MAR)" else "Dropout at random (MCAR)",
              reps, mean(cA), mean(cB)))
  cat(sprintf("  mean estimate %.4f (true %.2f); empirical SD %.4f vs formula %.4f (ratio %.3f)\n",
              mean(est), d, sd(est), se_formula, sd(est) / se_formula))
  cat(sprintf("  power with the formula SE: %.3f (Monte Carlo SE %.3f)\n",
              pw, sqrt(pw * (1 - pw) / reps)))
  if (!is.null(cmp))
    cat(sprintf("  full MMRM on the first %d trials: largest |difference| between estimates %.4f; mean MMRM SE %.4f vs formula %.4f\n",
                nrow(cmp), max(abs(cmp[, "seq"] - cmp[, "mmrm"])), mean(cmp[, "se"]), se_formula))
}
run(4000, mar = FALSE, full = 20)
run(4000, mar = TRUE, full = 20)

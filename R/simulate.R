cs_sim_trial <- function(n, p, mu, sigma) {
  comp <- as.logical(rbinom(n, 1, p))
  y <- rep(NA_real_, n)
  y[comp] <- rnorm(sum(comp), mu, sigma)
  list(comp = comp, y = y)
}

cs_simulate <- function(design, n, nsim = 10000,
                        dropout = c("random", "appendixC"), model = NULL,
                        n_cal = 1e5,
                        null = FALSE, muB = NULL, s = 0, conf_level = 0.95,
                        comparators = character(0), M = 20,
                        mmrm_engine = c("nlme", "mmrm"), progress = NULL,
                        seed = NULL) {
  dropout <- match.arg(dropout)
  mmrm_engine <- match.arg(mmrm_engine)
  if (length(comparators))
    comparators <- match.arg(comparators, c("mmrm", "pmi"), several.ok = TRUE)
  if (length(comparators) && dropout != "appendixC")
    stop("Comparators need visit-level data: use dropout = 'appendixC'.")
  if (design$efficacy != "continuous" || design$margin_d > 0 ||
      design$margin_e > 0)
    stop("Simulation currently supports continuous efficacy, superiority form.")
  if (is.null(seed)) seed <- sample.int(1e8, 1)
  set.seed(seed)
  pA <- design$pA; pB <- if (null) design$pA else design$pB
  delta <- if (null) 0 else design$delta
  sigma <- design$sigma
  mB0 <- if (is.null(muB)) 0 else muB
  muA0 <- mB0 + delta
  nA <- n; nB <- ceiling(design$r * n)
  crit <- cs_critical(design)
  w <- design$weights; zq <- qnorm(1 - (1 - conf_level) / 2)

  policy_effect <- hyp_effect <- NA
  if (dropout == "appendixC") {
    if (!inherits(model, "cs_appc_model"))
      stop("Supply 'model' = cs_appc_model(...) for dropout = 'appendixC'.")
    if (is.null(muB)) stop("Supply 'muB' for dropout = 'appendixC'.")
    calB <- cs_calibrate_appc(model, "B", pB, mB0, n_cal = n_cal, seed = seed)
    calA <- cs_calibrate_appc(model, "A", pA, muA0, n_cal = n_cal,
                              seed = seed + 7)
    tru <- list(pA = pA, pB = pB, Delta = delta,
                sigma = sqrt((calA$completers_sd^2 + calB$completers_sd^2) / 2),
                muA = muA0, muB = mB0)
    policy_effect <- cs_appc_policy_effect(model, calA, calB,
                                           n = max(2e4, n_cal / 5),
                                           seed = seed + 11)
    hyp_effect <- calA$final_mean - calB$final_mean
    set.seed(seed + 1)
  } else {
    tru <- list(pA = pA, pB = pB, Delta = delta, sigma = sigma,
                muA = muA0, muB = mB0)
  }
  psi_true <- .psi_fun(tru$pA, tru$pB, tru$Delta, tru$sigma, w,
                       design$r)$psi
  theta_true <- (tru$pA * tru$muA + (1 - tru$pA) * s) -
    (tru$pB * tru$muB + (1 - tru$pB) * s)

  res <- matrix(NA_real_, nsim, 8,
                dimnames = list(NULL, c("Zd", "Ze", "reject", "rej_d",
                                        "rej_e", "psi_cover", "theta_cover",
                                        "theta_reject")))
  za <- qnorm(1 - design$alpha)
  cres <- lapply(stats::setNames(comparators, comparators), function(x)
    matrix(NA_real_, nsim, 3, dimnames = list(NULL, c("est", "se", "df"))))
  for (i in seq_len(nsim)) {
    if (!is.null(progress)) progress(i, nsim)
    if (dropout == "random") {
      a <- cs_sim_trial(nA, pA, muA0, sigma)
      b <- cs_sim_trial(nB, pB, mB0, sigma)
    } else {
      trial <- cs_appc_trial(model, calA, calB, nA, nB, post = FALSE)
      tr <- trial$patients
      for (cm in comparators) {
        f <- if (cm == "mmrm")
          cs_mmrm(trial$long, model$k, engine = mmrm_engine)
        else cs_pmi(trial$long, model$k, M = M)
        cres[[cm]][i, ] <- c(f$estimate, f$se, f$df)
      }
      a <- list(comp = tr$completed[tr$arm == "A"], y = tr$y_final[tr$arm == "A"])
      b <- list(comp = tr$completed[tr$arm == "B"], y = tr$y_final[tr$arm == "B"])
    }
    mA <- sum(a$comp); mB <- sum(b$comp)
    if (mA < 2 || mB < 2) next
    pAh <- mA / nA; pBh <- mB / nB; pbar <- (mA + mB) / (nA + nB)
    Zd <- (pAh - pBh) / sqrt(pbar * (1 - pbar) * (1 / nA + 1 / nB))
    yA <- a$y[a$comp]; yB <- b$y[b$comp]
    Dh <- mean(yA) - mean(yB)
    dfe <- mA + mB - 2
    sp <- sqrt(((mA - 1) * stats::var(yA) + (mB - 1) * stats::var(yB)) / dfe)
    seD <- sp * sqrt(1 / mA + 1 / mB)
    Ze <- qnorm(pt(Dh / seD, dfe))
    gate <- !design$contingent ||
      .gate_open(Zd, design$c, if (is.null(design$gate_ties)) "stop" else design$gate_ties)
    if (design$method == "wz") {
      ex <- (w[1] * Zd + w[2] * Ze) > crit$k
    } else {
      ex <- pnorm(Zd, lower.tail = FALSE) * pnorm(Ze, lower.tail = FALSE) < crit$k
    }
    est <- list(pA = pAh, pB = pBh, Delta = Dh, seDelta = seD, sigma = sp,
                df = dfe, nA = nA, nB = nB, r = nB / nA)
    ps <- .psi_fun(pAh, pBh, Dh, sp, w, nB / nA)$psi
    pse <- .psi_se(est, w, "information")
    th_cov <- th_rej <- NA
    if (!is.null(muB)) {
      vA <- c(yA, rep(s, nA - mA)); vB <- c(yB, rep(s, nB - mB))
      th <- mean(vA) - mean(vB)
      sth <- sqrt(stats::var(vA) / nA + stats::var(vB) / nB)
      th_cov <- abs(th - theta_true) <= zq * sth
      th_rej <- th / sth > za
    }
    res[i, ] <- c(Zd, Ze, gate && ex, Zd > za, Ze > za,
                  abs(ps - psi_true) <= zq * pse, th_cov, th_rej)
  }
  ok <- !is.na(res[, "reject"])
  r <- res[ok, , drop = FALSE]
  mc <- function(x) c(mean(x), sqrt(mean(x) * (1 - mean(x)) / length(x)))
  out <- list(
    nsim = sum(ok), n = c(A = nA, B = nB), null = null, dropout = dropout,
    reject = mc(r[, "reject"]),
    reject_completion = mean(r[, "rej_d"]),
    reject_completers = mean(r[, "rej_e"]),
    cor_Zd_Ze = stats::cor(r[, "Zd"], r[, "Ze"]),
    psi_true = psi_true, psi_coverage = mean(r[, "psi_cover"]),
    theta_true = if (!is.null(muB)) theta_true else NA,
    theta_coverage = if (!is.null(muB)) mean(r[, "theta_cover"]) else NA,
    theta_reject = if (!is.null(muB)) mean(r[, "theta_reject"]) else NA,
    calculated_power = if (null) NA else cs_power(design, n)$power,
    truth = tru, treatment_policy_effect = policy_effect,
    hypothetical_effect = hyp_effect,
    comparators = if (length(comparators)) do.call(rbind, lapply(
      comparators, function(cm) cbind(method = cm, .summarise_comparator(
        cres[[cm]][, "est"], cres[[cm]][, "se"], cres[[cm]][, "df"],
        policy_effect, hyp_effect, design$alpha, conf_level)))) else NULL)
  class(out) <- "cs_sim"
  out
}

print.cs_sim <- function(x, ...) {
  cat(sprintf("Simulation: %d trials, nA = %d, nB = %d, %s dropout, %s\n",
              x$nsim, x$n["A"], x$n["B"], x$dropout,
              if (x$null) "under H0" else "under H1"))
  lab <- if (x$null) "Type I error" else "Power"
  cat(sprintf("  %s (composite): %.4f (MC SE %.4f)", lab, x$reject[1],
              x$reject[2]))
  if (!x$null) cat(sprintf("; calculated %.4f", x$calculated_power))
  cat("\n")
  cat(sprintf("  Component rejection: completion %.4f, completers %.4f\n",
              x$reject_completion, x$reject_completers))
  cat(sprintf("  cor(Z(d), Z(e)) = %.4f\n", x$cor_Zd_Ze))
  cat(sprintf("  psi: true %.4f, CI coverage %.4f\n", x$psi_true,
              x$psi_coverage))
  if (!is.na(x$treatment_policy_effect))
    cat(sprintf("  Treatment policy effect implied by the model: %.4f\n",
                x$treatment_policy_effect))
  if (!is.na(x$hypothetical_effect))
    cat(sprintf("  Hypothetical (on-treatment) effect at final visit: %.4f\n",
                x$hypothetical_effect))
  if (!is.null(x$comparators)) {
    cat("  Comparator methods (judged against the treatment policy effect):\n")
    print(format(x$comparators, digits = 3), row.names = FALSE)
  }
  if (!is.na(x$theta_true))
    cat(sprintf("  theta: true %.4f, CI coverage %.4f, rejection %.4f\n",
                x$theta_true, x$theta_coverage, x$theta_reject))
  invisible(x)
}


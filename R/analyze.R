cs_analyze <- function(data, arm, active, completed, outcome,
                       covariates = NULL, strata = NULL,
                       weights = c(0.5, 0.5), method = c("wz", "fisher"),
                       contingent = FALSE, c = 0,
                       calibration = c("recalibrated", "usual"),
                       alpha = 0.025, s = 0,
                       psi_type = c("information", "clinical"),
                       conf_level = 0.95, bootstrap = 0, seed = NULL,
                       completion_model = c("default", "logistic"),
                       gate_ties = c("stop", "pass")) {
  method <- match.arg(method); calibration <- match.arg(calibration)
  psi_type <- match.arg(psi_type)
  completion_model <- match.arg(completion_model)
  gate_ties <- match.arg(gate_ties)
  miss <- setdiff(c(arm, completed, outcome, covariates, strata), names(data))
  if (length(miss)) stop("Columns not found: ", paste(miss, collapse = ", "))
  df <- data.frame(A = as.numeric(data[[arm]] == active),
                   comp = as.logical(data[[completed]]),
                   y = data[[outcome]])
  if (anyNA(df$comp)) stop("Completion status must be known for every patient.")
  if (length(unique(df$A)) != 2) stop("'arm' must have exactly two groups.")
  for (v in covariates) df[[v]] <- data[[v]]
  if (!is.null(strata)) df$.strata <- factor(data[[strata]])
  n_miss <- sum(df$comp & is.na(df$y))

  aligned <- identical(weights, "aligned")
  if (aligned) {
    if (method != "wz")
      stop("Estimand-aligned weights apply to the weighted Z (method = \"wz\").")
    if (contingent)
      stop("Estimand-aligned weights are for the unconditional weighted Z; ",
           "the contingent test uses fixed weights.")
  }
  # dummy planning values: only the test-structure fields are used
  des <- cs_design(0.5, 0.5, delta = 0, sigma = 1,
                   weights = if (aligned) c(0.5, 0.5) else weights,
                   method = method, contingent = contingent, c = c,
                   calibration = calibration, alpha = alpha)
  crit <- cs_critical(des)
  est <- .cs_estimates(df, covariates, !is.null(strata), s)
  if (completion_model == "logistic")
    est$Zd <- .logit_zd(df, covariates, !is.null(strata))

  Zd <- est$Zd; Ze <- est$Ze
  pd <- pnorm(Zd, lower.tail = FALSE); pe <- pnorm(Ze, lower.tail = FALSE)
  w <- des$weights
  if (aligned) {
    # Plug-in weights from this trial: w(d) Z(d) + w(e) Z(e) is then, up to the
    # t-to-normal transformation of Z(e), the estimated theta difference.
    pbar <- mean(df$comp)
    se_d <- sqrt(pbar * (1 - pbar) * (1 / est$nA + 1 / est$nB))
    w <- cs_aligned_weights(est$pA, est$pB, est$muA, est$muB, se_d,
                            est$seDelta, s)$weights
  }
  gate_open <- !contingent || .gate_open(Zd, c, gate_ties)
  if (method == "wz") {
    stat <- w[1] * Zd + w[2] * Ze
    p_uncond <- pnorm(stat, lower.tail = FALSE)
    exceeds <- stat > crit$k
  } else {
    stat <- pd * pe
    p_uncond <- stat * (1 - log(stat))
    exceeds <- stat < crit$k
  }
  p_contingent <- if (!contingent) p_uncond else if (!gate_open) 1 else
    .null_size(des, stat)
  reject <- gate_open && exceeds

  zq <- qnorm(1 - (1 - conf_level) / 2)
  psi <- .psi_fun(est$pA, est$pB, est$Delta, est$sigma, w, est$r, psi_type)$psi
  psi_se <- .psi_se(est, w, psi_type)
  out <- list(
    n = c(A = est$nA, B = est$nB),
    completion = c(A = est$pA, B = est$pB),
    completers_effect = c(estimate = est$Delta, se = est$seDelta,
                          sigma = est$sigma),
    components = c(Z_d = Zd, Z_e = Ze, P_d = pd, P_e = pe),
    test = list(method = method, weights = w, aligned = aligned, statistic = stat,
                critical = crit$k, scale = crit$scale, contingent = contingent,
                c = c, gate_ties = gate_ties, gate_open = gate_open,
                calibration = calibration,
                p_unconditional = p_uncond, p_value = p_contingent,
                reject = reject, alpha = alpha,
                completion_model = completion_model),
    psi = c(estimate = psi, se = psi_se, lower = psi - zq * psi_se,
            upper = psi + zq * psi_se),
    psi_type = psi_type,
    theta = c(estimate = est$theta, se = est$se_theta,
              lower = est$theta - zq * est$se_theta,
              upper = est$theta + zq * est$se_theta,
              p_one_sided = pnorm(est$theta / est$se_theta, lower.tail = FALSE)),
    s = s, conf_level = conf_level, n_completers_missing_outcome = n_miss)

  if (bootstrap > 0) {
    if (!is.null(seed)) set.seed(seed)
    idxA <- which(df$A == 1); idxB <- which(df$A == 0)
    bs <- t(replicate(bootstrap, {
      ii <- c(sample(idxA, replace = TRUE), sample(idxB, replace = TRUE))
      e <- tryCatch(.cs_estimates(df[ii, ], covariates, !is.null(strata), s,
                                  tests = FALSE),
                    error = function(err) NULL)
      if (is.null(e)) c(NA, NA) else
        c(.psi_fun(e$pA, e$pB, e$Delta, e$sigma, w, e$r, psi_type)$psi, e$theta)
    }))
    pr <- c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2)
    out$bootstrap <- list(
      B = bootstrap,
      psi = quantile(bs[, 1], pr, na.rm = TRUE),
      theta = quantile(bs[, 2], pr, na.rm = TRUE))
  }
  class(out) <- "cs_analysis"
  out
}

.cs_estimates <- function(df, covariates, stratified, s, tests = TRUE) {
  nA <- sum(df$A == 1); nB <- sum(df$A == 0)
  pA <- mean(df$comp[df$A == 1]); pB <- mean(df$comp[df$A == 0])
  rhs <- paste(c("A", covariates, if (stratified) ".strata"), collapse = " + ")
  # completers' efficacy
  cdf <- df[df$comp & !is.na(df$y), , drop = FALSE]
  cdf <- cdf[complete.cases(cdf[, c("y", "A", covariates), drop = FALSE]), ,
             drop = FALSE]
  if (stratified) cdf$.strata <- droplevels(cdf$.strata)
  rhs_e <- if (stratified && nlevels(cdf$.strata) < 2)
    paste(c("A", covariates), collapse = " + ") else rhs
  fit <- lm(as.formula(paste("y ~", rhs_e)), data = cdf)
  Delta <- unname(coef(fit)["A"]); seDelta <- sqrt(vcov(fit)["A", "A"])
  dfe <- df.residual(fit); sig <- sigma(fit)
  # theta: dropouts scored s
  tdf <- df[!(df$comp & is.na(df$y)), , drop = FALSE]
  tdf$v <- ifelse(tdf$comp, tdf$y, s)
  tdf <- tdf[complete.cases(tdf[, c("v", "A", covariates), drop = FALSE]), ,
             drop = FALSE]
  if (stratified) tdf$.strata <- droplevels(tdf$.strata)
  rhs_t <- if (stratified && nlevels(tdf$.strata) < 2)
    paste(c("A", covariates), collapse = " + ") else rhs
  tfit <- lm(as.formula(paste("v ~", rhs_t)), data = tdf)
  out <- list(nA = nA, nB = nB, r = nB / nA, pA = pA, pB = pB,
              muA = mean(cdf$y[cdf$A == 1]), muB = mean(cdf$y[cdf$A == 0]),
              Delta = Delta, seDelta = seDelta, sigma = sig, df = dfe,
              theta = unname(coef(tfit)["A"]),
              se_theta = .hc2_se(tfit, "A"))
  if (tests) {
    out$Zd <- if (stratified) .cmh_z(df) else {
      pbar <- mean(df$comp)
      (pA - pB) / sqrt(pbar * (1 - pbar) * (1 / nA + 1 / nB))
    }
    tval <- Delta / seDelta
    out$Ze <- qnorm(pt(tval, dfe))
  }
  out
}

cs_theta_sensitivity <- function(data, arm, active, completed, outcome,
                                 s = c(0, -5, -10), covariates = NULL,
                                 strata = NULL, conf_level = 0.95) {
  miss <- setdiff(c(arm, completed, outcome, covariates, strata), names(data))
  if (length(miss)) stop("Columns not found: ", paste(miss, collapse = ", "))
  df <- data.frame(A = as.numeric(data[[arm]] == active),
                   comp = as.logical(data[[completed]]), y = data[[outcome]])
  if (anyNA(df$comp)) stop("Completion status must be known for every patient.")
  for (v in covariates) df[[v]] <- data[[v]]
  if (!is.null(strata)) df$.strata <- factor(data[[strata]])
  zq <- qnorm(1 - (1 - conf_level) / 2)
  ok <- !(df$comp & is.na(df$y))
  yA <- df$y[ok & df$comp & df$A == 1]; yB <- df$y[ok & df$comp & df$A == 0]
  pA <- mean(df$comp[ok & df$A == 1]); pB <- mean(df$comp[ok & df$A == 0])
  out <- do.call(rbind, lapply(s, function(si) {
    e <- .cs_estimates(df, covariates, !is.null(strata), si, tests = FALSE)
    data.frame(s = si, theta = e$theta, se = e$se_theta,
               lower = e$theta - zq * e$se_theta,
               upper = e$theta + zq * e$se_theta,
               p_one_sided = pnorm(e$theta / e$se_theta, lower.tail = FALSE),
               completion_part = (pA - pB) * ((mean(yA) + mean(yB)) / 2 - si),
               completers_part = (pA + pB) / 2 * (mean(yA) - mean(yB)))
  }))
  attr(out, "adjusted") <- !is.null(covariates) || !is.null(strata)
  out
}

# HC2 robust standard error of one coefficient. The per-patient benefit
# (completers' improvement, dropouts scored s) has variance
# p sigma^2 + p (1 - p) (mu - s)^2, which differs between arms, so the
# model-based lm variance is wrong under unequal allocation. Without
# covariates HC2 reduces to sqrt(var_A / n_A + var_B / n_B).
.hc2_se <- function(fit, term) {
  X <- stats::model.matrix(fit)[, !is.na(coef(fit)), drop = FALSE]
  h <- stats::hatvalues(fit)
  u <- stats::residuals(fit) / sqrt(1 - h)
  XtXi <- summary(fit)$cov.unscaled
  V <- XtXi %*% crossprod(X * u) %*% XtXi
  sqrt(V[term, term])
}

# Cochran-Mantel-Haenszel Z for completion, oriented so positive favours A
.cmh_z <- function(df) {
  num <- 0; v <- 0
  for (lv in levels(df$.strata)) {
    d <- df[df$.strata == lv, ]
    N <- nrow(d); if (N < 2) next
    n1 <- sum(d$A == 1); n0 <- N - n1; x <- sum(d$comp); x1 <- sum(d$comp & d$A == 1)
    num <- num + x1 - n1 * x / N
    v <- v + n1 * n0 * x * (N - x) / (N^2 * (N - 1))
  }
  num / sqrt(v)
}

# delta-method SE of psi-hat, treating (pA, pB, Delta, sigma) as independent
.psi_se <- function(est, w, type) {
  th <- c(est$pA, est$pB, est$Delta, est$sigma)
  V <- c(est$pA * (1 - est$pA) / est$nA, est$pB * (1 - est$pB) / est$nB,
         est$seDelta^2, est$sigma^2 / (2 * est$df))
  f <- function(x) .psi_fun(x[1], x[2], x[3], x[4], w, est$r, type)$psi
  g <- vapply(seq_along(th), function(j) {
    h <- 1e-6 * max(1, abs(th[j]))
    e <- replace(numeric(4), j, h)
    (f(th + e) - f(th - e)) / (2 * h)
  }, numeric(1))
  sqrt(sum(g^2 * V))
}

print.cs_analysis <- function(x, digits = 4, ...) {
  f <- function(v) formatC(v, digits = digits, format = "fg")
  cat("Composite-strategy analysis\n")
  cat(sprintf("  Randomized: A = %d, B = %d\n", x$n["A"], x$n["B"]))
  cat(sprintf("  Completion: A = %s, B = %s;  Z(d) = %s (P = %s)\n",
              f(x$completion["A"]), f(x$completion["B"]),
              f(x$components["Z_d"]), f(x$components["P_d"])))
  cat(sprintf("  Completers' effect: %s (SE %s);  Z(e) = %s (P = %s)\n",
              f(x$completers_effect["estimate"]), f(x$completers_effect["se"]),
              f(x$components["Z_e"]), f(x$components["P_e"])))
  t <- x$test
  cat(sprintf("  Composite test (%s, w = %s/%s%s): statistic %s, P = %s -> %s\n",
              if (t$method != "wz") "Fisher" else if (isTRUE(t$aligned))
                "weighted Z, estimand-aligned" else "weighted Z",
              f(t$weights[1]), f(t$weights[2]),
              if (t$contingent) sprintf(", contingent c = %g, %s", t$c,
                                        t$calibration) else "",
              f(t$statistic), f(t$p_value),
              if (t$reject) "reject H0" else "do not reject H0"))
  if (t$contingent && !t$gate_open)
    cat("    (gate closed: A observed worse on completion; testing stopped)\n")
  cl <- round(100 * x$conf_level)
  cat(sprintf("  psi (%s) = %s, %d%% CI %s to %s\n", x$psi_type,
              f(x$psi["estimate"]), cl, f(x$psi["lower"]), f(x$psi["upper"])))
  cat(sprintf("  theta (s = %g) = %s, %d%% CI %s to %s\n", x$s,
              f(x$theta["estimate"]), cl, f(x$theta["lower"]),
              f(x$theta["upper"])))
  if (!is.null(x$bootstrap))
    cat(sprintf("  Bootstrap (B = %d) CIs: psi %s to %s; theta %s to %s\n",
                x$bootstrap$B, f(x$bootstrap$psi[1]), f(x$bootstrap$psi[2]),
                f(x$bootstrap$theta[1]), f(x$bootstrap$theta[2])))
  if (x$n_completers_missing_outcome > 0)
    cat(sprintf("  Note: %d completer(s) with missing final outcome excluded.\n",
                x$n_completers_missing_outcome))
  invisible(x)
}


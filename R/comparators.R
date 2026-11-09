NULL

cs_mmrm <- function(long, k = max(long$visit), active = "A",
                    engine = c("nlme", "mmrm")) {
  engine <- match.arg(engine)
  dat <- data.frame(id = factor(long$id), visit = as.integer(long$visit),
                    visitf = factor(long$visit, levels = seq_len(k)),
                    A = as.numeric(long$arm == active), y = long$y)
  dat <- dat[order(dat$id, dat$visit), ]
  nm <- paste0("visitf", k, ":A")
  fail <- list(estimate = NA_real_, se = NA_real_, df = NA_real_,
               converged = FALSE)
  if (engine == "mmrm") {
    if (!requireNamespace("mmrm", quietly = TRUE))
      stop("Install the 'mmrm' package or use engine = 'nlme'.")
    fit <- tryCatch(mmrm::mmrm(y ~ 0 + visitf + visitf:A + us(visitf | id),
                               data = dat),
                    error = function(e) NULL)
    if (is.null(fit)) return(fail)
    ct <- summary(fit)$coefficients
    return(list(estimate = ct[nm, "Estimate"], se = ct[nm, "Std. Error"],
                df = ct[nm, "df"], converged = TRUE))
  }
  fit <- tryCatch(
    nlme::gls(y ~ 0 + visitf + visitf:A, data = dat,
              correlation = nlme::corSymm(form = ~ visit | id),
              weights = nlme::varIdent(form = ~ 1 | visitf),
              method = "REML", na.action = stats::na.omit),
    error = function(e) NULL)
  if (is.null(fit)) return(fail)
  b <- stats::coef(fit); V <- stats::vcov(fit)
  list(estimate = unname(b[nm]), se = sqrt(V[nm, nm]),
       df = nlevels(dat$id) - 2 * k, converged = TRUE)
}

cs_pmi <- function(long, k = max(long$visit), active = "A", reference = "B",
                   M = 20) {
  ids <- sort(unique(long$id))
  W <- matrix(NA_real_, length(ids), k)
  W[cbind(match(long$id, ids), long$visit)] <- long$y
  arm <- long$arm[match(ids, long$id)]
  if (anyNA(W[, 1])) stop("Visit 1 must be observed for every patient.")
  isref <- arm == reference
  est <- vv <- numeric(M)
  for (m in seq_len(M)) {
    Wm <- W
    for (j in 2:k) {
      miss <- is.na(Wm[, j])
      if (!any(miss)) next
      fitrows <- isref & !is.na(W[, j])          # observed reference patients
      X <- cbind(1, Wm[fitrows, 1:(j - 1), drop = FALSE])
      yj <- W[fitrows, j]
      if (nrow(X) <= ncol(X) + 1)
        stop("Too few observed reference patients at visit ", j, ".")
      XtX_inv <- solve(crossprod(X))
      bhat <- XtX_inv %*% crossprod(X, yj)
      dfr <- nrow(X) - ncol(X)
      s2 <- sum((yj - X %*% bhat)^2) / stats::rchisq(1, dfr)
      bdraw <- bhat + t(chol(s2 * XtX_inv)) %*% stats::rnorm(ncol(X))
      Xm <- cbind(1, Wm[miss, 1:(j - 1), drop = FALSE])
      Wm[miss, j] <- as.vector(Xm %*% bdraw) + stats::rnorm(sum(miss), 0, sqrt(s2))
    }
    yA <- Wm[arm == active, k]; yB <- Wm[arm == reference, k]
    est[m] <- mean(yA) - mean(yB)
    vv[m] <- stats::var(yA) / length(yA) + stats::var(yB) / length(yB)
  }
  .rubin(est, vv, n_complete = length(ids) - 2)
}

# Rubin's rules with Barnard-Rubin small-sample degrees of freedom
.rubin <- function(est, vv, n_complete) {
  M <- length(est)
  qbar <- mean(est); ubar <- mean(vv); b <- stats::var(est)
  tv <- ubar + (1 + 1 / M) * b
  lam <- (1 + 1 / M) * b / tv
  df_old <- if (lam > 0) (M - 1) / lam^2 else Inf
  df_obs <- (n_complete + 1) / (n_complete + 3) * n_complete * (1 - lam)
  df <- if (is.finite(df_old)) 1 / (1 / df_old + 1 / df_obs) else df_obs
  list(estimate = qbar, se = sqrt(tv), df = df, converged = TRUE,
       fraction_missing_info = lam)
}

# Summarise comparator results over simulated trials
.summarise_comparator <- function(est, se, df, target_policy, target_hyp,
                                  alpha, conf_level) {
  ok <- is.finite(est) & is.finite(se)
  e <- est[ok]; s <- se[ok]; d <- pmax(df[ok], 1)
  q <- stats::qt(1 - (1 - conf_level) / 2, d)
  data.frame(
    n_ok = sum(ok),
    mean_estimate = mean(e),
    bias_policy = mean(e) - target_policy,
    bias_hypothetical = mean(e) - target_hyp,
    empirical_se = stats::sd(e),
    mean_model_se = mean(s),
    coverage_policy = mean(abs(e - target_policy) <= q * s),
    coverage_hypothetical = mean(abs(e - target_hyp) <= q * s),
    rejection = mean(stats::pt(e / s, d, lower.tail = FALSE) < alpha))
}


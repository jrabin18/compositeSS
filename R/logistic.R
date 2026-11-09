# Covariate-adjusted completion statistic (Shih and Quan 2001, model 3a):
# Wald Z for treatment in a logistic regression of completion on treatment,
# baseline covariates and (optionally) stratum. Positive values favour A.
.logit_zd <- function(df, covariates, stratified) {
  rhs <- paste(c("A", covariates, if (stratified) ".strata"), collapse = " + ")
  d <- df[complete.cases(df[, c("comp", "A", covariates), drop = FALSE]), ,
          drop = FALSE]
  if (stratified) d$.strata <- droplevels(d$.strata)
  fit <- suppressWarnings(stats::glm(stats::as.formula(paste("comp ~", rhs)),
                                     family = stats::binomial(), data = d))
  unname(summary(fit)$coefficients["A", "z value"])
}

cs_sim_sq2001 <- function(nsim = 20000, n = 150, weights = c(0.3, 0.7), c = 0,
                          alpha = 0.025, errors = c("chisq3", "normal"),
                          seed = 1) {
  errors <- match.arg(errors)
  set.seed(seed)
  d <- cs_design(0.5, 0.5, delta = 0, sigma = 1, weights = weights,
                 alpha = alpha)
  dc <- cs_design(0.5, 0.5, delta = 0, sigma = 1, weights = weights,
                  alpha = alpha, contingent = TRUE, c = c)
  k <- cs_critical(d)$k; kc <- cs_critical(dc)$k; w <- d$weights
  A <- rep(1:0, each = n)
  Z <- matrix(NA_real_, nsim, 2, dimnames = list(NULL, c("Z_d", "Z_e")))
  for (i in seq_len(nsim)) {
    g <- stats::rnorm(2 * n)
    comp <- stats::runif(2 * n) < stats::plogis(0.4 + 0.8 * g)
    e <- if (errors == "chisq3") (stats::rchisq(2 * n, 3) - 3) / sqrt(6) else
      stats::rnorm(2 * n)
    y <- 10 + 4 * g + 5 * e
    fd <- suppressWarnings(stats::glm(comp ~ A + g, family = stats::binomial()))
    fe <- stats::lm(y ~ A + g, subset = comp)
    te <- summary(fe)$coefficients["A", "t value"]
    Z[i, ] <- c(summary(fd)$coefficients["A", "z value"],
                stats::qnorm(stats::pt(te, fe$df.residual)))
  }
  S <- w[1] * Z[, 1] + w[2] * Z[, 2]
  list(nsim = nsim, n_per_arm = n, errors = errors,
       cor_Zd_Ze = stats::cor(Z[, 1], Z[, 2]),
       mean_sd = rbind(Z_d = c(mean(Z[, 1]), stats::sd(Z[, 1])),
                       Z_e = c(mean(Z[, 2]), stats::sd(Z[, 2]))),
       size_weighted_z = mean(S > k),
       size_contingent = mean(.gate_open(Z[, 1], c) & S > kc),
       critical = c(weighted_z = k, contingent = kc),
       mc_se = sqrt(alpha * (1 - alpha) / nsim))
}


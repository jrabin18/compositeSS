cs_appc_model <- function(k = 6, shape = c("linear", "quadratic", "flat"),
                          cov = c("cs", "ris"), sd = 18, rho = 0.5,
                          sd_int = 12, sd_slope = 10, cor_is = 0,
                          sd_resid = 10, beta = 0, lambda = 0,
                          depends_on = c("current", "previous"),
                          post_sd = NULL) {
  if (is.numeric(shape)) {
    if (length(shape) != k || abs(shape[k] - 1) > 1e-12)
      stop("A numeric 'shape' must have length k and last element 1.")
    shp <- shape
  } else {
    shape <- match.arg(shape)
    t <- seq_len(k) / k
    shp <- switch(shape, linear = t, quadratic = 2 * t - t^2,
                  flat = rep(1, k))
  }
  cov <- match.arg(cov); depends_on <- match.arg(depends_on)
  if (k < 2) stop("'k' must be at least 2.")
  if (rho < 0 || rho >= 1) stop("'rho' must be in [0, 1).")
  two <- function(x, nm) {
    if (length(x) == 1) x <- c(x, x)
    if (length(x) != 2) stop(sprintf("'%s' must have length 1 or 2.", nm))
    stats::setNames(x, c("A", "B"))
  }
  lam_vec <- function(l) {
    l <- if (length(l) == 1) c(0, rep(l, k - 1)) else l
    if (length(l) != k || any(l < 0 | l >= 1))
      stop("'lambda' must be a probability or a vector of k probabilities.")
    l[1] <- 0
    l
  }
  lam <- if (is.list(lambda)) list(A = lam_vec(lambda$A), B = lam_vec(lambda$B))
         else list(A = lam_vec(lambda), B = lam_vec(lambda))
  sdv <- two(sd, "sd")
  if (is.null(post_sd))
    post_sd <- if (cov == "cs") mean(sdv) * sqrt(1 - rho) else sd_resid
  structure(list(k = k, shape = shp, cov = cov, sd = sdv, rho = rho,
                 sd_int = sd_int, sd_slope = sd_slope, cor_is = cor_is,
                 sd_resid = sd_resid, beta = two(beta, "beta"), lambda = lam,
                 depends_on = depends_on, post_sd = post_sd),
            class = "cs_appc_model")
}

# Common random numbers for one arm of n patients
.appc_noise <- function(n, k) {
  list(z1 = stats::rnorm(n), z2 = stats::rnorm(n),
       E = matrix(stats::rnorm(n * k), n, k),
       U1 = matrix(stats::runif(n * k), n, k),
       U2 = matrix(stats::runif(n * k), n, k),
       P = matrix(stats::rnorm(n * k), n, k))
}

# Deviations from the mean trajectory (n x k)
.appc_dev <- function(model, arm, nz) {
  k <- model$k; n <- length(nz$z1)
  if (model$cov == "cs") {
    s <- model$sd[[arm]]
    s * (sqrt(model$rho) * nz$z1 + sqrt(1 - model$rho) * nz$E)
  } else {
    t <- seq_len(k) / k
    b0 <- model$sd_int * nz$z1
    b1 <- model$sd_slope * (model$cor_is * nz$z1 +
                              sqrt(1 - model$cor_is^2) * nz$z2)
    b0 + outer(b1, t) + model$sd_resid * nz$E
  }
}

# Core generator for one arm given final-visit mean Fm and intercept alpha
.appc_arm <- function(model, arm, Fm, alpha, nz) {
  k <- model$k; n <- length(nz$z1)
  mu <- Fm * model$shape
  Y <- sweep(.appc_dev(model, arm, nz), 2, mu, "+")
  beta <- model$beta[[arm]]; lam <- model$lambda[[arm]]
  alive <- rep(TRUE, n); dv <- rep(NA_integer_, n); reason <- rep(NA_character_, n)
  for (j in 2:k) {
    ycur <- if (model$depends_on == "current") Y[, j] else Y[, j - 1]
    ae <- alive & nz$U1[, j] < lam[j]
    lor <- alive & !ae & nz$U2[, j] < stats::plogis(alpha + beta * ycur)
    dv[ae | lor] <- j
    reason[ae] <- "unrelated"
    reason[lor] <- if (beta == 0) "unrelated" else "lack_of_response"
    alive <- alive & !(ae | lor)
  }
  obs <- Y
  for (j in 2:k) obs[!is.na(dv) & dv <= j, j] <- NA
  list(Y = Y, obs = obs, mu = mu, comp = alive, drop_visit = dv,
       reason = reason, nz = nz)
}

# Part 3: outcomes after dropout (treatment policy trajectories)
.appc_post <- function(model, g, ref_mu) {
  # The same rules apply in both arms; ref_mu (the comparator's mean trajectory)
  # only supplies the increments for dropout unrelated to outcome.
  k <- model$k; Ytp <- g$obs; tt <- seq_len(k)
  for (d in 2:k) {
    idx <- which(!g$comp & g$drop_visit == d)
    if (!length(idx)) next
    last <- d - 1; later <- d:k
    ylast <- g$obs[idx, last]
    fl <- model$post_sd * g$nz$P[idx, later, drop = FALSE]
    un <- g$reason[idx] == "unrelated"
    if (any(un)) {
      inc <- ref_mu[later] - ref_mu[last]
      Ytp[idx[un], later] <- ylast[un] +
        matrix(inc, sum(un), length(later), byrow = TRUE) + fl[un, , drop = FALSE]
    }
    if (any(!un)) {
      yy <- g$obs[idx[!un], 1:last, drop = FALSE]
      slope <- as.vector(yy %*% tt[1:last]) / sum(tt[1:last]^2)  # own trend from baseline
      Ytp[idx[!un], later] <- ylast[!un] + outer(slope, later - last) +
        fl[!un, , drop = FALSE]
    }
  }
  Ytp
}

cs_calibrate_appc <- function(model, arm = c("A", "B"), p, completers_mean,
                              n_cal = 1e5, seed = 1, tol = 1e-3,
                              max_iter = 50) {
  arm <- match.arg(arm)
  .check_prob(p, "p")
  old <- if (exists(".Random.seed", envir = globalenv()))
    get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  nz <- .appc_noise(n_cal, model$k)
  comp_rate <- function(a, Fm) mean(.appc_arm(model, arm, Fm, a, nz)$comp)
  find_alpha <- function(Fm) {
    lo <- -40; hi <- 40
    if (comp_rate(lo, Fm) < p)
      stop("Outcome-unrelated dropout (lambda) alone already gives completion ",
           "below the target; reduce lambda.")
    if (comp_rate(hi, Fm) > p) stop("Cannot reach the target completion rate.")
    uniroot(function(a) comp_rate(a, Fm) - p, c(lo, hi), tol = 1e-8)$root
  }
  Fm <- completers_mean
  for (it in seq_len(max_iter)) {
    a <- find_alpha(Fm)
    g <- .appc_arm(model, arm, Fm, a, nz)
    cm <- mean(g$Y[g$comp, model$k])
    if (abs(cm - completers_mean) < tol) break
    Fm <- Fm + (completers_mean - cm)
  }
  g <- .appc_arm(model, arm, Fm, a, nz)
  list(arm = arm, alpha = a, final_mean = Fm, mu = Fm * model$shape,
       completion = mean(g$comp),
       completers_mean = mean(g$Y[g$comp, model$k]),
       completers_sd = stats::sd(g$Y[g$comp, model$k]),
       iterations = it)
}

cs_appc_trial <- function(model, calA, calB, nA, nB, post = TRUE) {
  gB <- .appc_arm(model, "B", calB$final_mean, calB$alpha,
                  .appc_noise(nB, model$k))
  gA <- .appc_arm(model, "A", calA$final_mean, calA$alpha,
                  .appc_noise(nA, model$k))
  k <- model$k
  mk <- function(g, arm, n, ref) {
    tp <- if (post) .appc_post(model, g, calB$mu)[, k] else NA
    data.frame(arm = arm, completed = g$comp, drop_visit = g$drop_visit,
               reason = g$reason,
               y_final = ifelse(g$comp, g$Y[, k], NA_real_),
               y_tp_final = tp)
  }
  pts <- rbind(mk(gA, "A", nA, FALSE), mk(gB, "B", nB, TRUE))
  pts$id <- seq_len(nrow(pts))
  obs <- rbind(gA$obs, gB$obs)
  long <- data.frame(id = rep(pts$id, each = k), arm = rep(pts$arm, each = k),
                     visit = rep(seq_len(k), nrow(pts)), y = as.vector(t(obs)))
  long <- long[!is.na(long$y), ]
  list(patients = pts[, c("id", setdiff(names(pts), "id"))], long = long)
}

cs_appc_policy_effect <- function(model, calA, calB, n = 1e5, seed = 99) {
  old <- if (exists(".Random.seed", envir = globalenv()))
    get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  # Common random numbers: both arms use the same noise, so the effect is
  # exactly zero when the arms are identical and Monte Carlo error is reduced.
  nz <- .appc_noise(n, model$k); k <- model$k
  gA <- .appc_arm(model, "A", calA$final_mean, calA$alpha, nz)
  gB <- .appc_arm(model, "B", calB$final_mean, calB$alpha, nz)
  mean(.appc_post(model, gA, calB$mu)[, k]) - mean(.appc_post(model, gB, calB$mu)[, k])
}

cs_appc_plot <- function(model, calA, calB, n = 2e4, seed = 5) {
  old <- if (exists(".Random.seed", envir = globalenv()))
    get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  k <- model$k
  gA <- .appc_arm(model, "A", calA$final_mean, calA$alpha, .appc_noise(n, k))
  gB <- .appc_arm(model, "B", calB$final_mean, calB$alpha, .appc_noise(n, k))
  tpA <- .appc_post(model, gA, calB$mu)
  tpB <- .appc_post(model, gB, calB$mu)
  res <- data.frame(
    visit = rep(seq_len(k), 2), arm = rep(c("A", "B"), each = k),
    on_treatment = c(calA$mu, calB$mu),
    observed = c(colMeans(gA$obs, na.rm = TRUE), colMeans(gB$obs, na.rm = TRUE)),
    policy = c(colMeans(tpA), colMeans(tpB)),
    remaining = c(colMeans(!is.na(gA$obs)), colMeans(!is.na(gB$obs))))
  op <- graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  on.exit(graphics::par(op), add = TRUE)
  yl <- range(c(0, res$on_treatment, res$observed, res$policy))
  plot(c(0, k), yl, type = "n", xlab = "Visit", ylab = "Mean improvement",
       main = "Trajectories")
  for (a in c("A", "B")) {
    d <- res[res$arm == a, ]; col <- if (a == "A") 2 else 4
    lines(c(0, d$visit), c(0, d$on_treatment), col = col, lwd = 2)
    lines(c(0, d$visit), c(0, d$observed), col = col, lwd = 2, lty = 2)
    lines(c(0, d$visit), c(0, d$policy), col = col, lwd = 2, lty = 3)
  }
  legend("topleft", bty = "n", cex = 0.75,
         legend = c("A", "B", "on treatment", "observed (remaining)",
                    "treatment policy"),
         col = c(2, 4, 1, 1, 1), lty = c(1, 1, 1, 2, 3), lwd = 2)
  plot(c(0, k), c(0, 1), type = "n", xlab = "Visit",
       ylab = "Proportion on treatment", main = "Retention")
  for (a in c("A", "B")) {
    d <- res[res$arm == a, ]
    lines(c(0, d$visit), c(1, d$remaining), type = "s",
          col = if (a == "A") 2 else 4, lwd = 2)
  }
  invisible(res)
}


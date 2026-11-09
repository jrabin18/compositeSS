cs_assurance <- function(design, n, prior_pA = NULL, prior_pB = NULL,
                         prior_delta = NULL, nsim = 2000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  dr <- function(pr, fixed) {
    if (is.null(pr)) return(rep(fixed, nsim))
    if (length(pr) != 2 || any(pr <= 0))
      stop("Beta prior must be c(a, b) with a, b > 0.")
    rbeta(nsim, pr[1], pr[2])
  }
  pA <- dr(prior_pA, design$pA); pB <- dr(prior_pB, design$pB)
  if (inherits(design, "theta_design")) {
    eff0 <- design$muA - design$muB
    eff <- if (is.null(prior_delta)) rep(eff0, nsim) else
      rnorm(nsim, prior_delta[1], prior_delta[2])
    pw <- vapply(seq_len(nsim), function(i) {
      d <- design; d$pA <- pA[i]; d$pB <- pB[i]; d$muA <- d$muB + eff[i]
      theta_power(d, n)
    }, numeric(1))
    p0 <- theta_power(design, n)
  } else {
    eff <- if (is.null(prior_delta)) rep(design$delta, nsim) else
      rnorm(nsim, prior_delta[1], prior_delta[2])
    crit <- cs_critical(design)
    pw <- vapply(seq_len(nsim), function(i) {
      d <- .update_design(design, pA = pA[i], pB = pB[i], delta = eff[i])
      .power_one(d, cs_components(d, n), crit$k)
    }, numeric(1))
    p0 <- cs_power(design, n)$power
  }
  list(assurance = mean(pw), power_at_planning_values = p0,
       draws = data.frame(pA = pA, pB = pB, effect = eff, power = pw))
}

cs_assurance_sample_size <- function(design, target = 0.70, prior_pA = NULL,
                                     prior_pB = NULL, prior_delta = NULL,
                                     nsim = 1000, seed = 2026, n_max = 2e4) {
  if (!is.null(seed)) set.seed(seed)
  dr <- function(pr, fixed) {
    if (is.null(pr)) return(rep(fixed, nsim))
    if (length(pr) != 2 || any(pr <= 0))
      stop("Beta prior must be c(a, b) with a, b > 0.")
    rbeta(nsim, pr[1], pr[2])
  }
  pA <- dr(prior_pA, design$pA); pB <- dr(prior_pB, design$pB)
  is_theta <- inherits(design, "theta_design")
  eff0 <- if (is_theta) design$muA - design$muB else design$delta
  eff <- if (is.null(prior_delta)) rep(eff0, nsim) else
    rnorm(nsim, prior_delta[1], prior_delta[2])
  # the same prior draws at every n (common random numbers), so that the
  # Monte Carlo assurance is monotone in n and the search is well defined
  if (is_theta) {
    ds <- lapply(seq_len(nsim), function(i) {
      d <- design; d$pA <- pA[i]; d$pB <- pB[i]; d$muA <- d$muB + eff[i]; d })
    af <- function(n) mean(vapply(ds, theta_power, numeric(1), n = n))
  } else {
    crit <- cs_critical(design)
    ds <- lapply(seq_len(nsim), function(i)
      .update_design(design, pA = pA[i], pB = pB[i], delta = eff[i]))
    af <- function(n) mean(vapply(ds, function(d)
      .power_one(d, cs_components(d, n), crit$k), numeric(1)))
  }
  a_max <- af(n_max)
  if (a_max < target)
    stop(sprintf(paste0("Target assurance %.2f not reached by n_max = %d ",
                        "(assurance there %.3f). Assurance is bounded by the ",
                        "prior probability that A is truly better."),
                 target, as.integer(n_max), a_max))
  n <- .search_n(af, target, n_max)
  list(nA = n, nB = ceiling(design$r * n), total = n + ceiling(design$r * n),
       assurance = af(n), target = target, nsim = nsim,
       assurance_upper_bound = a_max)
}

cs_beta_prior <- function(mean, sd) {
  .check_prob(mean, "mean")
  k <- mean * (1 - mean) / sd^2 - 1
  if (sd <= 0 || k <= 0)
    stop("'sd' must be positive and smaller than sqrt(mean * (1 - mean)).")
  c(a = mean * k, b = (1 - mean) * k)
}

cs_power_grid <- function(design, n, pA = design$pA, pB = design$pB,
                          delta = design$delta) {
  if (is.null(delta)) delta <- NA_real_
  g <- expand.grid(pA = pA, pB = pB, delta = delta)
  crit <- cs_critical(design)
  g$power <- vapply(seq_len(nrow(g)), function(i) {
    d <- .update_design(design, pA = g$pA[i], pB = g$pB[i])
    if (!is.na(g$delta[i])) d$delta <- g$delta[i]
    .power_one(d, cs_components(d, n), crit$k)
  }, numeric(1))
  g
}

cs_power_curve <- function(designs, n = seq(20, 400, by = 5), target = 0.80,
                           ...) {
  if (inherits(designs, c("cs_design", "theta_design")))
    designs <- list(design = designs)
  if (is.null(names(designs))) names(designs) <- paste("Design", seq_along(designs))
  vals <- lapply(designs, function(d)
    if (inherits(d, "theta_design")) theta_power(d, n) else cs_power(d, n)$power)
  plot(range(n), c(0, 1), type = "n", xlab = "Patients per arm (A)",
       ylab = "Power", ...)
  for (i in seq_along(vals)) lines(n, vals[[i]], lty = i, col = i, lwd = 2)
  abline(h = target, col = "grey50")
  legend("bottomright", legend = names(designs), lty = seq_along(vals),
         col = seq_along(vals), lwd = 2, bty = "n", cex = 0.8)
  invisible(data.frame(n = n, do.call(cbind, vals)))
}

cs_estimand_statement <- function(design,
                                  population = "patients as defined by the trial's eligibility criteria",
                                  treatment = "the active treatment (A) versus the comparator (B), each as randomized",
                                  endpoint = "improvement in the primary efficacy score at the final visit") {
  ice <- paste0("Treatment discontinuation (including switching antipsychotic or ",
                "adding rescue medication not defined as part of the regimen) is ",
                "handled by the composite strategy: it is part of the outcome, ",
                "and no outcomes after discontinuation are modelled.")
  if (inherits(design, "theta_design")) {
    var <- sprintf(paste0("Benefit per randomized patient: %s for patients who ",
                          "complete the trial on assigned treatment; a score of %g ",
                          "for patients who discontinue."), endpoint, design$s)
    sm <- sprintf(paste0("Difference in mean benefit per randomized patient ",
                         "(theta, in clinical units)%s."),
                  if (design$margin > 0) sprintf("; non-inferiority margin %g",
                                                 design$margin) else "")
  } else {
    w <- design$weights / sum(design$weights)
    var <- sprintf(paste0("Two components: completion of the trial on assigned ",
                          "treatment, and %s among completers."),
                   if (design$efficacy == "continuous") endpoint else
                     "response among completers")
    sm <- sprintf(paste0("Composite of the difference in completion and the ",
                         "difference among completers, combined by %s with ",
                         "prespecified relative weights %.2f (completion) and ",
                         "%.2f (efficacy)%s; summarized by the standardized ",
                         "composite effect size psi with its confidence interval, ",
                         "reported with both components."),
                  if (design$method == "wz") "the weighted Z" else "Fisher's method",
                  w[1], w[2],
                  if (design$contingent) sprintf(paste0("; efficacy counts only if ",
                    "A is not worse on completion (gate Z(d) >= %g, %s critical value)"),
                    design$c, design$calibration) else "")
  }
  txt <- paste0("Population: ", population, ".\n",
                "Treatment: ", treatment, ".\n",
                "Variable: ", var, "\n",
                "Intercurrent events: ", ice, "\n",
                "Population-level summary: ", sm)
  cat(txt, "\n")
  invisible(txt)
}

cs_compare_designs <- function(pA, pB, delta, sigma, muB, s = 0,
                               alpha = 0.025, power = 0.80, gamma = "zero",
                               c = 0, equal_contingent = FALSE) {
  opts <- list(
    "Weighted Z, equal weights (0.5 / 0.5)" = list(w = c(.5, .5)),
    "Weighted Z, 0.7 / 0.3" = list(w = c(.7, .3)),
    "Weighted Z, 0.3 / 0.7" = list(w = c(.3, .7)),
    "Weighted Z, 0.3 / 0.7, contingent, recalibrated" =
      list(w = c(.3, .7), cont = TRUE, cal = "recalibrated"),
    "Weighted Z, 0.3 / 0.7, contingent, usual critical value" =
      list(w = c(.3, .7), cont = TRUE, cal = "usual"),
    "Fisher (Shih-Quan), equal weights" = list(w = c(.5, .5), m = "fisher"),
    "Fisher, equal weights, contingent, recalibrated" =
      list(w = c(.5, .5), m = "fisher", cont = TRUE, cal = "recalibrated"),
    "Completion rates only" = list(special = "completion"),
    "Completers analysis only" = list(w = c(0, 1)),
    "Expected benefit per randomized patient" = list(special = "theta"))
  if (equal_contingent)   # the design specified in Box B1 of the proposal
    opts <- append(opts, list(
      "Weighted Z, equal weights, contingent, recalibrated" =
        list(w = c(.5, .5), cont = TRUE, cal = "recalibrated"),
      "Weighted Z, equal weights, contingent, usual critical value" =
        list(w = c(.5, .5), cont = TRUE, cal = "usual")), after = 5)
  one <- function(o) {
    if (identical(o$special, "completion")) {
      r <- completion_sample_size(pA, pB, alpha = alpha, power = power)
      return(c(r$nA, NA))
    }
    if (identical(o$special, "theta")) {
      r <- theta_sample_size(theta_design(pA, pB, muB + delta, muB, sigma,
                                          s = s, alpha = alpha), power)
      return(c(r$nA, r$power))
    }
    d <- cs_design(pA, pB, delta = delta, sigma = sigma, weights = o$w,
                   method = if (is.null(o$m)) "wz" else o$m,
                   contingent = isTRUE(o$cont), c = c,
                   calibration = if (is.null(o$cal)) "recalibrated" else o$cal,
                   alpha = alpha, gamma = gamma)
    r <- cs_sample_size(d, power)
    c(r$nA, r$power)
  }
  res <- t(vapply(opts, one, numeric(2)))
  data.frame(design = names(opts), n_per_arm = as.numeric(res[, 1]),
             power = res[, 2], row.names = NULL)
}

cs_table6 <- function(gamma = "zero", equal_contingent = FALSE) {
  p <- cs_compare_designs(0.65, 0.50, 6, 18, muB = 15, gamma = gamma,
                          equal_contingent = equal_contingent)
  a <- cs_compare_designs(0.72, 0.64, 3, 18, muB = 18, gamma = gamma,
                          equal_contingent = equal_contingent)
  data.frame(design = p$design, placebo_controlled = p$n_per_arm,
             active_comparator = a$n_per_arm, check.names = FALSE)
}

run_calculator <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE))
    stop("Install the 'shiny' package to use the interface.")
  shiny::runApp(system.file("app", package = "compositeSS"), ...)
}

cs_app <- function(...) run_calculator(...)


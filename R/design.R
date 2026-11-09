cs_design <- function(pA, pB, delta = NULL, sigma = NULL, r = 1,
                      weights = c(0.5, 0.5),
                      method = c("wz", "fisher"),
                      contingent = FALSE, c = 0,
                      calibration = c("recalibrated", "usual"),
                      alpha = 0.025, gamma = "compute",
                      efficacy = c("continuous", "binary"),
                      respA = NULL, respB = NULL,
                      margin_d = 0, margin_e = 0,
                      gate_ties = c("stop", "pass")) {
  method <- match.arg(method)
  calibration <- match.arg(calibration)
  gate_ties <- match.arg(gate_ties)
  efficacy <- match.arg(efficacy)
  .check_prob(pA, "pA"); .check_prob(pB, "pB")
  if (efficacy == "continuous") {
    if (is.null(delta) || is.null(sigma)) stop("Supply 'delta' and 'sigma'.")
    if (sigma <= 0) stop("'sigma' must be positive.")
  } else {
    if (is.null(respA) || is.null(respB))
      stop("Supply 'respA' and 'respB' for a binary efficacy outcome.")
    .check_prob(respA, "respA"); .check_prob(respB, "respB")
  }
  aligned <- identical(weights, "aligned")
  if (aligned) {
    if (efficacy != "binary")
      stop("weights = \"aligned\" needs efficacy = \"binary\": with a continuous ",
           "outcome there is no response probability to align to. Use ",
           "theta_design() for the direct test of theta instead.")
    if (method != "wz")
      stop("Estimand-aligned weights apply to the weighted Z (method = \"wz\").")
    if (contingent)
      stop("Estimand-aligned weights are for the unconditional weighted Z; ",
           "the contingent test uses fixed weights.")
    if (margin_d != 0 || margin_e != 0)
      stop("Estimand-aligned weights are implemented for superiority only.")
    weights <- .aligned_planning_weights(pA, pB, respA, respB, r)
  }
  if (!is.numeric(weights) || length(weights) != 2 || any(weights < 0) ||
      sum(weights) == 0)
    stop("'weights' must be two non-negative numbers, not both zero, ",
         "or \"aligned\".")
  weights <- weights / sqrt(sum(weights^2))
  if (method == "fisher" && abs(weights[1] - weights[2]) > 1e-8)
    stop("Fisher's method is implemented with equal weights only; use method = 'wz'.")
  if (method == "fisher" && (margin_d != 0 || margin_e != 0))
    stop("Non-inferiority margins are implemented for method = 'wz' only.")
  if (!(identical(gamma, "compute") || identical(gamma, "zero") ||
        (is.numeric(gamma) && length(gamma) == 1)))
    stop("'gamma' must be 'compute', 'zero' or a number.")
  if (margin_d < 0 || margin_e < 0) stop("Margins must be non-negative.")
  if (r <= 0) stop("'r' must be positive.")
  if (alpha <= 0 || alpha >= 0.5) stop("'alpha' must be in (0, 0.5).")
  if (!contingent) c <- -Inf
  structure(list(pA = pA, pB = pB, delta = delta, sigma = sigma, r = r,
                 weights = weights, method = method, contingent = contingent,
                 c = c, calibration = calibration, alpha = alpha,
                 gamma = gamma, efficacy = efficacy, respA = respA,
                 respB = respB, margin_d = margin_d, margin_e = margin_e,
                 gate_ties = gate_ties, aligned = aligned),
            class = "cs_design")
}

# Estimand-aligned weights at the planning stage (binary efficacy).
# Delta theta = sbar (pA - pB) + pbar (sA - sB), with sbar and pbar the simple
# averages over arms. Z(d) = (pA - pB) / se_d and Z(e) = (sA - sB) / se_e, so
# weights proportional to (sbar se_d, pbar se_e) make w(d) Z(d) + w(e) Z(e)
# proportional to Delta theta. se_d and se_e are the null standard errors used
# in cs_components(); their ratio does not depend on n, so n = 1 is used.
.aligned_planning_weights <- function(pA, pB, respA, respB, r) {
  pw <- (pA + r * pB) / (1 + r)
  se_d <- sqrt(pw * (1 - pw) * (1 + 1 / r))
  qbar <- (pA * respA + r * pB * respB) / (pA + r * pB)
  se_e <- sqrt(qbar * (1 - qbar) * (1 / pA + 1 / (r * pB)))
  c((respA + respB) / 2 * se_d, (pA + pB) / 2 * se_e)
}

cs_aligned_weights <- function(pA, pB, muA, muB, se_d, se_e, s = 0) {
  for (nm in c("pA", "pB")) .check_prob(get(nm), nm)
  if (se_d <= 0 || se_e <= 0) stop("'se_d' and 'se_e' must be positive.")
  m_bar <- (muA + muB) / 2 - s
  p_bar <- (pA + pB) / 2
  if (m_bar <= 0)
    warning("The completers' mean outcome is not above the dropout score s, so ",
            "keeping patients in the trial lowers theta; the completion weight ",
            "is negative.")
  w <- c(m_bar * se_d, p_bar * se_e)
  list(weights = w / sqrt(sum(w^2)), c_sq = p_bar / m_bar,
       m_bar = m_bar, p_bar = p_bar)
}

print.cs_design <- function(x, ...) {
  cat("Composite-strategy design (Shih-Quan test)\n")
  cat(sprintf("  Completion: pA = %.3f, pB = %.3f; allocation nB/nA = %g\n",
              x$pA, x$pB, x$r))
  if (x$efficacy == "continuous")
    cat(sprintf("  Completers' efficacy: delta = %g, sigma = %g (d = %.3f)\n",
                x$delta, x$sigma, x$delta / x$sigma))
  else
    cat(sprintf("  Completers' response: %.3f vs %.3f\n", x$respA, x$respB))
  cat(sprintf("  Combination: %s, weights w(d) = %.3f, w(e) = %.3f%s\n",
              if (x$method == "wz") "weighted Z" else "Fisher",
              x$weights[1], x$weights[2],
              if (isTRUE(x$aligned)) " (estimand-aligned)" else ""))
  if (x$contingent)
    cat(sprintf("  Contingent: stop if Z(d) %s %g; %s critical value\n",
                if (identical(x$gate_ties, "pass")) "<" else "<=",
                x$c, x$calibration))
  if (x$margin_d > 0 || x$margin_e > 0)
    cat(sprintf("  Non-inferiority margins: completion %g, efficacy %g\n",
                x$margin_d, x$margin_e))
  cat(sprintf("  One-sided alpha = %g\n", x$alpha))
  invisible(x)
}

# Gate of the contingent test. Z(d) is discrete (ties between the numbers of
# completers give Z(d) = 0 exactly), so the rule at Z(d) = c must be stated.
# "stop" (default): testing continues only if Z(d) > c, which keeps the size
# at or below alpha; "pass": continues if Z(d) >= c, which lets ties through
# and gives a size slightly above alpha in finite samples.
.gate_open <- function(zd, c, ties = "stop") {
  if (!is.finite(c)) return(rep(TRUE, length(zd)))
  if (identical(ties, "pass")) zd >= c - 1e-12 else zd > c + 1e-12
}

.check_prob <- function(p, name) {
  if (!is.numeric(p) || length(p) != 1 || p <= 0 || p >= 1)
    stop(sprintf("'%s' must be a single number strictly between 0 and 1.", name))
}

.update_design <- function(design, ...) {
  a <- list(...)
  for (nm in names(a)) design[[nm]] <- a[[nm]]
  design
}


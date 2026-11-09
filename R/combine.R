cs_combine <- function(est_d = NULL, se_d = NULL, est_e = NULL, se_e = NULL,
                       p_d = NULL, p_e = NULL,
                       method = c("wz", "fisher", "simes"),
                       weights = c(0.5, 0.5), c_sq = NULL,
                       contingent = FALSE, c = 0,
                       calibration = c("recalibrated", "usual"),
                       alpha = 0.025, gate_ties = c("stop", "pass")) {
  method <- match.arg(method); calibration <- match.arg(calibration)
  gate_ties <- match.arg(gate_ties)
  zfun <- function(est, se, p, nm) {
    if (!is.null(est) && !is.null(se)) return(est / se)
    if (!is.null(p)) return(stats::qnorm(p, lower.tail = FALSE))
    stop("Give est_", nm, " and se_", nm, ", or p_", nm, ".")
  }
  z_d <- zfun(est_d, se_d, p_d, "d"); z_e <- zfun(est_e, se_e, p_e, "e")
  pd <- stats::pnorm(z_d, lower.tail = FALSE)
  pe <- stats::pnorm(z_e, lower.tail = FALSE)
  if (!is.null(c_sq)) {
    if (is.null(se_d) || is.null(se_e))
      stop("The Shih-Quan weight c_sq needs se_d and se_e.")
    weights <- c(se_d, c_sq * se_e)   # T_G = (est_d + c est_e) / sqrt(se_d^2 + c^2 se_e^2)
  }
  if (method == "simes" && contingent)
    stop("The contingent rule is implemented for 'wz' and 'fisher'.")
  if (method == "fisher") weights <- c(0.5, 0.5)
  des <- cs_design(0.5, 0.5, delta = 0, sigma = 1, weights = weights,
                   method = if (method == "fisher") "fisher" else "wz",
                   contingent = contingent, c = c, calibration = calibration,
                   alpha = alpha)
  gate_open <- !contingent || .gate_open(z_d, c, gate_ties)
  if (method == "simes") {
    stat <- NA_real_; crit <- NA_real_
    p_joint <- min(2 * min(pd, pe), max(pd, pe))
    adj <- c(d = if (pd <= pe) min(2 * pd, pe) else pd,
             e = if (pe < pd) min(2 * pe, pd) else pe)
    reject <- p_joint <= alpha
    reject_ind <- adj <= alpha
  } else {
    cr <- cs_critical(des); crit <- cr$k
    w <- des$weights
    if (method == "wz") {
      stat <- w[1] * z_d + w[2] * z_e
      p_uncond <- stats::pnorm(stat, lower.tail = FALSE)
      exceeds <- stat > crit
    } else {
      stat <- pd * pe
      p_uncond <- stat * (1 - log(stat))
      exceeds <- stat < crit
    }
    p_joint <- if (!contingent) p_uncond else if (!gate_open) 1 else
      .null_size(des, stat)
    reject <- gate_open && exceeds
    reject_ind <- c(d = reject && pd <= alpha, e = reject && pe <= alpha)
  }
  structure(list(method = method, z = c(d = z_d, e = z_e), p = c(d = pd, e = pe),
                 weights = if (method == "wz") des$weights else NULL,
                 statistic = stat, critical = crit, p_joint = p_joint,
                 reject_joint = reject, reject_individual = reject_ind,
                 contingent = contingent, gate_open = gate_open, alpha = alpha),
            class = "cs_combined")
}

print.cs_combined <- function(x, digits = 4, ...) {
  f <- function(v) formatC(v, digits = digits, format = "fg")
  lab <- c(wz = "weighted Z", fisher = "Fisher (Weiler) T_Q",
           simes = "Simes-Hommel T_B")[x$method]
  cat(sprintf("Composite test from summary statistics: %s%s\n", lab,
              if (x$contingent) " (contingent)" else ""))
  cat(sprintf("  Completion: Z = %s, one-sided P = %s\n", f(x$z["d"]), f(x$p["d"])))
  cat(sprintf("  Completers: Z = %s, one-sided P = %s\n", f(x$z["e"]), f(x$p["e"])))
  if (!is.null(x$weights))
    cat(sprintf("  Normalised weights: %s / %s\n", f(x$weights[1]), f(x$weights[2])))
  if (!is.na(x$statistic))
    cat(sprintf("  Statistic %s (critical value %s)\n", f(x$statistic), f(x$critical)))
  cat(sprintf("  Joint P = %s -> %s at one-sided alpha = %g\n", f(x$p_joint),
              if (x$reject_joint) "reject H0" else "do not reject H0", x$alpha))
  cat(sprintf("  Step-down: completion %s, completers %s\n",
              if (x$reject_individual["d"]) "rejected" else "not rejected",
              if (x$reject_individual["e"]) "rejected" else "not rejected"))
  invisible(x)
}


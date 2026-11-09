cs_psi <- function(design, type = c("information", "clinical")) {
  type <- match.arg(type)
  if (design$efficacy != "continuous")
    stop("psi is defined for a continuous efficacy outcome.")
  .psi_fun(design$pA, design$pB, design$delta, design$sigma,
           design$weights, design$r, type)
}

.psi_fun <- function(pA, pB, delta, sigma, weights, r = 1,
                     type = "information") {
  pbar <- (pA + r * pB) / (1 + r)
  h <- (pA - pB) / sqrt(pbar * (1 - pbar))
  d <- delta / sigma
  w <- weights / sqrt(sum(weights^2))
  # information factor: E[Z(e)] / sqrt(n_eff) = d sqrt(pi_w), with
  # n_eff = nA nB / (nA + nB) and pi_w = (1/nA + 1/nB) / (1/(nA pA) + 1/(nB pB)),
  # a weighted harmonic mean of the completion rates (the harmonic mean when
  # r = 1). With this factor psi * sqrt(n_eff) is exactly the weighted Z.
  pw <- (1 + 1 / r) / (1 / pA + 1 / (r * pB))
  f <- if (type == "information") sqrt(pw) else 1
  list(psi = w[1] * h + w[2] * d * f, h = h, d = d, pbar = pbar, pi_w = pw)
}


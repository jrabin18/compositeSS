theta_design <- function(pA, pB, muA, muB, sigma, s = 0, r = 1,
                         alpha = 0.025, margin = 0) {
  .check_prob(pA, "pA"); .check_prob(pB, "pB")
  if (sigma <= 0) stop("'sigma' must be positive.")
  if (margin < 0) stop("'margin' must be non-negative.")
  structure(list(pA = pA, pB = pB, muA = muA, muB = muB, sigma = sigma,
                 s = s, r = r, alpha = alpha, margin = margin),
            class = "theta_design")
}

print.theta_design <- function(x, ...) {
  cat("Expected benefit per randomized patient (composite strategy)\n")
  cat(sprintf("  theta = %.4g (dropout score s = %g)\n", theta_value(x), x$s))
  dc <- theta_decompose(x)
  cat(sprintf("  = completion part %.4g + completers part %.4g\n",
              dc$completion, dc$completers))
  invisible(x)
}

theta_value <- function(design) {
  with(design, (pA * muA + (1 - pA) * s) - (pB * muB + (1 - pB) * s))
}

theta_decompose <- function(design) {
  with(design, {
    comp <- (pA - pB) * ((muA + muB) / 2 - s)
    eff <- ((pA + pB) / 2) * (muA - muB)
    list(completion = comp, completers = eff, theta = comp + eff)
  })
}

.theta_V <- function(p, mu, sigma, s) p * sigma^2 + p * (1 - p) * (mu - s)^2

theta_power <- function(design, n) {
  with(design, {
    VA <- .theta_V(pA, muA, sigma, s); VB <- .theta_V(pB, muB, sigma, s)
    se <- sqrt(VA / n + VB / (r * n))
    pnorm(qnorm(1 - alpha) - (theta_value(design) + margin) / se,
          lower.tail = FALSE)
  })
}

theta_sample_size <- function(design, power = 0.80) {
  with(design, {
    th <- theta_value(design) + margin
    if (th <= 0) stop("theta (+ margin) must be positive under the alternative.")
    VA <- .theta_V(pA, muA, sigma, s); VB <- .theta_V(pB, muB, sigma, s)
    nA <- ceiling((qnorm(1 - alpha) + qnorm(power))^2 * (VA + VB / r) / th^2)
    list(nA = nA, nB = ceiling(r * nA), total = nA + ceiling(r * nA),
         theta = theta_value(design), VA = VA, VB = VB,
         power = theta_power(design, nA))
  })
}

completion_sample_size <- function(pA, pB, r = 1, alpha = 0.025, power = 0.80) {
  pbar <- (pA + r * pB) / (1 + r)
  num <- qnorm(1 - alpha) * sqrt(pbar * (1 - pbar) * (1 + 1 / r)) +
    qnorm(power) * sqrt(pA * (1 - pA) + pB * (1 - pB) / r)
  nA <- ceiling(num^2 / (pA - pB)^2)
  list(nA = nA, nB = ceiling(r * nA), total = nA + ceiling(r * nA))
}


cs_critical <- function(design) {
  alpha <- design$alpha
  if (design$method == "wz") {
    k0 <- qnorm(1 - alpha)
    if (!design$contingent || design$calibration == "usual" ||
        design$weights[2] == 0) {
      k <- k0
    } else {
      k <- uniroot(function(x) .null_size(design, x) - alpha, lower = k0 - 8, upper = k0, tol = 1e-10)$root
    }
    list(k = k, scale = "Z", size = .null_size(design, k))
  } else {
    t0 <- exp(-qchisq(1 - alpha, df = 4) / 2)
    if (!design$contingent || design$calibration == "usual") {
      t <- t0
    } else {
      t <- uniroot(function(x) .null_size(design, x) - alpha, lower = t0, upper = min(0.5, 50 * t0), tol = 1e-12)$root
    }
    list(k = t, scale = "product", size = .null_size(design, t))
  }
}

# Boundary b(z; k): value Z(e) must exceed when Z(d) = z (vectorised in z)
.boundary <- function(z, k, design) {
  if (design$method == "wz") {
    wd <- design$weights[1]; we <- design$weights[2]
    (k - wd * z) / we   # weights are normalised so sqrt(wd^2 + we^2) = 1
  } else {
    pd <- pnorm(z, lower.tail = FALSE)
    q <- k / pd
    out <- rep(-Inf, length(z))
    ok <- q < 1
    out[ok] <- qnorm(q[ok], lower.tail = FALSE)
    out
  }
}

# Size under the composite null: Z(d), Z(e) independent N(0, 1)
.null_size <- function(design, k) {
  c0 <- design$c
  if (design$method == "wz" && design$weights[2] == 0)
    return(pnorm(max(k, c0), lower.tail = FALSE))
  if (design$method == "wz" && !is.finite(c0))
    return(pnorm(k, lower.tail = FALSE))
  lo <- if (is.finite(c0)) c0 else -12
  f <- function(z) dnorm(z) * pnorm(.boundary(z, k, design), lower.tail = FALSE)
  integrate(f, lo, 12, rel.tol = 1e-10, abs.tol = 0,
            subdivisions = 500L)$value
}


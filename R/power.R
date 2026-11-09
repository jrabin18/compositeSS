cs_power <- function(design, n) {
  crit <- cs_critical(design)
  za <- qnorm(1 - design$alpha)
  res <- lapply(n, function(ni) {
    cmp <- cs_components(design, ni)
    data.frame(
      nA = ni, nB = design$r * ni,
      power = .power_one(design, cmp, crit$k),
      power_completion = pnorm((za - cmp$eta_d) / sqrt(cmp$v_d),
                               lower.tail = FALSE),
      power_completers = pnorm(za - cmp$eta_e, lower.tail = FALSE))
  })
  out <- do.call(rbind, res)
  attr(out, "critical") <- crit
  class(out) <- c("cs_power", "data.frame")
  out
}

print.cs_power <- function(x, digits = 4, ...) {
  crit <- attr(x, "critical")
  cat(sprintf("Critical value: %s = %.5g (size %.4f)\n",
              if (crit$scale == "Z") "k" else "product threshold",
              crit$k, crit$size))
  print(as.data.frame(x), digits = digits, row.names = FALSE)
  invisible(x)
}

.power_one <- function(design, cmp, k) {
  wd <- design$weights[1]; we <- design$weights[2]
  sv <- sqrt(cmp$v_d)
  g <- cmp$gamma
  if (design$method == "wz" && we == 0) {
    thr <- max(k, design$c)
    return(pnorm((thr - cmp$eta_d) / sv, lower.tail = FALSE))
  }
  if (design$method == "wz" && !is.finite(design$c)) {
    mu <- wd * cmp$eta_d + we * cmp$eta_e
    tau <- sqrt(wd^2 * cmp$v_d + we^2 + 2 * wd * we * g)
    return(pnorm((k - mu) / tau, lower.tail = FALSE))
  }
  rho2 <- g^2 / cmp$v_d
  if (rho2 >= 1) stop("Implausible gamma: |gamma| must be < sqrt(v(d)).")
  f <- function(z) {
    me <- cmp$eta_e + g * (z - cmp$eta_d) / cmp$v_d
    dnorm(z, cmp$eta_d, sv) *
      pnorm((.boundary(z, k, design) - me) / sqrt(1 - rho2), lower.tail = FALSE)
  }
  lo <- cmp$eta_d - 12 * sv
  if (is.finite(design$c)) lo <- max(lo, design$c)
  hi <- cmp$eta_d + 12 * sv
  if (lo >= hi) return(0)
  integrate(f, lo, hi, rel.tol = 1e-9, abs.tol = 0,
            subdivisions = 500L)$value
}

cs_sample_size <- function(design, power = 0.80, n_max = 1e5) {
  crit <- cs_critical(design)
  pw <- function(n) .power_one(design, cs_components(design, n), crit$k)
  n <- .search_n(pw, power, n_max)
  structure(list(nA = n, nB = ceiling(design$r * n),
                 total = n + ceiling(design$r * n), power = pw(n),
                 target = power, critical = crit, design = design),
            class = "cs_ss")
}

print.cs_ss <- function(x, ...) {
  cat(sprintf("Required sample size: nA = %d, nB = %d (total %d)\n",
              x$nA, x$nB, x$total))
  cat(sprintf("Power at this size: %.4f (target %.2f)\n", x$power, x$target))
  invisible(x)
}

# smallest integer n >= 2 with f(n) >= target (f assumed increasing)
.search_n <- function(f, target, n_max) {
  lo <- 2L
  if (f(lo) >= target) return(lo)
  hi <- 4L
  while (f(hi) < target) {
    lo <- hi; hi <- hi * 2L
    if (hi > n_max) stop("Target power not reached below n_max.")
  }
  while (hi - lo > 1L) {
    mid <- (lo + hi) %/% 2L
    if (f(mid) >= target) hi <- mid else lo <- mid
  }
  hi
}

cs_opt_weights <- function(design, n) {
  cmp <- cs_components(design, n)
  w <- pmax(c(cmp$eta_d, cmp$eta_e), 0)
  if (sum(w) == 0) stop("Neither component favours A under these planning values.")
  w <- w / sum(w)
  d_opt <- .update_design(design, weights = w / sqrt(sum(w^2)))
  list(weights = w,
       power_design = cs_power(design, n)$power,
       power_optimal = cs_power(d_opt, n)$power)
}


# Other composite-strategy analyses for a continuous endpoint, used as
# comparators: trimmed means (Permutt and Li 2017), worst-rank analysis
# with all dropouts tied at the worst rank (Lachin 1999), and the win odds
# for a hierarchy in which earlier dropout is worse than later dropout and
# completers are compared on improvement (Brunner et al. 2021).

.alt_prep <- function(arm, active, completed, y) {
  A <- as.numeric(arm == active)
  if (length(unique(A)) != 2) stop("'arm' must have exactly two groups.")
  comp <- as.logical(completed)
  if (anyNA(comp)) stop("Completion status must be known for every patient.")
  if (any(comp & is.na(y))) stop("Completers with a missing outcome: impute or remove them first.")
  list(A = A, comp = comp, y = y)
}

# Mann-Whitney comparison of scores (higher = better), oriented so that
# positive z favours A. theta = P(A better) + P(tie)/2.
.wmw <- function(sA, sB) {
  nA <- length(sA); nB <- length(sB); N <- nA + nB
  r <- rank(c(sA, sB))
  U <- sum(r[seq_len(nA)]) - nA * (nA + 1) / 2
  theta <- U / (nA * nB)
  tt <- table(r)
  v0 <- nA * nB / 12 * ((N + 1) - sum(tt^3 - tt) / (N * (N - 1)))
  z <- (U - nA * nB / 2) / sqrt(v0)
  # placement values for a variance under the alternative
  rA <- rank(sA); rB <- rank(sB)
  pA <- (r[seq_len(nA)] - rA) / nB          # proportion of B that each A beats (ties 1/2)
  pB <- (r[nA + seq_len(nB)] - rB) / nA     # proportion of A that each B beats
  v1 <- stats::var(pA) / nA + stats::var(pB) / nB
  list(theta = theta, z = z, se = sqrt(v1))
}

cs_worst_rank <- function(arm, active, completed, y, alpha = 0.025) {
  d <- .alt_prep(arm, active, completed, y)
  s <- ifelse(d$comp, d$y, -Inf)
  w <- .wmw(s[d$A == 1], s[d$A == 0])
  list(theta = w$theta, z = w$z, p = stats::pnorm(w$z, lower.tail = FALSE),
       reject = w$z > stats::qnorm(1 - alpha))
}

cs_win_odds <- function(arm, active, completed, y, drop_visit, alpha = 0.025,
                        conf_level = 0.95) {
  d <- .alt_prep(arm, active, completed, y)
  if (any(!d$comp & is.na(drop_visit))) stop("'drop_visit' is needed for every dropout.")
  # dropouts ordered by visit of dropout (earlier = worse), all below completers
  ymin <- min(d$y[d$comp], 0)
  s <- ifelse(d$comp, d$y - ymin + max(drop_visit, na.rm = TRUE) + 1, drop_visit)
  w <- .wmw(s[d$A == 1], s[d$A == 0])
  wo <- w$theta / (1 - w$theta)
  se_log <- w$se / (w$theta * (1 - w$theta))
  zq <- stats::qnorm(1 - (1 - conf_level) / 2)
  list(win_odds = wo, lower = exp(log(wo) - zq * se_log),
       upper = exp(log(wo) + zq * se_log), theta = w$theta, z = w$z,
       p = stats::pnorm(w$z, lower.tail = FALSE),
       reject = w$z > stats::qnorm(1 - alpha))
}

cs_trimmed_mean <- function(arm, active, completed, y, nperm = 1000,
                            alpha = 0.025, seed = NULL) {
  d <- .alt_prep(arm, active, completed, y)
  s <- ifelse(d$comp, d$y, -Inf)
  stat <- function(A) {
    k <- max(1 - mean(d$comp[A == 1]), 1 - mean(d$comp[A == 0]))
    tm <- function(x) {
      m <- length(x) - round(k * length(x))
      if (m < 1) return(NA_real_)
      mean(sort(x, decreasing = TRUE)[seq_len(m)])
    }
    c(tm(s[A == 1]) - tm(s[A == 0]), k)
  }
  obs <- stat(d$A)
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = globalenv()))
      get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv())
            else assign(".Random.seed", old, envir = globalenv()))
    set.seed(seed)
  }
  perm <- vapply(seq_len(nperm), function(i) stat(sample(d$A))[1], numeric(1))
  p <- (1 + sum(perm >= obs[1], na.rm = TRUE)) / (1 + sum(!is.na(perm)))
  list(estimate = obs[1], trim = obs[2], p = p, reject = p <= alpha,
       nperm = nperm)
}

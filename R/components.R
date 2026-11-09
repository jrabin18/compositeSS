cs_components <- function(design, n) {
  with(design, {
    nA <- n; nB <- r * n
    pbar <- (nA * pA + nB * pB) / (nA + nB)
    if (margin_d > 0) {
      se0 <- sqrt(pA * (1 - pA) / nA + pB * (1 - pB) / nB)
      eta_d <- (pA - pB + margin_d) / se0
      v_d <- 1
    } else {
      se0 <- sqrt(pbar * (1 - pbar) * (1 / nA + 1 / nB))
      eta_d <- (pA - pB) / se0
      v_d <- (pA * (1 - pA) / nA + pB * (1 - pB) / nB) / se0^2
    }
    d_std <- .std_efficacy(design)
    mA <- nA * pA; mB <- nB * pB
    eta_e <- d_std / sqrt(1 / mA + 1 / mB)
    g <- if (identical(gamma, "compute")) {
      .gamma_delta(d_std, mA, mB, pA, pB, se0)
    } else if (identical(gamma, "zero")) 0 else gamma
    list(eta_d = eta_d, v_d = v_d, eta_e = eta_e, gamma = g,
         nA = nA, nB = nB, se0 = se0, d_std = d_std)
  })
}

# standardized per-completer efficacy effect (incl. NI margin)
.std_efficacy <- function(design) {
  with(design, {
    if (efficacy == "continuous") {
      (delta + margin_e) / sigma
    } else {
      qbar <- (pA * respA + r * pB * respB) / (pA + r * pB)
      (respA - respB + margin_e) / sqrt(qbar * (1 - qbar))
    }
  })
}

# Delta-method covariance of Z(d) and Z(e) under the alternative.
# Conditioning on the numbers of completers mA, mB, E[Z(e) | m] = d / sqrt(S),
# S = 1/mA + 1/mB. Expanding about the expected numbers and using
# Cov(mA, Z(d)) = pA(1-pA)/se0, Cov(mB, Z(d)) = -pB(1-pB)/se0 gives
# gamma = (d/2) S^{-3/2} [pA(1-pA)/mA^2 - pB(1-pB)/mB^2] / se0,
# which does not depend on n for fixed allocation.
.gamma_delta <- function(d_std, mA, mB, pA, pB, se0) {
  S <- 1 / mA + 1 / mB
  0.5 * d_std * S^(-1.5) * (pA * (1 - pA) / mA^2 - pB * (1 - pB) / mB^2) / se0
}

cs_gamma <- function(design) {
  cs_components(.update_design(design, gamma = "compute"), 100)$gamma
}


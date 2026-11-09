sim_binary <- function(n, pA, pB, sA, sB, seed) {
  set.seed(seed)
  arm <- rep(c("A", "B"), each = n)
  comp <- rbinom(2 * n, 1, ifelse(arm == "A", pA, pB))
  resp <- ifelse(comp == 1, rbinom(2 * n, 1, ifelse(arm == "A", sA, sB)), NA)
  data.frame(arm = arm, comp = comp, resp = resp)
}

test_that("planning weights are proportional to (sbar se_d, pbar se_e)", {
  d <- cs_design(0.65, 0.55, efficacy = "binary", respA = 0.726, respB = 0.58,
                 weights = "aligned")
  expect_true(d$aligned)
  pw <- 0.6
  se_d <- sqrt(pw * (1 - pw) * 2)
  qbar <- (0.65 * 0.726 + 0.55 * 0.58) / 1.2
  se_e <- sqrt(qbar * (1 - qbar) * (1 / 0.65 + 1 / 0.55))
  w <- c((0.726 + 0.58) / 2 * se_d, 0.6 * se_e)
  expect_equal(d$weights, w / sqrt(sum(w^2)), tolerance = 1e-12)
})

test_that("planned drift of the aligned weighted Z is Delta theta over its SE", {
  d <- cs_design(0.65, 0.55, efficacy = "binary", respA = 0.726, respB = 0.58,
                 weights = "aligned")
  n <- 200
  cmp <- cs_components(d, n)
  drift <- sum(d$weights * c(cmp$eta_d, cmp$eta_e))
  sbar <- (0.726 + 0.58) / 2
  se_d <- cmp$se0
  qbar <- (0.65 * 0.726 + 0.55 * 0.58) / 1.2
  se_e <- sqrt(qbar * (1 - qbar) * (1 / (0.65 * n) + 1 / (0.55 * n)))
  dtheta <- 0.65 * 0.726 - 0.55 * 0.58
  expect_equal(drift, dtheta / sqrt(sbar^2 * se_d^2 + 0.6^2 * se_e^2),
               tolerance = 1e-10)
})

test_that("aligned sample size is close to the direct test of Delta theta", {
  # Section 6 reference design: d = 0.40 with sP = 0.58 gives sA = 0.726;
  # the direct two-proportion test needs 213 per arm.
  d <- cs_design(0.65, 0.55, efficacy = "binary", respA = 0.726, respB = 0.58,
                 weights = "aligned")
  n <- cs_sample_size(d, power = 0.90)$nA
  expect_lt(abs(n - 213) / 213, 0.05)
})

test_that("aligned weights are refused where they do not apply", {
  expect_error(cs_design(0.65, 0.55, delta = 6, sigma = 18,
                         weights = "aligned"), "binary")
  expect_error(cs_design(0.65, 0.55, efficacy = "binary", respA = 0.6,
                         respB = 0.5, weights = "aligned", method = "fisher"),
               "weighted Z")
  expect_error(cs_design(0.65, 0.55, efficacy = "binary", respA = 0.6,
                         respB = 0.5, weights = "aligned", contingent = TRUE),
               "contingent")
  expect_error(cs_design(0.65, 0.55, efficacy = "binary", respA = 0.6,
                         respB = 0.5, weights = "aligned", margin_d = 0.05),
               "superiority")
  expect_error(cs_design(0.65, 0.55, delta = 6, sigma = 18,
                         weights = "equal"), "aligned")
})

test_that("cs_analyze with aligned weights tracks the estimated theta", {
  for (seed in 1:5) {
    dat <- sim_binary(300, 0.65, 0.55, 0.48, 0.60, seed)
    a <- cs_analyze(dat, "arm", "A", "comp", "resp", weights = "aligned")
    expect_true(a$test$aligned)
    pA <- a$completion[["A"]]; pB <- a$completion[["B"]]
    sA <- mean(dat$resp[dat$arm == "A" & dat$comp == 1])
    sB <- mean(dat$resp[dat$arm == "B" & dat$comp == 1])
    dtheta <- pA * sA - pB * sB
    # same sign as the estimated difference in completing and responding
    expect_equal(sign(a$test$statistic), sign(dtheta))
    # and close to the direct test of theta
    z_theta <- a$theta[["estimate"]] / a$theta[["se"]]
    expect_lt(abs(a$test$statistic - z_theta), 0.15)
  }
})

test_that("equal weights can favour A when theta favours B; aligned do not", {
  # Paper, Table 2: retention 0.65 vs 0.55 favours A, completers' response
  # 0.48 vs 0.60 favours B, theta 0.312 vs 0.330 favours B; 1,500 per arm,
  # observed values equal to their expectations.
  n <- 1500
  se_d <- sqrt(0.6 * 0.4 * 2 / n)
  q <- (0.65 * 0.48 + 0.55 * 0.60) / 1.2
  se_e <- sqrt(q * (1 - q) * (1 / (0.65 * n) + 1 / (0.55 * n)))
  eq <- cs_combine(est_d = 0.10, se_d = se_d, est_e = -0.12, se_e = se_e)
  aw <- cs_aligned_weights(0.65, 0.55, 0.48, 0.60, se_d, se_e)
  al <- cs_combine(est_d = 0.10, se_d = se_d, est_e = -0.12, se_e = se_e,
                   weights = pmax(aw$weights, 0))
  expect_gt(eq$statistic, 0)
  expect_lt(al$statistic, 0)
})

test_that("cs_analyze refuses aligned weights with Fisher or the contingent rule", {
  dat <- sim_binary(100, 0.65, 0.55, 0.5, 0.4, 3)
  expect_error(cs_analyze(dat, "arm", "A", "comp", "resp", weights = "aligned",
                          method = "fisher"), "weighted Z")
  expect_error(cs_analyze(dat, "arm", "A", "comp", "resp", weights = "aligned",
                          contingent = TRUE), "contingent")
})

test_that("cs_aligned_weights matches cs_combine through c_sq", {
  aw <- cs_aligned_weights(0.65, 0.55, 0.48, 0.60, se_d = 0.049, se_e = 0.056)
  expect_equal(sqrt(sum(aw$weights^2)), 1)
  expect_equal(aw$c_sq, aw$p_bar / aw$m_bar)
  est_d <- 0.10; est_e <- -0.12
  viaw <- cs_combine(est_d = est_d, se_d = 0.049, est_e = est_e, se_e = 0.056,
                     weights = aw$weights)
  viac <- cs_combine(est_d = est_d, se_d = 0.049, est_e = est_e, se_e = 0.056,
                     c_sq = aw$c_sq)
  expect_equal(viaw$statistic, viac$statistic, tolerance = 1e-10)
  # numerator is the estimated Delta theta
  dtheta <- aw$m_bar * est_d + aw$p_bar * est_e
  expect_equal(viaw$statistic,
               dtheta / sqrt((aw$m_bar * 0.049)^2 + (aw$p_bar * 0.056)^2),
               tolerance = 1e-10)
  expect_warning(cs_aligned_weights(0.6, 0.5, -1, -2, 0.05, 0.5), "negative")
})

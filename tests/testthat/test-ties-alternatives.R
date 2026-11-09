test_that("gate: ties at Z(d) = c stop by default and pass on request", {
  r1 <- cs_combine(est_d = 0, se_d = 1, est_e = 3, se_e = 1, contingent = TRUE)
  r2 <- cs_combine(est_d = 0, se_d = 1, est_e = 3, se_e = 1, contingent = TRUE,
                   gate_ties = "pass")
  expect_false(r1$gate_open); expect_false(r1$reject_joint); expect_equal(r1$p_joint, 1)
  expect_true(r2$gate_open); expect_true(r2$reject_joint)
  # equal numbers of completers give Z(d) = 0 exactly in cs_analyze
  set.seed(9)
  arm <- rep(c("A", "B"), each = 100)
  comp <- rep(c(rep(TRUE, 60), rep(FALSE, 40)), 2)
  y <- ifelse(comp, rnorm(200, ifelse(arm == "A", 30, 15), 18), NA)
  dat <- data.frame(arm, comp, y)
  a <- cs_analyze(dat, "arm", "A", "comp", "y", contingent = TRUE)
  expect_equal(unname(a$components["Z_d"]), 0)
  expect_false(a$test$gate_open); expect_false(a$test$reject)
  a2 <- cs_analyze(dat, "arm", "A", "comp", "y", contingent = TRUE, gate_ties = "pass")
  expect_true(a2$test$gate_open); expect_true(a2$test$reject)
})

test_that("strict gate keeps the contingent test at or below alpha", {
  d <- cs_design(.575, .575, delta = 6, sigma = 18, weights = c(.3, .7),
                 contingent = TRUE)
  expect_equal(d$gate_ties, "stop")
  s <- cs_simulate(d, 50, nsim = 40000, null = TRUE, seed = 11)
  dp <- cs_design(.575, .575, delta = 6, sigma = 18, weights = c(.3, .7),
                  contingent = TRUE, gate_ties = "pass")
  sp <- cs_simulate(dp, 50, nsim = 40000, null = TRUE, seed = 11)
  expect_lt(s$reject[1], sp$reject[1])         # same trials: ties make the difference
  expect_lt(s$reject[1], 0.025 + 2 * s$reject[2])
})

test_that("worst-rank equals the Wilcoxon test with dropouts tied worst", {
  set.seed(3)
  arm <- rep(c("A", "B"), c(70, 80))
  comp <- rbinom(150, 1, ifelse(arm == "A", .7, .5)) == 1
  y <- ifelse(comp, rnorm(150, ifelse(arm == "A", 20, 16), 18), NA)
  w <- cs_worst_rank(arm, "A", comp, y)
  s <- ifelse(comp, y, -1e9)
  ref <- suppressWarnings(wilcox.test(s[arm == "A"], s[arm == "B"], alternative = "greater",
                                      exact = FALSE, correct = FALSE))
  expect_equal(w$p, ref$p.value, tolerance = 1e-10)
  expect_equal(w$theta, mean(outer(s[arm == "A"], s[arm == "B"], ">") +
                               0.5 * outer(s[arm == "A"], s[arm == "B"], "==")))
})

test_that("win odds orders earlier dropout below later dropout below completers", {
  arm <- c("A", "A", "B", "B")
  comp <- c(TRUE, FALSE, FALSE, FALSE)
  y <- c(5, NA, NA, NA)
  v <- c(NA, 3, 1, 3)
  wo <- cs_win_odds(arm, "A", comp, y, v)
  # A beats both B (completer) ; A's visit-3 dropout beats B's visit-1, ties visit-3
  expect_equal(wo$theta, (2 + 1 + 0.5) / 4)
  expect_equal(wo$win_odds, wo$theta / (1 - wo$theta))
})

test_that("trimmed mean with no dropout is the difference in means", {
  set.seed(4)
  arm <- rep(c("A", "B"), each = 40)
  y <- rnorm(80, ifelse(arm == "A", 3, 0))
  tm <- cs_trimmed_mean(arm, "A", rep(TRUE, 80), y, nperm = 200, seed = 1)
  expect_equal(tm$trim, 0)
  expect_equal(tm$estimate, mean(y[arm == "A"]) - mean(y[arm == "B"]))
  expect_true(tm$p > 0 && tm$p <= 1)
  # the global random state is left unchanged
  set.seed(5); a <- runif(1); set.seed(5)
  invisible(cs_trimmed_mean(arm, "A", rep(TRUE, 80), y, nperm = 50, seed = 2))
  expect_equal(runif(1), a)
})
test_that("treatment policy effect is zero when the arms are identical", {
  m <- cs_appc_model(k = 6, shape = "quadratic", beta = -0.08, lambda = 0.02)
  cA <- cs_calibrate_appc(m, "A", p = .65, completers_mean = 15, n_cal = 2e4)
  cB <- cs_calibrate_appc(m, "B", p = .65, completers_mean = 15, n_cal = 2e4)
  expect_equal(cs_appc_policy_effect(m, cA, cB, n = 2e4), 0)
  m2 <- cs_appc_model(k = 6, shape = "quadratic", beta = -0.08, lambda = 0.02,
                      depends_on = "previous")
  cA2 <- cs_calibrate_appc(m2, "A", p = .65, completers_mean = 15, n_cal = 2e4)
  expect_equal(cs_appc_policy_effect(m2, cA2, cA2, n = 2e4), 0)
})

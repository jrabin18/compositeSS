test_that("cs_analyze runs and is consistent with the components", {
  set.seed(42)
  n <- 150
  arm <- rep(c("A", "B"), each = n)
  comp <- c(rbinom(n, 1, .65), rbinom(n, 1, .5)) == 1
  y <- ifelse(comp, rnorm(2 * n, ifelse(arm == "A", 21, 15), 18), NA)
  dat <- data.frame(arm, comp, y, base = rnorm(2 * n, 90, 10),
                    site = sample(letters[1:5], 2 * n, TRUE))
  a <- cs_analyze(dat, "arm", "A", "comp", "y")
  expect_s3_class(a, "cs_analysis")
  w <- a$test$weights
  expect_equal(unname(a$test$statistic),
               unname(w[1] * a$components["Z_d"] + w[2] * a$components["Z_e"]))
  a2 <- cs_analyze(dat, "arm", "A", "comp", "y", covariates = "base",
                   strata = "site", contingent = TRUE, bootstrap = 50, seed = 1)
  expect_true(is.finite(a2$psi["se"]))
  expect_length(a2$bootstrap$psi, 2)
})

test_that("simulated power is close to calculated power", {
  d <- cs_design(.65, .5, delta = 6, sigma = 18)
  s <- cs_simulate(d, n = 102, nsim = 2000, seed = 7, muB = 15)
  expect_lt(abs(s$reject[1] - s$calculated_power), 0.04)
  s0 <- cs_simulate(d, n = 102, nsim = 2000, null = TRUE, seed = 8)
  expect_lt(abs(s0$reject[1] - 0.025), 0.015)
})

test_that("Appendix C model calibrates completion and completers' mean", {
  m <- cs_appc_model(k = 6, shape = "quadratic", sd = 18, rho = 0.5,
                     beta = -0.08, lambda = 0.02)
  ca <- cs_calibrate_appc(m, "A", p = 0.65, completers_mean = 21, n_cal = 2e4)
  cb <- cs_calibrate_appc(m, "B", p = 0.50, completers_mean = 15, n_cal = 2e4)
  expect_equal(ca$completion, 0.65, tolerance = 0.005)
  expect_equal(ca$completers_mean, 21, tolerance = 0.01)
  expect_equal(cb$completers_mean, 15, tolerance = 0.01)
  expect_lt(cb$final_mean, 15)   # informative dropout: completers better than all
  tr <- cs_appc_trial(m, ca, cb, 100, 100)
  expect_equal(nrow(tr$patients), 200)
  expect_true(all(is.finite(tr$patients$y_tp_final)))
  expect_true(all(c("unrelated", "lack_of_response") %in% tr$patients$reason))
  d <- cs_design(.65, .5, delta = 6, sigma = 18)
  s <- cs_simulate(d, 102, nsim = 200, dropout = "appendixC", model = m,
                   muB = 15, seed = 3, n_cal = 2e4)
  expect_s3_class(s, "cs_sim")
  expect_true(is.finite(s$treatment_policy_effect))
})


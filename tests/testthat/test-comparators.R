test_that("MMRM and placebo-based MI run on Appendix C data", {
  set.seed(11)
  m <- cs_appc_model(k = 5, shape = "linear", sd = 18, rho = 0.5,
                     beta = -0.08, lambda = 0.03)
  ca <- cs_calibrate_appc(m, "A", 0.65, 21, n_cal = 2e4)
  cb <- cs_calibrate_appc(m, "B", 0.50, 15, n_cal = 2e4)
  tr <- cs_appc_trial(m, ca, cb, 150, 150, post = FALSE)
  f1 <- cs_mmrm(tr$long, 5)
  expect_true(f1$converged && is.finite(f1$estimate) && f1$se > 0)
  f2 <- cs_pmi(tr$long, 5, M = 10)
  expect_true(is.finite(f2$estimate) && f2$se > 0 && f2$df > 1)
})

test_that("cs_simulate reports comparator summaries", {
  m <- cs_appc_model(k = 4, sd = 18, rho = 0.5, beta = -0.05, lambda = 0.02)
  d <- cs_design(.65, .5, delta = 6, sigma = 18)
  s <- cs_simulate(d, 80, nsim = 20, dropout = "appendixC", model = m,
                   muB = 15, n_cal = 1e4, comparators = c("mmrm", "pmi"),
                   M = 5, seed = 4)
  expect_equal(nrow(s$comparators), 2)
  expect_true(all(c("bias_policy", "coverage_policy") %in% names(s$comparators)))
})

test_that("Rubin's rules reduce to the complete-data result with no missing info", {
  r <- compositeSS:::.rubin(rep(2, 5), rep(0.25, 5), n_complete = 100)
  expect_equal(r$estimate, 2); expect_equal(r$se, 0.5)
})


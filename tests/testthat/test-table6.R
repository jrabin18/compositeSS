test_that("Table 6 of the proposal is reproduced", {
  tab <- cs_table6()
  expect_equal(tab$placebo_controlled,
               c(102, 110, 129, 122, 129, 108, 105, 170, 250, 110))
  expect_equal(tab$active_comparator,
               c(330, 350, 419, 399, 422, 347, 337, 533, 834, 366))
})

test_that("gamma is about -0.02 and does not depend on n", {
  d <- cs_design(0.70, 0.50, delta = 0.3, sigma = 1)
  expect_equal(cs_gamma(d), -0.0195, tolerance = 1e-3)
  expect_equal(cs_components(d, 50)$gamma, cs_components(d, 5000)$gamma)
})

test_that("recalibrated contingent critical values give size alpha", {
  for (m in c("wz", "fisher")) {
    w <- if (m == "wz") c(.3, .7) else c(.5, .5)
    d <- cs_design(.65, .5, delta = 6, sigma = 18, weights = w, method = m,
                   contingent = TRUE, c = 0)
    expect_equal(cs_critical(d)$size, 0.025, tolerance = 1e-6)
    du <- cs_design(.65, .5, delta = 6, sigma = 18, weights = w, method = m,
                    contingent = TRUE, c = 0, calibration = "usual")
    expect_lt(cs_critical(du)$size, 0.025)
  }
})

test_that("theta decomposes exactly", {
  td <- theta_design(.65, .5, 21, 15, 18, s = 2)
  dc <- theta_decompose(td)
  expect_equal(dc$theta, theta_value(td))
})

test_that("contingent power with c = -Inf equals closed form", {
  d1 <- cs_design(.65, .5, delta = 6, sigma = 18, weights = c(.3, .7))
  d2 <- d1; d2$contingent <- TRUE; d2$c <- -Inf
  expect_equal(compositeSS:::.power_one(d2, cs_components(d2, 120), qnorm(.975)),
               cs_power(d1, 120)$power, tolerance = 1e-6)
})


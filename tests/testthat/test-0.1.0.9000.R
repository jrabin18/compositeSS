sim_dat <- function(nA, nB, seed) {
  set.seed(seed)
  arm <- rep(c("A", "B"), c(nA, nB))
  comp <- rbinom(nA + nB, 1, ifelse(arm == "A", .65, .5)) == 1
  y <- ifelse(comp, rnorm(nA + nB, ifelse(arm == "A", 21, 15), 18), NA)
  data.frame(arm, comp, y, base = rnorm(nA + nB))
}

test_that("unadjusted theta SE is the unpooled (Welch-type) SE", {
  dat <- sim_dat(60, 240, 1)
  a <- cs_analyze(dat, "arm", "A", "comp", "y")
  v <- ifelse(dat$comp, dat$y, 0)
  se <- sqrt(var(v[dat$arm == "A"]) / 60 + var(v[dat$arm == "B"]) / 240)
  expect_equal(unname(a$theta["se"]), se, tolerance = 1e-10)
  expect_equal(unname(a$theta["estimate"]),
               mean(v[dat$arm == "A"]) - mean(v[dat$arm == "B"]))
})

test_that("psi * sqrt(n_eff) equals the weighted Z (normal approximation)", {
  dat <- sim_dat(100, 200, 2)
  a <- cs_analyze(dat, "arm", "A", "comp", "y", weights = c(.3, .7))
  z_e_normal <- a$completers_effect["estimate"] / a$completers_effect["se"]
  zw <- sum(a$test$weights * c(a$components["Z_d"], z_e_normal))
  expect_equal(unname(a$psi["estimate"]) * sqrt(100 * 200 / 300), unname(zw),
               tolerance = 1e-10)
  # and at the planning level, with r = 2
  d <- cs_design(.65, .5, delta = 6, sigma = 18, r = 2, gamma = "zero")
  cmp <- cs_components(d, 120)
  expect_equal(cs_psi(d)$psi * sqrt(120 * 240 / 360),
               sum(d$weights * c(cmp$eta_d, cmp$eta_e)), tolerance = 1e-10)
})

test_that("Box B1 design rows: equal weights, contingent", {
  tab <- cs_table6(equal_contingent = TRUE)
  i <- grep("^Weighted Z, equal weights, contingent", tab$design)
  expect_equal(tab$placebo_controlled[i], c(102, 103))
  expect_equal(tab$active_comparator[i], c(329, 331))
  expect_equal(nrow(cs_table6()), 10)          # default unchanged
})

test_that("cs_beta_prior reproduces the mean and SD", {
  ab <- cs_beta_prior(0.65, 0.05)
  m <- ab[1] / sum(ab); v <- ab[1] * ab[2] / (sum(ab)^2 * (sum(ab) + 1))
  expect_equal(unname(m), 0.65); expect_equal(unname(sqrt(v)), 0.05)
  expect_error(cs_beta_prior(0.5, 0.6))
})

test_that("assurance sample size is the smallest n reaching the target", {
  d <- cs_design(0.65, 0.50, delta = 6, sigma = 18)
  pa <- cs_beta_prior(.65, .05); pb <- cs_beta_prior(.5, .05)
  r <- cs_assurance_sample_size(d, 0.70, pa, pb, c(6, 2), nsim = 300, seed = 4)
  a_n <- cs_assurance(d, r$nA, pa, pb, c(6, 2), nsim = 300, seed = 4)$assurance
  a_m <- cs_assurance(d, r$nA - 1, pa, pb, c(6, 2), nsim = 300, seed = 4)$assurance
  expect_gte(a_n, 0.70); expect_lt(a_m, 0.70)
  expect_error(cs_assurance_sample_size(d, 0.99, pa, pb, c(6, 6), nsim = 100,
                                        n_max = 2000), "not reached")
  td <- theta_design(.65, .5, 21, 15, 18)
  rt <- cs_assurance_sample_size(td, 0.70, pa, pb, c(6, 2), nsim = 200)
  expect_true(rt$assurance >= 0.70)
})

test_that("theta sensitivity matches cs_analyze at each s", {
  dat <- sim_dat(120, 120, 3)
  ts <- cs_theta_sensitivity(dat, "arm", "A", "comp", "y", s = c(0, -10))
  for (k in 1:2) {
    a <- cs_analyze(dat, "arm", "A", "comp", "y", s = ts$s[k])
    expect_equal(ts$theta[k], unname(a$theta["estimate"]))
    expect_equal(ts$se[k], unname(a$theta["se"]))
  }
  expect_equal(ts$theta, ts$completion_part + ts$completers_part)
  tsa <- cs_theta_sensitivity(dat, "arm", "A", "comp", "y", s = 0, covariates = "base")
  expect_true(attr(tsa, "adjusted"))
})

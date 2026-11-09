# Worked examples of Shih and Quan (1997), Stat Med 16:1225-1239 (alpha = 0.05).
# Two misprints in the paper are corrected: Example 2's T_Q uses
# 0.0062 x 0.4556 (printed "0.062 x 0.04556"), and Example 3's optimal c is
# 0.4 (Table III prints 0.04; the text gives 0.4, which reproduces T_G = 2.08).

test_that("Example 1 (Heyting et al., antidepressant vs placebo)", {
  b <- 30 / 123 - 26 / 133                 # placebo minus active dropout
  sb <- b / qnorm(0.172, lower.tail = FALSE)
  q <- cs_combine(p_d = 0.172, p_e = 0.042, method = "fisher", alpha = 0.05)
  expect_equal(round(q$p_joint, 3), 0.043)
  expect_true(q$reject_joint)
  expect_equal(unname(q$reject_individual), c(FALSE, TRUE))
  s <- cs_combine(p_d = 0.172, p_e = 0.042, method = "simes", alpha = 0.05)
  expect_false(s$reject_joint)
  g <- cs_combine(est_d = b, se_d = sb, est_e = 1.9, se_e = 1.1, c_sq = 0.09,
                  alpha = 0.05)
  expect_equal(round(g$statistic, 2), 1.97)
  expect_equal(round(g$p_joint, 3), 0.024)
})

test_that("Example 2 (Gould, two active treatments)", {
  q <- cs_combine(p_d = 0.0062, p_e = 0.4556, method = "fisher", alpha = 0.05)
  expect_equal(round(q$p_joint, 4), 0.0194)
  expect_equal(unname(q$reject_individual), c(TRUE, FALSE))
  s <- cs_combine(p_d = 0.0062, p_e = 0.4556, method = "simes", alpha = 0.05)
  expect_true(s$reject_joint)
  expect_equal(unname(s$reject_individual), c(TRUE, FALSE))
  b <- 72 / 138 - 49 / 131; sb <- b / qnorm(0.0062, lower.tail = FALSE)
  l <- -log(0.965); sl <- l / qnorm(0.4556, lower.tail = FALSE)
  g <- cs_combine(est_d = b, se_d = sb, est_e = l, se_e = sl, c_sq = 0.00894,
                  alpha = 0.05)
  expect_equal(round(g$statistic, 2), 2.50, tolerance = 0.011)  # paper 2.51
  expect_equal(round(g$p_joint, 3), 0.006)
})

test_that("Example 3 (two active treatments, BPRS, baseline-adjusted)", {
  args <- list(est_d = 0.4890, se_d = 0.3750, est_e = 1.8991, se_e = 1.1717,
               alpha = 0.05)
  q <- do.call(cs_combine, c(args, method = "fisher"))
  expect_equal(round(q$p, 3), c(d = 0.096, e = 0.053))
  expect_equal(round(q$p_joint, 3), 0.032)
  expect_equal(unname(q$reject_individual), c(FALSE, FALSE))
  s <- do.call(cs_combine, c(args, method = "simes"))
  expect_false(s$reject_joint)
  g <- do.call(cs_combine, c(args, c_sq = 0.4))
  expect_equal(round(g$statistic, 2), 2.08)
  expect_equal(round(g$p_joint, 3), 0.019)
  e <- do.call(cs_combine, args)            # package default: equal weights
  expect_equal(round(e$statistic, 2), 2.07)
})

test_that("cs_combine agrees with cs_critical for the contingent test", {
  r <- cs_combine(p_d = 0.2, p_e = 0.01, weights = c(0.3, 0.7),
                  contingent = TRUE)
  d <- cs_design(0.5, 0.5, delta = 0, sigma = 1, weights = c(0.3, 0.7),
                 contingent = TRUE)
  expect_equal(r$critical, cs_critical(d)$k)
  expect_false(cs_combine(p_d = 0.7, p_e = 0.001, contingent = TRUE)$reject_joint)
})


test_that("no correlation between visits gives the completers-only variance", {
  des <- cs_mmrm_design(0.4, 0.65, 0.55, structure = "ar1", rho = 0)
  expect_equal(des$phiA, 1 / 0.65)
  expect_equal(des$phiB, 1 / 0.55)
})

test_that("almost no dropout gives the complete-data size", {
  des <- cs_mmrm_design(0.4, 0.999999, 0.999999)
  expect_equal(des$phiA, 1, tolerance = 1e-5)
  n <- cs_mmrm_sample_size(des, power = 0.90)$nA
  expect_equal(n, ceiling((qnorm(0.975) + qnorm(0.90))^2 * 2 / 0.16))
})

test_that("near-perfect correlation needs only the first visit", {
  des <- cs_mmrm_design(0.4, 0.65, 0.55, structure = "cs", rho = 0.999999,
                        pattern = "even")
  r1 <- 1 - 0.35 / 6
  expect_equal(des$phiA, 1 / r1, tolerance = 1e-4)
})

test_that("the NEWMEDS design lies between complete data and completers only", {
  des <- cs_mmrm_design(0.4, 0.65, 0.55)
  ss <- cs_mmrm_sample_size(des, power = 0.90)
  full <- cs_mmrm_sample_size(cs_mmrm_design(0.4, 0.999999, 0.999999))$nA
  comp <- cs_mmrm_sample_size(cs_mmrm_design(0.4, 0.65, 0.55, structure = "ar1",
                                             rho = 0))$nA
  expect_gt(ss$nA, full)
  expect_lt(ss$nA, comp)
  expect_equal(ss$total, 346)
  expect_gte(ss$power, 0.90)
})

test_that("earlier dropout needs more patients than later dropout", {
  n <- function(pt) cs_mmrm_sample_size(cs_mmrm_design(0.4, 0.65, 0.55,
                                                       pattern = pt))$total
  expect_gt(n("early"), n("even"))
  expect_gt(n("even"), n("late"))
})

test_that("variance factors reproduce Table 2 of Lu, Luo and Chen (2008)", {
  # AR(1) with correlation p1J between first and last visit, retention falling
  # linearly from 1 at visit 1 to 1 - aJ at visit J. Published values as
  # transcribed in the tests of the longpower package (power.mmrm.ar1).
  pub <- matrix(c(
    1.11111111111111, 1.11, 1.10111111111111, 1.08333333333333, 1.05666666666667, 1.02111111111111,
    1.25, 1.2475, 1.2275, 1.1875, 1.1275, 1.0475,
    1.42857142857143, 1.42428571428571, 1.39, 1.32142857142857, 1.21857142857143, 1.08142857142857,
    1.66666666666667, 1.66, 1.60666666666667, 1.5, 1.34, 1.12666666666667,
    1.11111111111111, 1.10050206679506, 1.08280448732802, 1.06283003239702, 1.03996686151349, 1.01408550239776,
    1.25, 1.2247445108099, 1.18392908471514, 1.13886229048648, 1.08796962333115, 1.03089517940421,
    1.42857142857143, 1.38254161826624, 1.31065330422472, 1.23318260782795, 1.14702199890211, 1.0514263712482,
    1.66666666666667, 1.59010470647863, 1.47488822527398, 1.35405269950478, 1.221989488974, 1.07728040119987),
    8, 6, byrow = TRUE)
  i <- 0
  for (J in c(2, 4)) for (aJ in c(0.1, 0.2, 0.3, 0.4)) {
    i <- i + 1
    for (k in 1:6) {
      p1J <- c(0, 0.1, 0.3, 0.5, 0.7, 0.9)[k]
      rho <- if (p1J == 0) 0 else p1J^(1 / (J - 1))
      R <- rho^abs(outer(seq_len(J), seq_len(J), "-"))
      expect_equal(compositeSS:::.mmrm_phi(R, seq(1, 1 - aJ, length = J)),
                   pub[i, k], tolerance = 1e-8)
    }
  }
})

test_that("sample size reproduces the longpower example of Iddi and Donohue (2022)", {
  # R Journal 14(1), Section 3.1: power.mmrm() with four visits, retention
  # 1, 0.9, 0.8, 0.7, exchangeable correlation 0.25, delta 0.5, two-sided
  # alpha 0.05, power 0.80. Published n per group: 86.99175.
  R <- matrix(0.25, 4, 4); diag(R) <- 1
  phi <- compositeSS:::.mmrm_phi(R, c(1, 0.9, 0.8, 0.7))
  n <- (qnorm(0.975) + qnorm(0.80))^2 * 2 * phi / 0.5^2
  expect_equal(n, 86.99175, tolerance = 1e-7)
})

test_that("invalid inputs are refused", {
  expect_error(cs_mmrm_design(0.4, 0.65, 0.55, visits = 8), "6 weekly visits")
  expect_error(cs_mmrm_design(0.4, 0.65, 0.55, structure = "ar1"), "rho")
  expect_error(cs_mmrm_design(0, 0.65, 0.55), "non-zero")
})

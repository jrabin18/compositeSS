test_that("adjusted Z(d) and Z(e) are independent under H0 (Shih and Quan 2001)", {
  set.seed(2001)
  nsim <- 400; n <- 150
  z <- matrix(NA_real_, nsim, 2)
  for (i in seq_len(nsim)) {
    A <- rep(1:0, each = n); g <- rnorm(2 * n)
    comp <- runif(2 * n) < plogis(0.4 + 0.8 * g)
    y <- 10 + 4 * g + 5 * (rchisq(2 * n, 3) - 3) / sqrt(6)
    dat <- data.frame(arm = ifelse(A == 1, "A", "B"), comp = comp,
                      y = ifelse(comp, y, NA), g = g)
    a <- cs_analyze(dat, "arm", "A", "comp", "y", covariates = "g",
                    completion_model = "logistic")
    z[i, ] <- a$components[c("Z_d", "Z_e")]
  }
  expect_lt(abs(cor(z[, 1], z[, 2])), 0.15)
  expect_lt(abs(mean(z[, 1])), 0.2); expect_lt(abs(sd(z[, 1]) - 1), 0.12)
  expect_lt(abs(mean(z[, 2])), 0.2); expect_lt(abs(sd(z[, 2]) - 1), 0.12)
})

test_that("cs_sim_sq2001 gives size close to alpha", {
  r <- cs_sim_sq2001(nsim = 1500, seed = 3)
  expect_lt(abs(r$cor_Zd_Ze), 0.1)
  expect_lt(abs(r$size_weighted_z - 0.025), 0.015)
  expect_lt(abs(r$size_contingent - 0.025), 0.015)
})


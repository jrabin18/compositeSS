# Head-to-head simulation of Section 7.3 / Table 2B of the manuscript.
# Usage: Rscript headtohead.R <scenario> <nsim> <nperm> [seed offset]
#   scenario: random, mnar_weak, mnar_strong, mar, null_mnar
# The null scenario uses one calibration for both arms, so the arms are identical.
# Results are written to <scenario>[_offset].rds in the working directory.

args <- commandArgs(TRUE); id <- args[1]; nsim <- as.integer(args[2]); nperm <- as.integer(args[3])
library(compositeSS)
sc <- list(random = list(beta = 0, dep = "current"), mnar_weak = list(beta = -0.04, dep = "current"),
           mnar_strong = list(beta = -0.12, dep = "current"), mar = list(beta = -0.08, dep = "previous"),
           null_mnar = list(beta = -0.08, dep = "current"))[[id]]
isnull <- id == "null_mnar"
n <- 102; pA <- .65; pB <- if (isnull) .65 else .50; muB <- 15; muA <- if (isnull) 15 else 21
m <- cs_appc_model(k = 6, shape = "quadratic", sd = 18, rho = 0.5, beta = sc$beta, lambda = 0.02, depends_on = sc$dep)
seed <- 2026
calB <- cs_calibrate_appc(m, "B", pB, muB, n_cal = 1e5, seed = seed)
calA <- if (isnull) { ca <- calB; ca$arm <- "A"; ca } else cs_calibrate_appc(m, "A", pA, muA, n_cal = 1e5, seed = seed + 7)
set.seed(seed + 1 + as.integer(if (length(args) > 3) args[4] else 0))
kz <- qnorm(.975)
k37 <- cs_critical(cs_design(.65, .5, delta = 6, sigma = 18, weights = c(.3, .7), contingent = TRUE))$k
kf <- cs_critical(cs_design(.65, .5, delta = 6, sigma = 18, method = "fisher"))$k
w37 <- c(.3, .7) / sqrt(.58)
res <- matrix(NA_real_, nsim, 12, dimnames = list(NULL, c("wz55", "wz37c", "fisher", "theta", "worst", "winodds", "trimmed",
                                                          "zd", "ze", "sd_compl", "tm_est", "wo")))
for (i in seq_len(nsim)) {
  tr <- cs_appc_trial(m, calA, calB, n, n, post = FALSE)$patients
  A <- tr$arm == "A"; cp <- tr$completed; y <- tr$y_final
  mA <- sum(cp & A); mB <- sum(cp & !A); pa <- mA / n; pb <- mB / n; pbar <- (mA + mB) / (2 * n)
  zd <- (pa - pb) / sqrt(pbar * (1 - pbar) * 2 / n)
  yA <- y[cp & A]; yB <- y[cp & !A]; dfe <- mA + mB - 2
  sp <- sqrt(((mA - 1) * var(yA) + (mB - 1) * var(yB)) / dfe)
  ze <- qnorm(pt((mean(yA) - mean(yB)) / (sp * sqrt(1 / mA + 1 / mB)), dfe))
  vA <- ifelse(cp[A], y[A], 0); vB <- ifelse(cp[!A], y[!A], 0)
  zt <- (mean(vA) - mean(vB)) / sqrt(var(vA) / n + var(vB) / n)
  wr <- cs_worst_rank(tr$arm, "A", cp, y)
  wo <- cs_win_odds(tr$arm, "A", cp, y, tr$drop_visit)
  tm <- cs_trimmed_mean(tr$arm, "A", cp, y, nperm = nperm)
  res[i, ] <- c((zd + ze) / sqrt(2) > kz,
                zd > 0 && (w37[1] * zd + w37[2] * ze) > k37,
                pnorm(zd, lower.tail = FALSE) * pnorm(ze, lower.tail = FALSE) < kf,
                zt > kz, wr$reject, wo$reject, tm$reject, zd, ze, sp, tm$estimate, wo$win_odds)
}
saveRDS(list(res = res, id = id, calA = calA, calB = calB, model = m, n = n), sprintf("%s%s.rds", id, if (length(args) > 3) paste0("_", args[4]) else ""))
cat("done", id, "\n")

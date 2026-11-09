# Comparisons of the planning formulas with published values.
# Run from the package root: Rscript inst/validation/published_comparisons.R
suppressMessages(devtools::load_all(".", quiet = TRUE))

# 1. MMRM variance factors: Lu, Luo and Chen (2008), Table 2, as transcribed in
#    the tests of the longpower package. AR(1) with correlation p1J between
#    the first and last visits; retention falls linearly from 1 to 1 - aJ.
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
ours <- pub * NA; i <- 0
for (J in c(2, 4)) for (aJ in c(0.1, 0.2, 0.3, 0.4)) {
  i <- i + 1
  for (k in 1:6) {
    p1J <- c(0, 0.1, 0.3, 0.5, 0.7, 0.9)[k]
    rho <- if (p1J == 0) 0 else p1J^(1 / (J - 1))
    ours[i, k] <- compositeSS:::.mmrm_phi(rho^abs(outer(1:J, 1:J, "-")),
                                          seq(1, 1 - aJ, length = J))
  }
}
cat(sprintf("Lu, Luo and Chen (2008) Table 2: 48 variance factors, largest difference %.1e\n",
            max(abs(ours - pub))))

# 1b. MMRM sample size: Iddi and Donohue (2022), R Journal 14(1), Section 3.1,
#     the worked example for power.mmrm() in the longpower package and its
#     Shiny app. Four visits, retention 1, 0.9, 0.8, 0.7, exchangeable
#     correlation 0.25, effect size 0.5, two-sided alpha 0.05, power 0.80.
#     Published: 86.99175 per group.
R4 <- matrix(0.25, 4, 4); diag(R4) <- 1
phi4 <- compositeSS:::.mmrm_phi(R4, c(1, 0.9, 0.8, 0.7))
n4 <- (qnorm(0.975) + qnorm(0.80))^2 * 2 * phi4 / 0.5^2
cat(sprintf("Iddi and Donohue (2022) power.mmrm example: %.5f per group (published 86.99175)\n", n4))

# 2. Tied worst-rank Wilcoxon power: Matsouaka and Betensky (2015), Table 2,
#    analytical power, 50 per arm, two-sided alpha 0.05. Survival q2 in the
#    active arm, control survival q2^HR, outcome N(0, 1) vs N(sqrt(2) Delta, 1),
#    so the completers' effect size is d = sqrt(2) Delta. Our formula is the
#    Noether-type approximation of Section 6.1 (null variance with ties).
wr_power <- function(pA, pB, d, n, alpha = 0.05) {
  pi_ <- pA * (1 - pB) + pA * pB * pnorm(d / sqrt(2)) + 0.5 * (1 - pA) * (1 - pB)
  t <- 1 - (pA + pB) / 2
  z <- (pi_ - 0.5) / sqrt((1 - t^3) / (6 * n))
  pnorm(z - qnorm(1 - alpha / 2)) + pnorm(-z - qnorm(1 - alpha / 2))
}
mb <- rbind(
  data.frame(HR = 1.0, q2 = 0.6, pub = c(.05, .06, .08, .12, .17, .23, .31)),
  data.frame(HR = 1.0, q2 = 0.8, pub = c(.05, .07, .14, .26, .40, .57, .71)),
  data.frame(HR = 1.4, q2 = 0.6, pub = c(.17, .23, .29, .37, .45, .53, .60)),
  data.frame(HR = 1.4, q2 = 0.8, pub = c(.09, .17, .28, .43, .59, .73, .84)),
  data.frame(HR = 2.0, q2 = 0.6, pub = c(.61, .67, .73, .79, .83, .87, .90)),
  data.frame(HR = 2.0, q2 = 0.8, pub = c(.28, .42, .56, .70, .81, .89, .94)))
mb$Delta <- rep(seq(0, 0.6, 0.1), 6)
mb$ours <- with(mb, round(wr_power(q2, q2^HR, sqrt(2) * Delta, 50), 2))
cat(sprintf("Matsouaka and Betensky (2015) Table 2: 42 tied worst-rank powers, largest difference %.2f, mean %.3f\n",
            max(abs(mb$ours - mb$pub)), mean(abs(mb$ours - mb$pub))))
print(mb[abs(mb$ours - mb$pub) > 0.005, ], row.names = FALSE)

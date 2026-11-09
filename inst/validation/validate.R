# Validation report for compositeSS: observed values, not pass/fail flags.
# Run: Rscript inst/validation/validate.R  (about 3 minutes)
library(compositeSS)
cat("compositeSS", as.character(packageVersion("compositeSS")), "-", R.version.string, "\n\n")

cat("1. Table 6 (gamma = 0, as published) with the Box B1 rows, and with gamma computed\n")
t0 <- cs_table6(equal_contingent = TRUE)
t1 <- cs_table6("compute", equal_contingent = TRUE)
print(data.frame(design = t0$design, placebo = t0$placebo_controlled,
                 placebo_gamma = t1$placebo_controlled, active = t0$active_comparator,
                 active_gamma = t1$active_comparator), row.names = FALSE)

cat("\n2. gamma (closed form)\n")
for (sc in list(c(.70, .50, .3, 1), c(.65, .50, 6, 18), c(.72, .64, 3, 18))) {
  d <- cs_design(sc[1], sc[2], delta = sc[3], sigma = sc[4])
  cat(sprintf("  pA %.2f pB %.2f d %.3f: gamma %.5f (n = 50: %.5f, n = 5000: %.5f)\n",
              sc[1], sc[2], sc[3] / sc[4], cs_gamma(d), cs_components(d, 50)$gamma,
              cs_components(d, 5000)$gamma))
}

cat("\n3. Critical values and exact size\n")
for (o in list(list("wz", c(.3, .7), "recalibrated"), list("wz", c(.3, .7), "usual"),
               list("wz", c(.5, .5), "recalibrated"), list("fisher", c(.5, .5), "recalibrated"),
               list("fisher", c(.5, .5), "usual"))) {
  d <- cs_design(.65, .5, delta = 6, sigma = 18, method = o[[1]], weights = o[[2]],
                 contingent = TRUE, calibration = o[[3]])
  cr <- cs_critical(d)
  cat(sprintf("  %-6s %.1f/%.1f contingent, %-12s k = %.6g (%s scale), size %.7f\n",
              o[[1]], o[[2]][1] / sum(o[[2]]), o[[2]][2] / sum(o[[2]]), o[[3]], cr$k,
              cr$scale, cr$size))
}

cat("\n4. Simulated vs calculated power and type I error (random dropout, 100,000 trials each)\n")
runs <- list(
  list("WZ 0.5/0.5", list(weights = c(.5, .5)), 102),
  list("WZ 0.3/0.7 contingent recal.", list(weights = c(.3, .7), contingent = TRUE), 122),
  list("WZ 0.5/0.5 contingent recal.", list(weights = c(.5, .5), contingent = TRUE), 102),
  list("Fisher", list(method = "fisher"), 108),
  list("Fisher contingent recal.", list(method = "fisher", contingent = TRUE), 105))
for (r in runs) {
  d <- do.call(cs_design, c(list(.65, .5, delta = 6, sigma = 18), r[[2]]))
  s1 <- cs_simulate(d, r[[3]], nsim = 1e5, seed = 11, muB = 15)
  s0 <- cs_simulate(d, r[[3]], nsim = 1e5, seed = 12, null = TRUE)
  cat(sprintf("  %-30s n = %3d  power sim %.4f (MC SE %.4f) calc %.4f | type I %.4f (MC SE %.4f) | psi cov %.4f theta cov %.4f\n",
              r[[1]], r[[3]], s1$reject[1], s1$reject[2], s1$calculated_power,
              s0$reject[1], s0$reject[2], s1$psi_coverage, s1$theta_coverage))
}

cat("\n5. Coverage of 95% CIs for psi and theta, unequal allocation (20,000 trials each)\n")
for (r in c(0.5, 2)) {
  d <- cs_design(.65, .5, delta = 6, sigma = 18, r = r)
  nA <- if (r == 2) 60 else 240
  s <- cs_simulate(d, nA, nsim = 2e4, seed = 21, muB = 15)
  cat(sprintf("  nA = %d, nB = %d: psi coverage %.4f, theta coverage %.4f (MC SE %.4f)\n",
              s$n["A"], s$n["B"], s$psi_coverage, s$theta_coverage,
              sqrt(.95 * .05 / s$nsim)))
}

cat("\n6. Assurance, placebo-controlled scenario (priors: completion SD 0.05, delta ~ N(6, 2))\n")
d <- cs_design(.65, .5, delta = 6, sigma = 18)
pa <- cs_beta_prior(.65, .05); pb <- cs_beta_prior(.50, .05)
a <- cs_assurance(d, 102, pa, pb, c(6, 2), nsim = 4000, seed = 31)
cat(sprintf("  n = 102: power at planning values %.4f, assurance %.4f\n",
            a$power_at_planning_values, a$assurance))
n70 <- cs_assurance_sample_size(d, 0.70, pa, pb, c(6, 2), nsim = 2000, seed = 31)
n80 <- tryCatch(cs_assurance_sample_size(d, 0.80, pa, pb, c(6, 2), nsim = 2000, seed = 31),
                error = function(e) conditionMessage(e))
cat(sprintf("  n for assurance 0.70: %d per arm\n", n70$nA))
cat("  n for assurance 0.80:", if (is.list(n80)) n80$nA else n80, "\n")

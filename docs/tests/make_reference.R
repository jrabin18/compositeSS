# Reference values for the JavaScript port (assets/compositess.js).
# Computes results with the compositeSS R package source and writes
# tests/reference.json, which tests/validate.html compares against the
# JavaScript results. Also writes the synthetic example data set
# (example_trial.csv, no real patient data) used by the Analyze tab.
#
# Usage, from the repository root: Rscript docs/tests/make_reference.R
# (the package source is two folders up; another path can be given as the
# first argument)

args <- commandArgs(trailingOnly = TRUE)
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
pkg <- if (length(args)) args[1] else file.path(here, "..", "..")
for (f in list.files(file.path(pkg, "R"), pattern = "\\.R$", full.names = TRUE))
  source(f, local = globalenv())

js <- function(x) {
  if (is.null(x)) return("null")
  if (is.list(x)) {
    nm <- names(x)
    body <- vapply(seq_along(x), function(i) js(x[[i]]), "")
    if (is.null(nm)) return(paste0("[", paste(body, collapse = ","), "]"))
    return(paste0("{", paste0('"', nm, '":', body, collapse = ","), "}"))
  }
  one <- function(v) {
    if (is.character(v)) return(paste0('"', gsub('"', '\\\\"', gsub("\n", "\\\\n", v)), '"'))
    if (is.logical(v)) return(if (is.na(v)) "null" else tolower(as.character(v)))
    if (is.na(v) || !is.finite(v)) return(if (!is.na(v) && v == Inf) '"Inf"' else if (!is.na(v) && v == -Inf) '"-Inf"' else "null")
    formatC(v, digits = 17, format = "g")
  }
  if (length(x) == 1) one(x) else paste0("[", paste(vapply(x, one, ""), collapse = ","), "]")
}

out <- list()

# ---- distribution functions
zs <- c(-8, -3.2, -1.5, -0.3, 0, 0.2, 0.9, 1.96, 2.7, 4.5, 7)
ps <- c(1e-10, 1e-4, 0.001, 0.025, 0.2, 0.5, 0.8, 0.975, 0.999)
tt <- expand.grid(t = c(-4, -1.2, 0.3, 2.1, 6), df = c(3, 17, 98, 400))
out$dist <- list(
  z = zs, pnorm = pnorm(zs), pnorm_upper = pnorm(zs, lower.tail = FALSE),
  p = ps, qnorm = qnorm(ps),
  t = tt$t, df = tt$df, pt = pt(tt$t, tt$df), pt_z = qnorm(pt(tt$t, tt$df)),
  alpha = c(0.025, 0.05, 0.01), qchisq4 = qchisq(1 - c(0.025, 0.05, 0.01), 4))

# ---- composite designs
cases <- list()
add <- function(args, n_fixed = 80, power = 0.8) {
  d <- do.call(cs_design, args)
  cr <- cs_critical(d)
  ss <- tryCatch(cs_sample_size(d, power), error = function(e) NULL)
  pw <- cs_power(d, n_fixed)
  res <- list(args = args, power_target = power, n_fixed = n_fixed,
              weights = d$weights, k = cr$k, size = cr$size,
              gamma = cs_components(d, n_fixed)$gamma,
              nA = if (is.null(ss)) NA else ss$nA,
              power_at_nA = if (is.null(ss)) NA else ss$power,
              power = pw$power, power_completion = pw$power_completion,
              power_completers = pw$power_completers)
  if (d$method == "wz" && !d$aligned) {
    ow <- tryCatch(cs_opt_weights(d, n_fixed), error = function(e) NULL)
    if (!is.null(ow)) res$opt <- list(weights = ow$weights, power_optimal = ow$power_optimal)
  }
  if (d$efficacy == "continuous") {
    res$psi_info <- cs_psi(d)$psi; res$psi_clin <- cs_psi(d, "clinical")$psi
  }
  res$statement <- capture.output(invisible(cs_estimand_statement(d)))
  cases[[length(cases) + 1]] <<- res
}
base <- list(pA = 0.65, pB = 0.50, delta = 6, sigma = 18)
for (method in c("wz", "fisher")) for (cont in c(FALSE, TRUE))
  for (cal in c("recalibrated", "usual")) for (gm in c("compute", "zero")) {
    if (!cont && cal == "usual") next
    a <- c(base, list(method = method, contingent = cont, calibration = cal, gamma = gm))
    if (method == "wz") a$weights <- c(0.3, 0.7)
    add(a)
  }
add(c(base, list(weights = c(0.7, 0.3), r = 2, alpha = 0.05)))
add(c(base, list(weights = c(0, 1))))
add(c(base, list(weights = c(1, 0))))
add(c(base, list(weights = c(0.5, 0.5), contingent = TRUE, c = 0.5)))
add(c(base, list(weights = c(0.5, 0.5), margin_d = 0.05, margin_e = 2)))
add(list(pA = 0.72, pB = 0.64, delta = 3, sigma = 18, weights = c(0.3, 0.7), gamma = 1e-3))
add(list(pA = 0.6, pB = 0.62, delta = 5, sigma = 15), power = 0.9)
add(list(pA = 0.65, pB = 0.55, efficacy = "binary", respA = 0.55, respB = 0.40))
add(list(pA = 0.65, pB = 0.55, efficacy = "binary", respA = 0.55, respB = 0.40, weights = "aligned"))
add(list(pA = 0.55, pB = 0.45, efficacy = "binary", respA = 0.50, respB = 0.42, weights = "aligned", r = 0.5))
add(list(pA = 0.65, pB = 0.55, efficacy = "binary", respA = 0.55, respB = 0.40, method = "fisher"))
add(list(pA = 0.65, pB = 0.55, efficacy = "binary", respA = 0.55, respB = 0.40, contingent = TRUE, weights = c(0.4, 0.6)))
out$cs <- cases

# ---- theta
th <- list()
for (a in list(list(pA = 0.65, pB = 0.5, muA = 21, muB = 15, sigma = 18),
               list(pA = 0.65, pB = 0.5, muA = 21, muB = 15, sigma = 18, s = -10, r = 2),
               list(pA = 0.72, pB = 0.64, muA = 21, muB = 18, sigma = 18, margin = 2, alpha = 0.05))) {
  d <- do.call(theta_design, a)
  ss <- theta_sample_size(d, 0.8)
  dc <- theta_decompose(d)
  th[[length(th) + 1]] <- list(args = a, nA = ss$nA, power_at_nA = ss$power,
    power_120 = theta_power(d, 120), theta = dc$theta, completion = dc$completion,
    completers = dc$completers, statement = capture.output(invisible(cs_estimand_statement(d))))
}
out$theta <- th
out$completion <- list(nA = completion_sample_size(0.65, 0.5)$nA,
                       nA_r2 = completion_sample_size(0.65, 0.5, r = 2, alpha = 0.05, power = 0.9)$nA)

# ---- MMRM planning
mm <- list()
for (a in list(list(d = 0.4, pA = 0.65, pB = 0.55),
               list(d = 0.35, pA = 0.7, pB = 0.6, structure = "ar1", rho = 0.8, visits = 8, pattern = "early"),
               list(d = 0.3, pA = 0.75, pB = 0.6, structure = "cs", rho = 0.5, visits = 4, pattern = "late", r = 2, alpha = 0.05))) {
  d <- do.call(cs_mmrm_design, a)
  ss <- cs_mmrm_sample_size(d, 0.9)
  mm[[length(mm) + 1]] <- list(args = a, phiA = d$phiA, phiB = d$phiB, nA = ss$nA, power = ss$power)
}
out$mmrm <- mm

# ---- comparisons
cmp <- cs_compare_designs(0.65, 0.5, 6, 18, muB = 15, gamma = "compute")
out$compare <- list(n = cmp$n_per_arm, power = cmp$power)
t6 <- cs_table6()
out$table6 <- list(placebo = t6$placebo_controlled, active = t6$active_comparator)

# ---- power grid
d <- cs_design(0.65, 0.5, delta = 6, sigma = 18, weights = c(0.3, 0.7), contingent = TRUE)
g <- cs_power_grid(d, 100, pA = c(0.5, 0.7), pB = c(0.4, 0.6), delta = c(3, 8))
out$grid <- list(power = g$power)

# ---- analysis of a synthetic trial
set.seed(20261109)
n <- 300
arm <- rep(c("drug", "placebo"), each = n / 2)
site <- sample(c("S1", "S2", "S3", "S4"), n, replace = TRUE)
baseline <- round(rnorm(n, 95, 12))
comp <- rbinom(n, 1, ifelse(arm == "drug", 0.68, 0.55))
imp <- round(ifelse(arm == "drug", 20, 14) + 0.25 * (baseline - 95) + rnorm(n, 0, 17), 1)
imp[comp == 0] <- NA
ex <- data.frame(patient = sprintf("P%03d", 1:n), arm, site, baseline, completed = comp, improvement = imp)
write.csv(ex, file.path(here, "..", "example_trial.csv"), row.names = FALSE, na = "")
an <- list()
spec <- list(
  list(),
  list(weights = "aligned"),
  list(method = "fisher"),
  list(weights = c(0.3, 0.7), contingent = TRUE, c = 0),
  list(method = "fisher", contingent = TRUE),
  list(covariates = "baseline"),
  list(covariates = "baseline", strata = "site", weights = "aligned", s = -5),
  list(strata = "site", psi_type = "clinical", s = -10))
for (sp in spec) {
  a <- do.call(cs_analyze, c(list(data = ex, arm = "arm", active = "drug", completed = "completed",
                                  outcome = "improvement"), sp))
  an[[length(an) + 1]] <- list(spec = sp, Zd = a$components[["Z_d"]], Ze = a$components[["Z_e"]],
    weights = a$test$weights, statistic = a$test$statistic, critical = a$test$critical,
    p_value = a$test$p_value, reject = a$test$reject,
    effect = a$completers_effect[["estimate"]], se = a$completers_effect[["se"]],
    psi = a$psi[["estimate"]], psi_se = a$psi[["se"]],
    theta = a$theta[["estimate"]], theta_se = a$theta[["se"]])
}
out$analysis <- an

writeLines(js(out), file.path(here, "reference.json"))
cat("Wrote reference.json and example_trial.csv\n")

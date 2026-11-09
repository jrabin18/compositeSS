#' compositeSS: design and sample size for estimands that include dropout
#'
#' See the README and ?cs_design, ?theta_design, ?cs_simulate, ?run_calculator.
#' @keywords internal
"_PACKAGE"

#' Analyze a completed trial under the composite strategy
#'
#' Computes, from patient-level data, the completion statistic Z(d), the
#' completers' efficacy statistic Z(e), the weighted (or Fisher) composite
#' test with optional contingent rule, the standardized composite effect size
#' psi with delta-method (and optionally bootstrap) confidence interval, and
#' the expected benefit per randomized patient theta. The standard error of
#' theta is HC2 robust, because the per-patient benefit has a different
#' variance in each arm; the model-based variance undercovers when allocation
#' is unequal.
#'
#' Completion must be defined in the protocol. Following the proposal:
#' patients who stop treatment but continue to be assessed, or who switch or
#' add rescue medication, are non-completers; patients who complete treatment
#' but miss the final assessment are completers with a missing efficacy value
#' (handled within the completers analysis under MAR, i.e. dropped from the
#' efficacy model, and excluded from theta with a message).
#'
#' @param data Data frame with one row per randomized patient.
#' @param arm Name of the treatment-arm column.
#' @param active Value of \code{arm} identifying the active treatment (A).
#' @param completed Name of a logical / 0-1 column: completed the trial.
#' @param outcome Name of the improvement column (positive = better). Used
#'   for completers only.
#' @param covariates Optional character vector of baseline covariates (e.g.
#'   baseline PANSS) for the efficacy model and theta.
#' @param strata Optional name of a stratification column (site, country). If
#'   given, Z(d) is the Cochran-Mantel-Haenszel statistic and the stratum
#'   enters the efficacy model as a factor.
#' @param weights As in \code{\link{cs_design}}. With \code{"aligned"} the
#'   weights are computed from the trial's own data (completion rates,
#'   completers' mean outcomes minus \code{s}, and the standard errors behind
#'   Z(d) and Z(e)), so the weighted Z tests theta. For a 0/1 response outcome
#'   with \code{s = 0} this is the difference in the probability of completing
#'   and responding. Not available with \code{method = "fisher"} or the
#'   contingent rule.
#' @param method,contingent,c,calibration,alpha,gate_ties As in
#'   \code{\link{cs_design}}.
#' @param s Dropout score for theta.
#' @param psi_type \code{"information"} or \code{"clinical"}; see
#'   \code{\link{cs_psi}}.
#' @param conf_level Confidence level for two-sided intervals.
#' @param bootstrap Number of bootstrap resamples (0 = none). Resampling is
#'   within arm.
#' @param seed Optional seed for the bootstrap.
#' @param completion_model How Z(d) is formed: \code{"default"} (pooled
#'   two-proportion test, or Cochran-Mantel-Haenszel if \code{strata} is
#'   given) or \code{"logistic"} (Wald Z for treatment in a logistic
#'   regression of completion on treatment, \code{covariates} and
#'   \code{strata}; Shih and Quan 2001, model 3a).
#' @return An object of class \code{cs_analysis}.
#' @name cs_analyze
#' @usage cs_analyze(data, arm, active, completed, outcome, covariates = NULL, strata = NULL, weights = c(0.5, 0.5), method = c("wz", "fisher"), contingent = FALSE, c = 0, calibration = c("recalibrated", "usual"), alpha = 0.025, s = 0, psi_type = c("information", "clinical"), conf_level = 0.95, bootstrap = 0, seed = NULL, completion_model = c("default", "logistic"), gate_ties = c("stop", "pass"))
NULL

#' Simulation model of Appendix C (informative dropout)
#'
#' Specifies the three-part model of Appendix C of the proposal:
#' \enumerate{
#'   \item \strong{Outcome trajectories} while on assigned treatment:
#'     \eqn{Y_i \sim N(\mu_i, \Sigma_i)} at post-baseline visits 1..k, with a
#'     linear, quadratic (diminishing gains) or flat mean trajectory, and
#'     \eqn{\Sigma_i} either compound symmetry (variance \eqn{\sigma_i^2},
#'     within-patient correlation \eqn{\rho}) or random intercept and slope.
#'   \item \strong{Dropout}: dropout unrelated to outcome with probability
#'     \eqn{\lambda_{ij}} in visit interval j (\eqn{\lambda_{i1} = 0}, so every
#'     patient has a first post-baseline assessment), and/or informative dropout
#'     with \eqn{logit(p_{ijm}) = \alpha_i + \beta_i Y_{ijm}}; mixtures of the two
#'     are allowed. \eqn{\alpha_i} is calibrated to give the planned completion
#'     rate (always calibrated; with beta = 0 it adds outcome-unrelated dropout
#'     on top of lambda).
#'   \item \strong{Outcomes after dropout}, used only to define the treatment
#'     policy effect against which comparator methods are judged. In the active
#'     arm (A), a patient who drops out for reasons unrelated to outcome
#'     (adverse effects) follows the mean incremental trajectory of the
#'     reference arm (B) from the time of dropout; a patient who drops out for
#'     lack of response continues his or her own trend up to dropout, with random
#'     fluctuation. In the reference arm, dropout does not change the outcome
#'     the patient would have had (the part-1 trajectory).
#' }
#'
#' \strong{Sign convention.} As elsewhere in the package, \eqn{Y} is
#' \emph{improvement} (e.g. baseline minus current PANSS), so that positive
#' values are better. Dropout for lack of response therefore corresponds to
#' \eqn{\beta_i < 0} (less improvement, more dropout). If you work with raw
#' change in PANSS (negative = better), reverse the sign of beta.
#'
#' The mean trajectory is calibrated, together with \eqn{\alpha_i}, so that
#' the completers' mean at the final visit (and hence the completers effect)
#' matches the planning values; see \code{\link{cs_calibrate_appc}}.
#'
#' @param k Number of post-baseline visits (usually 4 to 8).
#' @param shape Mean trajectory: \code{"linear"} (\eqn{\mu_j = F j/k}),
#'   \code{"quadratic"} (diminishing gains, \eqn{\mu_j = F(2t - t^2)},
#'   \eqn{t = j/k}), \code{"flat"} (\eqn{\mu_j = F}), or a numeric vector of
#'   length k with last element 1 giving \eqn{\mu_j / F}. F, the final-visit
#'   mean on treatment, is set by calibration.
#' @param cov \code{"cs"} (compound symmetry) or \code{"ris"} (random
#'   intercept and slope).
#' @param sd Per-visit SD for compound symmetry; length 1 or 2 (A, B).
#' @param rho Within-patient correlation for compound symmetry (0.3-0.7).
#' @param sd_int,sd_slope,cor_is,sd_resid Random-intercept SD, random-slope
#'   SD (slope per whole trial, i.e. over t = j/k from 0 to 1), their
#'   correlation, and residual SD, for \code{cov = "ris"}.
#' @param beta Dependence of dropout on the outcome, length 1 or 2 (A, B);
#'   0 = no informative dropout. Negative = dropout for lack of response.
#' @param lambda Probability of outcome-unrelated dropout per visit interval:
#'   a scalar (same for intervals 2..k), a vector of length k (first element
#'   forced to 0), or a list with elements \code{A} and \code{B}.
#' @param depends_on \code{"current"} (default; dropout at visit j depends on
#'   the unobserved \eqn{Y_{ij}}: missing not at random) or \code{"previous"}
#'   (depends on the last observed value: missing at random).
#' @param post_sd SD of the random fluctuation added to post-dropout
#'   outcomes; default is the within-patient residual SD.
#' @return An object of class \code{cs_appc_model}.
#' @examples
#' m <- cs_appc_model(k = 6, shape = "quadratic", sd = 18, rho = 0.5,
#'                    beta = -0.05, lambda = 0.02)
#' @name cs_appc_model
#' @usage cs_appc_model(k = 6, shape = c("linear", "quadratic", "flat"), cov = c("cs", "ris"), sd = 18, rho = 0.5, sd_int = 12, sd_slope = 10, cor_is = 0, sd_resid = 10, beta = 0, lambda = 0, depends_on = c("current", "previous"), post_sd = NULL)
NULL

#' Calibrate the Appendix C model for one arm
#'
#' Finds the dropout intercept \eqn{\alpha_i} and the final-visit mean on
#' treatment F such that, under the model, the completion rate equals
#' \code{p} and the mean improvement among completers at the final visit
#' equals \code{completers_mean}. Uses common random numbers and alternates a
#' root search for \eqn{\alpha_i} with a fixed-point update of F.
#'
#' @param model A \code{cs_appc_model}.
#' @param arm \code{"A"} or \code{"B"}.
#' @param p Target completion rate.
#' @param completers_mean Target mean improvement among completers at the
#'   final visit.
#' @param n_cal Size of the calibration sample.
#' @param seed Seed for the calibration sample.
#' @param tol Tolerance on the completers' mean.
#' @param max_iter Maximum number of fixed-point iterations.
#' @return A list with \code{alpha}, \code{final_mean} (F), the mean
#'   trajectory \code{mu}, and the attained \code{completion},
#'   \code{completers_mean} and \code{completers_sd}.
#' @name cs_calibrate_appc
#' @usage cs_calibrate_appc(model, arm = c("A", "B"), p, completers_mean, n_cal = 1e5, seed = 1, tol = 1e-3, max_iter = 50)
NULL

#' Generate one trial from the calibrated Appendix C model
#'
#' @param model A \code{cs_appc_model}.
#' @param calA,calB Calibrations from \code{\link{cs_calibrate_appc}} for
#'   arms A and B.
#' @param nA,nB Numbers randomized.
#' @param post Logical: also generate outcomes after dropout (part 3).
#' @return A list with \code{patients} (one row per patient: arm, completed,
#'   drop_visit, reason, y_final observed for completers, and y_tp_final, the
#'   treatment-policy final outcome if \code{post = TRUE}) and \code{long}
#'   (observed visit-level data for MMRM or multiple imputation comparators).
#' @name cs_appc_trial
#' @usage cs_appc_trial(model, calA, calB, nA, nB, post = TRUE)
NULL

#' Treatment policy effect implied by a calibrated Appendix C model
#'
#' Large-sample value of the difference in mean final outcome when outcomes
#' after dropout are generated by part 3 of the model; the target against
#' which the comparator methods (MMRM, placebo-based multiple imputation) are
#' judged. The same rules for outcomes after dropout apply in both arms, and
#' both arms are generated from common random numbers, so the effect is zero
#' when the two arms are identical.
#'
#' @param model A \code{cs_appc_model}.
#' @param calA,calB Calibrations from \code{\link{cs_calibrate_appc}} for
#'   arms A and B.
#' @param n Large-sample size per arm.
#' @param seed Seed.
#' @return Numeric scalar.
#' @name cs_appc_policy_effect
#' @usage cs_appc_policy_effect(model, calA, calB, n = 1e5, seed = 99)
NULL

#' Plot a calibrated Appendix C model
#'
#' Two panels: (left) mean improvement by visit on treatment (part 1, solid),
#' among patients still on treatment at each visit (dashed) and under the
#' treatment policy (part 3, dotted); (right) proportion still on treatment by
#' visit.
#'
#' @param model A \code{cs_appc_model}.
#' @param calA,calB Calibrations from \code{\link{cs_calibrate_appc}} for
#'   arms A and B.
#' @param n Patients per arm used to compute the curves.
#' @param seed Seed.
#' @return Invisibly, a data frame of the plotted values.
#' @name cs_appc_plot
#' @usage cs_appc_plot(model, calA, calB, n = 2e4, seed = 5)
NULL

#' Comparator analyses: MMRM and placebo-based multiple imputation
#'
#' Analyses of the visit-level data that do not treat dropout as part of the
#' outcome. They are the comparators of Aim 4 and are judged, in
#' \code{\link{cs_simulate}}, against the treatment policy effect defined by
#' part 3 of the Appendix C model.
#'
#' \describe{
#'   \item{\code{cs_mmrm()}}{Mixed model for repeated measures: visit-by-arm
#'     means, unstructured within-patient covariance, REML, fitted to all
#'     observed visits. The estimate is the difference between arms at the
#'     final visit. Under MAR this targets the effect had patients stayed on
#'     treatment (a hypothetical estimand), not the treatment policy effect.
#'     Engine \code{"nlme"} uses \code{nlme::gls} (always available);
#'     \code{"mmrm"} uses the \pkg{mmrm} package if installed (faster,
#'     Kenward-Roger degrees of freedom).}
#'   \item{\code{cs_pmi()}}{Placebo-based (copy-reference) multiple imputation
#'     by sequential regression (Ratitch, O'Kelly and Tosiello 2013). For each
#'     visit j = 2..k in turn, a regression of Y_j on Y_1..Y_{j-1} is fitted
#'     in the reference arm only; its parameters are drawn from their
#'     posterior and used to impute the missing values at visit j in both
#'     arms, conditional on each patient's own observed (or already imputed)
#'     history. Each completed data set is analyzed by the difference in final
#'     means, and the results are combined by Rubin's rules.}
#' }
#'
#' @param long Visit-level data with columns \code{id}, \code{arm},
#'   \code{visit} (1..k) and \code{y}, as returned in \code{$long} by
#'   \code{\link{cs_appc_trial}}. Dropout must be monotone and visit 1
#'   observed for every patient (as in the Appendix C model).
#' @param k Number of post-baseline visits.
#' @param active Label of the active arm.
#' @param engine For \code{cs_mmrm}: \code{"nlme"} or \code{"mmrm"}.
#' @return A list with \code{estimate}, \code{se}, \code{df}, and
#'   \code{converged}.
#' @name comparators
NULL

#' @rdname comparators
#' @name cs_mmrm
#' @usage cs_mmrm(long, k = max(long$visit), active = "A", engine = c("nlme", "mmrm"))
NULL

#' @rdname comparators
#' @param M Number of imputations.
#' @param reference Label of the reference (placebo or comparator) arm.
#' @name cs_pmi
#' @usage cs_pmi(long, k = max(long$visit), active = "A", reference = "B", M = 20)
NULL

#' Distribution of the two component statistics under the design alternative
#'
#' Returns the means and variance of the completion statistic Z(d) and the
#' completers' efficacy statistic Z(e) for nA patients on A and nB = r * nA on
#' B, together with their covariance gamma.
#'
#' Z(d) = (pA_hat - pB_hat) / sqrt(pbar (1 - pbar) (1/nA + 1/nB)) has mean
#' eta(d) and variance v(d); Z(e) has mean
#' eta(e) = delta / (sigma sqrt(1/(nA pA) + 1/(nB pB))) and unit variance.
#' For non-inferiority (margin_d > 0) Z(d) uses the unpooled standard error,
#' so v(d) = 1.
#'
#' @param design A \code{cs_design} object.
#' @param n Number of patients on A (nA).
#' @return A list with elements \code{eta_d}, \code{v_d}, \code{eta_e},
#'   \code{gamma}, \code{nA}, \code{nB}, \code{se0}, \code{d_std}.
#' @name cs_components
#' @usage cs_components(design, n)
NULL

#' Covariance of the component statistics under the alternative
#'
#' Convenience wrapper returning the delta-method gamma for a design
#' (independent of n for a fixed allocation ratio).
#'
#' @param design A \code{cs_design} object.
#' @return Numeric scalar.
#' @examples
#' cs_gamma(cs_design(pA = 0.70, pB = 0.50, delta = 0.3, sigma = 1))  # about -0.02
#' @name cs_gamma
#' @usage cs_gamma(design)
NULL

#' Critical value of the composite test
#'
#' For the weighted Z the critical value k applies to Z(w). For Fisher's
#' method it is returned on the scale of the product P(d) P(e): H0 is rejected
#' if the product is below k (equivalently if the Fisher combined P value
#' P(d)P(e)[1 - ln(P(d)P(e))] is below alpha when the test is not contingent).
#'
#' For a contingent design with \code{calibration = "recalibrated"}, k(c)
#' solves \eqn{\int_c^\infty \phi(z)\{1 - \Phi[b(z; k)]\} dz = \alpha}, so that
#' the gated test has size exactly alpha under the composite null hypothesis.
#'
#' @param design A \code{cs_design} object.
#' @return A list with \code{k}, \code{scale} ("Z" or "product") and the
#'   attained \code{size}.
#' @name cs_critical
#' @usage cs_critical(design)
NULL

#' Specify a composite-test design
#'
#' Collects the planning values and design choices for the Shih-Quan
#' composite test of completion and completers' efficacy, with optional
#' unequal weights, contingent (gated) testing and non-inferiority margins.
#'
#' @param pA,pB Planned probabilities of completing the trial on the active
#'   treatment (A) and the comparator (B).
#' @param delta Planned difference in mean improvement between completers on A
#'   and completers on B (clinical units, e.g. PANSS points; positive favours A).
#'   Ignored when \code{efficacy = "binary"}.
#' @param sigma Common within-arm standard deviation of improvement among
#'   completers. Ignored when \code{efficacy = "binary"}.
#' @param r Allocation ratio nB / nA (default 1).
#' @param weights Numeric length 2: weights for completion, w(d), and for
#'   completers' efficacy, w(e). Rescaled internally; only the ratio matters.
#'   Or \code{"aligned"} (binary efficacy, weighted Z, no contingent rule or
#'   margins): estimand-aligned weights proportional to
#'   \eqn{(\bar{s}\,se_d, \bar{p}\,se_e)}, with \eqn{\bar{s}} and
#'   \eqn{\bar{p}} the average response among completers and the average
#'   completion over arms. With these weights the combined statistic tracks
#'   the difference in the probability of completing and responding,
#'   \eqn{\Delta\theta = \bar{s}(p_A - p_B) + \bar{p}(s_A - s_B)}, so it tests
#'   that difference rather than the global null of no difference in either
#'   component. Fixed weights, including equal weights, test the global null
#'   only and can favour A when \eqn{\Delta\theta \le 0}. See
#'   \code{\link{cs_aligned_weights}}.
#' @param method \code{"wz"} (weighted Z, default) or \code{"fisher"} (Fisher's
#'   method as used by Shih and Quan; equal weights only).
#' @param contingent Logical. If \code{TRUE}, testing stops (H0 not rejected)
#'   when the completion statistic Z(d) is below \code{c}.
#' @param c Gate threshold for contingent testing (default 0: stop whenever A
#'   is observed to be worse on completion).
#' @param calibration \code{"recalibrated"} (critical value chosen so the
#'   contingent test has size exactly alpha under the composite null) or
#'   \code{"usual"} (unconditional critical value; conservative).
#' @param alpha One-sided significance level (per comparison).
#' @param gamma \code{"compute"} (delta-method covariance of Z(d) and Z(e)
#'   under the alternative), \code{"zero"}, or a numeric value.
#' @param efficacy \code{"continuous"} (default) or \code{"binary"} (response
#'   among completers, e.g. >= 20\% PANSS reduction).
#' @param respA,respB Planned response proportions among completers when
#'   \code{efficacy = "binary"}.
#' @param margin_d,margin_e Non-inferiority margins (positive numbers) for the
#'   completion difference and the completers' efficacy difference. Both 0
#'   (default) gives the superiority form.
#' @param gate_ties How a completion statistic exactly equal to \code{c} is
#'   treated. Z(d) is discrete: when the numbers of completers are equal it is
#'   exactly 0. \code{"stop"} (default) continues testing only if Z(d) > c,
#'   which keeps the size at or below alpha (2.42\% to 2.46\% for weights
#'   0.3/0.7, c = 0, 50 to 500 per arm); \code{"pass"} continues if
#'   Z(d) >= c, which lets ties through and gives a size slightly above alpha
#'   in finite samples (2.59\% at 50 per arm, 2.51\% at 500).
#' @return An object of class \code{cs_design}.
#' @examples
#' d <- cs_design(pA = 0.65, pB = 0.50, delta = 6, sigma = 18)
#' cs_sample_size(d)
#' # Estimand-aligned weights for a responder outcome
#' da <- cs_design(pA = 0.65, pB = 0.55, efficacy = "binary",
#'                 respA = 0.726, respB = 0.58, weights = "aligned")
#' cs_sample_size(da, power = 0.90)
#' @name cs_design
#' @usage cs_design(pA, pB, delta = NULL, sigma = NULL, r = 1, weights = c(0.5, 0.5), method = c("wz", "fisher"), contingent = FALSE, c = 0, calibration = c("recalibrated", "usual"), alpha = 0.025, gamma = "compute", efficacy = c("continuous", "binary"), respA = NULL, respB = NULL, margin_d = 0, margin_e = 0, gate_ties = c("stop", "pass"))
NULL

#' Sample size for the MMRM comparator
#'
#' Sample size for a mixed model for repeated measures (MMRM) of the
#' continuous outcome at the final visit, fitted to all observed visits. The
#' MMRM estimates a different quantity from the composite strategy: the
#' hypothetical effect had every patient stayed on treatment, assuming that
#' dropout is missing at random. It is offered as a reference for the analysis
#' most trials would otherwise use, not as a test of the composite estimand.
#'
#' With monotone dropout, the variance of the final-visit mean in an arm is
#' \eqn{(\sigma^2/n)\sum_j \tau_j / r_j}, where \eqn{r_j} is the share of the arm
#' still observed at visit j and \eqn{\tau_j} the share of the final-visit
#' variance first explained by visit j (Lu, Luo and Chen 2008). With no
#' correlation between visits this is the completers-only variance; with no
#' dropout it is the complete-data variance.
#'
#' @param d Standardized effect at the final visit for the hypothetical
#'   estimand (difference in means divided by the final-visit SD). This is not
#'   in general the completers' effect size used by \code{\link{cs_design}}.
#' @param pA,pB Probabilities of completing the trial on A and B.
#' @param visits Number of post-baseline visits, equally spaced.
#' @param structure Correlation between visits: \code{"newmeds"} (the matrix
#'   for PANSS total change at weeks 1 to 6, pooled within arm and trial over
#'   the 16 six-week placebo-controlled NEWMEDS trials; 6 visits only),
#'   \code{"ar1"} or \code{"cs"} (compound symmetry), with \code{rho}.
#' @param rho Correlation for \code{"ar1"} (between adjacent visits) or
#'   \code{"cs"} (between any two visits).
#' @param pattern When dropout happens: \code{"even"} (the same share at each
#'   visit, close to the NEWMEDS pattern), \code{"early"} or \code{"late"}.
#' @param r Allocation ratio nB / nA.
#' @param alpha One-sided significance level.
#' @param design A \code{cs_mmrm_design} object.
#' @param nA Number of patients on A.
#' @param power Target power.
#' @return \code{cs_mmrm_design} returns the design with the variance
#'   inflation factors \code{phiA} and \code{phiB}; \code{cs_mmrm_power} the
#'   power at \code{nA}; \code{cs_mmrm_sample_size} a list with \code{nA},
#'   \code{nB}, \code{total} and the power achieved.
#' @references Lu K, Luo X, Chen PY. Sample size estimation for repeated
#'   measures analysis in randomized clinical trials with missing data.
#'   Int J Biostat 2008; 4(1): Article 9.
#' @examples
#' m <- cs_mmrm_design(d = 0.40, pA = 0.65, pB = 0.55)
#' cs_mmrm_sample_size(m, power = 0.90)
#' @name cs_mmrm_design
#' @usage cs_mmrm_design(d, pA, pB, visits = 6, structure = c("newmeds", "ar1", "cs"), rho = NULL, pattern = c("even", "early", "late"), r = 1, alpha = 0.025)
NULL

#' @rdname cs_mmrm_design
#' @name cs_mmrm_power
#' @usage cs_mmrm_power(design, nA)
NULL

#' @rdname cs_mmrm_design
#' @name cs_mmrm_sample_size
#' @usage cs_mmrm_sample_size(design, power = 0.90)
NULL

#' Estimand-aligned weights for the weighted Z
#'
#' The composite estimand is the difference in expected benefit per
#' randomized patient, theta, which decomposes exactly as
#' \deqn{\Delta\theta = \bar{m}(p_A - p_B) + \bar{p}(\mu_A - \mu_B),}
#' where \eqn{\bar{m}} is the average over arms of the completers' mean
#' outcome minus the dropout score \code{s}, and \eqn{\bar{p}} the average
#' completion. With a 0/1 response outcome and \code{s = 0}, \eqn{\mu} is the
#' response probability among completers and theta is the probability of
#' completing and responding. Because Z(d) = (pA - pB) / se_d and
#' Z(e) = (muA - muB) / se_e, weights proportional to
#' \eqn{(\bar{m}\,se_d, \bar{p}\,se_e)} make the weighted Z proportional to the
#' estimated \eqn{\Delta\theta}: it then tests the estimand, and is
#' asymptotically equivalent to the direct test of theta. Equal or other
#' fixed weights test only the global null of no difference in either
#' component.
#'
#' In plain terms: keeping one more patient in the trial is worth
#' \eqn{\bar{m}} (for a response outcome, the share of completers who
#' respond), and a gain in the completers' outcome counts only for the
#' \eqn{\bar{p}} share of patients who complete.
#'
#' The weights depend on the data, so at analysis they are computed from the
#' trial (\code{cs_analyze(weights = "aligned")}); at the design stage
#' \code{cs_design(weights = "aligned")} computes them from the planning
#' values. Use this function with summary results, for example to supply
#' \code{weights} or \code{c_sq} to \code{\link{cs_combine}}.
#'
#' @param pA,pB Completion proportions on A and B.
#' @param muA,muB Completers' mean outcome on A and B (for a response
#'   outcome, the response proportions among completers).
#' @param se_d,se_e Standard errors of the completion difference and of the
#'   completers' difference used to form Z(d) and Z(e).
#' @param s Dropout score (default 0).
#' @return A list with \code{weights} (normalised to unit length),
#'   \code{c_sq} (the equivalent Shih-Quan weight on the estimate scale,
#'   \eqn{\bar{p}/\bar{m}}), \code{m_bar} and \code{p_bar}.
#' @examples
#' cs_aligned_weights(pA = 0.65, pB = 0.55, muA = 0.48, muB = 0.60,
#'                    se_d = 0.049, se_e = 0.056)
#' @seealso \code{\link{cs_design}}, \code{\link{cs_analyze}},
#'   \code{\link{cs_combine}}
#' @name cs_aligned_weights
#' @usage cs_aligned_weights(pA, pB, muA, muB, se_d, se_e, s = 0)
NULL

#' Assurance: power averaged over prior distributions for the planning values
#'
#' The probability of success averaged over Beta priors for the completion
#' rates (e.g. from the dropout-prediction formulae) and a normal prior for
#' the completers' effect (e.g. from the NEWMEDS trials). The trial analysis
#' itself remains frequentist and prespecified.
#'
#' @param design A \code{cs_design} or \code{theta_design}.
#' @param n Patients on A.
#' @param prior_pA,prior_pB Beta shape parameters \code{c(a, b)}; NULL keeps
#'   the design value fixed.
#' @param prior_delta Normal \code{c(mean, sd)} for delta (cs_design) or for
#'   muA - muB holding muB fixed (theta_design); NULL keeps it fixed.
#' @param nsim Number of prior draws.
#' @param seed Optional seed.
#' @return A list with \code{assurance}, \code{power_at_planning_values} and
#'   the prior draws.
#' @name cs_assurance
#' @usage cs_assurance(design, n, prior_pA = NULL, prior_pB = NULL, prior_delta = NULL, nsim = 2000, seed = NULL)
NULL

#' Power when the true values differ from those assumed
#'
#' Evaluates the power of a design at a fixed sample size over a grid of true
#' completion rates and completers' effects.
#'
#' @param design A \code{cs_design}.
#' @param n Patients on A (typically the planned size).
#' @param pA,pB,delta Vectors of true values (default: the design values).
#' @return A data frame with one row per combination.
#' @name cs_power_grid
#' @usage cs_power_grid(design, n, pA = design$pA, pB = design$pB, delta = design$delta)
NULL

#' Power curves
#'
#' Plots power against patients per arm (nA) for one or more designs.
#'
#' @param designs A named list of \code{cs_design} and/or \code{theta_design}
#'   objects.
#' @param n Vector of nA values.
#' @param target Horizontal reference line (default 0.80).
#' @param ... Passed to \code{plot}.
#' @return Invisibly, a data frame of the plotted values.
#' @name cs_power_curve
#' @usage cs_power_curve(designs, n = seq(20, 400, by = 5), target = 0.80, ...)
NULL

#' Statement of the estimand in ICH E9(R1) terms
#'
#' @param design A \code{cs_design} or \code{theta_design}.
#' @param population Text describing the population.
#' @param treatment Text describing the treatment conditions compared.
#' @param endpoint Text for the efficacy variable (e.g. "change from baseline
#'   in PANSS total score at week 6").
#' @return A character string (also printed).
#' @name cs_estimand_statement
#' @usage cs_estimand_statement(design, population = "patients as defined by the trial's eligibility criteria", treatment = "the active treatment (A) versus the comparator (B), each as randomized", endpoint = "improvement in the primary efficacy score at the final visit")
NULL

#' Compare the designs of Table 6 for any planning scenario
#'
#' Sample size per arm (nA, 1:1 allocation) for each of the designs compared
#' in Table 6 of the proposal: weighted Z with three weightings, contingent
#' weighted Z (recalibrated and usual critical value), Fisher's method
#' (unconditional and contingent), completion rates only, completers only, and
#' the expected benefit per randomized patient.
#'
#' @param pA,pB Completion probabilities.
#' @param delta,sigma Completers' difference and SD.
#' @param muB Mean improvement among completers on B (muA = muB + delta);
#'   used only for the expected-benefit row.
#' @param s Dropout score for the expected-benefit row.
#' @param alpha One-sided level.
#' @param power Target power.
#' @param gamma Passed to \code{\link{cs_design}}.
#' @param c Gate for the contingent designs.
#' @param equal_contingent If \code{TRUE}, also include the equally weighted
#'   contingent weighted Z (recalibrated and usual critical value), the
#'   design specified in Box B1 of the proposal.
#' @return A data frame with columns \code{design}, \code{n_per_arm} and
#'   \code{power}.
#' @name cs_compare_designs
#' @usage cs_compare_designs(pA, pB, delta, sigma, muB, s = 0, alpha = 0.025, power = 0.80, gamma = "zero", c = 0, equal_contingent = FALSE)
NULL

#' Reproduce Table 6 of the proposal
#'
#' Placebo-controlled (completion 65\% vs 50\%, completers' difference 6
#' PANSS points, SD 18, completers' means 21 vs 15) and active-comparator
#' (72\% vs 64\%, 3 points, SD 18, 21 vs 18) scenarios, one-sided
#' alpha = 0.025, 80\% power, 1:1 allocation, gamma = 0 as in the proposal.
#'
#' @param gamma Passed to \code{\link{cs_design}}.
#' @param equal_contingent Add the rows for the Box B1 design (see
#'   \code{\link{cs_compare_designs}}).
#' @return A data frame.
#' @name cs_table6
#' @usage cs_table6(gamma = "zero", equal_contingent = FALSE)
NULL

#' Launch the graphical interface
#'
#' Opens the Shiny application shipped in \code{inst/app}: sample size and
#' power, comparison of designs, sensitivity and assurance, simulation under
#' the Appendix C model with MMRM and placebo-based multiple imputation
#' comparators, and analysis of a trial data set. Requires \pkg{shiny}.
#'
#' @param ... Passed to \code{shiny::runApp} (e.g. \code{launch.browser}).
#' @name run_calculator
#' @usage run_calculator(...)
NULL

#' @rdname run_calculator
#' @name cs_app
#' @usage cs_app(...)
NULL

#' Power of the composite test
#'
#' Power of the (weighted or Fisher) Shih-Quan composite test, optionally
#' contingent, at nA = \code{n} patients on A and nB = r * n on B. For the
#' non-contingent weighted Z the closed form
#' \eqn{1 - \Phi\{[z_{1-\alpha} - \mu(w)] / \tau(w)\}} is used. Otherwise power
#' is the one-dimensional integral
#' \deqn{\int_c^\infty \phi([z - \eta_d]/\sqrt{v_d})/\sqrt{v_d}\,
#'   \{1 - \Phi([b(z; k) - \eta_{e|z}] / \sqrt{1 - \rho^2})\} dz,}
#' with \eqn{\eta_{e|z} = \eta_e + \gamma (z - \eta_d)/v_d} and
#' \eqn{\rho = \gamma / \sqrt{v_d}}.
#'
#' @param design A \code{cs_design} object.
#' @param n Number of patients on A (may be a vector).
#' @return An object of class \code{cs_power}: a data frame with the power of
#'   the composite test and of each component claim (completion alone and
#'   completers' efficacy alone, each at level alpha).
#' @examples
#' d <- cs_design(0.65, 0.50, delta = 6, sigma = 18)
#' cs_power(d, n = c(80, 102, 120))
#' @name cs_power
#' @usage cs_power(design, n)
NULL

#' Sample size for the composite test
#'
#' Smallest nA (with nB = ceiling(r * nA)) for which the composite test
#' reaches the target power. No inflation for anticipated dropout is needed:
#' dropout is part of the endpoint.
#'
#' @param design A \code{cs_design} object.
#' @param power Target power (default 0.80).
#' @param n_max Upper limit for the search.
#' @return An object of class \code{cs_ss}.
#' @examples
#' cs_sample_size(cs_design(0.65, 0.50, delta = 6, sigma = 18))  # 102 per arm
#' @name cs_sample_size
#' @usage cs_sample_size(design, power = 0.80, n_max = 1e5)
NULL

#' Power-maximizing weights
#'
#' For the weighted Z, power is greatest when the weights are proportional to
#' eta(d) and eta(e). Reported so investigators can see what their clinically
#' chosen weights cost in power; the analysis weights must still be fixed in
#' the protocol. These weights maximize power against the global null of no
#' difference in either component. They are not estimand-aligned: when the
#' composite estimand is the target, use \code{weights = "aligned"} in
#' \code{\link{cs_design}} and \code{\link{cs_analyze}}.
#'
#' @param design A \code{cs_design} object.
#' @param n Number of patients on A at which to evaluate.
#' @return A list with the optimal (normalised) weights and the power under the
#'   design's weights and under the optimal weights.
#' @name cs_opt_weights
#' @usage cs_opt_weights(design, n)
NULL

#' Standardized composite effect size psi
#'
#' Population-level summary of the weighted composite estimand,
#' \deqn{\psi = [w_d h + w_e d \sqrt{\pi_w}] / \sqrt{w_d^2 + w_e^2},}
#' with \eqn{h = (\pi_A - \pi_B)/\sqrt{\bar\pi(1-\bar\pi)}},
#' \eqn{d = \Delta/\sigma}, \eqn{\bar\pi = (\pi_A + r\pi_B)/(1 + r)} and
#' \eqn{\pi_w = (1 + 1/r)/(1/\pi_A + 1/(r\pi_B))}, a weighted harmonic mean
#' of the completion rates (the harmonic mean when r = 1). Then
#' \eqn{\psi\sqrt{n_{eff}}}, with \eqn{n_{eff} = n_A n_B/(n_A + n_B)}, equals
#' the expected weighted Z exactly, and a test of psi = 0 is equivalent to the
#' weighted Z test. (Versions before 0.1.0.9000 used \eqn{\sqrt{\bar\pi}},
#' for which the relation held only approximately.)
#' \code{type = "clinical"} drops the \eqn{\sqrt{\pi_w}} factor,
#' giving purely clinical weights (this corresponds to a weighted Z with
#' completion-adjusted weights; the protocol must say which is intended).
#'
#' @param design A \code{cs_design} with \code{efficacy = "continuous"}.
#' @param type \code{"information"} (default) or \code{"clinical"}.
#' @return A list with \code{psi}, \code{h}, \code{d}, \code{pbar} and
#'   \code{pi_w}.
#' @examples
#' cs_psi(cs_design(0.65, 0.50, delta = 6, sigma = 18))
#' @name cs_psi
#' @usage cs_psi(design, type = c("information", "clinical"))
NULL

#' Simulate one arm with dropout unrelated to outcome
#'
#' Completion is Bernoulli(p) independently of the outcome, and completers'
#' improvement is N(mu, sigma). This is the quick check used in Appendix E.
#' For trajectories and informative dropout use the Appendix C model
#' (\code{\link{cs_appc_model}}, \code{\link{cs_appc_trial}}).
#'
#' @param n Number randomized to the arm.
#' @param p Completion probability.
#' @param mu Mean improvement among completers.
#' @param sigma SD of improvement among completers.
#' @return A list with logical \code{comp} and numeric \code{y} (NA for
#'   dropouts).
#' @name cs_sim_trial
#' @usage cs_sim_trial(n, p, mu, sigma)
NULL

#' Validate a design by simulation
#'
#' Simulates \code{nsim} trials at nA = \code{n}, nB = r * n and reports the
#' empirical power (or type I error when \code{null = TRUE}) of the composite
#' test and its components, the empirical correlation of Z(d) and Z(e), and
#' the coverage of the delta-method confidence interval for psi (and for
#' theta if \code{muB} is given).
#'
#' Under the null, pB is set to pA and the completers' difference to 0.
#'
#' @param design A \code{cs_design} (continuous efficacy, superiority).
#' @param n Patients on A.
#' @param nsim Number of simulated trials (proposal: 10,000 for power,
#'   20,000 under the null).
#' @param dropout \code{"random"} (completion and efficacy independent, as in
#'   Appendix E) or \code{"appendixC"} (longitudinal trajectories with
#'   outcome-unrelated and/or informative dropout, Appendix C).
#' @param model A \code{\link{cs_appc_model}} (required for
#'   \code{dropout = "appendixC"}). Each arm is calibrated so that its
#'   completion rate and completers' mean match the planning values.
#' @param n_cal Calibration sample size for the Appendix C model.
#' @param null Logical: simulate under the composite null.
#' @param muB Mean improvement among completers on B (muA = muB + delta).
#'   Needed for theta, and required for \code{"appendixC"} because informative
#'   dropout depends on the level of the outcome. Default NULL skips theta.
#' @param s Dropout score for theta.
#' @param conf_level Level of the two-sided intervals whose coverage is checked.
#' @param comparators Character vector, any of \code{"mmrm"} and
#'   \code{"pmi"} (placebo-based multiple imputation). Only with
#'   \code{dropout = "appendixC"}. Each is summarised by bias, empirical and
#'   model-based SE, coverage relative to the treatment policy effect (part 3
#'   of Appendix C) and to the hypothetical on-treatment effect, and rejection
#'   rate. See \code{\link{comparators}}.
#' @param M Number of imputations for \code{"pmi"}.
#' @param mmrm_engine \code{"nlme"} or \code{"mmrm"}; see \code{\link{cs_mmrm}}.
#' @param progress Optional function called as \code{progress(i, nsim)} after
#'   each simulated trial (used by the graphical interface).
#' @param seed Optional seed.
#' @return An object of class \code{cs_sim}.
#' @name cs_simulate
#' @usage cs_simulate(design, n, nsim = 10000, dropout = c("random", "appendixC"), model = NULL, n_cal = 1e5, null = FALSE, muB = NULL, s = 0, conf_level = 0.95, comparators = character(0), M = 20, mmrm_engine = c("nlme", "mmrm"), progress = NULL, seed = NULL)
NULL

#' Expected benefit per randomized patient: design
#'
#' Each completer contributes the observed improvement and each dropout a
#' prespecified score \code{s} (e.g. 0, no improvement). The estimand is
#' \deqn{\theta = [\pi_A\mu_A + (1-\pi_A)s] - [\pi_B\mu_B + (1-\pi_B)s].}
#'
#' @param pA,pB Planned completion probabilities.
#' @param muA,muB Planned mean improvement among completers on A and B
#'   (note: means, not just their difference, are needed).
#' @param sigma Within-arm SD of improvement among completers.
#' @param s Score assigned to dropouts (default 0).
#' @param r Allocation ratio nB / nA.
#' @param alpha One-sided level.
#' @param margin Non-inferiority margin in clinical units (default 0 =
#'   superiority). The hypothesis tested is theta > -margin.
#' @return An object of class \code{theta_design}.
#' @name theta_design
#' @usage theta_design(pA, pB, muA, muB, sigma, s = 0, r = 1, alpha = 0.025, margin = 0)
NULL

#' Value of theta for a design
#' @param design A \code{theta_design}.
#' @name theta_value
#' @usage theta_value(design)
NULL

#' Exact decomposition of theta into the two Shih-Quan components
#'
#' \eqn{\theta = (\pi_A - \pi_B)(\bar\mu - s) + \bar\pi(\mu_A - \mu_B)}, where
#' \eqn{\bar\mu} and \eqn{\bar\pi} are averages over the two arms.
#' @param design A \code{theta_design}.
#' @return A list with \code{completion}, \code{completers} and \code{theta}.
#' @name theta_decompose
#' @usage theta_decompose(design)
NULL

#' Power for the expected benefit per randomized patient
#' @param design A \code{theta_design}.
#' @param n Patients on A (vector allowed).
#' @name theta_power
#' @usage theta_power(design, n)
NULL

#' Sample size for the expected benefit per randomized patient
#'
#' Closed form \eqn{n_A = [z_{1-\alpha} + z_{1-\beta}]^2 (V_A + V_B/r) /
#' (\theta + M)^2}, with \eqn{V = \pi\sigma^2 + \pi(1-\pi)(\mu - s)^2}.
#'
#' @param design A \code{theta_design}.
#' @param power Target power.
#' @return A list with \code{nA}, \code{nB}, \code{total}, \code{theta},
#'   \code{VA}, \code{VB}.
#' @examples
#' theta_sample_size(theta_design(0.65, 0.50, muA = 21, muB = 15, sigma = 18))
#' @name theta_sample_size
#' @usage theta_sample_size(design, power = 0.80)
NULL

#' Sample size for a test of completion rates alone
#'
#' Standard (Fleiss) two-proportion formula, one-sided.
#' @param pA,pB Completion probabilities.
#' @param r Allocation ratio.
#' @param alpha One-sided level.
#' @param power Target power.
#' @name completion_sample_size
#' @usage completion_sample_size(pA, pB, r = 1, alpha = 0.025, power = 0.80)
NULL

#' Check the independence of adjusted component statistics (Shih and Quan 2001)
#'
#' Simulates trials under the composite null in which a baseline covariate
#' affects both completion (logistic model) and outcome, with skewed
#' (centred chi-square, 3 df) or normal errors. Z(d) is the Wald Z from a
#' logistic regression of completion on treatment and the covariate, and Z(e)
#' the ANCOVA statistic among completers, as in Shih and Quan (2001,
#' Statistica Sinica 11:53-61), who proved their asymptotic independence.
#' Returns their correlation and the type I error of the weighted Z and the
#' contingent test with the package's critical values.
#'
#' @param nsim Number of simulated trials (about a minute for 20,000).
#' @param n Patients per arm.
#' @param weights Weights for completion and efficacy.
#' @param c Gate for the contingent test.
#' @param alpha One-sided level.
#' @param errors \code{"chisq3"} (skewed) or \code{"normal"}.
#' @param seed Seed.
#' @return A list with the correlation, means and SDs of Z(d) and Z(e), the
#'   simulated sizes, the critical values used and the Monte Carlo SE.
#' @examples
#' \donttest{cs_sim_sq2001(nsim = 2000)}
#' @name cs_sim_sq2001
#' @usage cs_sim_sq2001(nsim = 20000, n = 150, weights = c(0.3, 0.7), c = 0, alpha = 0.025, errors = c("chisq3", "normal"), seed = 1)
NULL

#' Composite test from summary statistics
#'
#' Runs the Shih-Quan composite test when only summary results are
#' available: the estimate and standard error (or the one-sided P value) for
#' the completion comparison and for the completers' comparison, each
#' oriented so that positive values favour the active treatment. Methods are
#' the weighted Z (including Shih and Quan's \eqn{T_G = (\hat\beta + c
#' \hat\lambda) / \sqrt{se_\beta^2 + c^2 se_\lambda^2}} via \code{c_sq}),
#' Fisher's method (Weiler's \eqn{T_Q}, combined P = q(1 - ln q)), and the
#' Simes-Hommel closed procedure (\eqn{T_B}), each followed by the step-down
#' tests of the individual hypotheses. The contingent rule is available for
#' the weighted Z and Fisher's method.
#'
#' @param est_d,se_d Estimate and SE for completion (e.g. difference in
#'   dropout proportions, or log odds ratio of completion).
#' @param est_e,se_e Estimate and SE for the completers' comparison.
#' @param p_d,p_e One-sided P values, used if the estimate and SE are not
#'   given.
#' @param method \code{"wz"}, \code{"fisher"} or \code{"simes"}.
#' @param weights Weights for the weighted Z (ignored if \code{c_sq} is given).
#' @param c_sq Shih and Quan's weight c on the estimate scale. When the
#'   estimates are differences in completion and in completers' response
#'   proportions, \code{c_sq} equal to the \code{c_sq} returned by
#'   \code{\link{cs_aligned_weights}} gives the estimand-aligned test.
#' @param contingent,c,calibration,gate_ties As in \code{\link{cs_design}}.
#' @param alpha One-sided level (Shih and Quan used 0.05).
#' @return An object of class \code{cs_combined}.
#' @references Shih WJ, Quan H. Testing for treatment differences with
#'   dropouts present in clinical trials: a composite approach. Stat Med
#'   1997; 16:1225-1239.
#' @examples
#' # Shih and Quan (1997), Example 3: two active treatments in acute
#' # schizophrenia, BPRS at week 6, adjusted for baseline BPRS.
#' # Positive values favour the arm labelled as active here.
#' cs_combine(est_d = 0.4890, se_d = 0.3750, est_e = 1.8991, se_e = 1.1717,
#'            method = "fisher", alpha = 0.05)            # joint P = 0.032
#' cs_combine(est_d = 0.4890, se_d = 0.3750, est_e = 1.8991, se_e = 1.1717,
#'            c_sq = 0.4, alpha = 0.05)                   # T_G = 2.08, P = 0.019
#' cs_combine(est_d = 0.4890, se_d = 0.3750, est_e = 1.8991, se_e = 1.1717,
#'            method = "simes", alpha = 0.05)             # not significant
#' @seealso \code{\link{cs_aligned_weights}}
#' @name cs_combine
#' @usage cs_combine(est_d = NULL, se_d = NULL, est_e = NULL, se_e = NULL, p_d = NULL, p_e = NULL, method = c("wz", "fisher", "simes"), weights = c(0.5, 0.5), c_sq = NULL, contingent = FALSE, c = 0, calibration = c("recalibrated", "usual"), alpha = 0.025, gate_ties = c("stop", "pass"))
NULL

#' Sample size for a target assurance
#'
#' Smallest nA at which the assurance (power averaged over the priors, as in
#' \code{\link{cs_assurance}}) reaches \code{target}. The same prior draws are
#' used at every n (common random numbers), so the Monte Carlo assurance is
#' monotone in n. Assurance is bounded above by the prior probability that A
#' is truly better; if \code{target} exceeds the assurance at \code{n_max} the
#' function stops and reports that value.
#'
#' @param design A \code{cs_design} or \code{theta_design}.
#' @param target Target assurance.
#' @param prior_pA,prior_pB,prior_delta Priors, as in \code{\link{cs_assurance}}.
#' @param nsim Number of prior draws (reused at every n).
#' @param seed Seed for the prior draws (NULL leaves the RNG alone).
#' @param n_max Upper limit of the search.
#' @return A list with \code{nA}, \code{nB}, \code{total}, the achieved
#'   \code{assurance}, and \code{assurance_upper_bound} (the assurance at
#'   \code{n_max}).
#' @examples
#' d <- cs_design(0.65, 0.50, delta = 6, sigma = 18)
#' cs_assurance_sample_size(d, target = 0.70, prior_pA = cs_beta_prior(0.65, 0.05),
#'   prior_pB = cs_beta_prior(0.50, 0.05), prior_delta = c(6, 2), nsim = 200)
#' @name cs_assurance_sample_size
#' @usage cs_assurance_sample_size(design, target = 0.70, prior_pA = NULL, prior_pB = NULL, prior_delta = NULL, nsim = 1000, seed = 2026, n_max = 2e4)
NULL

#' Beta prior for a completion rate from a mean and SD
#'
#' Method-of-moments Beta parameters, for example from the completion rate
#' predicted by the dropout formulae (Rabinowitz et al. 2009) and its standard
#' error, for use in \code{\link{cs_assurance}}.
#'
#' @param mean Prior mean, strictly between 0 and 1.
#' @param sd Prior SD, smaller than \code{sqrt(mean * (1 - mean))}.
#' @return Named vector \code{c(a, b)}.
#' @examples
#' cs_beta_prior(0.65, 0.05)
#' @name cs_beta_prior
#' @usage cs_beta_prior(mean, sd)
NULL

#' Expected benefit per randomized patient over a range of dropout scores
#'
#' Estimates theta for each dropout score in \code{s}, with robust (HC2)
#' standard errors and confidence intervals, so that the influence of the
#' valuation of dropout can be shown, as the proposal requires. The
#' completion and completers parts of the exact decomposition are computed
#' from the unadjusted estimates; with covariates or strata theta is adjusted
#' and need not equal their sum.
#'
#' @param data,arm,active,completed,outcome,covariates,strata,conf_level As in
#'   \code{\link{cs_analyze}}.
#' @param s Vector of dropout scores.
#' @return A data frame with one row per score.
#' @name cs_theta_sensitivity
#' @usage cs_theta_sensitivity(data, arm, active, completed, outcome, s = c(0, -5, -10), covariates = NULL, strata = NULL, conf_level = 0.95)
NULL

#' Other composite-strategy analyses: worst-rank, win odds, trimmed means
#'
#' Comparators for the composite estimands of this package, each treating
#' dropout as part of the outcome without a model for unobserved outcomes.
#' \code{cs_worst_rank}: Wilcoxon-Mann-Whitney test in which every dropout is
#' tied at the worst rank and completers are ranked by improvement (Lachin
#' 1999, untimed form). \code{cs_win_odds}: the same test for a hierarchy in
#' which earlier dropout is worse than later dropout and completers rank above
#' all dropouts by improvement; the win odds theta/(1 - theta) is reported with
#' a placement-value confidence interval (Brunner et al. 2021).
#' \code{cs_trimmed_mean}: dropouts are assigned the worst value, the same
#' proportion (the larger dropout rate of the two arms) is trimmed from the
#' bottom of each arm, and the difference in trimmed means is tested by
#' permutation (Permutt and Li 2017). All tests are one-sided, positive
#' favouring the active arm.
#'
#' @param arm Vector of arm labels (two groups).
#' @param active Label of the active arm.
#' @param completed Logical or 0/1 completion indicator.
#' @param y Improvement (positive = better); ignored for dropouts.
#' @param drop_visit Visit of dropout (1, 2, ...) for dropouts; ignored for
#'   completers.
#' @param alpha One-sided level.
#' @param conf_level Confidence level for the win odds.
#' @param nperm Number of permutations.
#' @param seed Optional seed; the global random state is restored afterwards.
#' @return A list with the estimate, the one-sided P value and the decision.
#' @examples
#' set.seed(1)
#' arm <- rep(c("A", "B"), each = 60)
#' comp <- rbinom(120, 1, ifelse(arm == "A", 0.65, 0.5)) == 1
#' y <- ifelse(comp, rnorm(120, ifelse(arm == "A", 21, 15), 18), NA)
#' visit <- ifelse(comp, NA, sample(1:5, 120, TRUE))
#' cs_worst_rank(arm, "A", comp, y)
#' cs_win_odds(arm, "A", comp, y, visit)
#' cs_trimmed_mean(arm, "A", comp, y, nperm = 200, seed = 1)
#' @name composite_alternatives
NULL

#' @rdname composite_alternatives
#' @name cs_worst_rank
#' @usage cs_worst_rank(arm, active, completed, y, alpha = 0.025)
NULL

#' @rdname composite_alternatives
#' @name cs_win_odds
#' @usage cs_win_odds(arm, active, completed, y, drop_visit, alpha = 0.025, conf_level = 0.95)
NULL

#' @rdname composite_alternatives
#' @name cs_trimmed_mean
#' @usage cs_trimmed_mean(arm, active, completed, y, nperm = 1000, alpha = 0.025, seed = NULL)
NULL

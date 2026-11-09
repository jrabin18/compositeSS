# compositeSS (development)

* MMRM comparator sample size: `cs_mmrm_design()`, `cs_mmrm_power()` and `cs_mmrm_sample_size()` give the size for an MMRM of the final-visit outcome under dropout at random (hypothetical estimand), from the number of visits, the correlation between visits (the NEWMEDS week 1 to 6 matrix by default, or AR(1) or compound symmetry) and when dropout happens. For d = 0.40 at the final visit and completion 0.65 vs 0.55 it gives 346 patients at 90% power, against 442 for completers only. `inst/validation/mmrm_sample_size_check.R` checks the power by simulation with `cs_mmrm()`.
* Example and test text no longer name medications.

* Estimand-aligned weights for the weighted Z. `cs_design(weights = "aligned")` (binary efficacy) sets the weights from the planning values, and `cs_analyze(weights = "aligned")` computes them from the trial's data, so the combined statistic tests the difference in theta (for a response outcome, the probability of completing and responding) rather than the global null of no difference in either component. New `cs_aligned_weights()` gives the weights and the equivalent Shih-Quan `c_sq` from summary results for `cs_combine()`. For completion 0.65 vs 0.55 and completer response 0.726 vs 0.58, the aligned design needs 213 per arm at 90% power, the same as the direct test of theta. The calculator app has the option in both tabs. `cs_opt_weights()` is now documented as maximizing power against the global null, not as estimand-aligned.

* `cs_appc_policy_effect()` and `cs_appc_plot()`: outcomes after dropout now follow the same rules in both arms (the comparator previously kept its unobserved on-treatment trajectory), and `cs_appc_policy_effect()` uses common random numbers, so the treatment policy effect is exactly zero when the arms are identical. Earlier versions gave nonzero values under the null (for example -3.4 with dropout depending on the current outcome).
* Documented all arguments of the Appendix C functions; removed a duplicate local definition in `cs_critical()`. `R CMD check` is clean apart from the note for the optional `mmrm` package when it is not installed.

# compositeSS 0.1.0.9000

## Contingent test: rule for ties at the gate

* New argument `gate_ties` in `cs_design()`, `cs_analyze()` and `cs_combine()`,
  also used by `cs_simulate()`. The completion statistic is discrete: with
  equal numbers of completers Z(d) = 0 exactly. `"stop"` (new default) continues
  testing only if Z(d) > c; `"pass"` (previous behaviour) continues if
  Z(d) >= c. For weights 0.3/0.7 and c = 0, size with `"pass"` was 2.59%,
  2.54% and 2.51% at 50, 122 and 500 per arm; with `"stop"` it is 2.42% to
  2.46% (2,000,000 null trials). Planning calculations are unaffected.

## Other composite analyses

* `cs_worst_rank()`: Wilcoxon-Mann-Whitney test with all dropouts tied worst.
* `cs_win_odds()`: win odds for a hierarchy (earlier dropout worse than later
  dropout, completers above all dropouts by improvement), with confidence
  interval from placement values.
* `cs_trimmed_mean()`: Permutt-Li trimmed means with a permutation test.

## Bug fix

* `cs_analyze()`, `cs_theta_sensitivity()`: the standard error of theta is now
  HC2 robust. The model-based `lm` variance assumed equal per-patient
  variance in the two arms, which does not hold for the benefit variable
  (V = p sigma^2 + p(1 - p)(mu - s)^2). With 60 vs 240 patients per arm the
  old 95% interval covered 91.5%; with 240 vs 60, 96.6%. Results with equal
  allocation are essentially unchanged.

## Changed

* `cs_psi()` and the psi reported by `cs_analyze()` / `cs_simulate()`: the
  information factor is now sqrt(pi_w), pi_w a weighted harmonic mean of the
  completion rates (harmonic mean for 1:1), instead of sqrt(pbar). With this
  factor psi * sqrt(nA nB / (nA + nB)) equals the weighted Z exactly. The
  change in psi is below 1% for the Table 6 scenarios.

## New

* `cs_assurance_sample_size()`: smallest n for a target assurance, with
  common random numbers across n.
* `cs_beta_prior()`: Beta parameters from a predicted completion rate and
  its SE.
* `cs_theta_sensitivity()`: theta, robust SE and CI over a range of dropout
  scores s.
* `cs_compare_designs()` / `cs_table6()` gain `equal_contingent`, adding the
  equally weighted contingent design of Box B1. Defaults reproduce Table 6
  as published.

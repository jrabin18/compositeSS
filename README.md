# compositeSS

Design, sample size and analysis for randomized trials in which **dropout is part of the outcome**: the composite strategy of the ICH E9(R1) addendum.

- **Web calculator** (runs in the browser; no installation): https://jrabin18.github.io/compositeSS/
- **R package**: this repository

## What the package does

- **Composite test of completion and efficacy among completers.** This is the Shih–Quan test.
  - It combines the two by weighted Z or Fisher.
  - The weights can be fixed or estimand-aligned. Aligned weights make the combined test a test of the composite estimand.
  - It has an optional contingent rule with a recalibrated critical value, and non-inferiority margins.
  - Functions: `cs_design()`, `cs_sample_size()`, `cs_power()`, `cs_analyze()`, `cs_combine()`.
- **Expected benefit per randomized patient (θ) and its direct test.** Functions: `theta_design()`, `theta_sample_size()`.
- **Standardized composite effect size ψ** with a confidence interval. Function: `cs_psi()`.
- **Sample size for an MMRM**, as a comparator for the hypothetical estimand. Functions: `cs_mmrm_design()`, `cs_mmrm_sample_size()`.
- **Assurance, sensitivity grids and simulation.** Simulation can use dropout unrelated to outcome or a longitudinal dropout model with MMRM and placebo-based multiple imputation as comparators. Functions: `cs_assurance()`, `cs_power_grid()`, `cs_simulate()`.
- **Other composite analyses.** Functions: `cs_worst_rank()` (worst-rank), `cs_win_odds()` (win odds), `cs_trimmed_mean()` (trimmed mean).
- **A Shiny interface.** Function: `run_calculator()`.

## Installation

```r
# install.packages("remotes")
remotes::install_github("jrabin18/compositeSS")
```

The package needs R 4.0 or later and the `nlme` package. `shiny` (for the interface) and `mmrm` (as an MMRM engine) are optional.

## Examples

```r
library(compositeSS)

# Completion 0.65 vs 0.50, completers' mean improvement 21 vs 15 (SD 18),
# equal weights, one-sided alpha 0.025, power 0.80: 102 per arm
d <- cs_design(pA = 0.65, pB = 0.50, delta = 6, sigma = 18)
cs_sample_size(d, power = 0.80)

# Response among completers, with estimand-aligned weights
# (tests the difference in the probability of completing and responding)
da <- cs_design(pA = 0.65, pB = 0.55, efficacy = "binary",
                respA = 0.55, respB = 0.40, weights = "aligned")
cs_sample_size(da, power = 0.90)

# Expected benefit per randomized patient, dropouts scored 0
cs_estimand_statement(theta_design(0.65, 0.50, muA = 21, muB = 15, sigma = 18))

# MMRM comparator (hypothetical estimand, missing at random)
cs_mmrm_sample_size(cs_mmrm_design(d = 0.40, pA = 0.65, pB = 0.55), power = 0.90)

# Calculator in the browser
run_calculator()
```

## Web calculator

`docs/` holds a JavaScript version of the calculator, served by GitHub Pages. It ports the package's computations step by step.

- `docs/tests/validate.html` recomputes values produced by this package in the browser and compares them with R (285 checks; sample sizes agree exactly).
- `docs/tests/make_reference.R` regenerates those values.
- See `docs/README.md` for details.

## Validation

- `tests/testthat` holds the unit tests.
- `inst/validation` compares the planning formulas with published values: Lu, Luo and Chen (2008); the longpower worked example of Iddi and Donohue (2022); and Matsouaka and Betensky (2015).

## Citation

Please cite the paper describing the methods (manuscript submitted) and this repository.

## Licence

MIT. See `LICENSE.md`.

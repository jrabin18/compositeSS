# compositeSS calculator (web)

This folder is the browser version of the compositeSS calculator. GitHub Pages serves it from the `docs/` folder of the repository. It is static HTML and JavaScript: everything runs in the visitor's browser, and data files loaded on the Analyze tab are never uploaded.

- `index.html`, `assets/`: the page. `assets/compositess.js` ports the package's computations step by step; `assets/app.js` is the interface.
- `tests/make_reference.R`: computes reference values with the package source in this repository and writes `tests/reference.json`. It also writes `example_trial.csv`, a simulated data set (300 patients, no real data).
- `tests/validate.html`: recomputes the reference values in the browser and compares them with R (285 checks).
  - Sample sizes agree exactly, and probabilities to at least six significant digits.
  - Monte Carlo results (assurance, simulation, bootstrap) use a different random number generator, so they agree with R only up to simulation error.

After changing the package, regenerate the reference values from the repository root with `Rscript docs/tests/make_reference.R`. Then open `docs/tests/validate.html` from a local web server, for example `python -m http.server --directory docs`, and confirm that every check agrees.

The Appendix C longitudinal simulation, the worst-rank, win-odds and trimmed-mean analyses, the logistic completion statistic and assurance-based sample size are only in the R package.

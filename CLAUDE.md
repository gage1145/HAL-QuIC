# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

HAL-QuIC is an R project for modeling RT-QuIC (Real-Time Quaking-Induced Conversion) kinetics — a prion detection assay. It fits a double sigmoidal model to time-series fluorescence data to characterize the growth and decay phases of the reaction.

The core model is:

```
f(t) = S1 / (1 + exp(a1 * (b1 - t))) + S2 / (1 + exp(a2 * (b2 - t)))
```

where S1/S2 are asymptotes, a1/a2 are steepness parameters, and b1/b2 are midpoint inflection time points.

## Commands

### Run the main model pipeline
```r
# In R or RStudio
source("model.R")
```

### Run the kinetics analysis
```r
source("kinetics.R")
```

### Restore R package environment
```r
renv::restore()
```

### Install a new package (then snapshot)
```r
install.packages("pkg")
renv::snapshot()
```

### Install quicR (GitHub-only package) and snapshot
```r
renv::install("gage1145/quicR")
renv::snapshot()
```

## Architecture

### Data flow

1. Raw `.parquet` data is loaded from `data/data.parquet`
2. `normalize()` (`R/normalize.R`) smooths RFU values with a rolling mean, normalizes to the 8th time point, then computes 1st and 2nd derivatives
3. `calculate_metrics()` (from the `quicR` package) computes kinetic metrics (MPR, AUC, time-to-threshold)
4. `estimate_params()` (`R/estimate_params.R`) extracts landmark-based starting values for the NLS fitter (peak, inflection times, decay characteristics)
5. `fit_model()` (`R/fit_model.R`) attempts a double sigmoidal NLS fit using `algorithm = "port"` with bounded parameters. If the double fit fails, it falls back to single sigmoidal, then re-attempts double using the single-fit coefficients as starting values. Each `nls()` call has a 10-second timeout.
6. Results (model coefficients, fit quality metrics, derived parameters `Kapp`, `tlag`, `tlagprime`) are saved to `data/results.parquet`

### Key files

| File | Purpose |
|------|---------|
| `kinetics.R` | Secondary analysis on raw `.xlsx` files from `raw/kinetics/`, parsed via `quicR::get_quic()`. Produces `figures/kinetics/` plots including parameter vs. seed dilution analysis |
| `R/normalize.R` | `normalize()` function |
| `R/fit_model.R` | `fit_model()` function |
| `R/estimate_params.R` | `estimate_params()` function |
| `R/main_theme.R` | Shared `ggplot2` theme (`main_theme`) for publication-quality figures |
| `model.qmd` | Quarto document explaining the model, rendered to `docs/` |

### Key derived parameters

- `Kapp = 1 / b1` — apparent rate constant
- `tlag = b1 - (S1 / (2 * a1))` — lag time (x-intercept of the tangent at inflection)
- `tlagprime = b1 - a1 / 2` — alternative lag time estimate

### Dependencies

R 4.6.0 with `renv` for package management. Key packages: `tidyverse`, `arrow`, `zoo`, `modelr`, `investr`, `quicR` (custom package providing `get_quic()` and `calculate_metrics()`), `ggpubr`, `latex2exp`, `patchwork`.

The `quicR` package is a domain-specific package for QuIC assay data and is not on CRAN. It is installed using `renv::install("gage1145/quicR")`.

### Output

- `figures/` — PNG plots from `model.R`
- `figures/kinetics/` — PNG plots from `kinetics.R`
- `data/results.parquet` — fitted model results
- `docs/` — rendered Quarto site (`_quarto.yml` sets `output-dir: docs`)

### Website publishing

The Quarto website is **not rendered locally**. It is built and deployed via GitHub Actions (`.github/workflows/publish.yml`) on every push to the `website` branch.

The workflow:
1. Installs R (version must match `renv.lock`, currently 4.6.1) and Quarto
2. Restores R packages via `renv::restore()`
3. Runs `quarto publish` targeting the `gh-pages` branch, which triggers `pre-render: model.R` from `_quarto.yml` before rendering the `.qmd` files

Key constraints:
- The R version in the workflow (`r-version`) must match the version in `renv.lock` — mismatches cause package install failures (e.g., MASS requires R >= 4.4.0)
- `quicR` is a GitHub-only package (`gage1145/quicR`) and must be present in `renv.lock`; after installing it locally, run `renv::snapshot()` and commit the updated lockfile
- `data/data.parquet` must be committed to the repo since `model.R` reads it at render time

# cdmeta: Confidence-Distribution-Based Inference for Random-Effects Meta-Analysis

[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![GitHub development version](https://img.shields.io/badge/devel-1.1--1-purple.svg)](https://github.com/guido-s/cdmeta)

<!-- Activate these badges after the first CRAN release.
[![CRAN status](https://www.r-pkg.org/badges/version/cdmeta)](https://CRAN.R-project.org/package=cdmeta)
[![Monthly downloads](https://cranlogs.r-pkg.org/badges/cdmeta)](https://cran.r-project.org/package=cdmeta)
[![Total downloads](https://cranlogs.r-pkg.org/badges/grand-total/cdmeta)](https://cran.r-project.org/package=cdmeta)
-->

Development repository for the R package **cdmeta**.

## Overview

**cdmeta** implements confidence-distribution propagation for frequentist inference in random-effects meta-analysis. Its main purpose is to construct prediction intervals for the true effect in a future study while propagating uncertainty through the random-effects hierarchy.

For each Monte Carlo replication, the method:

1. draws the between-study variance, $\tau^2$, from a confidence distribution based on the exact distribution of Cochran's $Q$;
2. draws the average effect, $\mu$, from its conditional confidence distribution given the sampled $\tau^2$; and
3. generates the true effect in a future study, $\theta_{\mathrm{new}}$, from the random-effects model.

Empirical quantiles of the resulting Monte Carlo samples provide:

- a prediction interval for $\theta_{\mathrm{new}}$;
- confidence intervals for $\mu$, $\tau^2$, $\tau$, and $I^2$; and
- Monte Carlo point summaries and draws for further analysis.

The package also provides forest plots, predictive-distribution plots, optional effect-scale transformations, and example datasets for ratio measures analyzed on the logarithmic scale.

## Status

Version **1.1-1** is the current development release. The package is being prepared for submission to CRAN.

## Installation

Install the development version directly from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("guido-s/cdmeta")
```

Then load the package:

```r
library(cdmeta)
```

After the first CRAN release, the stable version will be installable with:

```r
install.packages("cdmeta")
```

## Quick start: intravenous iron for heart failure

The `hf_iron` dataset contains study-specific log risk ratios and standard errors from six randomized trials of intravenous iron therapy for patients with heart failure and iron deficiency.

```r
library(cdmeta)

data("hf_iron")

fit_hf <- cdmeta(
  y = hf_iron$yi,
  se = hf_iron$sei,
  B = 20000,
  seed = 3333,
  transf = exp,
  transf_name = "exp"
)

fit_hf
```

Draw a forest plot on the risk-ratio scale:

```r
forest(
  fit_hf,
  slab = hf_iron$study,
  at = log(c(0.25, 0.5, 1, 2, 4)),
  xlab = "Risk ratio",
  mark_summary_estimate = TRUE,
  mark_prediction_estimate = TRUE
)
```

Plot the Monte Carlo predictive distribution of the future true effect:

```r
plot(fit_hf)
```

For final analyses, a larger Monte Carlo sample may be desirable; the default is `B = 25000`.

## Second example: vitamin K antagonists for atrial fibrillation

The `af_vka` dataset contains direct comparisons of vitamin K antagonists versus placebo or control for prevention of stroke in patients with non-valvular atrial fibrillation. Effect estimates are log odds ratios.

```r
data("af_vka")

fit_af <- cdmeta(
  y = af_vka$yi,
  se = af_vka$sei,
  B = 20000,
  seed = 3333,
  transf = exp,
  transf_name = "exp"
)

fit_af

forest(
  fit_af,
  slab = af_vka$study,
  order = "weight",
  at = log(c(0.125, 0.25, 0.5, 1, 2)),
  xlab = "Odds ratio",
  main = "Vitamin K antagonists for atrial fibrillation",
  mark_summary_estimate = TRUE,
  mark_prediction_estimate = TRUE
)
```

## Methodological defaults

The defaults reproduce the confidence-distribution propagation method evaluated in the accompanying methodological article:

```r
i2_method = "typical_se2"
mu_dist = "normal"
```

The alternative settings

```r
i2_method = "mean_se2"
i2_method = "harmonic_mean_se2"
mu_dist = "t"
```

are available for sensitivity analyses. In particular, `mu_dist = "t"` is an optional extension and is not the primary method evaluated in the article.

## Supplying external heterogeneity draws

Advanced users may supply externally generated non-negative draws of $\tau^2$ through `tau2_samples`. In that case, `pimeta::pima()` is not called, and the number of supplied draws determines the Monte Carlo sample size.

```r
# tau2_draws should contain non-negative draws obtained from
# an appropriate external confidence-distribution procedure.

fit_external <- cdmeta(
  y = hf_iron$yi,
  se = hf_iron$sei,
  tau2_samples = tau2_draws,
  seed = 3333,
  i2_method = "typical_se2",
  mu_dist = "normal"
)
```

## Main functions

- `cdmeta()` fits the confidence-distribution propagation model and returns Monte Carlo summaries, intervals, and draws.
- `print.cdmeta()` prints point estimates, confidence intervals, and the prediction interval.
- `plot.cdmeta()` displays the predictive distribution of the future true effect.
- `forest_cdmeta()` produces a forest plot with the pooled effect, prediction interval, and heterogeneity summaries.

## Reproducibility

The accompanying methodological article is:

> Noma H, Schwarzer G. *Frequentist prediction intervals for random-effects meta-analysis via confidence-distribution propagation*.

The arXiv DOI and the permanent repository for the manuscript's simulation and application code will be added after public release.

## Citation

After installation, citation information can be obtained with:

```r
citation("cdmeta")
```

## Authors

- Hisashi Noma
- Guido Schwarzer

## Bug reports

Please report bugs and feature requests through the GitHub issue tracker:

<https://github.com/guido-s/cdmeta/issues>

Alternatively, contact the package maintainer, Hisashi Noma, at <noma@ism.ac.jp>.

## License

**cdmeta** is distributed under the GNU General Public License version 3 (GPL-3).

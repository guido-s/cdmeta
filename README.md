# cdmeta: Confidence-Distribution-Based Inference for Random-Effects Meta-Analysis

Official Git repository of R package **cdmeta**

[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![CRAN Version](https://www.r-pkg.org/badges/version/cdmeta)](https://cran.r-project.org/package=cdmeta)
[![GitHub develop](https://img.shields.io/badge/develop-0.1--0-purple)](https://img.shields.io/badge/develop-0.1--0-purple)
[![Monthly Downloads](https://cranlogs.r-pkg.org/badges/cdmeta)](https://cranlogs.r-pkg.org/badges/cdmeta)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/cdmeta)](https://cranlogs.r-pkg.org/badges/grand-total/cdmeta)

## Authors

Hisashi Noma,
Guido Schwarzer

## Description

**cdmeta** provides computational tools for confidence-distribution-propagation-based inference in random-effects meta-analysis. The package implements frequentist prediction intervals for the effect in a future study by propagating uncertainty in the between-study variance through the confidence distribution.

The main function, `cdmeta()`, provides interval estimation for the overall mean effect, between-study variance, between-study standard deviation, heterogeneity proportion, and the effect in a future study. The package also includes a forest plot method for `cdmeta` objects and example datasets for applications with log ratio measures.

## Installation

<!--
### Current stable [![CRAN Version](https://www.r-pkg.org/badges/version/cdmeta)](https://cran.r-project.org/package=cdmeta) release:
```r
install.packages("cdmeta")
```
-->

### Current [![GitHub develop](https://img.shields.io/badge/develop-0.1--0-purple)](https://img.shields.io/badge/develop-0.1--0-purple) release on GitHub:

Installation using R package [**remotes**](https://cran.r-project.org/package=remotes):

```r
install.packages("remotes")
remotes::install_github("guido-s/cdmeta")
```

## Example

```r
library(cdmeta)

data(hf_iron)

fit_hf <- cdmeta(
  y = hf_iron$yi,
  se = hf_iron$sei,
  B = 10000,
  seed = 11111,
  transf = exp,
  transf_name = "exp")

fit_hf

forest(
  fit_hf,
  slab = hf_iron$study,
  at = log(c(0.25, 0.5, 1, 2, 4)),
  xlab = "Odds ratio",
  mark_summary_estimate = TRUE,
  mark_prediction_estimate = TRUE)
```

## Bug Reports

```r
bug.report(package = "cdmeta")
```

The `bug.report()` function is not supported in RStudio. Please send an email to Hisashi Noma <noma@ism.ac.jp> if you use RStudio.

#' The 'cdmeta' package
#'
#' Computational tools for confidence-distribution-propagation-based inference
#' in random-effects meta-analysis. The package implements frequentist
#' prediction intervals for the effect in a future study by propagating
#' uncertainty in the between-study variance through the confidence distribution.
#' It also provides interval estimation for the overall mean effect,
#' between-study variance, between-study standard deviation, and heterogeneity
#' proportion. The implemented methods are described in Noma and Schwarzer
#' (2026).
#'
#' The main function is [cdmeta()], which performs unified confidence-
#' distribution-based inference for random-effects meta-analysis. Given sampled
#' values of the between-study variance \eqn{\tau^2}, the function conducts
#' conditional sampling of the overall mean effect \eqn{\mu} and predictive
#' sampling of the effect in a future study \eqn{\theta_{\mathrm{new}}}. The
#' function also provides inference for heterogeneity measures, including
#' \eqn{\tau^2}, \eqn{\tau}, and \eqn{I^2}.
#'
#' @references
#' Higgins, J. P. T., Thompson, S. G., and Spiegelhalter, D. J. (2009).
#' A re-evaluation of random-effects meta-analysis.
#' \emph{Journal of the Royal Statistical Society: Series A.}
#' \strong{172}(1): 137-159.
#' \doi{10.1111/j.1467-985X.2008.00552.x}
#'
#' Noma, H., and Schwarzer, G. (2026).
#' Frequentist prediction intervals for random-effects meta-analysis via
#' confidence-distribution propagation.
#' \emph{arXiv.}
#' \href{https://arxiv.org/abs/XXXX.XXXXX}{arXiv:XXXX.XXXXX}
#'
#' Partlett, C., and Riley, R. D. (2017).
#' Random effects meta-analysis: Coverage performance of 95 percent confidence
#' and prediction intervals following REML estimation.
#' \emph{Statistics in Medicine.}
#' \strong{36}(2): 301-317.
#' \doi{10.1002/sim.7140}
#'
#' Viechtbauer, W. (2010).
#' Conducting meta-analyses in R with the metafor package.
#' \emph{Journal of Statistical Software.}
#' \strong{36}(3): 1-48.
#' \doi{10.18637/jss.v036.i03}
#'
#' @seealso [cdmeta()], [pimeta::pima()], [metafor::rma()], [meta::metagen()]
#' 
#' @importFrom metafor forest
#' @export forest

"_PACKAGE"

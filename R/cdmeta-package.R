#' cdmeta: Confidence-Distribution-Based Inference for Random-Effects Meta-Analysis
#'
#' @description
#' Computational tools for confidence-distribution-propagation-based inference
#' in random-effects meta-analysis. The package provides prediction intervals
#' for the effect in a future study and interval estimation for the overall mean
#' effect and heterogeneity measures.
#'
#' @details
#' The main function is \code{\link{cdmeta}}. Forest plots are produced by the
#' exported S3 generic \code{\link{forest}}, which dispatches to
#' \code{forest.cdmeta()} for objects of class \code{"cdmeta"}.
#'
#' @references
#' Higgins, J. P. T., Thompson, S. G., and Spiegelhalter, D. J. (2009).
#' A re-evaluation of random-effects meta-analysis. \emph{Journal of the Royal
#' Statistical Society: Series A}, \strong{172}(1), 137--159.
#' \doi{10.1111/j.1467-985X.2008.00552.x}
#'
#' Noma, H., and Schwarzer, G. (2026). Frequentist prediction intervals for
#' random-effects meta-analysis via confidence-distribution propagation.
#'
#' Partlett, C., and Riley, R. D. (2017). Random effects meta-analysis: Coverage
#' performance of 95 percent confidence and prediction intervals following REML
#' estimation. \emph{Statistics in Medicine}, \strong{36}(2), 301--317.
#' \doi{10.1002/sim.7140}
#'
#' Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor
#' package. \emph{Journal of Statistical Software}, \strong{36}(3), 1--48.
#' \doi{10.18637/jss.v036.i03}
#'
#' @seealso \code{\link{cdmeta}}, \code{\link{forest}}
#' @keywords internal
#' @aliases cdmeta-package NULL
"_PACKAGE"

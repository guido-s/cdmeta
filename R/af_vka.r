#' Vitamin K antagonists for prevention of stroke in non-valvular atrial fibrillation
#'
#' A pairwise subset extracted from a network meta-analysis of antithrombotic
#' treatments for prevention of stroke in patients with non-valvular atrial
#' fibrillation. The dataset includes direct comparisons of vitamin K
#' antagonists versus placebo or control, with study-specific log odds ratios and
#' standard errors for the stroke outcome.
#'
#' @format A data frame with 6 rows and 11 variables:
#' \describe{
#' \item{study_id}{Study identifier.}
#' \item{study}{Trial name.}
#' \item{treatment}{Treatment group label.}
#' \item{control}{Control group label.}
#' \item{treatment_events}{Number of stroke events in the treatment group.}
#' \item{treatment_total}{Number of participants in the treatment group.}
#' \item{control_events}{Number of stroke events in the control group.}
#' \item{control_total}{Number of participants in the control group.}
#' \item{yi}{Study-specific log odds ratio comparing vitamin K antagonists with placebo or control.}
#' \item{sei}{Standard error of `yi`.}
#' \item{vi}{Sampling variance of `yi`, equal to `sei^2`.}
#' }
#'
#' @references
#' Dogliotti, A., Paolasso, E., and Giugliano, R. P. (2014).
#' Current and new oral antithrombotics in non-valvular atrial fibrillation: a
#' network meta-analysis of 79808 patients.
#' \emph{Heart.}
#' \strong{100}(5): 396--405.
#' \doi{10.1136/heartjnl-2013-304347}
#'
#' Viechtbauer, W., White, T., Noble, D., Senior, A., Hamilton, W., Schwarzer,
#' G., and Rover, C. (2026).
#' \emph{metadat: Meta-Analysis Datasets.}
#' R package.
#'
#' @examples
#' data(af_vka)
#'
#' fit_af <- cdmeta(
#'   y = af_vka$yi,
#'   se = af_vka$sei,
#'   B = 10000,
#'   seed = 11111,
#'   i2_method = "harmonic_mean_se2",
#'   mu_dist = "t",
#'   df = nrow(af_vka) - 1,
#'   transf = exp,
#'   transf_name = "exp"
#' )
#'
#' fit_af
#'
#' forest(
#'   fit_af,
#'   slab = af_vka$study,
#'   order = "weight",
#'   summary_stat = "median",
#'   at = log(c(0.125, 0.25, 0.5, 1, 2)),
#'   xlab = "Odds ratio",
#'   main = "Vitamin K antagonists for atrial fibrillation",
#'   mark_summary_estimate = TRUE,
#'   mark_prediction_estimate = TRUE
#' )
#'
#' @keywords datasets
"af_vka"

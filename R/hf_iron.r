#' Intravenous iron therapy for patients with heart failure and iron deficiency
#'
#' A dataset from a systematic review and meta-analysis of randomized trials
#' evaluating intravenous iron therapy for patients with heart failure and iron
#' deficiency. The outcome is represented by study-specific log risk ratios and
#' their standard errors.
#'
#' @format A data frame with 6 rows and 5 variables:
#' \describe{
#' \item{study_id}{Study identifier.}
#' \item{study}{Trial name.}
#' \item{yi}{Study-specific log risk ratio.}
#' \item{sei}{Standard error of `yi`.}
#' \item{vi}{Sampling variance of `yi`, equal to `sei^2`.}
#' }
#'
#' @references
#' Anker, S. D., Karakas, M., Mentz, R. J., Ponikowski, P., Butler, J., Khan,
#' M. S., Talha, K. M., Kalra, P. R., Hernandez, A. F., Mulder, H., Rockhold,
#' F. W., Placzek, M., Rover, C., Cleland, J. G. F., and Friede, T. (2025).
#' Systematic review and meta-analysis of intravenous iron therapy for patients
#' with heart failure and iron deficiency.
#' \emph{Nature Medicine.}
#' \strong{31}: 2640--2646.
#' \doi{10.1038/s41591-025-03671-1}
#'
#' Viechtbauer, W., White, T., Noble, D., Senior, A., Hamilton, W., Schwarzer,
#' G., and Rover, C. (2026).
#' \emph{metadat: Meta-Analysis Datasets.}
#' R package.
#'
#' @examples
#' data(hf_iron)
#'
#' fit_hf <- cdmeta(
#'   y = hf_iron$yi,
#'   se = hf_iron$sei,
#'   B = 10000,
#'   seed = 11111,
#'   transf = exp,
#'   transf_name = "exp"
#' )
#'
#' fit_hf
#'
#' forest(
#'   fit_hf,
#'   slab = hf_iron$study,
#'   at = log(c(0.25, 0.5, 1, 2, 4)),
#'   xlab = "Odds ratio",
#'   mark_summary_estimate = TRUE,
#'   mark_prediction_estimate = TRUE
#' )
#'
#' @keywords datasets
"hf_iron"

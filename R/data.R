#' Vitamin K antagonists for prevention of stroke in non-valvular atrial fibrillation
#'
#' @description
#' A pairwise subset extracted from a network meta-analysis of antithrombotic
#' treatments for prevention of stroke in patients with non-valvular atrial
#' fibrillation. The dataset includes direct comparisons of vitamin K antagonists
#' versus placebo or control, with study-specific log odds ratios and standard
#' errors for the stroke outcome.
#'
#' @format A data frame with 6 rows and 11 variables:
#' \describe{
#'   \item{study_id}{Study identifier.}
#'   \item{study}{Trial name.}
#'   \item{treatment}{Label for the treatment group.}
#'   \item{control}{Label for the control group.}
#'   \item{treatment_events}{Number of stroke events in the treatment group.}
#'   \item{treatment_total}{Number of participants in the treatment group.}
#'   \item{control_events}{Number of stroke events in the control group.}
#'   \item{control_total}{Number of participants in the control group.}
#'   \item{yi}{Study-specific log odds ratio comparing vitamin K antagonists with placebo or control.}
#'   \item{sei}{Standard error of \code{yi}.}
#'   \item{vi}{Sampling variance of \code{yi}, equal to \code{sei^2}.}
#' }
#'
#' @source
#' Dogliotti, A., Paolasso, E., and Giugliano, R. P. (2014). Current and new oral
#' antithrombotics in non-valvular atrial fibrillation: a network meta-analysis
#' of 79808 patients. \emph{Heart}, \strong{100}(5), 396--405.
#' \doi{10.1136/heartjnl-2013-304347}
#'
#' @docType data
#' @usage data(af_vka)
#' @keywords datasets
"af_vka"


#' Intravenous iron therapy for patients with heart failure and iron deficiency
#'
#' @description
#' A dataset from a systematic review and meta-analysis of randomized trials
#' evaluating intravenous iron therapy for patients with heart failure and iron
#' deficiency. The outcome is represented by study-specific log risk ratios and
#' their standard errors.
#'
#' @format A data frame with 6 rows and 5 variables:
#' \describe{
#'   \item{study_id}{Study identifier.}
#'   \item{study}{Trial name.}
#'   \item{yi}{Study-specific log risk ratio.}
#'   \item{sei}{Standard error of \code{yi}.}
#'   \item{vi}{Sampling variance of \code{yi}, equal to \code{sei^2}.}
#' }
#'
#' @source
#' Anker, S. D., et al. (2025). Systematic review and meta-analysis of
#' intravenous iron therapy for patients with heart failure and iron deficiency.
#' \emph{Nature Medicine}, \strong{31}, 2640--2646.
#' \doi{10.1038/s41591-025-03671-1}
#'
#' @docType data
#' @usage data(hf_iron)
#' @keywords datasets
"hf_iron"

\docType{data}
\name{hf_iron}
\alias{hf_iron}
\title{Intravenous iron therapy for patients with heart failure and iron deficiency}
\format{A data frame with 6 rows and 5 variables}
\usage{
data(hf_iron)
}
\description{
A dataset from a systematic review and meta-analysis of randomized trials evaluating intravenous iron therapy for patients with heart failure and iron deficiency. The outcome is represented by study-specific log risk ratios and their standard errors.
\itemize{
\item \code{study}: Trial name
\item \code{yi}: Study-specific log risk ratio
\item \code{sei}: Standard error of \code{yi}
\item \code{vi}: Sampling variance of \code{yi}, equal to \code{sei^2}
}
}
\references{
Anker, S. D., Karakas, M., Mentz, R. J., Ponikowski, P., Butler, J., Khan, M. S., Talha, K. M., Kalra, P. R., Hernandez, A. F., Mulder, H., Rockhold, F. W., Placzek, M., Rover, C., Cleland, J. G. F., and Friede, T. (2025).
Systematic review and meta-analysis of intravenous iron therapy for patients with heart failure and iron deficiency.
\emph{Nature Medicine.}
\strong{31}: 2640--2646.
\doi{10.1038/s41591-025-03671-1}

Viechtbauer, W., White, T., Noble, D., Senior, A., Hamilton, W., Schwarzer, G., and Rover, C. (2026).
\emph{metadat: Meta-Analysis Datasets.}
R package.
}
\keyword{datasets}

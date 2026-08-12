\docType{data}
\name{af_vka}
\alias{af_vka}
\title{Vitamin K antagonists for prevention of stroke in non-valvular atrial fibrillation}
\format{A data frame with 6 rows and 11 variables}
\usage{
data(af_vka)
}
\description{
A pairwise subset extracted from a network meta-analysis of antithrombotic treatments for prevention of stroke in patients with non-valvular atrial fibrillation. The dataset includes direct comparisons of vitamin K antagonists versus placebo or control, with study-specific log odds ratios and standard errors for the stroke outcome.
\itemize{
\item \code{study}: Trial name
\item \code{events_vka}: Number of stroke events in the vitamin K antagonist group
\item \code{n_vka}: Number of participants in the vitamin K antagonist group
\item \code{events_control}: Number of stroke events in the placebo or control group
\item \code{n_control}: Number of participants in the placebo or control group
\item \code{yi}: Study-specific log odds ratio comparing vitamin K antagonists with placebo or control
\item \code{sei}: Standard error of \code{yi}
\item \code{vi}: Sampling variance of \code{yi}, equal to \code{sei^2}
}
}
\references{
Dogliotti, A., Paolasso, E., and Giugliano, R. P. (2014).
Current and new oral antithrombotics in non-valvular atrial fibrillation: a network meta-analysis of 79808 patients.
\emph{Heart.}
\strong{100}(5): 396--405.
\doi{10.1136/heartjnl-2013-304347}

Viechtbauer, W., White, T., Noble, D., Senior, A., Hamilton, W., Schwarzer, G., and Rover, C. (2026).
\emph{metadat: Meta-Analysis Datasets.}
R package.
}
\keyword{datasets}

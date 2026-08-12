\name{forest_cdmeta}
\alias{forest_cdmeta}
\title{Forest plot for cdmeta objects}
\usage{
forest_cdmeta(
  x,
  slab = NULL,
  order = c("none", "increasing", "decreasing", "precision", "weight"),
  level = NULL,
  summary_stat = c("mean", "median"),
  refline = 0,
  xlab = "Effect size",
  main = NULL,
  atransf = NULL,
  at = NULL,
  alim = NULL,
  xlim = NULL,
  digits = 2,
  ci_digits = digits,
  weight_digits = 1,
  show_weights = TRUE,
  show_pi = TRUE,
  show_het = TRUE,
  header = TRUE,
  annotate = TRUE,
  cex = 0.85,
  psize = NULL,
  pch = 22,
  box_col = "black",
  box_bg = "white",
  ci_col = "black",
  summary_col = "black",
  summary_bg = "gray20",
  pi_col = "gray40",
  refline_col = "gray70",
  grid = TRUE,
  grid_col = "gray90",
  qtype = 8,
  symmetric_shapes = TRUE,
  mark_summary_estimate = FALSE,
  mark_prediction_estimate = FALSE,
  estimate_mark_col = "black",
  estimate_mark_lwd = 1,
  mar = c(4.5, 1.0, 3.0, 1.0),
  ...
)
}
\arguments{
\item{x}{An object of class \code{"cdmeta"}.}

\item{slab}{An optional character vector of study labels. If \code{NULL},
the names of \code{x$data$y} are used when available; otherwise generic study
labels are created. Use this argument to print study names inside the forest
plot.}

\item{order}{A character string specifying the order of studies in the plot.
The available options are \code{"none"}, \code{"increasing"}, 
\code{"decreasing"}, \code{"precision"}, and \code{"weight"}.}

\item{level}{The confidence level used for study-specific confidence intervals
and confidence-distribution-based summary intervals. If \code{NULL}, 
\code{1 - x$alpha} is used.}

\item{summary_stat}{A character string specifying whether Monte Carlo means
(\code{"mean"}) or medians (\code{"median"}) are used as point estimates for
the summary rows.}

\item{refline}{A numeric value indicating the location of the reference line.
For log odds ratios, log risk ratios, or log hazard ratios, this is usually 0.}

\item{xlab}{The label for the x-axis.}

\item{main}{An optional main title.}

\item{atransf}{An optional function used to transform axis labels and printed
estimates. For example, \code{atransf = exp} can be used when the analysis
scale is the log odds ratio scale. The plotting coordinates themselves remain
on the original analysis scale. If \code{atransf = NULL} and the \code{cdmeta}
object contains a transformation supplied through \code{cdmeta(transf = ...)},
that transformation is used automatically.}

\item{at}{An optional numeric vector specifying the locations of tick marks on
the original analysis scale. For example, when \code{atransf = exp}, tick marks
for odds ratios may be specified as \code{log(c(0.25, 0.5, 1, 2, 4))}.}

\item{alim}{An optional numeric vector of length 2 specifying the plotting
limits for the effect-size axis on the original analysis scale.}

\item{xlim}{An optional numeric vector of length 2 specifying the full plotting
limits including the left and right text columns.}

\item{digits}{The number of digits used for axis labels.}

\item{ci_digits}{The number of digits used for effect estimates and interval
limits.}

\item{weight_digits}{The number of digits used for study weights.}

\item{show_weights}{A logical value indicating whether study weights are shown.}

\item{show_pi}{A logical value indicating whether the prediction interval for
a future study effect is shown.}

\item{show_het}{A logical value indicating whether heterogeneity statistics are
shown.}

\item{header}{A logical value indicating whether column headers are shown.}

\item{annotate}{A logical value indicating whether estimates and intervals are
printed on the right side of the plot.}

\item{cex}{Character expansion factor.}

\item{psize}{Optional point sizes for study-specific estimates. If \code{NULL},
point sizes are proportional to inverse-variance weights based on
\eqn{1 / (se_i^2 + \hat{\tau}^2)}.}

\item{pch}{Plotting symbol for study-specific estimates.}

\item{box_col}{Color of the border of study-specific points.}

\item{box_bg}{Fill color of study-specific points.}

\item{ci_col}{Color of study-specific confidence intervals.}

\item{summary_col}{Border color of the summary diamond.}

\item{summary_bg}{Fill color of the summary diamond.}

\item{pi_col}{Color of the prediction interval.}

\item{refline_col}{Color of the reference line.}

\item{grid}{A logical value indicating whether vertical grid lines are drawn.}

\item{grid_col}{Color of the grid lines.}

\item{qtype}{The quantile type used in \code{stats::quantile()} for calculating
confidence-distribution-based intervals from Monte Carlo samples. Default is 8.}

\item{symmetric_shapes}{A logical value indicating whether the summary diamond
and the prediction-interval marker are drawn symmetrically around the midpoint
of the corresponding interval. Default is \code{TRUE}. This is useful for 
publication-quality plots when confidence-distribution-based intervals are 
asymmetric or when transformed labels such as odds ratios are displayed.}

\item{mark_summary_estimate}{A logical value indicating whether the actual
summary point estimate is marked inside the summary diamond. Default is
\code{FALSE}.}

\item{mark_prediction_estimate}{A logical value indicating whether the actual
predictive point estimate is marked on the prediction-interval row. Default is
\code{FALSE}.}

\item{estimate_mark_col}{Color of the marks for the actual summary and
predictive point estimates.}

\item{estimate_mark_lwd}{Line width of the marks for the actual summary and
predictive point estimates.}

\item{mar}{Plot margins passed to \code{graphics::par()}.}

\item{...}{Additional arguments passed to \code{graphics::plot()}.}
}
\value{
Invisibly returns a data frame containing the plotted study-specific estimates,
standard errors, confidence intervals, and weights. The returned data frame has
the following columns:
\itemize{
\item \code{study}: Study label.
\item \code{estimate}: Study-specific effect estimate.
\item \code{se}: Within-study standard error.
\item \code{ci.lb}: Lower limit of the study-specific confidence interval.
\item \code{ci.ub}: Upper limit of the study-specific confidence interval.
\item \code{weight}: Study weight in percent.
}
The returned object also has attributes \code{"overall"} and \code{"prediction"}
containing the plotted overall effect and prediction interval summaries. These
returned values are on the original analysis scale, even when \code{atransf} is
used for plotting.
}
\description{
Draws a forest plot from an object of class \code{"cdmeta"}. The plot shows
study-specific effect estimates and confidence intervals, the 
confidence-distribution-based summary estimate for the overall mean effect,
and optionally the prediction interval for the effect in a future study.
}
\details{
Study-specific confidence intervals are calculated using the normal 
approximation
\deqn{
y_i \pm z_{1-\alpha/2} se_i,
}
where \eqn{\alpha = 1 - level}. The summary estimate for the overall mean
effect \eqn{\mu} is obtained from the \code{cdmeta} object. The interval for
\eqn{\mu} is calculated from the Monte Carlo samples in \code{x$draws$mu} when
available. Similarly, the prediction interval is calculated from
\code{x$draws$theta_new} when available.

By default, the plot is drawn on the original analysis scale. For ratio
measures analyzed on the log scale, such as log odds ratios, log risk ratios,
or log hazard ratios, use \code{atransf = exp} to print axis labels and
estimates on the ratio scale. In that case, tick marks should generally be
specified on the log scale via the \code{at} argument.

If the \code{cdmeta} object was created with \code{transf = exp}, then
\code{forest_cdmeta()} automatically uses this transformation when 
\code{atransf = NULL}. Thus, for log odds ratios, the following two approaches
give the same displayed scale:
\preformatted{
fit <- cdmeta(y, se, transf = exp, transf_name = "exp")
forest_cdmeta(fit, at = log(c(0.25, 0.5, 1, 2, 4)), xlab = "Odds ratio")

fit <- cdmeta(y, se)
forest_cdmeta(fit, atransf = exp,
              at = log(c(0.25, 0.5, 1, 2, 4)), xlab = "Odds ratio")
}

When \code{symmetric_shapes = TRUE}, the summary diamond and the 
prediction-interval marker are drawn around the midpoint of their corresponding
intervals in the plotting coordinates. This keeps the graphical shapes visually
balanced, especially when confidence-distribution-based intervals are asymmetric.
The printed numerical annotations still show the actual point estimates and
interval limits. If \code{mark_summary_estimate = TRUE} or 
\code{mark_prediction_estimate = TRUE}, the actual point estimates are also
marked in the plot.
}
\examples{
data(hf_iron)

fit_hf <- cdmeta(
  y = hf_iron$yi,
  se = hf_iron$sei,
  B = 10000,
  seed = 11111,
  transf = exp,
  transf_name = "exp"
)

forest_cdmeta(
  fit_hf,
  slab = hf_iron$study,
  at = log(c(0.25, 0.5, 1, 2, 4)),
  xlab = "Risk ratio",
  mark_summary_estimate = TRUE,
  mark_prediction_estimate = TRUE
)

data(af_vka)

fit_af <- cdmeta(
  y = af_vka$yi,
  se = af_vka$sei,
  B = 10000,
  seed = 11111,
  i2_method = "harmonic_mean_se2",
  mu_dist = "t",
  df = nrow(af_vka) - 1,
  transf = exp,
  transf_name = "exp"
)

forest_cdmeta(
  fit_af,
  slab = af_vka$study,
  order = "weight",
  summary_stat = "median",
  at = log(c(0.125, 0.25, 0.5, 1, 2)),
  xlab = "Odds ratio",
  main = "Vitamin K antagonists for atrial fibrillation",
  mark_summary_estimate = TRUE,
  mark_prediction_estimate = TRUE
)
}
\seealso{
\code{\link{cdmeta}},
}

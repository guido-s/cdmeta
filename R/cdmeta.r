#' Confidence-distribution-based inference for random-effects meta-analysis
#'
#' Performs unified confidence-distribution-based inference for random-effects
#' meta-analysis. The function provides inference for the overall mean effect
#' \eqn{\mu}, heterogeneity variance \eqn{\tau^2}, heterogeneity standard
#' deviation \eqn{\tau}, heterogeneity proportion \eqn{I^2}, and the effect in
#' a future study \eqn{\theta_{\mathrm{new}}}.
#'
#' The function first obtains Monte Carlo samples of the between-study variance
#' \eqn{\tau^2}. These samples are obtained either from
#' `pimeta::pima(..., method = "boot")` or from a user-supplied vector
#' `tau2_samples`.
#'
#' Given sampled values of \eqn{\tau^2}, the function performs conditional
#' sampling of the overall mean effect. If `mu_dist = "normal"`, then
#' \deqn{\mu \mid \tau^2, y \sim N\{\hat{\mu}(\tau^2), V_{\mu}(\tau^2)\}.}
#' If `mu_dist = "t"`, then a \eqn{t} distribution with `df` degrees of freedom
#' is used instead. The predictive distribution for a future study effect is
#' then generated as
#' \deqn{\theta_{\mathrm{new}} \mid \mu, \tau^2 \sim N(\mu, \tau^2).}
#'
#' For \eqn{I^2}, the sampled \eqn{\tau^2} values are transformed using a
#' reference within-study variance. The reference variance is specified by
#' `i2_method`.
#'
#' If `transf` is supplied, all calculations are still performed on the original
#' analysis scale of `y`. The transformation is applied after Monte Carlo
#' sampling. Point estimates for `mu`, `theta_new`, and `mu_plugin` are shown on
#' the transformed scale, and interval estimates for `mu` and `theta_new` are
#' calculated from transformed Monte Carlo samples. Heterogeneity measures
#' \eqn{\tau^2}, \eqn{\tau}, and \eqn{I^2} are not transformed.
#'
#' For ratio measures analyzed on the log scale, such as log odds ratios, log
#' risk ratios, or log hazard ratios, `transf = exp` can be used to print the
#' main effect-scale summaries on the ratio scale.
#'
#' @param y A numeric vector of study-specific effect estimates (e.g., MD, SMD,
#'   log OR, log RR, or log HR).
#' @param se A numeric vector of within-study standard errors of `y`.
#' @param alpha The significance level for interval estimation. Default is 0.05;
#'   the `alpha/2`th and `1-alpha/2`th quantiles are outputted as the lower and
#'   upper limits of the interval estimates.
#' @param B The number of Monte Carlo samples. When `tau2_samples` is not
#'   supplied, this value is also passed to `pimeta::pima()` as the number of
#'   bootstrap samples. Default is 25000.
#' @param seed An optional numeric value that determines the random seed for
#'   reproducibility. Default is `NULL`.
#' @param parallel A logical value passed to `pimeta::pima()` specifying whether
#'   parallel computation is used. Default is `FALSE`.
#' @param tau2_samples An optional numeric vector of externally supplied
#'   \eqn{\tau^2} samples. If supplied, `pimeta::pima()` is not called.
#' @param i2_method A character string specifying the reference within-study
#'   variance used for calculating \eqn{I^2}. The available options are
#'   `"mean_se2"`, `"harmonic_mean_se2"`, and `"typical_se2"`. Default is
#'   `"mean_se2"`.
#' @param mu_dist A character string specifying the distribution used for
#'   conditional sampling of the overall mean effect \eqn{\mu}. The available
#'   options are `"normal"` and `"t"`. Default is `"normal"`.
#' @param df The degrees of freedom used when `mu_dist = "t"`. If `NULL`, the
#'   default value is \eqn{K - 1}, where \eqn{K} is the number of studies.
#' @param qtype The quantile type used in `stats::quantile()`. Default is 8.
#' @param transf An optional transformation function applied to effect-scale
#'   summaries. For example, `transf = exp` can be used when `y` is on the log
#'   odds ratio, log risk ratio, or log hazard ratio scale. The transformation is
#'   applied to `mu`, `theta_new`, and `mu_plugin`, but not to heterogeneity
#'   measures such as \eqn{\tau^2}, \eqn{\tau}, or \eqn{I^2}.
#' @param transf_name An optional character string giving the name of the
#'   transformation function. For example, `transf_name = "exp"`. If `NULL`, the
#'   name is inferred when possible.
#' @param ... Additional arguments passed to `pimeta::pima()` in `cdmeta()`. For
#'   the print method, additional arguments are currently not used.
#'
#' @return An object of class `"cdmeta"`. The main components are as follows.
#' \itemize{
#' \item `call`: The matched function call.
#' \item `K`: The number of studies.
#' \item `alpha`: The significance level used for interval estimation.
#' \item `B`: The number of Monte Carlo samples.
#' \item `i2_method`: The method used to define the reference within-study variance for \eqn{I^2}.
#' \item `mu_dist`: The distribution used for conditional sampling of \eqn{\mu}.
#' \item `df`: The degrees of freedom used when `mu_dist = "t"`.
#' \item `se2_reference`: The reference within-study variance used for calculating \eqn{I^2}.
#' \item `estimate`: A list of point estimates based on Monte Carlo means, including `mu`, `tau2`, `tau`, `theta_new`, and `I2`. It also includes the plug-in summaries `mu_plugin` and `var_mu_plugin`.
#' \item `median`: A list of confidence-distribution medians for `mu`, `tau2`, `tau`, `theta_new`, and `I2`.
#' \item `interval`: A list of interval estimates for `mu`, `tau2`, `tau`, `theta_new`, and `I2`.
#' \item `draws`: A list of Monte Carlo samples of `tau2`, `tau`, `muhat`, `var_mu`, `mu`, `theta_new`, and `I2`.
#' \item `data`: A list containing the input data `y`, `se`, and `var`.
#' \item `pimeta`: The output from `pimeta::pima()` when `tau2_samples` is not supplied; otherwise `NULL`.
#' \item `transformation`: A list documenting the transformation supplied through `transf`, if any.
#' \item `transformed`: A list of transformed summaries and transformed Monte Carlo samples for `mu` and `theta_new`. This component is `NULL` when `transf` is not supplied.
#' }
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
#' @seealso [forest.cdmeta()]
#' @export

cdmeta <- function(
  y,
  se,
  alpha = 0.05,
  B = 25000,
  seed = NULL,
  parallel = FALSE,
  tau2_samples = NULL,
  i2_method = c("mean_se2", "harmonic_mean_se2", "typical_se2"),
  mu_dist = c("normal", "t"),
  df = NULL,
  qtype = 8,
  transf = NULL,
  transf_name = NULL,
  ...) {
  # -----------------------------
  # Input checks
  # -----------------------------
  if (!is.numeric(y) || !is.numeric(se)) {
    stop("'y' and 'se' must be numeric vectors.")
  }
  if (length(y) != length(se)) {
    stop("'y' and 'se' must have the same length.")
  }
  if (length(y) < 2L) {
    stop("At least 2 studies are required.")
  }
  if (anyNA(y) || anyNA(se) || any(!is.finite(y)) || any(!is.finite(se))) {
    stop("'y' and 'se' must contain only finite, non-missing values.")
  }
  if (any(se <= 0)) {
    stop("All values of 'se' must be positive.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single number in (0, 1).")
  }
  if (!is.numeric(B) || length(B) != 1L || B <= 0) {
    stop("'B' must be a positive integer.")
  }
  if (!is.numeric(qtype) || length(qtype) != 1L) {
    stop("'qtype' must be a single numeric value.")
  }
  if (!is.null(transf) && !is.function(transf)) {
    stop("'transf' must be a function, for example exp.")
  }

  B <- as.integer(B)
  K <- length(y)
  v <- se^2
  i2_method <- match.arg(i2_method)
  mu_dist <- match.arg(mu_dist)

  if (is.null(df)) {
    df <- K - 1
  }
  if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df <= 0) {
    stop("'df' must be a single positive number.")
  }

  if (!is.null(transf) && is.null(transf_name)) {
    transf_name <- deparse(substitute(transf))
  }
  if (!is.null(transf_name)) {
    transf_name <- as.character(transf_name)[1L]
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # -----------------------------
  # tau^2 samples
  # -----------------------------
  pimeta_fit <- NULL

  if (!is.null(tau2_samples)) {
    tau2_draws <- as.numeric(tau2_samples)
    tau2_draws <- tau2_draws[is.finite(tau2_draws) & tau2_draws >= 0]

    if (length(tau2_draws) < 10L) {
      stop("The supplied 'tau2_samples' contains too few valid values.")
    }
  } else {
    if (!requireNamespace("pimeta", quietly = TRUE)) {
      stop("Package 'pimeta' is required when 'tau2_samples' is not supplied.")
    }

    pimeta_fit <- pimeta::pima(
      y = y,
      se = se,
      alpha = alpha,
      method = "boot",
      B = B,
      parallel = parallel,
      seed = seed,
      ...
    )

    if (is.null(pimeta_fit$rnd)) {
      stop("Could not extract tau^2 samples from pimeta output ('rnd').")
    }

    tau2_draws <- as.numeric(pimeta_fit$rnd)
    tau2_draws <- tau2_draws[is.finite(tau2_draws) & tau2_draws >= 0]

    if (length(tau2_draws) < 10L) {
      stop("No sufficient valid tau^2 samples were obtained from pimeta.")
    }
  }

  # Make exactly B draws
  if (length(tau2_draws) >= B) {
    tau2_draws <- tau2_draws[seq_len(B)]
  } else {
    tau2_draws <- sample(tau2_draws, size = B, replace = TRUE)
  }

  tau_draws <- sqrt(tau2_draws)

  # -----------------------------
  # I^2 denominator
  # -----------------------------
  se2_ref <- switch(
    i2_method,
    mean_se2 = mean(v),
    harmonic_mean_se2 = 1 / mean(1 / v),
    typical_se2 = {
      w0 <- 1 / v
      df_typical <- K - 1
      sum_w0 <- sum(w0)
      typical_v <- (df_typical * sum_w0) / (sum_w0^2 - sum(w0^2))
      typical_v
    }
  )

  # -----------------------------
  # Monte Carlo sampling
  # -----------------------------
  muhat_draws <- numeric(B)
  var_mu_draws <- numeric(B)
  mu_draws <- numeric(B)
  theta_new_draws <- numeric(B)

  for (b in seq_len(B)) {
    tau2_b <- tau2_draws[b]

    w_b <- 1 / (v + tau2_b)
    sum_w_b <- sum(w_b)

    muhat_b <- sum(w_b * y) / sum_w_b
    var_mu_b <- 1 / sum_w_b

    if (mu_dist == "normal") {
      mu_b <- stats::rnorm(1L, mean = muhat_b, sd = sqrt(var_mu_b))
    } else {
      mu_b <- muhat_b + sqrt(var_mu_b) * stats::rt(1L, df = df)
    }

    theta_b <- stats::rnorm(1L, mean = mu_b, sd = sqrt(tau2_b))

    muhat_draws[b] <- muhat_b
    var_mu_draws[b] <- var_mu_b
    mu_draws[b] <- mu_b
    theta_new_draws[b] <- theta_b
  }

  # -----------------------------
  # I^2 draws
  # -----------------------------
  I2_draws <- tau2_draws / (tau2_draws + se2_ref)
  I2_draws[I2_draws < 0] <- 0
  I2_draws[I2_draws > 1] <- 1

  # -----------------------------
  # Summaries
  # -----------------------------
  probs <- c(alpha / 2, 1 - alpha / 2)

  qfun <- function(x) {
    stats::quantile(x, probs = probs, names = FALSE, type = qtype, na.rm = TRUE)
  }

  tau2_plugin <- mean(tau2_draws)
  w_plugin <- 1 / (v + tau2_plugin)
  mu_plugin <- sum(w_plugin * y) / sum(w_plugin)
  var_mu_plugin <- 1 / sum(w_plugin)

  estimate <- list(
    mu = mean(mu_draws),
    tau2 = mean(tau2_draws),
    tau = mean(tau_draws),
    theta_new = mean(theta_new_draws),
    I2 = mean(I2_draws),
    mu_plugin = mu_plugin,
    var_mu_plugin = var_mu_plugin
  )

  median <- list(
    mu = stats::median(mu_draws),
    tau2 = stats::median(tau2_draws),
    tau = stats::median(tau_draws),
    theta_new = stats::median(theta_new_draws),
    I2 = stats::median(I2_draws)
  )

  interval <- list(
    mu = stats::setNames(qfun(mu_draws), c("lower", "upper")),
    tau2 = stats::setNames(qfun(tau2_draws), c("lower", "upper")),
    tau = stats::setNames(qfun(tau_draws), c("lower", "upper")),
    theta_new = stats::setNames(qfun(theta_new_draws), c("lower", "upper")),
    I2 = stats::setNames(qfun(I2_draws), c("lower", "upper"))
  )

  draws <- list(
    tau2 = tau2_draws,
    tau = tau_draws,
    muhat = muhat_draws,
    var_mu = var_mu_draws,
    mu = mu_draws,
    theta_new = theta_new_draws,
    I2 = I2_draws
  )

  # -----------------------------
  # Optional transformed summaries
  # -----------------------------
  transformed <- NULL

  if (!is.null(transf)) {
    safe_transf <- function(z) {
      zz <- transf(z)
      zz <- as.numeric(zz)
      if (length(zz) != length(z)) {
        stop("'transf' must return a vector with the same length as its input.")
      }
      zz
    }

    safe_transf_scalar <- function(z) {
      zz <- transf(z)
      zz <- as.numeric(zz)
      if (length(zz) < 1L) {
        stop("'transf' must return at least one numeric value.")
      }
      zz[1L]
    }

    mu_draws_t <- safe_transf(mu_draws)
    theta_new_draws_t <- safe_transf(theta_new_draws)

    qfun_t <- function(z) {
      z <- z[is.finite(z)]
      if (length(z) < 2L) {
        return(c(NA_real_, NA_real_))
      }
      stats::quantile(z, probs = probs, names = FALSE, type = qtype, na.rm = TRUE)
    }

    transformed <- list(
      estimate = list(
        mu = safe_transf_scalar(estimate$mu),
        theta_new = safe_transf_scalar(estimate$theta_new),
        mu_plugin = safe_transf_scalar(estimate$mu_plugin)
      ),
      median = list(
        mu = safe_transf_scalar(median$mu),
        theta_new = safe_transf_scalar(median$theta_new)
      ),
      interval = list(
        mu = stats::setNames(qfun_t(mu_draws_t), c("lower", "upper")),
        theta_new = stats::setNames(qfun_t(theta_new_draws_t), c("lower", "upper"))
      ),
      draws = list(
        mu = mu_draws_t,
        theta_new = theta_new_draws_t
      )
    )
  }

  out <- list(
    call = match.call(),
    K = K,
    alpha = alpha,
    B = B,
    i2_method = i2_method,
    mu_dist = mu_dist,
    df = df,
    se2_reference = se2_ref,
    estimate = estimate,
    median = median,
    interval = interval,
    draws = draws,
    data = list(
      y = y,
      se = se,
      var = v
    ),
    pimeta = pimeta_fit,
    transformation = list(
      transf = transf,
      transf_name = transf_name,
      applies_to = c("mu", "theta_new", "mu_plugin")
    ),
    transformed = transformed
  )

  class(out) <- "cdmeta"
  out
}


#' @param x An object of class `"cdmeta"`.
#' @param digits The number of digits to print.
#'
#' @return The print method returns the object invisibly.
#'
#' @rdname cdmeta
#' @method print cdmeta
#' @export

print.cdmeta <- function(
  x,
  digits = 4,
  transf = NULL,
  transf_name = NULL,
  qtype = 8,
  ...) {
  stopifnot(inherits(x, "cdmeta"))

  if (!is.null(transf) && !is.function(transf)) {
    stop("'transf' must be a function, for example exp.")
  }

  conf_level <- 100 * (1 - x$alpha)

  probs <- c(x$alpha / 2, 1 - x$alpha / 2)

  qfun <- function(z) {
    z <- z[is.finite(z)]
    if (length(z) < 2L) {
      return(c(NA_real_, NA_real_))
    }
    stats::quantile(z, probs = probs, names = FALSE, type = qtype, na.rm = TRUE)
  }

  use_transf <- FALSE
  transformed <- NULL
  tname <- NULL

  if (!is.null(transf)) {
    use_transf <- TRUE
    if (is.null(transf_name)) {
      transf_name <- deparse(substitute(transf))
    }
    tname <- as.character(transf_name)[1L]

    safe_transf_scalar <- function(z) as.numeric(transf(z))[1L]
    safe_transf <- function(z) as.numeric(transf(z))

    mu_draws_t <- safe_transf(x$draws$mu)
    theta_draws_t <- safe_transf(x$draws$theta_new)

    transformed <- list(
      estimate = list(
        mu = safe_transf_scalar(x$estimate$mu),
        theta_new = safe_transf_scalar(x$estimate$theta_new),
        mu_plugin = safe_transf_scalar(x$estimate$mu_plugin)
      ),
      interval = list(
        mu = stats::setNames(qfun(mu_draws_t), c("lower", "upper")),
        theta_new = stats::setNames(qfun(theta_draws_t), c("lower", "upper"))
      )
    )
  } else if (!is.null(x$transformed)) {
    use_transf <- TRUE
    transformed <- x$transformed
    tname <- x$transformation$transf_name
  }

  mu_est <- x$estimate$mu
  theta_est <- x$estimate$theta_new
  mu_plugin <- x$estimate$mu_plugin
  mu_ci <- x$interval$mu
  theta_ci <- x$interval$theta_new

  if (isTRUE(use_transf)) {
    mu_est <- transformed$estimate$mu
    theta_est <- transformed$estimate$theta_new
    mu_plugin <- transformed$estimate$mu_plugin
    mu_ci <- transformed$interval$mu
    theta_ci <- transformed$interval$theta_new
  }

  cat("Unified CD-based inference for random-effects meta-analysis\n")
  cat("----------------------------------------------------------\n")
  cat("Number of studies (K):", x$K, "\n")
  cat("Monte Carlo size (B):", x$B, "\n")
  cat("I^2 reference variance method:", x$i2_method, "\n")
  cat("mu distribution:", x$mu_dist, "\n")
  if (x$mu_dist == "t") {
    cat("degrees of freedom:", x$df, "\n")
  }
  if (isTRUE(use_transf)) {
    cat("Effect scale transformation:", tname, "\n")
    cat("Note: mu, theta_new, and mu_plugin are shown on the transformed scale.\n")
  }
  cat("\n")

  cat("Point estimates (Monte Carlo means)\n")
  cat("-----------------------------------\n")
  cat("mu        =", round(mu_est, digits), "\n")
  cat("tau^2     =", round(x$estimate$tau2, digits), "\n")
  cat("tau       =", round(x$estimate$tau, digits), "\n")
  cat("theta_new =", round(theta_est, digits), "\n")
  cat("I^2       =", round(x$estimate$I2, digits), "\n")
  cat("\n")

  cat("Additional plug-in summary\n")
  cat("--------------------------\n")
  cat("mu_plugin =", round(mu_plugin, digits), "\n")
  cat("Var(mu_plugin) =", round(x$estimate$var_mu_plugin, digits), "\n")
  cat("\n")

  cat(sprintf("%.1f%% interval estimates\n", conf_level))
  cat("--------------------------\n")
  cat("mu CI        = [",
      round(mu_ci[1], digits), ", ",
      round(mu_ci[2], digits), "]\n", sep = "")
  cat("tau^2 CI     = [",
      round(x$interval$tau2[1], digits), ", ",
      round(x$interval$tau2[2], digits), "]\n", sep = "")
  cat("tau CI       = [",
      round(x$interval$tau[1], digits), ", ",
      round(x$interval$tau[2], digits), "]\n", sep = "")
  cat("theta_new PI = [",
      round(theta_ci[1], digits), ", ",
      round(theta_ci[2], digits), "]\n", sep = "")
  cat("I^2 CI       = [",
      round(x$interval$I2[1], digits), ", ",
      round(x$interval$I2[2], digits), "]\n", sep = "")

  invisible(x)
}


plot_cdmeta <- function(
  x,
  show_interval = TRUE,
  breaks = 40,
  main = NULL,
  xlab = NULL,
  hist_col = "grey85",
  border_col = "white",
  density_lwd = 2,
  pi_col = "firebrick",
  pi_lwd = 3,
  rug = FALSE,
  ask_resize = TRUE,
  transf = NULL,
  transf_name = NULL,
  qtype = 8,
  ...
) {
  if (!inherits(x, "cdmeta")) {
    stop("'x' must be an object of class 'cdmeta'.")
  }

  if (!is.null(transf) && !is.function(transf)) {
    stop("'transf' must be a function, for example exp.")
  }

  if (is.null(transf) &&
      !is.null(x$transformation) &&
      !is.null(x$transformation$transf)) {
    transf <- x$transformation$transf
    transf_name <- x$transformation$transf_name
  }

  theta <- x$draws$theta_new
  theta_mean <- x$estimate$theta_new

  probs <- c(x$alpha / 2, 1 - x$alpha / 2)

  if (!is.null(transf)) {
    theta <- as.numeric(transf(theta))
    theta_mean <- as.numeric(transf(theta_mean))[1L]
    theta <- theta[is.finite(theta)]

    pi_lim <- stats::quantile(
      theta,
      probs = probs,
      names = FALSE,
      type = qtype,
      na.rm = TRUE
    )

    if (is.null(transf_name)) {
      transf_name <- deparse(substitute(transf))
    }
    if (is.null(main)) {
      main <- expression("Predictive distribution of transformed " * theta[new])
    }
    if (is.null(xlab)) {
      xlab <- paste0("Transformed theta_new")
    }
  } else {
    pi_lim <- x$interval$theta_new

    if (is.null(main)) {
      main <- expression("Predictive distribution of " * theta[new])
    }
    if (is.null(xlab)) {
      xlab <- expression(theta[new])
    }
  }

  tryCatch({
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)

    if (isTRUE(show_interval)) {
      graphics::layout(matrix(c(1, 2), nrow = 2), heights = c(4, 1))
      graphics::par(mar = c(2.5, 4.1, 2.6, 1.1))
    } else {
      graphics::par(mar = c(4.1, 4.1, 2.6, 1.1))
    }

    graphics::hist(
      theta,
      breaks = breaks,
      probability = TRUE,
      col = hist_col,
      border = border_col,
      main = main,
      xlab = if (show_interval) "" else xlab
    )

    dens <- stats::density(theta, ...)
    graphics::lines(dens, lwd = density_lwd)

    graphics::abline(v = pi_lim, col = pi_col, lwd = 2, lty = 2)
    graphics::abline(v = theta_mean, col = pi_col, lwd = 1.5, lty = 1)

    if (isTRUE(rug)) {
      graphics::rug(theta)
    }

    if (isTRUE(show_interval)) {
      graphics::par(mar = c(3.6, 4.1, 0.2, 1.1))
      xr <- range(theta, finite = TRUE)

      graphics::plot(
        NA, NA,
        xlim = xr,
        ylim = c(0, 1),
        xlab = xlab,
        ylab = "",
        yaxt = "n",
        bty = "n"
      )

      graphics::segments(
        pi_lim[1], 0.5, pi_lim[2], 0.5,
        lwd = pi_lwd, col = pi_col
      )

      graphics::points(
        theta_mean, 0.5,
        pch = 16, cex = 1.2, col = pi_col
      )

      graphics::text(
        theta_mean, 0.8,
        labels = paste0(round(100 * (1 - x$alpha)), "% PI"),
        col = pi_col
      )
    }

    invisible(x)
  }, error = function(e) {
    if (isTRUE(ask_resize)) {
      message("Plot skipped: increase the plot window size or plot to a file with png()/pdf().")
    }
    invisible(x)
  })
}

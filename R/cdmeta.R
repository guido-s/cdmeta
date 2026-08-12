# Internal utility functions -------------------------------------------------

.cdmeta_is_integerish <- function(x) {
  is.numeric(x) &&
    length(x) == 1L &&
    is.finite(x) &&
    x == floor(x)
}

.cdmeta_validate_qtype <- function(qtype) {
  if (!.cdmeta_is_integerish(qtype) || qtype < 1L || qtype > 9L) {
    stop("'qtype' must be a single integer from 1 to 9.", call. = FALSE)
  }
  as.integer(qtype)
}

.cdmeta_validate_transf_name <- function(transf_name) {
  if (is.null(transf_name)) {
    return(NULL)
  }

  if (!is.character(transf_name) ||
      length(transf_name) != 1L ||
      is.na(transf_name) ||
      !nzchar(transf_name)) {
    stop(
      "'transf_name' must be NULL or a single non-empty character string.",
      call. = FALSE
    )
  }

  transf_name
}

.cdmeta_apply_transf <- function(transf, x, argument = "'transf'") {
  out <- transf(x)

  if (!is.numeric(out)) {
    stop(argument, " must return numeric values.", call. = FALSE)
  }

  out <- as.numeric(out)

  if (length(out) != length(x)) {
    stop(
      argument,
      " must return a vector with the same length as its input.",
      call. = FALSE
    )
  }

  if (anyNA(out) || any(!is.finite(out))) {
    stop(argument, " returned non-finite or missing values.", call. = FALSE)
  }

  out
}

.cdmeta_apply_transf_scalar <- function(transf, x, argument = "'transf'") {
  out <- transf(x)

  if (!is.numeric(out)) {
    stop(argument, " must return a numeric value.", call. = FALSE)
  }

  out <- as.numeric(out)

  if (length(out) != 1L || is.na(out) || !is.finite(out)) {
    stop(
      argument,
      " must return one finite, non-missing numeric value for scalar input.",
      call. = FALSE
    )
  }

  out
}


# Main function -------------------------------------------------------------

cdmeta <- function(
  y,
  se,
  alpha = 0.05,
  B = 25000,
  seed = NULL,
  parallel = FALSE,
  tau2_samples = NULL,
  i2_method = c("typical_se2", "mean_se2", "harmonic_mean_se2"),
  mu_dist = c("normal", "t"),
  df = NULL,
  qtype = 8,
  transf = NULL,
  transf_name = NULL,
  ...
) {
  B_supplied <- !missing(B)

  # -----------------------------
  # Input checks
  # -----------------------------
  if (!is.numeric(y) || !is.null(dim(y)) ||
      !is.numeric(se) || !is.null(dim(se))) {
    stop("'y' and 'se' must be numeric vectors.", call. = FALSE)
  }
  if (length(y) != length(se)) {
    stop("'y' and 'se' must have the same length.", call. = FALSE)
  }
  if (length(y) < 2L) {
    stop("At least 2 studies are required.", call. = FALSE)
  }
  if (anyNA(y) || anyNA(se) || any(!is.finite(y)) || any(!is.finite(se))) {
    stop(
      "'y' and 'se' must contain only finite, non-missing values.",
      call. = FALSE
    )
  }
  if (any(se <= 0)) {
    stop("All values of 'se' must be positive.", call. = FALSE)
  }

  if (!is.numeric(alpha) ||
      length(alpha) != 1L ||
      !is.finite(alpha) ||
      alpha <= 0 ||
      alpha >= 1) {
    stop("'alpha' must be a single finite number in (0, 1).", call. = FALSE)
  }

  if (!.cdmeta_is_integerish(B) ||
      B < 10L ||
      B > .Machine$integer.max) {
    stop(
      "'B' must be a single integer from 10 to .Machine$integer.max.",
      call. = FALSE
    )
  }
  B <- as.integer(B)

  if (!is.null(seed)) {
    if (!.cdmeta_is_integerish(seed) ||
        seed < 0 ||
        seed > .Machine$integer.max) {
      stop(
        "'seed' must be NULL or a non-negative integer not exceeding ",
        ".Machine$integer.max.",
        call. = FALSE
      )
    }
    seed <- as.integer(seed)
  }

  parallel_is_valid <- identical(parallel, FALSE) ||
    (.cdmeta_is_integerish(parallel) &&
       parallel >= 1L &&
       parallel <= .Machine$integer.max)

  if (!parallel_is_valid) {
    stop(
      "'parallel' must be FALSE for single-threaded computation or a positive integer specifying the number of threads.",
      call. = FALSE
    )
  }
  if (!identical(parallel, FALSE)) {
    parallel <- as.integer(parallel)
  }

  qtype <- .cdmeta_validate_qtype(qtype)

  if (!is.null(transf) && !is.function(transf)) {
    stop("'transf' must be a function, for example exp.", call. = FALSE)
  }

  K <- length(y)
  v <- se^2
  i2_method <- match.arg(i2_method)
  mu_dist <- match.arg(mu_dist)

  if (mu_dist == "t") {
    if (is.null(df)) {
      df <- K - 1
    }
    if (!is.numeric(df) ||
        length(df) != 1L ||
        !is.finite(df) ||
        df <= 0) {
      stop(
        "'df' must be a single positive finite number when 'mu_dist = \"t\"'.",
        call. = FALSE
      )
    }
  } else {
    if (!is.null(df)) {
      warning("'df' is ignored when 'mu_dist = \"normal\"'.", call. = FALSE)
    }
    df <- NULL
  }

  if (!is.null(transf) && is.null(transf_name)) {
    transf_name <- paste(deparse(substitute(transf)), collapse = "")
  }
  transf_name <- .cdmeta_validate_transf_name(transf_name)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # -----------------------------
  # tau^2 samples
  # -----------------------------
  pimeta_fit <- NULL
  tau2_source <- NULL

  if (!is.null(tau2_samples)) {
    if (!is.numeric(tau2_samples) || !is.null(dim(tau2_samples))) {
      stop("'tau2_samples' must be a numeric vector.", call. = FALSE)
    }
    if (length(tau2_samples) < 10L) {
      stop(
        "'tau2_samples' must contain at least 10 values.",
        call. = FALSE
      )
    }
    if (anyNA(tau2_samples) || any(!is.finite(tau2_samples))) {
      stop(
        "'tau2_samples' must contain only finite, non-missing values.",
        call. = FALSE
      )
    }
    if (any(tau2_samples < 0)) {
      stop("All values of 'tau2_samples' must be non-negative.", call. = FALSE)
    }

    tau2_draws <- as.numeric(tau2_samples)

    if (B_supplied && B != length(tau2_draws)) {
      warning(
        "When 'tau2_samples' is supplied, its length determines the Monte Carlo sample size; 'B' is ignored.",
        call. = FALSE
      )
    }

    B <- length(tau2_draws)
    tau2_source <- "supplied"
  } else {
    if (!requireNamespace("pimeta", quietly = TRUE)) {
      stop(
        "Package 'pimeta' is required when 'tau2_samples' is not supplied.",
        call. = FALSE
      )
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
      stop(
        "The installed version of 'pimeta' did not return the required tau^2 draws in component 'rnd'.",
        call. = FALSE
      )
    }
    if (!is.numeric(pimeta_fit$rnd) || !is.null(dim(pimeta_fit$rnd))) {
      stop("'pimeta' returned tau^2 draws that are not a numeric vector.", call. = FALSE)
    }

    tau2_draws <- as.numeric(pimeta_fit$rnd)

    if (length(tau2_draws) != B) {
      stop(
        "'pimeta' returned ", length(tau2_draws),
        " tau^2 draws, but ", B, " were requested.",
        call. = FALSE
      )
    }
    if (anyNA(tau2_draws) || any(!is.finite(tau2_draws))) {
      stop(
        "'pimeta' returned non-finite or missing tau^2 draws.",
        call. = FALSE
      )
    }
    if (any(tau2_draws < 0)) {
      stop("'pimeta' returned negative tau^2 draws.", call. = FALSE)
    }

    tau2_source <- "pimeta"
  }

  if (B < 1000L) {
    warning(
      "Monte Carlo interval estimates may be unstable when fewer than 1000 draws are used.",
      call. = FALSE
    )
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
      sum_w0 <- sum(w0)
      denominator <- sum_w0^2 - sum(w0^2)

      if (!is.finite(denominator) || denominator <= 0) {
        stop(
          "The typical within-study variance could not be computed from the supplied standard errors.",
          call. = FALSE
        )
      }

      ((K - 1) * sum_w0) / denominator
    }
  )

  if (!is.finite(se2_ref) || se2_ref <= 0) {
    stop(
      "The reference within-study variance for I^2 must be finite and positive.",
      call. = FALSE
    )
  }

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
  I2_draws <- pmin(1, pmax(0, I2_draws))

  # -----------------------------
  # Summaries
  # -----------------------------
  probs <- c(alpha / 2, 1 - alpha / 2)

  qfun <- function(x) {
    stats::quantile(
      x,
      probs = probs,
      names = FALSE,
      type = qtype,
      na.rm = FALSE
    )
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
    mu_draws_t <- .cdmeta_apply_transf(transf, mu_draws)
    theta_new_draws_t <- .cdmeta_apply_transf(transf, theta_new_draws)

    qfun_t <- function(z) {
      stats::quantile(
        z,
        probs = probs,
        names = FALSE,
        type = qtype,
        na.rm = FALSE
      )
    }

    transformed <- list(
      estimate = list(
        mu = .cdmeta_apply_transf_scalar(transf, estimate$mu),
        theta_new = .cdmeta_apply_transf_scalar(transf, estimate$theta_new),
        mu_plugin = .cdmeta_apply_transf_scalar(transf, estimate$mu_plugin)
      ),
      median = list(
        mu = .cdmeta_apply_transf_scalar(transf, median$mu),
        theta_new = .cdmeta_apply_transf_scalar(transf, median$theta_new)
      ),
      interval = list(
        mu = stats::setNames(qfun_t(mu_draws_t), c("lower", "upper")),
        theta_new = stats::setNames(
          qfun_t(theta_new_draws_t),
          c("lower", "upper")
        )
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
    seed = seed,
    parallel = parallel,
    qtype = qtype,
    tau2_source = tau2_source,
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


# Print method --------------------------------------------------------------

print.cdmeta <- function(
  x,
  digits = 4,
  transf = NULL,
  transf_name = NULL,
  qtype = 8,
  ...
) {
  if (!inherits(x, "cdmeta")) {
    stop("'x' must be an object of class 'cdmeta'.", call. = FALSE)
  }

  if (!.cdmeta_is_integerish(digits) || digits < 0L || digits > 22L) {
    stop("'digits' must be a single integer from 0 to 22.", call. = FALSE)
  }
  digits <- as.integer(digits)
  qtype <- .cdmeta_validate_qtype(qtype)

  if (!is.null(transf) && !is.function(transf)) {
    stop("'transf' must be a function, for example exp.", call. = FALSE)
  }

  conf_level <- 100 * (1 - x$alpha)
  probs <- c(x$alpha / 2, 1 - x$alpha / 2)

  qfun <- function(z) {
    stats::quantile(
      z,
      probs = probs,
      names = FALSE,
      type = qtype,
      na.rm = FALSE
    )
  }

  use_transf <- FALSE
  transformed <- NULL
  tname <- NULL

  if (!is.null(transf)) {
    use_transf <- TRUE

    if (is.null(transf_name)) {
      transf_name <- paste(deparse(substitute(transf)), collapse = "")
    }
    tname <- .cdmeta_validate_transf_name(transf_name)

    mu_draws_t <- .cdmeta_apply_transf(transf, x$draws$mu)
    theta_draws_t <- .cdmeta_apply_transf(transf, x$draws$theta_new)

    transformed <- list(
      estimate = list(
        mu = .cdmeta_apply_transf_scalar(transf, x$estimate$mu),
        theta_new = .cdmeta_apply_transf_scalar(
          transf,
          x$estimate$theta_new
        ),
        mu_plugin = .cdmeta_apply_transf_scalar(
          transf,
          x$estimate$mu_plugin
        )
      ),
      interval = list(
        mu = stats::setNames(qfun(mu_draws_t), c("lower", "upper")),
        theta_new = stats::setNames(
          qfun(theta_draws_t),
          c("lower", "upper")
        )
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
  if (!is.null(x$tau2_source)) {
    cat("Source of tau^2 draws:", x$tau2_source, "\n")
  }
  cat("I^2 reference variance method:", x$i2_method, "\n")
  cat("mu distribution:", x$mu_dist, "\n")
  if (identical(x$mu_dist, "t")) {
    cat("degrees of freedom:", x$df, "\n")
  }
  if (isTRUE(use_transf)) {
    cat("Effect scale transformation:", tname, "\n")
    cat(
      "Note: mu, theta_new, and mu_plugin are shown on the transformed scale.\n"
    )
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
  cat(
    "mu CI        = [",
    round(mu_ci[1L], digits), ", ",
    round(mu_ci[2L], digits), "]\n",
    sep = ""
  )
  cat(
    "tau^2 CI     = [",
    round(x$interval$tau2[1L], digits), ", ",
    round(x$interval$tau2[2L], digits), "]\n",
    sep = ""
  )
  cat(
    "tau CI       = [",
    round(x$interval$tau[1L], digits), ", ",
    round(x$interval$tau[2L], digits), "]\n",
    sep = ""
  )
  cat(
    "theta_new PI = [",
    round(theta_ci[1L], digits), ", ",
    round(theta_ci[2L], digits), "]\n",
    sep = ""
  )
  cat(
    "I^2 CI       = [",
    round(x$interval$I2[1L], digits), ", ",
    round(x$interval$I2[2L], digits), "]\n",
    sep = ""
  )

  invisible(x)
}


# Plot method ---------------------------------------------------------------

plot.cdmeta <- function(
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
  transf = NULL,
  transf_name = NULL,
  qtype = 8,
  ...
) {
  if (!inherits(x, "cdmeta")) {
    stop("'x' must be an object of class 'cdmeta'.", call. = FALSE)
  }
  if (!is.logical(show_interval) ||
      length(show_interval) != 1L ||
      is.na(show_interval)) {
    stop("'show_interval' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(rug) || length(rug) != 1L || is.na(rug)) {
    stop("'rug' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(density_lwd) ||
      length(density_lwd) != 1L ||
      !is.finite(density_lwd) ||
      density_lwd <= 0) {
    stop("'density_lwd' must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(pi_lwd) ||
      length(pi_lwd) != 1L ||
      !is.finite(pi_lwd) ||
      pi_lwd <= 0) {
    stop("'pi_lwd' must be a single positive number.", call. = FALSE)
  }

  qtype <- .cdmeta_validate_qtype(qtype)

  if (!is.null(transf) && !is.function(transf)) {
    stop("'transf' must be a function, for example exp.", call. = FALSE)
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
    theta <- .cdmeta_apply_transf(transf, theta)
    theta_mean <- .cdmeta_apply_transf_scalar(transf, theta_mean)

    pi_lim <- stats::quantile(
      theta,
      probs = probs,
      names = FALSE,
      type = qtype,
      na.rm = FALSE
    )

    if (is.null(transf_name)) {
      transf_name <- paste(deparse(substitute(transf)), collapse = "")
    }
    transf_name <- .cdmeta_validate_transf_name(transf_name)

    if (is.null(main)) {
      main <- expression("Predictive distribution of transformed " * theta[new])
    }
    if (is.null(xlab)) {
      xlab <- paste0("Transformed theta_new (", transf_name, ")")
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

  if (length(theta) < 2L || anyNA(theta) || any(!is.finite(theta))) {
    stop(
      "The predictive draws must contain at least two finite, non-missing values.",
      call. = FALSE
    )
  }

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::layout(matrix(1L))
    graphics::par(oldpar)
  }, add = TRUE)

  if (isTRUE(show_interval)) {
    graphics::layout(matrix(c(1L, 2L), nrow = 2L), heights = c(4, 1))
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
    xlab = if (show_interval) "" else xlab,
    ...
  )

  if (length(unique(theta)) > 1L) {
    dens <- stats::density(theta)
    graphics::lines(dens, lwd = density_lwd)
  }

  graphics::abline(
    v = pi_lim,
    col = pi_col,
    lwd = max(1, pi_lwd * 2 / 3),
    lty = 2
  )
  graphics::abline(
    v = theta_mean,
    col = pi_col,
    lwd = max(1, pi_lwd / 2),
    lty = 1
  )

  if (isTRUE(rug)) {
    graphics::rug(theta)
  }

  if (isTRUE(show_interval)) {
    graphics::par(mar = c(3.6, 4.1, 0.2, 1.1))

    xr <- range(c(theta, pi_lim, theta_mean), finite = TRUE)
    if (diff(xr) == 0) {
      delta <- if (xr[1L] == 0) 1 else abs(xr[1L]) * 0.04
      xr <- xr + c(-delta, delta)
    }

    graphics::plot(
      NA_real_, NA_real_,
      xlim = xr,
      ylim = c(0, 1),
      xlab = xlab,
      ylab = "",
      yaxt = "n",
      bty = "n"
    )

    graphics::segments(
      pi_lim[1L], 0.5,
      pi_lim[2L], 0.5,
      lwd = pi_lwd,
      col = pi_col
    )

    graphics::points(
      theta_mean, 0.5,
      pch = 16,
      cex = 1.2,
      col = pi_col
    )

    graphics::text(
      theta_mean, 0.8,
      labels = paste0(round(100 * (1 - x$alpha)), "% PI"),
      col = pi_col
    )
  }

  invisible(x)
}

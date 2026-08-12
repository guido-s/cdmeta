forest_cdmeta <- function(
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
) {
  # -----------------------------
  # Input checks
  # -----------------------------
  if (!inherits(x, "cdmeta")) {
    stop("'x' must be an object of class 'cdmeta'.")
  }

  if (is.null(x$data$y) || is.null(x$data$se)) {
    stop("'x' must contain 'x$data$y' and 'x$data$se'.")
  }

  ynames <- names(x$data$y)
  y <- as.numeric(x$data$y)
  se <- as.numeric(x$data$se)

  if (length(y) != length(se)) {
    stop("'x$data$y' and 'x$data$se' must have the same length.")
  }
  if (length(y) < 2L) {
    stop("At least 2 studies are required.")
  }
  if (anyNA(y) || anyNA(se) || any(!is.finite(y)) || any(!is.finite(se))) {
    stop("'x$data$y' and 'x$data$se' must contain only finite, non-missing values.")
  }
  if (any(se <= 0)) {
    stop("All values of 'x$data$se' must be positive.")
  }

  K <- length(y)
  order <- match.arg(order)
  summary_stat <- match.arg(summary_stat)

  if (is.null(level)) {
    level <- 1 - x$alpha
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("'level' must be a single number in (0, 1).")
  }

  alpha_plot <- 1 - level
  zcrit <- stats::qnorm(1 - alpha_plot / 2)

  if (is.null(slab)) {
    if (!is.null(ynames) && length(ynames) == K && all(nzchar(ynames))) {
      slab <- ynames
    } else {
      slab <- paste("Study", seq_len(K))
    }
  }

  slab <- as.character(slab)
  if (length(slab) != K) {
    stop("'slab' must have the same length as the number of studies.")
  }

  if (is.null(atransf) &&
      !is.null(x$transformation) &&
      !is.null(x$transformation$transf)) {
    atransf <- x$transformation$transf
  }

  if (!is.null(atransf) && !is.function(atransf)) {
    stop("'atransf' must be a function, for example exp.")
  }

  # -----------------------------
  # Study-specific quantities
  # -----------------------------
  ci_lb <- y - zcrit * se
  ci_ub <- y + zcrit * se

  tau2_w <- x$estimate$tau2
  if (is.null(tau2_w) || length(tau2_w) != 1L || !is.finite(tau2_w) || tau2_w < 0) {
    tau2_w <- 0
  }

  weights <- 1 / (se^2 + tau2_w)
  weights_pct <- 100 * weights / sum(weights)

  ord <- switch(
    order,
    none = seq_len(K),
    increasing = base::order(y),
    decreasing = base::order(y, decreasing = TRUE),
    precision = base::order(se),
    weight = base::order(weights_pct, decreasing = TRUE)
  )

  y <- y[ord]
  se <- se[ord]
  ci_lb <- ci_lb[ord]
  ci_ub <- ci_ub[ord]
  slab <- slab[ord]
  weights_pct <- weights_pct[ord]

  if (is.null(psize)) {
    psize <- 0.7 + 1.2 * sqrt(weights_pct / max(weights_pct))
  } else {
    if (!is.numeric(psize)) {
      stop("'psize' must be numeric.")
    }
    if (length(psize) == 1L) {
      psize <- rep(psize, K)
    }
    if (length(psize) != K) {
      stop("'psize' must have length 1 or the same length as the number of studies.")
    }
    psize <- psize[ord]
  }

  # -----------------------------
  # CD-based summary quantities
  # -----------------------------
  qfun <- function(z) {
    z <- as.numeric(z)
    z <- z[is.finite(z)]
    if (length(z) < 2L) {
      return(c(NA_real_, NA_real_))
    }
    stats::quantile(
      z,
      probs = c(alpha_plot / 2, 1 - alpha_plot / 2),
      names = FALSE,
      type = qtype,
      na.rm = TRUE
    )
  }

  get_point <- function(name, allow_median = TRUE) {
    if (isTRUE(allow_median) &&
        summary_stat == "median" &&
        !is.null(x$median) &&
        !is.null(x$median[[name]])) {
      return(as.numeric(x$median[[name]])[1L])
    }

    if (!is.null(x$estimate) &&
        !is.null(x$estimate[[name]])) {
      return(as.numeric(x$estimate[[name]])[1L])
    }

    if (!is.null(x$draws) &&
        !is.null(x$draws[[name]])) {
      z <- as.numeric(x$draws[[name]])
      z <- z[is.finite(z)]
      if (length(z) > 0L) {
        if (isTRUE(allow_median) && summary_stat == "median") {
          return(stats::median(z))
        } else {
          return(mean(z))
        }
      }
    }

    NA_real_
  }

  get_interval <- function(name) {
    if (!is.null(x$draws) &&
        !is.null(x$draws[[name]])) {
      z <- as.numeric(x$draws[[name]])
      z <- z[is.finite(z)]
      if (length(z) > 1L) {
        return(qfun(z))
      }
    }

    if (!is.null(x$interval) &&
        !is.null(x$interval[[name]])) {
      return(as.numeric(x$interval[[name]])[1:2])
    }

    c(NA_real_, NA_real_)
  }

  # summary_stat can affect mu and theta_new
  mu_est <- get_point("mu", allow_median = TRUE)
  mu_ci <- get_interval("mu")

  theta_est <- get_point("theta_new", allow_median = TRUE)
  theta_ci <- get_interval("theta_new")

  # Heterogeneity summaries should match print.cdmeta():
  # always use Monte Carlo means from x$estimate, not medians.
  tau2_est <- get_point("tau2", allow_median = FALSE)
  tau2_ci <- get_interval("tau2")

  tau_est <- get_point("tau", allow_median = FALSE)
  tau_ci <- get_interval("tau")

  I2_est <- get_point("I2", allow_median = FALSE)
  I2_ci <- get_interval("I2")

  if (!is.finite(mu_est) || any(!is.finite(mu_ci))) {
    stop("Could not extract the summary estimate and interval for 'mu' from 'x'.")
  }

  # -----------------------------
  # Formatting helpers
  # -----------------------------
  format_number <- function(z, ndigits) {
    z <- as.numeric(z)
    out <- rep("", length(z))
    ok <- is.finite(z)
    out[ok] <- formatC(z[ok], digits = ndigits, format = "f")
    out
  }

  display_value <- function(z) {
    if (is.null(atransf)) {
      return(z)
    }
    atransf(z)
  }

  fmt_eff <- function(z, ndigits = ci_digits) {
    format_number(display_value(z), ndigits)
  }

  fmt_ci <- function(est, lb, ub) {
    paste0(
      fmt_eff(est), " [",
      fmt_eff(lb), ", ",
      fmt_eff(ub), "]"
    )
  }

  fmt_pct <- function(z, ndigits = 1) {
    format_number(100 * z, ndigits)
  }

  # -----------------------------
  # Axis limits and tick marks
  # -----------------------------
  range_values <- c(ci_lb, ci_ub, mu_ci, refline)
  if (isTRUE(show_pi)) {
    range_values <- c(range_values, theta_ci)
  }

  range_values <- range_values[is.finite(range_values)]

  if (length(range_values) == 0L) {
    stop("No finite values are available for plotting.")
  }

  if (is.null(alim)) {
    alim <- range(range_values)
    if (!all(is.finite(alim)) || alim[1] == alim[2]) {
      center <- alim[1]
      if (!is.finite(center)) {
        center <- 0
      }
      alim <- center + c(-1, 1)
    } else {
      pad <- 0.06 * diff(alim)
      alim <- alim + c(-pad, pad)
    }
  } else {
    if (!is.numeric(alim) || length(alim) != 2L || any(!is.finite(alim))) {
      stop("'alim' must be a finite numeric vector of length 2.")
    }
    alim <- sort(alim)
  }

  effect_width <- diff(alim)
  if (!is.finite(effect_width) || effect_width <= 0) {
    effect_width <- 1
  }

  if (is.null(at)) {
    at <- pretty(alim, n = 5)
  }
  at <- at[is.finite(at) & at >= alim[1] & at <= alim[2]]

  axis_labels <- fmt_eff(at, digits)

  if (is.null(xlim)) {
    left_space <- 0.85 * effect_width

    if (isTRUE(annotate) && isTRUE(show_weights)) {
      right_space <- 1.35 * effect_width
    } else if (isTRUE(annotate)) {
      right_space <- 0.95 * effect_width
    } else if (isTRUE(show_weights)) {
      right_space <- 0.80 * effect_width
    } else {
      right_space <- 0.25 * effect_width
    }

    xlim <- c(alim[1] - left_space, alim[2] + right_space)
  } else {
    if (!is.numeric(xlim) || length(xlim) != 2L || any(!is.finite(xlim))) {
      stop("'xlim' must be a finite numeric vector of length 2.")
    }
    xlim <- sort(xlim)
  }

  plot_width <- diff(xlim)
  x_study <- xlim[1] + 0.01 * plot_width
  x_annot <- alim[2] + 0.08 * effect_width
  x_weight <- alim[2] + 0.98 * effect_width

  # -----------------------------
  # Y positions
  # -----------------------------
  y_shift <- if (isTRUE(show_pi)) 4 else 3
  y_study <- rev(seq_len(K)) + y_shift
  y_overall <- if (isTRUE(show_pi)) 3 else 2
  y_pi <- 2
  y_het <- 0.7
  y_header <- max(y_study) + 0.7

  y_bottom <- if (isTRUE(show_het)) 0.1 else 0.8
  ylim <- c(y_bottom, y_header + 0.5)

  # -----------------------------
  # Plot helper
  # -----------------------------
  draw_interval <- function(lb, ub, yy, col = "black", lwd = 1, lty = 1,
                            cap = 0.08) {
    if (!is.finite(lb) || !is.finite(ub)) {
      return(invisible(NULL))
    }

    if (ub < alim[1]) {
      graphics::arrows(
        alim[1] + 0.05 * effect_width, yy,
        alim[1], yy,
        length = 0.06, angle = 30, code = 2,
        col = col, lwd = lwd, lty = lty
      )
      return(invisible(NULL))
    }

    if (lb > alim[2]) {
      graphics::arrows(
        alim[2] - 0.05 * effect_width, yy,
        alim[2], yy,
        length = 0.06, angle = 30, code = 2,
        col = col, lwd = lwd, lty = lty
      )
      return(invisible(NULL))
    }

    lb_clip <- max(lb, alim[1])
    ub_clip <- min(ub, alim[2])

    graphics::segments(lb_clip, yy, ub_clip, yy, col = col, lwd = lwd, lty = lty)

    if (lb >= alim[1]) {
      graphics::segments(lb, yy - cap, lb, yy + cap, col = col, lwd = lwd, lty = lty)
    } else {
      graphics::arrows(
        alim[1] + 0.05 * effect_width, yy,
        alim[1], yy,
        length = 0.06, angle = 30, code = 2,
        col = col, lwd = lwd, lty = lty
      )
    }

    if (ub <= alim[2]) {
      graphics::segments(ub, yy - cap, ub, yy + cap, col = col, lwd = lwd, lty = lty)
    } else {
      graphics::arrows(
        alim[2] - 0.05 * effect_width, yy,
        alim[2], yy,
        length = 0.06, angle = 30, code = 2,
        col = col, lwd = lwd, lty = lty
      )
    }

    invisible(NULL)
  }

  clip_x <- function(z) {
    min(max(z, alim[1]), alim[2])
  }

  midpoint <- function(lb, ub) {
    if (!is.finite(lb) || !is.finite(ub)) {
      return(NA_real_)
    }
    0.5 * (lb + ub)
  }

  # -----------------------------
  # Draw plot
  # -----------------------------
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)

  graphics::par(mar = mar, xpd = NA)

  graphics::plot(
    NA, NA,
    xlim = xlim,
    ylim = ylim,
    axes = FALSE,
    xlab = "",
    ylab = "",
    type = "n",
    bty = "n",
    ...
  )

  if (!is.null(main)) {
    graphics::title(main = main)
  }

  y_grid_low <- if (isTRUE(show_het)) 1.1 else y_bottom
  y_grid_high <- y_header - 0.35

  if (isTRUE(grid) && length(at) > 0L) {
    for (aa in at) {
      graphics::segments(
        aa, y_grid_low, aa, y_grid_high,
        col = grid_col, lty = 3, lwd = 0.8
      )
    }
  }

  if (!is.null(refline) && length(refline) == 1L && is.finite(refline)) {
    graphics::segments(
      refline, y_grid_low, refline, y_grid_high,
      col = refline_col, lty = 2, lwd = 1
    )
  }

  if (isTRUE(header)) {
    graphics::text(
      x_study, y_header, "Study",
      adj = 0, font = 2, cex = cex
    )

    if (isTRUE(annotate)) {
      graphics::text(
        x_annot, y_header,
        paste0("Estimate [", round(100 * level, 1), "% CI]"),
        adj = 0, font = 2, cex = cex
      )
    }

    if (isTRUE(show_weights)) {
      graphics::text(
        x_weight, y_header, "Weight",
        adj = 0.5, font = 2, cex = cex
      )
    }

    graphics::segments(
      xlim[1], y_header - 0.35,
      xlim[2], y_header - 0.35,
      lwd = 1
    )
  }

  # -----------------------------
  # Study rows
  # -----------------------------
  for (i in seq_len(K)) {
    graphics::text(
      x_study, y_study[i], slab[i],
      adj = 0, cex = cex
    )

    draw_interval(
      ci_lb[i], ci_ub[i], y_study[i],
      col = ci_col, lwd = 1, lty = 1
    )

    if (y[i] >= alim[1] && y[i] <= alim[2]) {
      graphics::points(
        y[i], y_study[i],
        pch = pch,
        col = box_col,
        bg = box_bg,
        cex = psize[i]
      )
    }

    if (isTRUE(annotate)) {
      graphics::text(
        x_annot, y_study[i],
        fmt_ci(y[i], ci_lb[i], ci_ub[i]),
        adj = 0, cex = cex
      )
    }

    if (isTRUE(show_weights)) {
      graphics::text(
        x_weight, y_study[i],
        paste0(format_number(weights_pct[i], weight_digits), "%"),
        adj = 0.5, cex = cex
      )
    }
  }

  # Separator before summary
  graphics::segments(
    xlim[1], y_overall + 0.55,
    xlim[2], y_overall + 0.55,
    lwd = 0.8
  )

  # -----------------------------
  # Overall-effect diamond
  # -----------------------------
  graphics::text(
    x_study, y_overall, "Overall effect",
    adj = 0, font = 2, cex = cex
  )

  diamond_h <- 0.28

  mu_lb_clip <- max(mu_ci[1], alim[1])
  mu_ub_clip <- min(mu_ci[2], alim[2])

  if (isTRUE(symmetric_shapes)) {
    mu_center <- midpoint(mu_ci[1], mu_ci[2])
  } else {
    mu_center <- mu_est
  }

  if (!is.finite(mu_center)) {
    mu_center <- mu_est
  }

  mu_center_clip <- clip_x(mu_center)

  graphics::polygon(
    x = c(mu_lb_clip, mu_center_clip, mu_ub_clip, mu_center_clip),
    y = c(y_overall, y_overall + diamond_h, y_overall, y_overall - diamond_h),
    border = summary_col,
    col = summary_bg
  )

  if (mu_ci[1] < alim[1]) {
    graphics::arrows(
      alim[1] + 0.05 * effect_width, y_overall,
      alim[1], y_overall,
      length = 0.06, angle = 30, code = 2,
      col = summary_col, lwd = 1.2
    )
  }

  if (mu_ci[2] > alim[2]) {
    graphics::arrows(
      alim[2] - 0.05 * effect_width, y_overall,
      alim[2], y_overall,
      length = 0.06, angle = 30, code = 2,
      col = summary_col, lwd = 1.2
    )
  }

  if (isTRUE(mark_summary_estimate) &&
      is.finite(mu_est) &&
      mu_est >= alim[1] &&
      mu_est <= alim[2]) {
    graphics::segments(
      mu_est,
      y_overall - diamond_h * 0.8,
      mu_est,
      y_overall + diamond_h * 0.8,
      col = estimate_mark_col,
      lwd = estimate_mark_lwd
    )
  }

  if (isTRUE(annotate)) {
    graphics::text(
      x_annot, y_overall,
      fmt_ci(mu_est, mu_ci[1], mu_ci[2]),
      adj = 0, font = 2, cex = cex
    )
  }

  # -----------------------------
  # Prediction interval
  # -----------------------------
  if (isTRUE(show_pi)) {
    graphics::text(
      x_study, y_pi, "Prediction interval",
      adj = 0, cex = cex
    )

    draw_interval(
      theta_ci[1], theta_ci[2], y_pi,
      col = pi_col, lwd = 2, lty = 2, cap = 0.10
    )

    if (isTRUE(symmetric_shapes)) {
      theta_center <- midpoint(theta_ci[1], theta_ci[2])
    } else {
      theta_center <- theta_est
    }

    if (!is.finite(theta_center)) {
      theta_center <- theta_est
    }

    if (is.finite(theta_center) &&
        theta_center >= alim[1] &&
        theta_center <= alim[2]) {
      graphics::points(
        theta_center, y_pi,
        pch = 18, col = pi_col, cex = 1.2
      )
    }

    if (isTRUE(mark_prediction_estimate) &&
        is.finite(theta_est) &&
        theta_est >= alim[1] &&
        theta_est <= alim[2]) {
      graphics::segments(
        theta_est,
        y_pi - 0.16,
        theta_est,
        y_pi + 0.16,
        col = estimate_mark_col,
        lwd = estimate_mark_lwd
      )
    }

    if (isTRUE(annotate)) {
      graphics::text(
        x_annot, y_pi,
        fmt_ci(theta_est, theta_ci[1], theta_ci[2]),
        adj = 0, cex = cex
      )
    }
  }

  # -----------------------------
  # Heterogeneity information
  # -----------------------------
  if (isTRUE(show_het)) {
    het_text <- paste0(
      "Heterogeneity: ",
      "tau^2 = ", format_number(tau2_est, ci_digits),
      " [", format_number(tau2_ci[1], ci_digits), ", ",
      format_number(tau2_ci[2], ci_digits), "]; ",
      "tau = ", format_number(tau_est, ci_digits),
      " [", format_number(tau_ci[1], ci_digits), ", ",
      format_number(tau_ci[2], ci_digits), "]; ",
      "I^2 = ", fmt_pct(I2_est, 1), "%",
      " [", fmt_pct(I2_ci[1], 1), "%, ",
      fmt_pct(I2_ci[2], 1), "%]"
    )

    graphics::text(
      x_study, y_het, het_text,
      adj = 0, cex = 0.9 * cex
    )
  }

  graphics::axis(
    side = 1,
    at = at,
    labels = axis_labels,
    cex.axis = cex
  )
  graphics::mtext(xlab, side = 1, line = 2.6, cex = cex)

  plotted <- data.frame(
    study = slab,
    estimate = y,
    se = se,
    ci.lb = ci_lb,
    ci.ub = ci_ub,
    weight = weights_pct,
    stringsAsFactors = FALSE
  )

  attr(plotted, "overall") <- c(
    estimate = mu_est,
    ci.lb = mu_ci[1],
    ci.ub = mu_ci[2]
  )

  attr(plotted, "prediction") <- c(
    estimate = theta_est,
    pi.lb = theta_ci[1],
    pi.ub = theta_ci[2]
  )

  attr(plotted, "heterogeneity") <- c(
    tau2 = tau2_est,
    tau2.lb = tau2_ci[1],
    tau2.ub = tau2_ci[2],
    tau = tau_est,
    tau.lb = tau_ci[1],
    tau.ub = tau_ci[2],
    I2 = I2_est,
    I2.lb = I2_ci[1],
    I2.ub = I2_ci[2]
  )

  invisible(plotted)
}

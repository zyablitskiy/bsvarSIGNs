
#' @title Plots the median and an interval between two specified percentiles 
#' for a sequence of \code{K} random variables
#'
#' @description Plots the median and an interval between two specified percentiles 
#' for a sequence of \code{K} random variables based on the \code{S} posterior 
#' draws provided for each of them.
#' 
#' @param draws a \code{K x S} matrix with \code{S} posterior draws of 
#' \code{K} random variables, or a \code{K x S x N} array with \code{N} such matrices
#' @param probability a number from interval \code{(0,1)} denoting the probability 
#' content of the plotted interval. The interval stretches from the 
#' \code{0.5 * (1 - probability)} to \code{1 - 0.5 * (1 - probability)} percentile 
#' of the posterior distribution.
#' @param col a colour of the plot
#' @param ylim the range of the \code{y} axis
#' @param ylab the label of the \code{y} axis
#' @param xlab the label of the \code{x} axis
#' @param start_at an integer to denote the beginning of the \code{x} axis range
#' @param add a logical value. If \code{TRUE} the current ribbon plot is added 
#' to an existing plot
#' @param ... other graphical parameters to be passed to \code{base::plot}
#' 
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @examples
#' data(us_fiscal_lsuw)                                               # upload data
#' set.seed(123)                                                      # set seed
#' specification  = specify_bsvar$new(us_fiscal_lsuw)                 # specify model
#' 
#' burn_in        = estimate(specification, 10)                       # run the burn-in
#' posterior      = estimate(burn_in, 20, thin = 1)                   # estimate the model
#' irf            = compute_impulse_responses(posterior, horizon = 4) # impulse responses
#' plot_ribbon(irf[1,1,,])
#' 
#' @export
plot_ribbon = function(
    draws,              # K x S x N
    probability = 0.9,
    col         = "#ff69b4",
    ylim,
    ylab,
    xlab,
    start_at    = 0,
    dates       = NULL,  # New parameter for dates
    add         = FALSE,
    ...
) {
  
  stopifnot("Argument draws must be a matrix or an array." = any(class(draws) == "matrix") || any(class(draws) == "array"))
  stopifnot("Argument probability must be a number from interval (0,1)." = is.numeric(probability) & length(probability) == 1 & probability > 0 & probability < 1)
  stopifnot("Argument start_at must be a matrix or an integer number." = start_at %% 1 == 0)
  stopifnot("Argument add must be a logical value." = is.logical(add) & length(add) == 1)
  
  mat_or_array_tmp  = length(dim(draws))
  
  K               = dim(draws)[1]
  S               = dim(draws)[2]
  if ( mat_or_array_tmp == 2 ) {
    N               = 1
    draws_tmp       = array(NA, c(K, S, 1)) 
    draws_tmp[,,1]  = draws
    draws           = draws_tmp
  } else if ( mat_or_array_tmp == 3 ) {
    N               = dim(draws)[3]
  }
  
  # Compute characteristics of draws
  draws_median    = apply(draws, c(1,3), stats::median) # K x N
  # draws_median    = apply(draws, c(1,3), base::mean) # K x N
  draws_lb        = apply(draws, c(1,3), stats::quantile, probs = 0.5 * (1 - probability)) # K x N  
  draws_ub        = apply(draws, c(1,3), stats::quantile, probs = 1 - 0.5 * (1 - probability)) # K x N
  draws_range     = range(draws_lb, draws_ub)
  
  # Create col_ribbon
  col_ribbon_rgb  = grDevices::col2rgb(col)
  col_ribbon      = grDevices::rgb(col_ribbon_rgb[1], col_ribbon_rgb[2], col_ribbon_rgb[3], 100, maxColorValue = 255)
  
  # Manage the arguments
  if ( missing(ylim) ) ylim = draws_range
  if ( missing(ylab) ) ylab = ""
  if ( missing(xlab) ) xlab = ""
  
  # Create x-axis values
  if (is.null(dates)) {
    x_vals = start_at:(K - 1 + start_at)
  } else {
    # If dates are provided, convert them to numeric for plotting
    # and use them directly
    if (length(dates) != K) {
      warning("Length of dates does not match number of periods. Using numeric indices.")
      x_vals = start_at:(K - 1 + start_at)
    } else {
      x_vals = as.numeric(dates)  # Convert to numeric for plotting
    }
  }
  
  if ( !add ) {
    if (is.null(dates)) {
      base::plot(
        x      = x_vals, 
        y      = draws_median[,1],
        type   = "n",
        ylim   = ylim,
        ylab   = ylab,
        xlab   = xlab, 
        ...
      )
    } else {
      # For date plotting, use the actual dates for x-axis
      base::plot(
        x      = dates, 
        y      = draws_median[,1],
        type   = "n",
        ylim   = ylim,
        ylab   = ylab,
        xlab   = xlab,
        xaxt = "n",
        ...
      )
    }
  }
  
  for (n in 1:N) {
    graphics::polygon(
      x      = c(x_vals, rev(x_vals)),
      y      = c(draws_lb[1:K, n], draws_ub[K:1, n]),
      col    = col_ribbon,
      border = col_ribbon
    )
    graphics::lines(
      x      = x_vals, 
      y      = draws_median[,n],
      type   = "l",
      lwd    = 2,
      col    = col
    )
  }
}







#' @title Plots structural shocks' conditional standard deviations
#'
#' @description Plots of structural shocks' conditional standard deviations 
#' including their median and percentiles.
#' 
#' @param x an object of class PosteriorSigma obtained using the
#' \code{compute_conditional_sd()} function containing posterior draws of 
#' conditional standard deviations of structural shocks.
#' @param probability a parameter determining the interval to be plotted. The 
#' interval stretches from the \code{0.5 * (1 - probability)} to 
#' \code{1 - 0.5 * (1 - probability)} percentile of the posterior distribution.
#' @param shock_names a vector of length \code{N} containing names of the structural shocks.
#' @param col a colour of the plot line and the ribbon
#' @param main an alternative main title for the plot
#' @param xlab an alternative x-axis label for the plot
#' @param mar.multi the default \code{mar} argument setting in \code{graphics::par}. Modify with care!
#' @param oma.multi the default \code{oma} argument setting in \code{graphics::par}. Modify with care!
#' @param ... additional arguments affecting the summary produced.
#' 
#' @method plot PosteriorSigma
#' 
#' @seealso \code{\link{compute_conditional_sd}}
#' 
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @examples
#' data(us_fiscal_lsuw)                                  # upload data
#' set.seed(123)                                         # set seed
#' specification  = specify_bsvar_sv$new(us_fiscal_lsuw) # specify model
#' burn_in        = estimate(specification, 5)           # run the burn-in
#' posterior      = estimate(burn_in, 5)                 # estimate the model
#' 
#' # compute structural shocks' conditional standard deviations
#' sigma          = compute_conditional_sd(posterior)
#' plot(sigma)                                            # plot conditional sds
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' us_fiscal_lsuw |>
#'   specify_bsvar_sv$new(p = 1) |>
#'   estimate(S = 5) |> 
#'   estimate(S = 5) |> 
#'   compute_conditional_sd() |>
#'   plot()
#' 
#' @export
plot.PosteriorSigma = function(
    x,
    probability = 0.9,
    shock_names,
    col = "#ff69b4",
    main,
    xlab,
    mar.multi = c(1, 4.6, 0, 2.1),
    oma.multi = c(6, 0, 5, 0),
    ...
) {

  if ( missing(main) ) main = "Strictural Shocks' Conditional Standard Deviations"
  if ( missing(xlab) ) xlab = "time"
  
  N = dim(x)[1]
  if ( missing(shock_names) ) shock_names = paste("shock", 1:N)
  
  oldpar <- graphics::par( 
    mfrow = c(N, 1),
    mar = mar.multi,
    oma = oma.multi
  )
  on.exit(graphics::par(oldpar))
  
  lw = which(as.numeric(dimnames(x)[[2]]) %% 1 == 0)
  ll = as.numeric(dimnames(x)[[2]])[lw]
  
  for (n in 1:N) {
    
    plot_ribbon(
      x[n,,],
      probability = probability,
      col         = col,
      main = "",
      ylab = shock_names[n],
      xlab = "",
      start_at    = 1,
      bty = "n",
      axes = FALSE,
      ...
    )
    
    graphics::axis(1, at = lw, labels = if (n == N) ll else FALSE)
    graphics::axis(2)
    
    graphics::abline(h = 1)
    
  } # END n loop
  
  graphics::mtext( # main title
    main,
    side = 3,
    line = 2,
    outer = TRUE
  )
  
  graphics::mtext( # x-axis label
    xlab,
    side = 1,
    line = 3,
    outer = TRUE
  )
  
  invisible(x)
} # END plot.PosteriorSigma




#' @title Plots fitted values of dependent variables
#'
#' @description Plots of fitted values of dependent variables including their 
#' median and percentiles.
#' 
#' @param x an object of class PosteriorFitted obtained using the
#' \code{compute_fitted_values()} function containing posterior draws of 
#' fitted values of dependent variables.
#' @param probability a parameter determining the interval to be plotted. The 
#' interval stretches from the \code{0.5 * (1 - probability)} to 
#' \code{1 - 0.5 * (1 - probability)} percentile of the posterior distribution.
#' @param col a colour of the plot line and the ribbon
#' @param main an alternative main title for the plot
#' @param xlab an alternative x-axis label for the plot
#' @param mar.multi the default \code{mar} argument setting in \code{graphics::par}. Modify with care!
#' @param oma.multi the default \code{oma} argument setting in \code{graphics::par}. Modify with care!
#' @param ... additional arguments affecting the summary produced.
#' 
#' @method plot PosteriorFitted
#' 
#' @seealso \code{\link{compute_fitted_values}}
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @examples
#' data(us_fiscal_lsuw)                                  # upload data
#' set.seed(123)                                         # set seed
#' specification  = specify_bsvar$new(us_fiscal_lsuw)    # specify model
#' burn_in        = estimate(specification, 10)          # run the burn-in
#' posterior      = estimate(burn_in, 20, thin = 1)      # estimate the model
#' 
#' # compute fitted values
#' fitted         = compute_fitted_values(posterior)
#' plot(fitted)                                          # plot fitted values
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' us_fiscal_lsuw |>
#'   specify_bsvar$new() |>
#'   estimate(S = 10) |> 
#'   estimate(S = 20, thin = 1) |> 
#'   compute_fitted_values() |>
#'   plot()
#' 
#' @export
plot.PosteriorFitted = function(
    x,
    probability = 0.9,
    col = "#ff69b4",
    main,
    xlab,
    mar.multi = c(1, 4.6, 0, 2.1),
    oma.multi = c(6, 0, 5, 0),
    ...
) {
  
  if ( missing(main) ) main = "Fitted Values of Dependent Variables"
  if ( missing(xlab) ) xlab = "time"
  
  N = dim(x)[1]
  
  oldpar <- graphics::par( 
    mfrow = c(N, 1),
    mar = mar.multi,
    oma = oma.multi
  )
  on.exit(graphics::par(oldpar))
  
  lw = which(as.numeric(dimnames(x)[[2]]) %% 1 == 0)
  ll = as.numeric(dimnames(x)[[2]])[lw]
  var_names = dimnames(x)[[1]]
  
  for (n in 1:N) {
    
    plot_ribbon(
      x[n,,],
      probability = probability,
      col         = col,
      main = "",
      ylab = var_names[n],
      xlab = "",
      start_at    = 1,
      bty = "n",
      axes = FALSE,
      ...
    )
    
    graphics::axis(1, at = lw, labels = if (n == N) ll else FALSE)
    graphics::axis(2)
    
    graphics::abline(h = 0)
    
  } # END n loop
  
  graphics::mtext( # main title
    main,
    side = 3,
    line = 2,
    outer = TRUE
  )
  
  graphics::mtext( # x-axis label
    xlab,
    side = 1,
    line = 3,
    outer = TRUE
  )
  
  invisible(x)
} # END plot.PosteriorFitted













#' @title Plots impulse responses
#'
#' @description Plots of of all variables to all shocks including their 
#' median and percentiles.
#' 
#' @param x an object of class PosteriorIR obtained using the
#' \code{compute_impulse_responses()} function containing posterior draws of 
#' impulse responses.
#' @param probability a parameter determining the interval to be plotted. The 
#' interval stretches from the \code{0.5 * (1 - probability)} to 
#' \code{1 - 0.5 * (1 - probability)} percentile of the posterior distribution.
#' @param shock_names a vector of length \code{N} containing names of the structural shocks.
#' @param col a colour of the plot line and the ribbon
#' @param main an alternative main title for the plot
#' @param xlab an alternative x-axis label for the plot
#' @param mar.multi the default \code{mar} argument setting in \code{graphics::par}. Modify with care!
#' @param oma.multi the default \code{oma} argument setting in \code{graphics::par}. Modify with care!
#' @param ... additional arguments affecting the summary produced.
#' 
#' @method plot PosteriorIR
#' 
#' @seealso \code{\link{compute_impulse_responses}}
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @examples
#' data(us_fiscal_lsuw)                                  # upload data
#' set.seed(123)                                         # set seed
#' specification  = specify_bsvar$new(us_fiscal_lsuw)    # specify model
#' burn_in        = estimate(specification, 10)          # run the burn-in
#' posterior      = estimate(burn_in, 20, thin = 1)      # estimate the model
#' 
#' # compute impulse responses
#' fitted         = compute_impulse_responses(posterior, horizon = 4)
#' plot(fitted)                                          # plot
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' us_fiscal_lsuw |>
#'   specify_bsvar$new() |>
#'   estimate(S = 10) |> 
#'   estimate(S = 20, thin = 1) |> 
#'   compute_impulse_responses(horizon = 4) |>
#'   plot()
#' 
#' @export
plot.PosteriorIR = function(
    x,
    probability = 0.9,
    shock_names,
    var_names,
    col = "#ff69b4",
    main,
    xlab,
    cex.axis = 1.2,
    cex.lab = 1.3,
    cex.main = 1.4,
    cex.shock = 0.9,
    y_label_offset = NULL,
    shock_label_offset = NULL,
    var_subset = NULL,
    shock_subset = NULL,
    ...
) {
  
  if (missing(main)) main = paste0("Impulse Responses (median and ", probability*100, "% CIs)")
  if (missing(xlab)) xlab = "horizon"
  
  # --- SUBSETTING (optional) ---
  if (!is.null(var_subset)) {
    x <- x[var_subset, , , , drop = FALSE]
    if (!missing(var_names)) var_names <- var_names[var_subset]
  }
  if (!is.null(shock_subset)) {
    x <- x[, shock_subset, , , drop = FALSE]
    if (!missing(shock_names)) shock_names <- shock_names[shock_subset]
  }
  
  n_vars   <- dim(x)[1]
  n_shocks <- dim(x)[2]
  
  if (n_vars == 0 || n_shocks == 0) stop("No IRFs to plot after subsetting.")
  
  # Default names
  if (missing(shock_names)) shock_names <- paste("shock", seq_len(n_shocks))
  if (missing(var_names))   var_names   <- paste("var",   seq_len(n_vars))
  
  # Pad or truncate names to match actual dimensions
  if (length(shock_names) < n_shocks) {
    shock_names <- c(shock_names, paste("shock", (length(shock_names)+1):n_shocks))
  } else if (length(shock_names) > n_shocks) {
    shock_names <- shock_names[seq_len(n_shocks)]
  }
  if (length(var_names) < n_vars) {
    var_names <- c(var_names, paste("var", (length(var_names)+1):n_vars))
  } else if (length(var_names) > n_vars) {
    var_names <- var_names[seq_len(n_vars)]
  }
  
  H = dim(x)[3] - 1
  
  # --- LAYOUT SETUP ---
  n_panels <- n_vars * n_shocks
  mat <- matrix(seq_len(n_panels), nrow = n_vars, byrow = TRUE)
  
  # First column slightly wider for y-labels
  col_widths <- rep(1, n_shocks)
  if (n_shocks > 0) col_widths[1] <- 1.3
  
  # Scale text down for large grids (use the larger of the two dimensions)
  max_dim <- max(n_vars, n_shocks)
  if (max_dim > 4) {
    scale <- max(0.5, 1 - (max_dim - 4) * 0.08)
    cex.axis  <- cex.axis * scale
    cex.lab   <- cex.lab * scale
    cex.shock <- cex.shock * scale
    cex.main  <- cex.main * scale
  }
  
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  
  graphics::par(omi = c(0.4, 0.1, 0.5, 0.1))
  graphics::layout(mat, widths = col_widths, heights = rep(1, n_vars))
  
  graphics::par(
    mai = c(0.35, 0.5, 0.2, 0.08),
    mgp = c(2.5, 0.6, 0),
    tcl = -0.25
  )
  
  var_line   <- if (is.null(y_label_offset)) 2.2 else y_label_offset
  shock_line <- if (is.null(shock_label_offset)) 0.3 else shock_label_offset
  
  for (n in seq_len(n_vars)) {
    for (i in seq_len(n_shocks)) {
      
      plot_ribbon(
        x[n, i, , ], probability = probability, col = col,
        main = "", 
        ylab = "",          
        xlab = "",          
        start_at = 0, bty = "n", axes = FALSE,
        cex.lab = cex.lab
      )
      
      graphics::axis(1, cex.axis = cex.axis)
      graphics::axis(2, cex.axis = cex.axis)
      graphics::abline(h = 0)
      
      # First column: variable names as y-labels
      if (i == 1) {
        graphics::mtext(var_names[n], side = 2, line = var_line,
                        cex = cex.lab, font = 1)
      }
      
      # Top row: shock names
      if (n == 1) {
        graphics::mtext(shock_names[i], side = 3, line = shock_line,
                        cex = cex.shock, font = 2)
      }
    }
  }
  
  graphics::mtext(main, side = 3, line = 1, outer = TRUE,
                  cex = cex.main, font = 2)
  
  if (!is.null(xlab)) {
    graphics::mtext(xlab, side = 1, line = 0.2, outer = TRUE,
                    cex = cex.lab * 0.9, font = 1)
  }
  
  invisible(x)
} # END plot.PosteriorIR




#' @title Plots structural shocks
#'
#' @description Plots of structural shocks including their median and percentiles.
#' 
#' @param x an object of class PosteriorShocks obtained using the
#' \code{compute_structural_shocks()} function containing posterior draws of 
#' structural shocks.
#' @param probability a parameter determining the interval to be plotted. The 
#' interval stretches from the \code{0.5 * (1 - probability)} to 
#' \code{1 - 0.5 * (1 - probability)} percentile of the posterior distribution.
#' @param col a colour of the plot line and the ribbon
#' @param main an alternative main title for the plot
#' @param xlab an alternative x-axis label for the plot
#' @param mar.multi the default \code{mar} argument setting in \code{graphics::par}. Modify with care!
#' @param oma.multi the default \code{oma} argument setting in \code{graphics::par}. Modify with care!
#' @param ... additional arguments affecting the summary produced.
#' 
#' @method plot PosteriorShocks
#' 
#' @seealso \code{\link{compute_structural_shocks}}
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @examples
#' data(us_fiscal_lsuw)                                  # upload data
#' set.seed(123)                                         # set seed
#' specification  = specify_bsvar$new(us_fiscal_lsuw)    # specify model
#' burn_in        = estimate(specification, 10)          # run the burn-in
#' posterior      = estimate(burn_in, 20, thin = 1)      # estimate the model
#' 
#' # compute structural shocks
#' shocks         = compute_structural_shocks(posterior)
#' plot(shocks)                                          # plot
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' us_fiscal_lsuw |>
#'   specify_bsvar$new() |>
#'   estimate(S = 10) |> 
#'   estimate(S = 20, thin = 1) |> 
#'   compute_structural_shocks() |>
#'   plot()
#' 
#' @export
plot.PosteriorShocks = function(
    x,
    spec,
    which_shocks = NULL,
    ncol = 1,
    probability = 0.9,
    col = "#ff69b4",
    main,
    shock_names,
    l = 6,
    xlab,
    start_date = NULL,
    end_date = NULL,
    cex.axis = 1.2,
    cex.lab = 1.4,
    cex.main = 1.6,
    legend.space = 0.3,
    ...
) {
  
  if ( missing(main) ) main = paste0("Structural Shocks (median and ", probability*100, "% CIs)")
  if ( missing(xlab) ) xlab = NULL
  
  N_all = dim(x)[1]
  K = dim(x)[2]
  
  if ( missing(shock_names) ) shock_names = paste("shock", 1:N_all)
  
  # --- which_shocks handling ---
  if (is.null(which_shocks)) {
    which_shocks <- 1:N_all
  } else if (is.logical(which_shocks)) {
    if (length(which_shocks) != N_all) stop("Logical which_shocks must have length N_all (", N_all, ").")
    which_shocks <- which(which_shocks)
  } else if (is.character(which_shocks)) {
    which_shocks <- match(which_shocks, shock_names)
    if (any(is.na(which_shocks))) stop("Some names in which_shocks not found in shock_names.")
  } else if (is.numeric(which_shocks)) {
    which_shocks <- as.integer(which_shocks)
    if (any(which_shocks < 1 | which_shocks > N_all)) stop("which_shocks indices out of bounds.")
  } else {
    stop("which_shocks must be NULL, logical, character, or integer.")
  }
  
  k <- length(which_shocks)
  if (ncol < 1) stop("ncol must be at least 1.")
  
  # Subset to selected shocks
  x <- x[which_shocks, , , drop = FALSE]
  shock_names <- shock_names[which_shocks]
  N <- k
  
  dates <- lubridate::floor_date(
    as.Date(lubridate::date_decimal(
      as.numeric(as.numeric(colnames(spec$data_matrices$Y))) + 0.01
    )), 
    unit = "months"
  )
  
  if (length(dates) != K) {
    warning("Length of dates does not match number of periods. Using numeric indices.")
    dates <- 1:K
  }
  
  # === DATE FILTER ===
  if (!is.null(start_date) && !inherits(start_date, "Date")) {
    start_date <- as.Date(start_date)
  }
  if (!is.null(end_date) && !inherits(end_date, "Date")) {
    end_date <- as.Date(end_date)
  }
  
  idx <- rep(TRUE, length(dates))
  if (!is.null(start_date)) idx <- idx & (dates >= start_date)
  if (!is.null(end_date))   idx <- idx & (dates <= end_date)
  
  dates <- dates[idx]
  K     <- length(dates)
  
  x <- x[,idx,, drop = FALSE]
  
  # --- multi-column layout ---
  nrow_plots <- ceiling(N / ncol)
  
  layout_matrix <- matrix(0, nrow = nrow_plots + 1, ncol = ncol)
  for (i in 1:N) {
    r <- ceiling(i / ncol)
    c <- (i - 1) %% ncol + 1
    layout_matrix[r, c] <- i
  }
  layout_matrix[nrow_plots + 1, ] <- N + 1   # bottom panel for outer title/xlab
  
  heights <- c(rep(1, nrow_plots), legend.space)
  graphics::layout(layout_matrix, heights = heights)
  
  oldpar <- graphics::par(mar = c(0, 5.5, 0, 2.1), oma = c(2, 0, 4, 0))
  on.exit({ graphics::par(oldpar); graphics::layout(1) })
  
  if (K <= 20) {
    tick_indices <- 1:K
  } else {
    tick_indices <- round(seq(1, K, length.out = l))
  }
  tick_labels <- format(dates[tick_indices], "%b %Y")
  
  # Bottom-most panel in EACH column gets the x-axis
  bottom_panels <- apply(layout_matrix[1:nrow_plots, , drop = FALSE], 2, function(col) {
    idx <- which(col > 0)
    if (length(idx) == 0) return(NA_integer_)
    col[max(idx)]
  })
  bottom_panels <- na.omit(bottom_panels)
  
  for (n in 1:N) {
    is_bottom <- n %in% bottom_panels
    
    graphics::par(mar = c(
      if (is_bottom) 3.5 else 0,
      5.5,
      0, 2.1
    ))
    
    plot_ribbon(
      x[n,,], probability = probability, col = col,
      main = "", ylab = "", xlab = "", start_at = 1,
      bty = "n", axes = FALSE, ...
    )
    
    graphics::mtext(shock_names[n], side = 2, line = 3.5,
                    cex = cex.lab, font = 1)
    
    if (is_bottom) {
      graphics::axis(1, at = tick_indices, labels = tick_labels,
                     las = 2, cex.axis = cex.axis)
    }
    graphics::axis(2, cex.axis = cex.axis)
    graphics::abline(h = 0)
  }
  
  # Bottom panel (spans all columns)
  graphics::par(mar = c(0, 0, 0, 0)); graphics::plot.new()
  
  graphics::mtext(main, side = 3, line = 1.5, outer = TRUE,
                  cex = cex.main, font = 2)
  
  if (!is.null(xlab)) {
    graphics::mtext(xlab, side = 1, line = 1, outer = TRUE,
                    cex = cex.lab * 0.9, font = 1)
  }
  
  invisible(x)
} # END plot.PosteriorShocks







#' @title Plots fitted values of dependent variables
#'
#' @description Plots of fitted values of dependent variables including their 
#' median and percentiles.
#' 
#' @param x an object of class Forecasts obtained using the
#' \code{forecast()} function containing posterior draws of 
#' fitted values of dependent variables.
#' @param probability a parameter determining the interval to be plotted. The 
#' interval stretches from the \code{0.5 * (1 - probability)} to 
#' \code{1 - 0.5 * (1 - probability)} percentile of the posterior distribution.
#' @param data_in_plot a fraction value in the range (0, 1) determining how many
#' of the last observations in the data should be plotted with the forecasts.
#' @param col a colour of the plot line and the ribbon
#' @param main an alternative main title for the plot
#' @param xlab an alternative x-axis label for the plot
#' @param mar.multi the default \code{mar} argument setting in \code{graphics::par}. Modify with care!
#' @param oma.multi the default \code{oma} argument setting in \code{graphics::par}. Modify with care!
#' @param ... additional arguments affecting the summary produced.
#' 
#' @method plot Forecasts
#' 
#' @seealso \code{\link{forecast}}
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @examples
#' data(us_fiscal_lsuw)                                  # upload data
#' set.seed(123)                                         # set seed
#' specification  = specify_bsvar$new(us_fiscal_lsuw)    # specify model
#' burn_in        = estimate(specification, 10)          # run the burn-in
#' posterior      = estimate(burn_in, 20, thin = 1)      # estimate the model
#' 
#' # compute forecasts
#' fore            = forecast(posterior, horizon = 4)
#' plot(fore)                                            # plot forecasts
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' us_fiscal_lsuw |>
#'   specify_bsvar$new() |>
#'   estimate(S = 10) |> 
#'   estimate(S = 20, thin = 1) |> 
#'   forecast(horizon = 4) |>
#'   plot()
#' 
#' @export
plot.Forecasts = function(
    x,
    probability = 0.9,
    data_in_plot = 1,
    col = "#ff69b4",
    main,
    xlab,
    mar.multi = c(1, 4.6, 0, 2.1),
    oma.multi = c(6, 0, 5, 0),
    ...
) {
  
  stopifnot("Argument data_in_plot must be a value within range (0,1]" = is.numeric(data_in_plot) & data_in_plot > 0 & data_in_plot <= 1)
  
  if ( missing(main) ) main = "Forecasting"
  if ( missing(xlab) ) xlab = "time"
  
  # create col_ribbon
  col_ribbon_rgb  = grDevices::col2rgb(col)
  col_ribbon      = grDevices::rgb(col_ribbon_rgb[1], col_ribbon_rgb[2], col_ribbon_rgb[3], 100, maxColorValue = 255)
  
  fore          = x$forecasts
  Y             = x$Y
  
  N             = dim(fore)[1]
  H             = dim(fore)[2]
  T             = dim(Y)[2]
  
  T_in_plot     = floor(data_in_plot * T)
  if (T_in_plot < 1) T_in_plot = 1
  obs_in_plot   = (T - T_in_plot + 1):T
  seq_in_plot   = 1:T_in_plot
  for_in_plot   = (T_in_plot + 1):(T_in_plot + H)
  
  oldpar <- graphics::par( 
    mfrow = c(N, 1),
    mar = mar.multi,
    oma = oma.multi
  )
  on.exit(graphics::par(oldpar))
  
  for (n in 1:N) {
    
    # compute forecasts characteristics
    fore_median   = apply(fore[n,,], 1, stats::median)
    fore_lb       = apply(fore[n,,], 1, stats::quantile, probs = 0.5 * (1 - probability)) # K x N  
    fore_ub       = apply(fore[n,,], 1, stats::quantile, probs = 1 - 0.5 * (1 - probability)) # K x N
    fore_range    = range(fore_lb, fore_ub, Y[n, obs_in_plot])
    
    base::plot(
      x      = c(seq_in_plot, for_in_plot),
      y      = c(Y[n, seq_in_plot], fore_median),
      type   = "n",
      ylim   = fore_range,
      main   = "",
      ylab   = paste("variable", n),
      xlab   = "",
      bty = "n",
      axes = FALSE,
      ...
    )

    graphics::polygon(
      x      = c(c(T_in_plot, for_in_plot), rev(c(T_in_plot, for_in_plot))),
      y      = c(c(Y[n, T], fore_lb), rev(c(Y[n, T], fore_ub))),
      col    = col_ribbon,
      border = col_ribbon
    )
    
    graphics::lines(
      x      = c(seq_in_plot, for_in_plot),
      y      = c(Y[n, obs_in_plot], fore_median),
      type   = "l",
      lwd    = 2,
      col    = col
    )
    
    graphics::axis(1, labels = if (n == N) TRUE else FALSE)
    graphics::axis(2)
    
    graphics::abline(h = 0)
    
  } # END n loop
  
  graphics::mtext( # main title
    main,
    side = 3,
    line = 2,
    outer = TRUE
  )
  
  graphics::mtext( # x-axis label
    xlab,
    side = 1,
    line = 3,
    outer = TRUE
  )
  
  invisible(x)
} # END plot.Forecasts








#' @title Plots forecast error variance decompositions
#'
#' @description Plots of the posterior means of the forecast error variance 
#' decompositions.
#' 
#' @param x an object of class PosteriorFEVD obtained using the
#' \code{compute_variance_decompositions()} function containing posterior draws of 
#' forecast error variance decompositions.
#' @param shock_names a vector of length \code{N} containing names of the structural shocks.
#' @param cols an \code{N}-vector with colours of the plot
#' @param main an alternative main title for the plot
#' @param xlab an alternative x-axis label for the plot
#' @param mar.multi the default \code{mar} argument setting in \code{graphics::par}. Modify with care!
#' @param oma.multi the default \code{oma} argument setting in \code{graphics::par}. Modify with care!
#' @param ... additional arguments affecting the summary produced.
#' 
#' @method plot PosteriorFEVD
#' 
#' @seealso \code{\link{compute_variance_decompositions}}
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @examples
#' data(us_fiscal_lsuw)                                  # upload data
#' set.seed(123)                                         # set seed
#' specification  = specify_bsvar$new(us_fiscal_lsuw)    # specify model
#' burn_in        = estimate(specification, 10)          # run the burn-in
#' posterior      = estimate(burn_in, 20, thin = 1)      # estimate the model
#' 
#' # compute forecast error variance decompositions
#' fevd           = compute_variance_decompositions(posterior, horizon = 4)
#' plot(fevd)
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' us_fiscal_lsuw |>
#'   specify_bsvar$new() |>
#'   estimate(S = 10) |> 
#'   estimate(S = 20, thin = 1) |> 
#'   compute_variance_decompositions(horizon = 4) |>
#'   plot()
#' 
#' @export
plot.PosteriorFEVD = function(
    x,
    shock_names,
    var_names,
    which_vars = NULL,
    ncol = 1,
    probability = NULL,
    cols,
    main,
    xlab,
    cex.axis = 1.2,
    cex.lab = 1.3,
    cex.main = 1.4,
    legend.ncol = NULL,
    cex.legend = 1.2,
    legend.space = 0.3,
    inset.space.y = 0.1,
    inset.space.x = 0.5,
    ...
) {
  if(!is.null(probability)) stopifnot("Argument probability must be a number from interval (0,1)." = is.numeric(probability) & length(probability) == 1 & probability > 0 & probability < 1)
  
  if ( missing(main) ) {
    if(!is.null(probability)) {
      main = paste0("Forecast Error Variance Decompositions (", probability*100, "% quantile)")
    } else {
      main = "Forecast Error Variance Decompositions (mean)"
    }
  }
  if ( missing(xlab) ) xlab = "horizon"
  
  N_all = dim(x)[1]
  if ( missing(shock_names) ) shock_names = paste("shock", 1:N_all)
  if ( missing(var_names) )   var_names   = paste("var", 1:N_all)
  
  # --- which_vars handling ---
  if (is.null(which_vars)) {
    which_vars <- 1:N_all
  } else if (is.logical(which_vars)) {
    if (length(which_vars) != N_all) stop("Logical which_vars must have length N_all (", N_all, ").")
    which_vars <- which(which_vars)
  } else if (is.character(which_vars)) {
    which_vars <- match(which_vars, var_names)
    if (any(is.na(which_vars))) stop("Some names in which_vars not found in var_names.")
  } else if (is.numeric(which_vars)) {
    which_vars <- as.integer(which_vars)
    if (any(which_vars < 1 | which_vars > N_all)) stop("which_vars indices out of bounds.")
  } else {
    stop("which_vars must be NULL, logical, character, or integer.")
  }
  
  k <- length(which_vars)
  if (ncol < 1) stop("ncol must be at least 1.")
  
  H = dim(x)[3] - 1
  
  if ( missing(cols) ) {
    fc   = grDevices::colorRampPalette(c("#ff69b4", "#ffd700"))
    cols = fc(N_all)
  }
  
  if (!is.null(probability)) {
    fevd = apply(x, 1:3, stats::quantile, probs = probability)
  } else {
    fevd = apply(x, 1:3, mean)
  }
  
  FEVD = list()
  for (n in 1:N_all) {
    FEVD[[n]] = rbind(rep(0, H + 1), apply(fevd[n,,], 2, cumsum))
  }
  
  FEVD_sub      <- FEVD[which_vars]
  var_names_sub <- var_names[which_vars]
  
  # --- multi-column layout ---
  nrow_plots <- ceiling(k / ncol)
  layout_matrix <- matrix(0, nrow = nrow_plots + 1, ncol = ncol)
  for (i in 1:k) {
    r <- ceiling(i / ncol)
    c <- (i - 1) %% ncol + 1
    layout_matrix[r, c] <- i
  }
  layout_matrix[nrow_plots + 1, ] <- k + 1   # legend panel
  
  heights <- c(rep(1, nrow_plots), legend.space)
  graphics::layout(layout_matrix, heights = heights)
  
  oldpar <- graphics::par(mar = c(0, 5.5, 0, 2.1), oma = c(2, 0, 4, 0))
  on.exit({ graphics::par(oldpar); graphics::layout(1) })
  
  # Bottom-most panel in EACH column gets the x-axis
  bottom_panels <- apply(layout_matrix[1:nrow_plots, , drop = FALSE], 2, function(col) {
    idx <- which(col > 0)
    if (length(idx) == 0) return(NA_integer_)
    col[max(idx)]
  })
  bottom_panels <- na.omit(bottom_panels)
  
  for (ii in 1:k) {
    is_bottom <- ii %in% bottom_panels
    
    graphics::par(mar = c(
      if (is_bottom) 3.5 else 0,
      5.5,
      0, 2.1
    ))
    
    graphics::plot(
      x = 0:H, y = FEVD_sub[[ii]][1,], type = "n",
      ylim = c(0, 100), 
      ylab = var_names_sub[ii],
      main = "", xlab = "", bty = "n", axes = FALSE,
      cex.lab = cex.lab, ...
    )
    
    for (i in 1:N_all) {
      graphics::polygon(
        c(0:H, H:0), 
        c(FEVD_sub[[ii]][i,], rev(FEVD_sub[[ii]][i + 1,])), 
        col = cols[i], border = cols[i]
      )
    }
    
    if (is_bottom) graphics::axis(1, cex.axis = cex.axis)
    graphics::axis(2, c(0, 50, 100), c(0, 50, 100), cex.axis = cex.axis)
  }
  
  # --- legend ---
  graphics::par(mar = c(0, 0, 0, 0)); graphics::plot.new()
  
  if (is.null(legend.ncol)) {
    legend.ncol <- if (N_all <= 4) N_all else min(4, ceiling(N_all/2))
  }
  
  leg_dim <- legend(
    x = 0.5, y = 0.5, plot = FALSE,
    legend = shock_names, fill = cols, border = cols,
    bty = "n", cex = cex.legend, ncol = legend.ncol, horiz = FALSE,
    title = "Shocks", title.adj = 0.5, text.font = 1,
    title.cex = cex.legend * 1.1
  )
  
  legend(
    x = inset.space.x - leg_dim$rect$w / 2,
    y = inset.space.y,
    xjust = 0, yjust = 0,
    legend = shock_names, fill = cols, border = cols,
    bty = "n", cex = cex.legend, ncol = legend.ncol, horiz = FALSE,
    title = "Shocks", title.adj = 0.5, text.font = 1,
    title.cex = cex.legend * 1.1
  )
  
  graphics::mtext(main, side = 3, line = 1.5, outer = TRUE, cex = cex.main, font = 2)
  if (!is.null(xlab)) {
    graphics::mtext(xlab, side = 1, line = 1, outer = TRUE, cex = cex.lab * 0.9, font = 1)
  }
  
  invisible(x)
} # END plot.PosteriorFEVD






#' @title Plots historical decompositions
#'
#' @description Plots of the posterior means of the historical decompositions.
#' 
#' @param x an object of class PosteriorHD obtained using the
#' \code{compute_historical_decompositions()} function containing posterior draws of 
#' historical decompositions.
#' @param shock_names a vector of length \code{N} containing names of the structural shocks.
#' @param cols an \code{N}-vector with colours of the plot
#' @param main an alternative main title for the plot
#' @param xlab an alternative x-axis label for the plot
#' @param mar.multi the default \code{mar} argument setting in \code{graphics::par}. Modify with care!
#' @param oma.multi the default \code{oma} argument setting in \code{graphics::par}. Modify with care!
#' @param ... additional arguments affecting the summary produced.
#' 
#' @method plot PosteriorHD
#' 
#' @seealso \code{\link{compute_historical_decompositions}}
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @examples
#' data(us_fiscal_lsuw)                                  # upload data
#' set.seed(123)                                         # set seed
#' specification  = specify_bsvar$new(us_fiscal_lsuw)    # specify model
#' burn_in        = estimate(specification, 10)          # run the burn-in
#' posterior      = estimate(burn_in, 20, thin = 1)      # estimate the model
#' 
#' # compute historical decompositions
#' fevd           = compute_historical_decompositions(posterior)
#' plot(fevd)                                            
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' us_fiscal_lsuw |>
#'   specify_bsvar$new() |>
#'   estimate(S = 10) |> 
#'   estimate(S = 20, thin = 1) |> 
#'   compute_historical_decompositions() |>
#'   plot()
#' 
#' @export

plot.PosteriorHD = function(
    x,
    spec,
    shock_names,
    var_names,
    which_vars = NULL,
    ncol = 1,
    probability = NULL,
    l = 6,
    cols,
    main,
    xlab,
    start_date = NULL,
    end_date = NULL,
    var_end_dates = NULL,
    include_deterministic = TRUE,
    det_name = "Deterministic",
    dummy_name = "Outlier dummy",
    det_col = "grey60",
    dummy_col = "grey20",
    legend.cex = 0.8,
    legend.ncol = NULL,
    cex.axis = 1.2,
    cex.lab = 1.3,
    cex.main = 1.4,
    cex.legend = 1.2,
    legend.space = 0.3,
    inset.space.y = 0.1,
    inset.space.x = 0.5,
    ...
) {
  
  if(!is.null(probability)) stopifnot("Argument probability must be a number from interval (0,1)." = is.numeric(probability) & length(probability) == 1 & probability > 0 & probability < 1)
  
  # --- pull custom attributes BEFORE any `[` subsetting, which drops them ---
  hd_det   = if (include_deterministic) attr(x, "deterministic") else NULL
  hd_dummy = if (include_deterministic) attr(x, "dummy") else NULL
  
  if ( missing(main) ) {
    if(!is.null(probability)) {
      main = paste0("Historical Decompositions (", probability*100, "% quantile)")
    } else {
      main = "Historical Decompositions (mean)"
    }
  }
  if ( missing(xlab) ) xlab = NULL
  
  dates <- lubridate::floor_date(
    as.Date(lubridate::date_decimal(
      as.numeric(as.numeric(colnames(spec$data_matrices$Y))) + 0.01
    )), 
    unit = "months"
  )
  
  # === DATE FILTER ===
  if (!is.null(start_date) && !inherits(start_date, "Date")) start_date <- as.Date(start_date)
  if (!is.null(end_date)   && !inherits(end_date, "Date"))   end_date   <- as.Date(end_date)
  
  idx <- rep(TRUE, length(dates))
  if (!is.null(start_date)) idx <- idx & (dates >= start_date)
  if (!is.null(end_date))   idx <- idx & (dates <= end_date)
  
  dates <- dates[idx]
  T <- length(dates)
  
  x <- x[,,idx,, drop = FALSE]
  if (!is.null(hd_det))   hd_det   = hd_det[, idx, , drop = FALSE]
  if (!is.null(hd_dummy)) hd_dummy = hd_dummy[, idx, , drop = FALSE]
  
  if (!is.null(probability)) {
    hd_mean = apply(x, 1:3, stats::quantile, probs = probability)
  } else {
    hd_mean = apply(x, 1:3, mean)
  }
  
  N_all = dim(hd_mean)[1]  # original number of response variables
  
  if ( missing(var_names) ) var_names = paste("var", 1:N_all)
  
  # === which_vars RESOLUTION (operates on ORIGINAL indices/names) ===
  if (is.null(which_vars)) {
    which_vars <- 1:N_all
  } else if (is.logical(which_vars)) {
    if (length(which_vars) != N_all) stop("Logical which_vars must have length N_all (", N_all, ").")
    which_vars <- which(which_vars)
  } else if (is.character(which_vars)) {
    which_vars <- match(which_vars, var_names)
    if (any(is.na(which_vars))) stop("Some names in which_vars not found in var_names.")
  } else if (is.numeric(which_vars)) {
    which_vars <- as.integer(which_vars)
    if (any(which_vars < 1 | which_vars > N_all)) stop("which_vars indices out of bounds.")
  } else {
    stop("which_vars must be NULL, logical, character, or integer.")
  }
  
  k <- length(which_vars)
  if (ncol < 1) stop("ncol must be at least 1.")
  
  # === CRITICAL: Subset hd_mean, hd_det, hd_dummy TOGETHER before any further processing ===
  hd_mean <- hd_mean[which_vars, , , drop = FALSE]
  if (!is.null(hd_det))   hd_det   <- hd_det[which_vars, , , drop = FALSE]
  if (!is.null(hd_dummy)) hd_dummy <- hd_dummy[which_vars, , , drop = FALSE]
  
  # From here on, N_resp refers to the SELECTED subset only
  N_resp <- k
  var_names <- var_names[which_vars]
  
  # --- helper: concatenate along the shock (2nd) dimension ---
  bind_shock = function(a, b) {
    da = dim(a); db = dim(b)
    out = array(0, dim = c(da[1], da[2] + db[2], da[3]))
    out[, 1:da[2], ] = a
    out[, (da[2] + 1):(da[2] + db[2]), ] = b
    out
  }
  
  extra_names = character(0)
  extra_cols  = character(0)
  
  if (!is.null(hd_det)) {
    if (!is.null(hd_dummy)) {
      baseline = hd_det - hd_dummy  # exclude the contemporaneous dummy piece
      if (!is.null(probability)) {
        baseline_mean = apply(baseline, 1:2, stats::quantile, probs = probability)
        dummy_mean    = apply(hd_dummy, 1:2, stats::quantile, probs = probability)
      } else {
        baseline_mean = apply(baseline, 1:2, mean)
        dummy_mean    = apply(hd_dummy, 1:2, mean)
      }
      hd_mean = bind_shock(hd_mean, array(baseline_mean, dim = c(N_resp, 1, T)))
      hd_mean = bind_shock(hd_mean, array(dummy_mean,    dim = c(N_resp, 1, T)))
      extra_names = c(det_name, dummy_name)
      extra_cols  = c(det_col, dummy_col)
    } else {
      if (!is.null(probability)) {
        det_mean = apply(hd_det, 1:2, stats::quantile, probs = probability)
      } else {
        det_mean = apply(hd_det, 1:2, mean)
      }
      hd_mean = bind_shock(hd_mean, array(det_mean, dim = c(N_resp, 1, T)))
      extra_names = det_name
      extra_cols  = det_col
    }
  }
  
  N_shock = dim(hd_mean)[2]  # structural shocks + deterministic/dummy
  
  if ( missing(shock_names) ) shock_names = paste("shock", 1:(N_shock - length(extra_names)))
  shock_names = c(shock_names, extra_names)
  
  # === PER-VARIABLE END DATES (now operates on the SELECTED subset) ===
  if (!is.null(var_end_dates)) {
    if (length(var_end_dates) != N_resp) {
      stop("Argument var_end_dates must be a vector of length N_resp (", N_resp, ").")
    }
    var_end_dates <- as.Date(var_end_dates)
    if (!is.null(names(var_end_dates))) {
      ord <- match(var_names, names(var_end_dates))
      if (any(is.na(ord))) stop("Names of var_end_dates must match var_names.")
      var_end_dates <- var_end_dates[ord]
    }
    for (n in 1:N_resp) {
      if (!is.na(var_end_dates[n])) {
        zero_idx <- which(dates > var_end_dates[n])
        if (length(zero_idx) > 0) {
          hd_mean[n, , zero_idx] <- 0
        }
      }
    }
  }
  
  if ( missing(cols) ) {
    fc   = grDevices::colorRampPalette(c("#ff69b4", "#ffd700"))
    cols = c(fc(N_shock - length(extra_names)), extra_cols)
  } else {
    stopifnot("cols must have length N_shock (structural + deterministic/dummy)." =
                length(cols) == N_shock)
  }
  
  y_hat = apply(hd_mean, c(1, 3), sum)
  
  # === STACKING (positive/negative contributions) ===
  lower = array(0, dim = c(N_resp, N_shock, T))
  upper = array(0, dim = c(N_resp, N_shock, T))
  
  for (n in 1:N_resp) {
    for (t in 1:T) {
      pos_cum = 0
      neg_cum = 0
      for (i in 1:N_shock) {
        contrib = hd_mean[n, i, t]
        if (contrib > 0) {
          lower[n, i, t] = pos_cum
          upper[n, i, t] = pos_cum + contrib
          pos_cum = pos_cum + contrib
        } else if (contrib < 0) {
          upper[n, i, t] = neg_cum
          lower[n, i, t] = neg_cum + contrib
          neg_cum = neg_cum + contrib
        }
      }
    }
  }
  
  # === MULTI-COLUMN LAYOUT ===
  nrow_plots <- ceiling(k / ncol)
  layout_matrix <- matrix(0, nrow = nrow_plots + 1, ncol = ncol)
  for (i in 1:k) {
    r <- ceiling(i / ncol)
    c <- (i - 1) %% ncol + 1
    layout_matrix[r, c] <- i
  }
  layout_matrix[nrow_plots + 1, ] <- k + 1  # legend panel
  
  heights <- c(rep(1, nrow_plots), legend.space)
  graphics::layout(layout_matrix, heights = heights)
  
  oldpar <- graphics::par(mar = c(0, 5.5, 0, 2.1), oma = c(2, 0, 4, 0))
  on.exit({ graphics::par(oldpar); graphics::layout(1) })
  
  if (T <= 20) {
    tick_indices <- 1:T
  } else {
    tick_indices <- round(seq(1, T, length.out = l))
  }
  tick_labels <- format(dates[tick_indices], "%b %Y")
  
  # Bottom-most panel in EACH column gets the x-axis
  bottom_panels <- apply(layout_matrix[1:nrow_plots, , drop = FALSE], 2, function(col) {
    idx <- which(col > 0)
    if (length(idx) == 0) return(NA_integer_)
    col[max(idx)]
  })
  bottom_panels <- na.omit(bottom_panels)
  
  for (ii in 1:k) {
    is_bottom <- ii %in% bottom_panels
    
    all_vals = c(0, lower[ii, , ], upper[ii, , ], y_hat[ii, ])
    range_n = range(all_vals, na.rm = TRUE)
    
    graphics::par(mar = c(
      if (is_bottom) 3.5 else 0,
      5.5, 0, 2.1
    ))
    
    graphics::plot(
      x = 1:T, y = y_hat[ii, ], type = "n",
      ylim = range_n, 
      ylab = var_names[ii],
      main = "", xlab = "", bty = "n", axes = FALSE,
      cex.lab = cex.lab, ...
    )
    
    for (i in 1:N_shock) {
      has_contrib <- (lower[ii, i, ] != upper[ii, i, ])
      if (any(has_contrib)) {
        t_idx <- which(has_contrib)
        graphics::rect(
          xleft   = t_idx - 0.5,
          ybottom = lower[ii, i, t_idx],
          xright  = t_idx + 0.5,
          ytop    = upper[ii, i, t_idx],
          col = cols[i], border = cols[i]
        )
      }
    }
    
    graphics::abline(h = 0, lty = 1, col = "black")
    
    if (is_bottom) {
      graphics::axis(1, at = tick_indices, labels = tick_labels,
                     las = 2, cex.axis = cex.axis)
    }
    graphics::axis(2, c(range_n[1], 0, range_n[2]), c(NA, 0, NA),
                   cex.axis = cex.axis)
  }
  
  # --- legend ---
  graphics::par(mar = c(0, 0, 0, 0)); graphics::plot.new()
  
  if (is.null(legend.ncol)) {
    legend.ncol <- if (N_shock <= 4) N_shock else min(4, ceiling(N_shock/2))
  }
  
  leg_dim <- legend(
    x = 0.5, y = 0.5, plot = FALSE,
    legend = shock_names, fill = cols, border = cols,
    bty = "n", cex = cex.legend, ncol = legend.ncol, horiz = FALSE,
    title = "Shocks", title.adj = 0.5, text.font = 1,
    title.cex = cex.legend * 1.1
  )
  
  legend(
    x = inset.space.x - leg_dim$rect$w / 2,
    y = inset.space.y,
    xjust = 0, yjust = 0,
    legend = shock_names, fill = cols, border = cols,
    bty = "n", cex = cex.legend, ncol = legend.ncol, horiz = FALSE,
    title = "Shocks", title.adj = 0.5, text.font = 1,
    title.cex = cex.legend * 1.1
  )
  
  graphics::mtext(main, side = 3, line = 1.5, outer = TRUE,
                  cex = cex.main, font = 2)
  if (!is.null(xlab)) {
    graphics::mtext(xlab, side = 1, line = 1, outer = TRUE,
                    cex = cex.lab * 0.9, font = 1)
  }
  
  invisible(x)
  
} # END plot.PosteriorHD

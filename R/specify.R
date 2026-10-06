
#' vector specifying one narrative restriction
#'
#' @description
#' The class narrative specifies a single narrative restriction.
#' 
#' @param start positive integer - the period in which the narrative starts (greater than the number of lags).
#' 
#' @param periods positive integer - the number of periods the narrative restriction lasts.
#' 
#' @param type character - the type of the narrative restriction (one of "S", "A", "B"), where:
#' "S" for restrictions on structural shocks;
#' "A" for type A restrictions on historical decomposition, i.e. if the absolute value of the historical decomposition
#' of shock to var is the greatest/least among all shocks;
#' "B" for type B restrictions on historical decomposition, i.e. if the absolute value of the historical decomposition
#' of shock to var is the greater/less than the sum of all other shocks.
#' 
#' @param sign integer - the sign of the narrative restriction (1 or -1).
#' 
#' @param shock positive integer - the index of the shock to which the narrative restriction applies.
#' 
#' @param var positive integer - the index of the variable to which the narrative restriction applies.
#' 
#' @return An object of class \code{narrative} specifying one narrative restriction.
#' 
#' @examples
#' # specify a narrative restriction
#' narrative       = specify_narrative(
#'                     start = 166, 
#'                     periods = 1, 
#'                     type = "S", 
#'                     sign = 1, 
#'                     shock = 1, 
#'                     var = 6
#'                   )
#' # use it to specify the model
#' specification   = specify_bsvarSIGN$new(monetary, sign_narrative = list(narrative))
#' 
#' @export
specify_narrative = function(start,
                             periods = 1,
                             type = "S",
                             sign = 1,
                             shock = 1,
                             var = NA) {
  
  if (start %% 1 != 0 || start <= 0) {
    stop("start must be a positive integer")
  }
  
  if (periods %% 1 != 0 || periods <= 0) {
    stop("periods must be a positive integer")
  }
  
  if (!(type %in% c("S", "A", "B", "C"))) {
    stop("type must be one of 'S', 'A', 'B', 'C'")
  }
  
  if (!(sign %in% c(-1, 1))) {
    stop("sign must be one of -1, 1")
  }
  
  if (!is.numeric(shock) ||
      length(shock) == 0 ||
      any(shock %% 1 != 0) ||
      any(shock <= 0)) {
    stop("shock must be a positive integer vector")
  }
  
  shock = unique(as.integer(shock))
  
  if (type != "C" && length(shock) != 1) {
    stop("For types 'S', 'A', and 'B', shock must have length 1. Use type 'C' for multiple shocks.")
  }
  
  if (!is.na(var)) {
    if (var %% 1 != 0 || var <= 0) {
      stop("var must be a positive integer")
    }
  }
  
  if (type == "C" && is.na(var)) {
    stop("For type 'C' narrative restrictions, var must be provided.")
  }
  
  narrative = list(
    start   = start,
    periods = periods,
    type    = type,
    sign    = sign,
    shock   = shock,
    var     = var
  )
  
  class(narrative) = "narrative"
  narrative
}

## helper function to add N extra columns for type "C" narratives
make_narrative_matrix = function(sign_narrative, p, N) {
  
  get_type = list(
    "S" = 1,
    "A" = 2,
    "B" = 3,
    "C" = 4
  )
  
  n_narratives = length(sign_narrative)
  
  if (n_narratives == 0) {
    narrative = matrix(0, 1, 6 + N)
    narrative[1, 1:6] = c(0, 1, 1, 1, 1, 1)
    return(narrative)
  }
  
  narrative = matrix(0, n_narratives, 6 + N)
  
  for (i in seq_len(n_narratives)) {
    
    nr = sign_narrative[[i]]
    
    if (!inherits(nr, "narrative")) {
      stop("Each element of sign_narrative must be of class 'narrative'.")
    }
    
    type_num = get_type[[nr$type]]
    
    if (is.null(type_num)) {
      stop("Unknown narrative type.")
    }
    
    narrative[i, 1] = type_num
    narrative[i, 2] = nr$sign
    
    # var is required for A, B, C.
    # For S it is not used, but set it to 1 if missing.
    narrative[i, 3] = ifelse(is.na(nr$var), 1, nr$var)
    
    narrative[i, 5] = nr$start - p
    narrative[i, 6] = nr$periods - 1
    
    shocks = unique(as.integer(nr$shock))
    
    if (any(shocks < 1 | shocks > N)) {
      stop("Narrative shock index is outside the valid range 1:N.")
    }
    
    if (nr$type == "C") {
      
      if (length(shocks) == 0) {
        stop("Type C narrative restriction requires at least one shock in the group.")
      }
      
      # The scalar shock column is unused for type C.
      narrative[i, 4] = 0
      
      # Group membership dummies.
      # Column 6 + j corresponds to shock j.
      narrative[i, 6 + shocks] = 1
      
    } else {
      
      if (length(shocks) != 1) {
        stop("For narrative types S, A, and B, shock must be scalar.")
      }
      
      narrative[i, 4] = shocks
      
      # Optional: also fill the group dummy column.
      # Not needed for S/A/B, but harmless.
      narrative[i, 6 + shocks] = 1
    }
  }
  
  narrative
}



#' Specify short-run elasticity restrictions
#'
#' @param shock integer vector: structural shock index.
#' @param num_var integer vector: numerator variable index.
#' @param den_var integer vector: denominator variable index.
#' @param lower numeric vector: lower bound for ratio.
#' @param upper numeric vector: upper bound for ratio.
#' @param horizon integer vector: horizon, 0 for impact.
#' @param use_abs logical vector: if TRUE, use abs(num)/abs(den).
#'
#' @return A matrix with columns:
#' shock, num_var, den_var, lower, upper, horizon, use_abs.
#'
#' @export
specify_elasticity = function(
    shock,
    num_var,
    den_var,
    lower = 0,
    upper = 1,
    horizon = 0,
    use_abs = FALSE
) {
  n = max(
    length(shock),
    length(num_var),
    length(den_var),
    length(lower),
    length(upper),
    length(horizon),
    length(use_abs)
  )
  
  shock   = rep_len(as.integer(shock), n)
  num_var = rep_len(as.integer(num_var), n)
  den_var = rep_len(as.integer(den_var), n)
  lower   = rep_len(as.numeric(lower), n)
  upper   = rep_len(as.numeric(upper), n)
  horizon = rep_len(as.integer(horizon), n)
  use_abs = rep_len(as.logical(use_abs), n)
  
  if (any(shock <= 0) || any(num_var <= 0) || any(den_var <= 0)) {
    stop("shock, num_var and den_var must be positive integers.")
  }
  
  if (any(shock %% 1 != 0) || any(num_var %% 1 != 0) || any(den_var %% 1 != 0)) {
    stop("shock, num_var and den_var must be integers.")
  }
  
  if (any(horizon < 0) || any(horizon %% 1 != 0)) {
    stop("horizon must be a non-negative integer.")
  }
  
  lower[is.na(lower)] = -Inf
  upper[is.na(upper)] =  Inf
  
  if (any(lower > upper)) {
    stop("lower bound cannot be larger than upper bound.")
  }
  
  
  out = cbind(
    shock   = shock,
    num_var = num_var,
    den_var = den_var,
    lower   = lower,
    upper   = upper,
    horizon = horizon,
    use_abs = as.numeric(use_abs)
  )
  
  
  class(out) = "elasticity"
  out
}



#' Specify response bound restrictions
#'
#' @param var integer vector: variable index.
#' @param shock integer vector: structural shock index.
#' @param horizon integer vector: IRF horizon, 0 for impact.
#' @param lower numeric vector: lower bound.
#' @param upper numeric vector: upper bound.
#' @param use_abs logical vector: if TRUE, restrict absolute value of response.
#'
#' @return A matrix with columns:
#' var, shock, horizon, lower, upper, use_abs.
#'
#' @export
specify_response_bound = function(
    var,
    shock,
    horizon = 0,
    lower = -Inf,
    upper =  Inf,
    use_abs = FALSE
) {
  n = max(
    length(var),
    length(shock),
    length(horizon),
    length(lower),
    length(upper),
    length(use_abs)
  )
  
  var     = rep_len(as.integer(var), n)
  shock   = rep_len(as.integer(shock), n)
  horizon = rep_len(as.integer(horizon), n)
  lower   = rep_len(as.numeric(lower), n)
  upper   = rep_len(as.numeric(upper), n)
  use_abs = rep_len(as.logical(use_abs), n)
  
  if (any(var <= 0) || any(shock <= 0)) {
    stop("var and shock must be positive integers.")
  }
  
  if (any(var %% 1 != 0) || any(shock %% 1 != 0) || any(horizon %% 1 != 0)) {
    stop("var, shock and horizon must be integers.")
  }
  
  if (any(horizon < 0)) {
    stop("horizon must be a non-negative integer.")
  }
  
  lower[is.na(lower)] = -Inf
  upper[is.na(upper)] =  Inf
  
  if (any(lower > upper)) {
    stop("lower bound cannot be larger than upper bound.")
  }
  
  if (any(use_abs & lower < 0)) {
    stop("If use_abs = TRUE, lower bound must be non-negative.")
  }
  
  out = cbind(
    var     = var,
    shock   = shock,
    horizon = horizon,
    lower   = lower,
    upper   = upper,
    use_abs = as.numeric(use_abs)
  )
  
  
  class(out) = "response_bound"
  out
}


# construct Z_j matrices
get_Z = function(sign_irf) {
  h = dim(sign_irf)[3]
  if (h >= 2) {
    test = sign_irf[, , 2:h]
    if (any(test[!is.na(test)] == 0)) {
      stop("Zero restrictions are not allowed for horizons >= 1")
    }
  }
  
  zero_irf = sign_irf[, , 1] == 0
  zero_irf[is.na(zero_irf)] = 0
  
  if (sum(zero_irf) == 0) {
    return(NULL)
  }
  
  N = dim(zero_irf)[2]
  
  Z = list()
  for (j in 1:N) {
    Z_j          = diag(zero_irf[, j])  
    nonzero_rows = rowSums(Z_j) > 0
    Z_j          = as.matrix(Z_j[nonzero_rows, ])
    
    if (dim(Z_j)[2] == 1) {
      Z_j = as.matrix(t(Z_j))
    }
    
    if (dim(Z_j)[1] > N-j) {
      stop("Too many zero restrictions for shock ", j)
    }
    
    Z[[j]]       = Z_j
  }
  
  Z
}


verify_elasticity = function(N, elasticity) {
  if (is.null(elasticity) || nrow(elasticity) == 0) {
    return(invisible(TRUE))
  }
  
  if (ncol(elasticity) != 7) {
    stop("elasticity matrix must have 7 columns: shock, num_var, den_var, lower, upper, horizon, use_abs.")
  }
  
  idx = elasticity[, 1:3, drop = FALSE]
  
  if (any(idx %% 1 != 0)) {
    stop("shock, num_var and den_var must be integers.")
  }
  
  if (any(idx < 1) || any(idx > N)) {
    stop("shock, num_var and den_var must be between 1 and N.")
  }
  
  if (any(elasticity[, 6] < 0)) {
    stop("elasticity horizon must be non-negative.")
  }
  
  if (any(elasticity[, 4] > elasticity[, 5])) {
    stop("elasticity lower bound cannot be larger than upper bound.")
  }
  
  invisible(TRUE)
}


verify_response_bounds = function(N, response_bounds) {
  if (is.null(response_bounds) || nrow(response_bounds) == 0) {
    return(invisible(TRUE))
  }
  
  if (ncol(response_bounds) != 6) {
    stop(
      paste0(
        "response_bounds matrix must have 6 columns: ",
        "var, shock, horizon, lower, upper, use_abs."
      )
    )
  }
  
  idx = response_bounds[, 1:2, drop = FALSE]
  
  if (any(idx %% 1 != 0)) {
    stop("var and shock in response_bounds must be integers.")
  }
  
  if (any(idx < 1) || any(idx > N)) {
    stop("var and shock in response_bounds must be between 1 and N.")
  }
  
  if (any(response_bounds[, 3] %% 1 != 0) || any(response_bounds[, 3] < 0)) {
    stop("horizon in response_bounds must be a non-negative integer.")
  }
  
  if (any(response_bounds[, 4] > response_bounds[, 5])) {
    stop("lower bound cannot be larger than upper bound in response_bounds.")
  }
  
  if (any(response_bounds[, 6] > 0.5 & response_bounds[, 4] < 0)) {
    stop("If use_abs = TRUE, lower bound must be non-negative.")
  }
  
  invisible(TRUE)
}


# verify if the matrix is NxN and has entries in {-1, 0, 1}
verify_traditional = function(N, A) {
  if (!(is.matrix(A) && all(dim(A) == c(N, N)))) {
    stop("Sign restriction matrix is not NxN.")
  }
  if (!(all(A %in% c(-1, 0, 1, NA)))) {
    stop("Sign restriction matrix has entries that are not in {-1, 0, 1, NA}.")
  }
}

# verify all restrictions
verify_all = function(N, sign_irf, sign_irf_cum, sign_narrative, sign_structural, elasticity, elasticity_cum, response_bounds, response_bounds_cum) {
  verify_traditional(N, sign_structural)
  if (any(sign_structural[!is.na(sign_structural)] == 0)) {
    stop("Zero restrictions are not allowed for sign_structural")
  }
  
  if (!is.list(sign_narrative)) {
    stop("sign_narrative must be a list")
  } else if (length(sign_narrative) > 0) {
    for (i in 1:length(sign_narrative)) {
      if (!inherits(sign_narrative[[i]], "narrative")) {
        stop("Each element of sign_narrative must be of class 'narrative'")
      }
    }
  }
  
  for (h in 1:dim(sign_irf)[3]) {
    verify_traditional(N, sign_irf[,,h])
  }
  
  for (h in 1:dim(sign_irf_cum)[3]) {
    verify_traditional(N, sign_irf_cum[,,h])
  }
  
  verify_response_bounds(N, response_bounds)
  verify_response_bounds(N, response_bounds_cum)
  verify_elasticity(N, elasticity)
  verify_elasticity(N, elasticity_cum)
}


# Minnesota prior of Normal-Inverse-Wishart form
# niw_prior = function(Y,
#                      p,
#                      non_stationary,
#                      lambda = 0.2) {
#   T = nrow(Y)
#   N = ncol(Y)
#   K = 1 + N * p
#   
#   B           = matrix(0, K, N)
#   B[1:N, 1:N] = diag(non_stationary)
#   
#   sigma2                  = sapply(1:N, \(i) summary(lm(Y[2:T, i] ~ Y[1:(T - 1), i]))$sigma^2)
#   V                       = matrix(0, K, K)
#   V[K, K]                 = 1e+6
#   V[1:(K - 1), 1:(K - 1)] = diag(lambda^2 * kronecker((1:p)^-2, sigma2^-1))
#   
#   S  = diag(sigma2)
#   nu = N + 2
#   
#   list(B = B, V = V, S = S, nu = nu)
# }

gamma_scale = function(mode, variance) {
  (2 * mode * variance) / (mode^2 + sqrt(mode^4 + 4 * variance * mode^2))
}

gamma_shape = function(mode, variance) {
  (mode^2 + sqrt(mode^4 + 4 * (mode^2) * variance) + 2 * variance) / (2 * variance)
}

igamma_scale = function(mode, variance) {
  (0.5 * mode * (sqrt(variance * (variance + 24 * mode^2)) - 5 * variance)) / (mode^2 - variance)
}

igamma_shape = function(mode, variance) {
  0.5 * (sqrt(variance * (variance + 24 * mode^2)) - 2 * mode^2 - 3 * variance ) / (mode^2 - variance)
}

#' R6 Class Representing PriorBSVAR
#'
#' @description
#' The class PriorBSVARSIGN presents a prior specification for the homoskedastic bsvar model.
#' 
#' @examples
#' # a prior for 5-variable example with one lag
#' data(optimism)
#' prior = specify_prior_bsvarSIGN$new(optimism, p = 1)
#' prior$A  # show autoregressive prior mean
#' 
#' @export
specify_prior_bsvarSIGN = R6::R6Class(
  "PriorBSVARSIGN",
  
  public = list(
    #' @field p a positive integer - the number of lags.
    p           = 1,
    
    #' @field hyper a \code{(N+4)xS} matrix of hyper-parameters \eqn{\mu, \delta, \lambda, \phi, \psi}.
    hyper      = matrix(),
    
    #' @field A a \code{NxK} normal prior mean matrix for the autoregressive 
    #' parameters.
    A          = matrix(),
    
    #' @field V a \code{KxK} matrix determining  the normal prior column-specific 
    #' covariance for the autoregressive parameters.
    V          = matrix(),
    
    #' @field S an \code{NxN} matrix determining the inverted-Wishart prior scale 
    #' of error terms covariance matrix.
    S          = matrix(),
    
    #' @field nu a positive scalar greater than \code{N+1} - the shape of the 
    #' inverted-Wishart prior for error terms covariance matrix.
    nu          = NA,
    
    #' @field data an \code{TxN} matrix of observations.
    data        = matrix(),
    
    #' @field Y an \code{NxT} matrix of dependent variables.
    Y           = matrix(),
    
    #' @field X an \code{KxT} matrix of independent variables.
    X           = matrix(),
    
    #' @field Ysoc an \code{NxN} matrix with the sum-of-coefficients dummy observations.
    Ysoc        = matrix(),
    
    #' @field Xsoc an \code{KxN} matrix with the sum-of-coefficients dummy observations.
    Xsoc        = matrix(),
    
    #' @field Ysur an \code{NxN} matrix with the single-unit-root dummy observations.
    Ysur        = matrix(),
    
    #' @field Xsur an \code{KxN} matrix with the single-unit-root dummy observations.
    Xsur        = matrix(),
    
    #' @field mu.scale a positive scalar - the shape of the gamma prior for \eqn{\mu}.
    mu.scale    = NA,
    
    #' @field mu.shape a positive scalar - the shape of the gamma prior for \eqn{\mu}.
    mu.shape    = NA,
    
    #' @field delta.scale a positive scalar - the shape of the gamma prior for \eqn{\delta}.
    delta.scale = NA,

    #' @field delta.shape a positive scalar - the shape of the gamma prior for \eqn{\delta}.
    delta.shape = NA,
    
    #' @field lambda.scale a positive scalar - the shape of the gamma prior for \eqn{\lambda}.
    lambda.scale = NA,
    
    #' @field lambda.shape a positive scalar - the shape of the gamma prior for \eqn{\lambda}.
    lambda.shape = NA,
    
    #' @field phi.scale a positive scalar - the scale of the gamma prior for \eqn{\phi}.
    phi.scale   = NA,
    
    #' @field phi.shape a positive scalar - the shape of the gamma prior for \eqn{\phi}.
    phi.shape   = NA,
    
    #' @field idx_dummy an integer vector of 0-based column positions in \code{X}
    #' holding the Pandemic-Priors time dummies.
    idx_dummy   = integer(0),
    
    #' @field psi.scale a positive scalar - the shape of the inverted gamma prior for \eqn{\psi}.
    psi.scale   = NA,
    
    #' @field psi.shape a positive scalar - the shape of the inverted gamma prior for \eqn{\psi}.
    psi.shape   = NA,
    
    #' @description
    #' Create a new prior specification PriorBSVAR.
    #' @param data the \code{TxN} data matrix of observations.
    #' @param p a positive integer - the autoregressive lag order of the SVAR model.
    #' @param exogenous a \code{Txd} matrix of exogenous variables.
    #' @param stationary an \code{N} logical vector - its element set to \code{FALSE} sets 
    #' the prior mean for the autoregressive parameters of the \code{N}th equation to the white noise process, 
    #' otherwise to random walk.
    #' @return A new prior specification PriorBSVARSIGN.
    #' @examples 
    #' # a prior for 5-variable example with one lag and stationary data
    #' data(optimism)
    #' prior = specify_prior_bsvarSIGN$new(optimism, p = 1)
    #' prior$B # show autoregressive prior mean
    #' 
    initialize = function(data, p, exogenous = NULL, stationary = rep(FALSE,  ncol(data)),
                          dummy_cols = NULL) {
      
      stopifnot("Argument p must be a positive integer number." = p > 0 & p %% 1 == 0)
      
      data_m  = bsvars::specify_data_matrices$new(data, p, exogenous)
      Y       = t(data_m$Y)
      X       = t(data_m$X)
      N       = ncol(Y)
      
      stopifnot("Argument stationary must be a logical vector of length equal to the number of columns in data." = length(stationary) == N & is.logical(stationary))
      
      d       = 0
      if (!is.null(exogenous)) {
        d     = ncol(exogenous)
      }
      T       = nrow(Y)
      K       = N * p + 1 + d
      
      B       = matrix(0, K, N)
      B[1:N,] = diag(!stationary)
      
      V       = diag(c(kronecker((1:p)^-2, rep(1, N)), rep(100, 1 + d)))
      
      s2.ols  = rep(NA, N)
      for (n in 1:N) {
        y     = as.matrix(Y[(p + 5 + 1):T, n])
        x     = matrix(1, T - p - 5, 1)
        for (i in 1:(p + 5)) {
          x   = cbind(x, Y[(p + 5 + 1):T - i, n])
        }
        s2.ols[n] = sum(((diag(T - p - 5) - x %*% solve(t(x) %*% x) %*% t(x)) %*% y)^2) / (T - p - 5)
      }
      
      hyper              = matrix(NA, N + 4, 1)
      hyper[1:4]         = c(1, 1, 0.2, 0.1)
      hyper[5:(N + 4),]  = s2.ols
      
      scale   = gamma_scale(1, 1)
      shape   = gamma_shape(1, 1)
      
      ybar    = colMeans(matrix(Y[1:p,], ncol = N))
      Ysoc    = diag(ybar)
      Ysur    = t(ybar)
      Xsoc    = cbind(kronecker(t(rep(1, p)), Ysoc), matrix(0, N, d + 1))
      Xsur    = cbind(kronecker(t(rep(1, p)), Ysur), 1, matrix(0, 1, d))
      
      idx_dummy = integer(0)
      if (!is.null(dummy_cols)) {
        stopifnot("dummy_cols must index columns of exogenous." =
                    d > 0 && all(dummy_cols %in% 1:d))
        idx_dummy = as.integer(N * p + dummy_cols)   # 0-based positions in X
      }
      
      # Ystar   = rbind(diag(ybar), ybar)
      # Xstar   = Ystar
      # if (p > 1) {
      #   for (i in 2:p) {
      #     Xstar = cbind(Xstar, Ystar)
      #   }
      # }
      # Xstar   = cbind(Xstar, c(rep(0, N), 1), matrix(0, N + 1, d))
      
      self$p             = p
      self$hyper         = hyper
      self$A             = t(B)
      self$V             = V
      self$S             = diag(N)
      self$nu            = N + 2
      self$Y             = t(Y)
      self$X             = t(X)
      self$Ysoc          = t(Ysoc)
      self$Xsoc          = t(Xsoc)
      self$Ysur          = t(Ysur)
      self$Xsur          = t(Xsur)
      self$mu.scale      = scale
      self$mu.shape      = shape
      self$delta.scale   = scale
      self$delta.shape   = shape
      self$lambda.scale  = gamma_scale(0.2, 0.4)
      self$lambda.shape  = gamma_shape(0.2, 0.4)
      self$phi.scale     = gamma_scale(0.1, 0.5)
      self$phi.shape     = gamma_shape(0.1, 0.5)
      self$idx_dummy     = idx_dummy
      self$psi.scale     = igamma_scale(0.02^2, 0.02^2)
      self$psi.shape     = igamma_shape(0.02^2, 0.02^2)
    }, # END initialize
    
    #' @description
    #' Returns the elements of the prior specification PriorBSVAR as a \code{list}.
    #' 
    #' @examples 
    #' # a prior for 5-variable example with four lags
    #' prior = specify_prior_bsvar$new(N = 5, p = 4)
    #' prior$get_prior() # show the prior as list
    #' 
    get_prior = function(){
      list(
        p            = self$p,
        hyper        = self$hyper,
        A            = self$A,
        V            = self$V,
        S            = self$S,
        nu           = self$nu,
        Ysoc         = self$Ysoc,
        Xsoc         = self$Xsoc,
        Ysur         = self$Ysur,
        Xsur         = self$Xsur,
        mu.scale     = self$mu.scale,
        mu.shape     = self$mu.shape,
        delta.scale  = self$delta.scale,
        delta.shape  = self$delta.shape,
        lambda.scale = self$lambda.scale,
        lambda.shape = self$lambda.shape,
        phi.scale    = self$phi.scale,
        phi.shape    = self$phi.shape,
        idx_dummy    = self$idx_dummy,
        psi.scale    = self$psi.scale,
        psi.shape    = self$psi.shape
      )
    }, # END get_prior
    
    #' @description
    #' Estimates hyper-parameters with adaptive Metropolis algorithm.
    #' 
    #' @param mu whether to estimate the hyper-parameter in the 
    #' sum-of-coefficients dummy prior.
    #' @param delta whether to estimate the hyper-parameter in the 
    #' single-unit-root dummy prior.
    #' @param lambda whether to estimate the hyper-parameter of the 
    #' shrinkage in the Minnesota prior.
    #' @param psi whether to estimate the hyper-parameter of the 
    #' variances in the Minnesota prior.
    #' @param S number of MCMC draws.
    #' @param burn_in number of burn-in draws.
    #' 
    #' @examples 
    #' # specify the model and set seed
    #' set.seed(123)
    #' data(optimism)
    #' prior = specify_prior_bsvarSIGN$new(optimism, p = 4)
    #' 
    #' # estimate hyper parameters with adaptive Metropolis algorithm
    #' prior$estimate_hyper(S = 10, psi = TRUE)
    #'
    #' # trace plot
    #' hyper = t(prior$hyper)
    #' colnames(hyper) = c("mu", "delta", "lambda", paste("psi", 1:5, sep = ""))
    #' plot.ts(hyper)
    #' 
    estimate_hyper = function(
      S = 10000, burn_in = S / 2,
      mu = FALSE, delta = FALSE, lambda = TRUE, phi = FALSE, psi = FALSE
      ) {
      
      model = c(mu, delta, lambda, phi, psi)
      
      if (all(!model)) {
        stop("At least one of the hyper-parameters must be estimated.")
      }
      
      hyper  = matrix(self$hyper[, ncol(self$hyper)])
      init   = narrow_hyper(model, hyper)
      prior  = self$get_prior()
      
      prior$B    = t(prior$A)
      prior$Ysoc = t(prior$Ysoc)
      prior$Xsoc = t(prior$Xsoc)
      prior$Ysur = t(prior$Ysur)
      prior$Xsur = t(prior$Xsur)
      
      result = stats::optim(
        init,
        \(x) -log_posterior_hyper(extend_hyper(hyper, model, matrix(x)), 
                                  model, t(self$Y), t(self$X), prior),
        method  = 'L-BFGS-B',
        lower   = rep(0, length(init)),
        upper   = init * 100,
        hessian = TRUE
        )

      mode       = extend_hyper(hyper, model, matrix(result$par))
      variance   = result$hessian

      if (length(init) == 1){
        variance = 1 / variance
      } else {
        e        = eigen(variance)
        variance = e$vectors %*% diag(as.vector(1 / abs(e$values))) %*% t(e$vectors)
      }
      
      self$hyper = sample_hyper(S, burn_in, mode, model, 
                                t(self$Y), t(self$X), variance, prior)
      self$hyper = self$hyper[, -(1:burn_in)]
    } # END estimate_hyper
    
  ) # END public
) # END specify_prior_bsvarSIGN


#' R6 Class Representing IdentificationBSVARSIGN
#'
#' @description
#' The class IdentificationBSVARSIGN presents the identifying restrictions for the Bayesian Structural VAR models with sign and narrative restrictions.
#'
#' @examples 
#' # recursive specification for a 5-variable system
#' specify_identification_bsvarSIGN$new(N = 5)
#' 
#' # specify sign restrictions of the first shock on the contemporaneous IRF
#' # + no effect on the first variable
#' # + positive effect on the second variable
#' sign_irf = matrix(c(0, 1, rep(NA, 23)), 5, 5)
#' specify_identification_bsvarSIGN$new(N = 5, sign_irf = sign_irf) 
#'
#' @export
specify_identification_bsvarSIGN = R6::R6Class(
  "IdentificationBSVARSIGN",
  
  public = list(
    
    #' @field VB a list of \code{N} matrices determining the unrestricted elements of matrix \eqn{B}.
    VB       = list(),
    #' @field sign_irf a \code{NxNxH} array of sign restrictions on the impulse response functions.
    sign_irf = array(),
    #' @field sign_irf_cum a \code{NxNxH} array of sign restrictions on the cumulative impulse response functions.
    sign_irf_cum = array(),
    #' @field sign_narrative a \code{ANYx6} matrix of narrative sign restrictions.
    sign_narrative  = matrix(),
    #' @field sign_structural a \code{NxN} matrix of sign restrictions on contemporaneous relations.
    sign_structural   = matrix(),
    #' @field elasticity a \code{ANYx7} matrix of elasticity restrictions.
    elasticity   = matrix(),
    #' @field elasticity_cum a \code{ANYx7} matrix of elasticity restrictions. Helpful when restricting level response to growth rate variable.
    elasticity_cum   = matrix(),
    #' @field response_bounds a \code{ANYx6} matrix of response bounds restrictions.
    response_bounds   = matrix(),
    #' @field response_bounds_cum a \code{ANYx6} matrix of response bounds restrictions. Helpful when restricting level response to growth rate variable.
    response_bounds_cum   = matrix(),
    #' @field max_tries a positive integer with the maximum number of iterations 
    #' for finding a rotation matrix \eqn{Q} that would satisfy sign restrictions.
    max_tries = Inf,
    
    #' @description
    #' Create new identifying restrictions IdentificationBSVARSIGN.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param sign_irf a \code{NxNxH} array - sign and zero restrictions 
    #' on the impulse response functions, ±1 for positive/negative sign restriction
    #' 0 for zero restrictions and NA for no restrictions,
    #' the \code{h}-th slice \code{NxN} matrix contains the
    #' restrictions on the \code{h-1} horizon.
    #' @param sign_irf_cum a \code{NxNxH} array - sign and zero restrictions 
    #' on the cumulative impulse response functions, ±1 for positive/negative sign restriction
    #' 0 for zero restrictions and NA for no restrictions,
    #' the \code{h}-th slice \code{NxN} matrix contains the
    #' restrictions on the \code{h-1} horizon. Helpful when restricting level response to growth rate variable.
    #' @param sign_narrative a list of objects of class "narrative" - narrative sign restrictions.
    #' @param sign_structural a \code{NxN} matrix with entries ±1 or NA - sign restrictions on the
    #' contemporaneous relations \code{B} between reduced-form errors \code{E} and
    #' structural shocks \code{U} where \code{BE=U}.
    #' @param elasticity an object of class "elasticity" - elasticity restrictions.
    #' @param elasticity_cum an object of class "elasticity" - elasticity restrictions. Helpful when restricting level response to growth rate variable.
    #' @param response_bounds an object of class "response_bounds" - response_bounds restrictions.
    #' @param response_bounds_cum an object of class "response_bounds" - response_bounds restrictions. Helpful when restricting level response to growth rate variable.
    #' @param max_tries a positive integer with the maximum number of iterations
    #' for finding a rotation matrix \eqn{Q} that would satisfy sign restrictions.
    #' @return Identifying restrictions IdentificationBSVARSIGN.
    initialize = function(N, sign_irf, sign_irf_cum, sign_narrative, sign_structural, elasticity, elasticity_cum, response_bounds, response_bounds_cum, max_tries = Inf) {
        
      missing_all   = TRUE
      if (missing(sign_irf)) {
        sign_irf = array(rep(NA, N^2), dim = c(N, N, 1))
      } else {
        missing_all = FALSE
      }
      if (missing(sign_irf_cum)) {
        sign_irf_cum = array(rep(NA, N^2), dim = c(N, N, 1))
      } else {
        missing_all = FALSE
      }
      if (missing(sign_narrative)) {
        sign_narrative = list()
      } else {
        missing_all = FALSE
      }
      if (missing(elasticity)) {
        elasticity = NULL
      } else {
        missing_all = FALSE
      }
      if (missing(elasticity_cum)) {
        elasticity_cum = NULL
      } else {
        missing_all = FALSE
      }
      if (missing(response_bounds)) {
        response_bounds = NULL
      } else {
        missing_all = FALSE
      }
      if (missing(response_bounds_cum)) {
        response_bounds_cum = NULL
      } else {
        missing_all = FALSE
      }
      if (missing(sign_structural)) {
        sign_structural = matrix(rep(NA, N^2), ncol = N, nrow = N)
        if (missing_all) {
          diag(sign_structural) = 1
        }
      }
      
      if (is.matrix(sign_irf)) {
        sign_irf = array(sign_irf, dim = c(dim(sign_irf), 1))
      }
      if (is.matrix(sign_irf_cum)) {
        sign_irf_cum = array(sign_irf_cum, dim = c(dim(sign_irf_cum), 1))
      }
      
      verify_all(N, sign_irf, sign_irf_cum, sign_narrative, sign_structural, elasticity, elasticity_cum, response_bounds, response_bounds_cum)
      
      B     = matrix(FALSE, N, N)
      B[lower.tri(B, diag = TRUE)] = TRUE
      
      self$VB = vector("list", N)
      for (n in 1:N) {
        self$VB[[n]] = matrix(diag(N)[B[n,],], ncol = N)
      }
      
      self$sign_irf       = sign_irf
      self$sign_irf_cum       = sign_irf_cum
      self$sign_narrative = sign_narrative
      self$sign_structural  = sign_structural
      self$max_tries      = max_tries
      self$elasticity     = elasticity
      self$elasticity_cum     = elasticity_cum
      self$response_bounds = response_bounds
      self$response_bounds_cum = response_bounds_cum
    }, # END initialize
    
    #' @description
    #' Returns the elements of the identification pattern IdentificationBSVARSIGN as a \code{list}.
    #'
    get_identification = function() {
      list(
        VB             = as.list(self$VB),
        sign_irf       = as.array(self$sign_irf),
        sign_irf_cum   = as.array(self$sign_irf_cum),
        sign_narrative = self$sign_narrative,
        sign_structural  = as.matrix(self$sign_structural),
        elasticity     = self$elasticity,
        elasticity_cum     = self$elasticity_cum,
        response_bounds = self$response_bounds,
        response_bounds_cum = self$response_bounds_cum,
        max_tries      = self$max_tries
        )
    }, # END get_identification
    
    #' @description
    #' Set new starting values StartingValuesBSVARSIGN.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param sign_irf a \code{NxNxH} array - sign and zero restrictions 
    #' on the impulse response functions, ±1 for positive/negative sign restriction
    #' 0 for zero restrictions and NA for no restrictions,
    #' the \code{h}-th slice \code{NxN} matrix contains the
    #' restrictions on the \code{h-1} horizon.
    #' @param sign_irf_cum a \code{NxNxH} array - sign and zero restrictions 
    #' on the cumulative impulse response functions, ±1 for positive/negative sign restriction
    #' 0 for zero restrictions and NA for no restrictions,
    #' the \code{h}-th slice \code{NxN} matrix contains the
    #' restrictions on the \code{h-1} horizon. Helpful when restricting level response to growth rate variable.
    #' @param sign_narrative a list of objects of class "narrative" - narrative sign restrictions.
    #' @param sign_structural a \code{NxN} matrix with entries ±1 or NA - sign restrictions on the
    #' contemporaneous relations \code{B} between reduced-form errors \code{E} and
    #' structural shocks \code{U} where \code{BE=U}.
    #' @param elasticity an object of class "elasticity" - elasticity restrictions.
    #' @param elasticity_cum an object of class "elasticity" - elasticity restrictions. Helpful when restricting level response to growth rate variable.
    #' @param response_bounds an object of class "response_bounds" - response_bounds restrictions.
    #' @param response_bounds_cum an object of class "response_bounds" - response_bounds restrictions. Helpful when restricting level response to growth rate variable.
    #' @param max_tries a positive integer with the maximum number of iterations
    #' for finding a rotation matrix \eqn{Q} that would satisfy sign restrictions
    set_identification = function(N, sign_irf, sign_narrative, sign_structural, elasticity, response_bounds) {
      B     = matrix(FALSE, N, N)
      B[lower.tri(B, diag = TRUE)] = TRUE
      
      self$VB          <- vector("list", N)
      for (n in 1:N) {
        self$VB[[n]]   <- matrix(diag(N)[B[n,],], ncol = N)
      }
      
      missing_all   = TRUE
      if (missing(sign_irf)) {
        sign_irf = array(rep(NA, N^2), dim = c(N, N, 1))
      } else {
        missing_all = FALSE
      }
      if (missing(sign_irf_cum)) {
        sign_irf_cum = array(rep(NA, N^2), dim = c(N, N, 1))
      } else {
        missing_all = FALSE
      }
      if (missing(sign_narrative)) {
        sign_narrative = list()
      } else {
        missing_all = FALSE
      }
      if (missing(elasticity)) {
        elasticity = NULL
      } else {
        missing_all = FALSE
      }
      if (missing(elasticity_cum)) {
        elasticity_cum = NULL
      } else {
        missing_all = FALSE
      }
      if (missing(response_bounds)) {
        response_bounds = NULL
      } else {
        missing_all = FALSE
      }
      if (missing(response_bounds_cum)) {
        response_bounds_cum = NULL
      } else {
        missing_all = FALSE
      }
      if (missing(sign_structural)) {
        sign_structural = matrix(rep(NA, N^2), ncol = N, nrow = N)
        if (missing_all) {
          diag(sign_structural) = 1
        }
      }

      
      if (is.matrix(sign_irf)) {
        sign_irf = array(sign_irf, dim = c(dim(sign_irf), 1))
      }
      if (is.matrix(sign_irf_cum)) {
        sign_irf_cum = array(sign_irf_cum, dim = c(dim(sign_irf_cum), 1))
      }
      
      verify_all(N, sign_irf, sign_irf_cum, sign_narrative, sign_structural, elasticity, elasticity_cum, response_bounds, response_bounds_cum)
      
      self$sign_irf       = sign_irf
      self$sign_irf_cum       = sign_irf_cum
      self$sign_narrative = sign_narrative
      self$sign_structural  = sign_structural
      self$elasticity     = elasticity
      self$elasticity_cum     = elasticity_cum
      self$response_bounds     = response_bounds
      self$response_bounds_cum     = response_bounds_cum
    } # END set_identification
  ) # END public
) # END specify_identification_bsvarSIGN



#' R6 Class representing the specification of the BSVARSIGN model
#'
#' @description
#' The class BSVARSIGN presents complete specification for the Bayesian Structural VAR model with sign and narrative restrictions.
#'
#' @seealso \code{\link{estimate.BSVARSIGN}}, \code{\link{specify_posterior_bsvarSIGN}}
#' 
#' @examples
#' # specify a model with the optimism data and 4 lags
#' 
#' data(optimism)
#' specification = specify_bsvarSIGN$new(
#'    data = optimism,
#'    p = 4
#' )
#'
#' @export
specify_bsvarSIGN = R6::R6Class(
  "BSVARSIGN",
  
  public = list(
    
    #' @field p a non-negative integer specifying the autoregressive lag order of the model.
    p                      = numeric(),
    
    #' @field identification an object IdentificationBSVARSIGN with the identifying restrictions.
    identification         = list(),
    
    #' @field prior an object PriorBSVARSIGN with the prior specification.
    prior                  = list(),
    
    #' @field data_matrices an object DataMatricesBSVARSIGN with the data matrices.
    data_matrices          = list(),
    
    #' @field starting_values an object StartingValuesBSVARSIGN with the starting values.
    starting_values        = list(),
    
    #' @description
    #' Create a new specification of the Bayesian Structural VAR model with sign and narrative restrictions BSVARSIGN.
    #' @param data a \code{(T+p)xN} matrix with time series data.
    #' @param p a positive integer providing model's autoregressive lag order.
    #' @param sign_irf a \code{NxNxH} array - sign and zero restrictions 
    #' on the impulse response functions, ±1 for positive/negative sign restriction
    #' 0 for zero restrictions and NA for no restrictions,
    #' the \code{h}-th slice \code{NxN} matrix contains the
    #' restrictions on the \code{h-1} horizon.
    #' @param sign_irf_cum a \code{NxNxH} array - sign and zero restrictions 
    #' on the cumulative impulse response functions, ±1 for positive/negative sign restriction
    #' 0 for zero restrictions and NA for no restrictions,
    #' the \code{h}-th slice \code{NxN} matrix contains the
    #' restrictions on the \code{h-1} horizon. Helpful when restricting level response to growth rate variable.
    #' @param sign_narrative a list of objects of class "narrative" - narrative sign restrictions.
    #' @param sign_structural a \code{NxN} matrix with entries ±1 or NA - sign restrictions on the
    #' contemporaneous relations \code{B} between reduced-form errors \code{E} and
    #' structural shocks \code{U} where \code{BE=U}.
    #' @param elasticity an object of class "elasticity" - elasticity restrictions.
    #' @param elasticity_cum an object of class "elasticity" - elasticity restrictions. Helpful when restricting level response to growth rate variable.
    #' @param response_bounds an object of class "response_bounds" - response_bounds restrictions.
    #' @param response_bounds_cum an object of class "response_bounds" - response_bounds restrictions. Helpful when restricting level response to growth rate variable.
    #' @param max_tries a positive integer with the maximum number of iterations
    #' for finding a rotation matrix \eqn{Q} that would satisfy sign restrictions
    #' @param exogenous a \code{(T+p)xd} matrix of exogenous variables.
    #' @param stationary an \code{N} logical vector - its element set to \code{FALSE} sets
    #' the prior mean for the autoregressive parameters of the \code{N}th equation to the white noise process,
    #' otherwise to random walk.
    #' @return A new complete specification for the Bayesian Structural VAR model BSVARSIGN.
    initialize = function(
    data,
    p = 1L,
    sign_irf,
    sign_irf_cum,
    sign_narrative,
    sign_structural,
    elasticity,
    elasticity_cum,
    response_bounds,
    response_bounds_cum,
    max_tries = Inf,
    exogenous = NULL,
    dummy_cols = NULL,
    stationary = rep(FALSE, ncol(data))
    ) {
      stopifnot("Argument p has to be a positive integer." = ((p %% 1) == 0 & p > 0))
      self$p        = p
      
      TT            = nrow(data)
      T             = TT - self$p
      N             = ncol(data)
      d             = 0
      if (!is.null(exogenous)) {
        d           = ncol(exogenous)
      }
      
      missing_all   = TRUE
      if (missing(sign_irf)) {
        sign_irf = array(rep(NA, N^2), dim = c(N, N, 1))
      } else {
        missing_all = FALSE
      }
      if (missing(sign_irf_cum)) {
        sign_irf_cum = array(rep(NA, N^2), dim = c(N, N, 1))
      } else {
        missing_all = FALSE
      }
      if (missing(sign_narrative)) {
        sign_narrative = list()
      } else {
        missing_all = FALSE
      }
      if (missing(elasticity)) {
        elasticity = NULL
      } else {
        missing_all = FALSE
      }
      if (missing(elasticity_cum)) {
        elasticity_cum = NULL
      } else {
        missing_all = FALSE
      }
      if (missing(response_bounds)) {
        response_bounds = NULL
      } else {
        missing_all = FALSE
      }
      if (missing(response_bounds_cum)) {
        response_bounds_cum = NULL
      } else {
        missing_all = FALSE
      }
      if (missing(sign_structural)) {
        sign_structural = matrix(rep(NA, N^2), ncol = N, nrow = N)
        if (missing_all) {
          diag(sign_structural) = 1
        }
      }
      
      if (is.matrix(sign_irf)) {
        sign_irf = array(sign_irf, dim = c(dim(sign_irf), 1))
      }
      if (is.matrix(sign_irf_cum)) {
        sign_irf_cum = array(sign_irf_cum, dim = c(dim(sign_irf_cum), 1))
      }
      
      verify_all(N, sign_irf, sign_irf_cum, sign_narrative, sign_structural, elasticity, elasticity_cum, response_bounds, response_bounds_cum)
      
      
      B                            = matrix(FALSE, N, N)
      B[lower.tri(B, diag = TRUE)] = TRUE
      
      self$data_matrices           = bsvars::specify_data_matrices$new(data, p, exogenous)
      self$identification          = specify_identification_bsvarSIGN$new(N,
                                                                          sign_irf,
                                                                          sign_irf_cum,
                                                                          sign_narrative,
                                                                          sign_structural,
                                                                          elasticity,
                                                                          elasticity_cum,
                                                                          response_bounds,
                                                                          response_bounds_cum,
                                                                          max_tries)
      self$prior                   = specify_prior_bsvarSIGN$new(data, p, exogenous,
                                                                 stationary, dummy_cols)
      self$starting_values         = bsvars::specify_starting_values_bsvar$new(N, self$p, d)
    }, # END initialize
    
    #' @description
    #' Returns the data matrices as the DataMatricesBSVAR object.
    #'
    #' @examples
    #' # specify a model with the optimism data and 4 lags
    #' 
    #' data(optimism)
    #' spec = specify_bsvarSIGN$new(
    #'    data = optimism,
    #'    p = 4
    #' )
    #' 
    #' # get the data matrices
    #' spec$get_data_matrices()
    #'
    get_data_matrices = function() {
      self$data_matrices$clone()
    }, # END get_data_matrices
    
    #' @description
    #' Returns the identifying restrictions as the IdentificationBSVARSIGN object.
    #'
    #' @examples 
    #' # specify a model with the optimism data and 4 lags
    #' data(optimism)
    #' spec = specify_bsvarSIGN$new(
    #'    data = optimism,
    #'    p = 4
    #' )
    #' 
    #' # get the identifying restrictions
    #' spec$get_identification()
    #'
    get_identification = function() {
      self$identification$clone()
    }, # END get_identification
    
    #' @description
    #' Returns the prior specification as the PriorBSVAR object.
    #'
    #' @examples 
    #' # specify a model with the optimism data and 4 lags
    #' 
    #' data(optimism)
    #' spec = specify_bsvarSIGN$new(
    #'    data = optimism,
    #'    p = 4
    #' )
    #' 
    #' # get the prior specification
    #' spec$get_prior()
    #'
    get_prior = function() {
      self$prior$clone()
    }, # END get_prior
    
    #' @description
    #' Returns the starting values as the StartingValuesBSVAR object.
    #'
    #' @examples 
    #' # specify a model with the optimism data and 4 lags
    #' 
    #' data(optimism)
    #' spec = specify_bsvarSIGN$new(
    #'    data = optimism,
    #'    p = 4
    #' )
    #' 
    #' # get the starting values
    #' spec$get_starting_values()
    #'
    get_starting_values = function() {
      self$starting_values$clone()
    } # END get_starting_values
  ) # END public
) # END specify_bsvarSIGN



#' R6 Class Representing PosteriorBSVARSIGN
#'
#' @description
#' The class PosteriorBSVARSIGN contains posterior output and the specification including
#' the last MCMC draw for the Bayesian Structural VAR model with sign and narrative restrictions.
#' Note that due to the thinning of the MCMC output the starting value in element \code{last_draw}
#' might not be equal to the last draw provided in element \code{posterior}.
#'
#' @seealso \code{\link{estimate.BSVARSIGN}}, \code{\link{specify_bsvarSIGN}}
#' 
#' @examples 
#' # This is a function that is used within estimate()
#' data(optimism)
#' specification  = specify_bsvarSIGN$new(optimism, p = 4)
#' set.seed(123)
#' posterior      = estimate(specification, 50)
#' class(posterior)
#'
#' @export
specify_posterior_bsvarSIGN = R6::R6Class(
  "PosteriorBSVARSIGN",
  
  private = list(
    normalised = TRUE
  ), # END private
  
  public = list(
    
    #' @field last_draw an object of class BSVARSIGN with the last draw of the current MCMC run as
    #' the starting value to be passed to the continuation of the MCMC estimation using \code{estimate()}.
    last_draw = list(),
    
    #' @field posterior a list containing Bayesian estimation output including:
    #' an \code{NxNxS} array \code{B}, an \code{NxKxS} array \code{A}, and a \code{5xS} matrix \code{hyper}.
    posterior = list(),
    
    #' @description
    #' Create a new posterior output PosteriorBSVARSIGN.
    #' @param specification_bsvarSIGN an object of class BSVARSIGN with the last draw of the current
    #' MCMC run as the starting value.
    #' @param posterior_bsvarSIGN a list containing Bayesian estimation output collected in elements
    #' an \code{NxNxS} array \code{B}, an \code{NxKxS} array \code{A}, and a \code{5xS} matrix \code{hyper}.
    #' @return A posterior output PosteriorBSVARSIGN.
    initialize = function(specification_bsvarSIGN, posterior_bsvarSIGN) {
      
      stopifnot("Argument specification_bsvarSIGN must be of class BSVARSIGN." = any(class(specification_bsvarSIGN) == "BSVARSIGN"))
      stopifnot("Argument posterior_bsvarSIGN must must contain MCMC output." = is.list(posterior_bsvarSIGN) & is.array(posterior_bsvarSIGN$B) & is.array(posterior_bsvarSIGN$A) & is.array(posterior_bsvarSIGN$hyper))
      
      self$last_draw    = specification_bsvarSIGN
      self$posterior    = posterior_bsvarSIGN
    }, # END initialize
    
    #' @description
    #' Returns a list containing Bayesian estimation output collected in elements
    #' an \code{NxNxS} array \code{B}, an \code{NxKxS} array \code{A}, and a \code{5xS} matrix \code{hyper}.
    #'
    #' @examples 
    #' data(optimism)
    #' specification  = specify_bsvarSIGN$new(optimism)
    #' set.seed(123)
    #' estimate       = estimate(specification, 50)
    #' estimate$get_posterior()
    #'
    get_posterior       = function(){
      self$posterior
    }, # END get_posterior
    
    #' @description
    #' Returns \code{TRUE} if the posterior has been normalised using \code{normalise_posterior()}
    #' and \code{FALSE} otherwise.
    #'
    #' @examples
    #' data(optimism)
    #' specification  = specify_bsvarSIGN$new(optimism)
    #' set.seed(123)
    #' estimate       = estimate(specification, 20)
    #' 
    #' # check normalisation status afterwards
    #' posterior$is_normalised()
    #'
    is_normalised      = function(){
      private$normalised
    } # END is_normalised
    
  ) # END public
) # END specify_posterior_bsvarSIGN


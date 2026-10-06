
# Given posterior draws with importance weights, sample with replacement
importance_sampling <- function(posterior) {
  
  # w <- posterior$posterior$w
  log_w <- as.vector(posterior$posterior$w)
  max_log_w <- max(log_w, na.rm = TRUE)
  w <- exp(log_w - max_log_w) # Safe exponentiation
  bad_w <- is.na(w) | is.infinite(w)
  w[bad_w] <- 0 # Safety catch
  
  if (sum(bad_w) > 0) {
    warning(paste0("There are exactly ", sum(bad_w), " of ", length(w)," NA or Inf importance weights."))
  }
  
  
  posterior$posterior$w <- NULL # remove weights
  
  # if (posterior$last_draw$identification$sign_narrative[1, 1] == 0) {
  #   return(posterior)
  # }
  
  if (sum(w) == 0) {
    stop("All draws were rejected or yielded zero weight. Try relaxing
          restrictions or increasing max_tries.")
  }
  
  ess <- sum(w)^2/sum(w^2)
  
  if (ess > 0.999*length(w)) {
    posterior$posterior$ess = ess
    return(posterior)
  }
  
  indices <- sample.int(length(w), length(w), replace = TRUE, prob = w)
  
  posterior$posterior$A      = posterior$posterior$A[, , indices, drop = F]
  posterior$posterior$B      = posterior$posterior$B[, , indices, drop = F]
  posterior$posterior$hyper  = posterior$posterior$hyper[, indices, drop = F]
  posterior$posterior$Q      = posterior$posterior$Q[, , indices, drop = F]
  posterior$posterior$Sigma  = posterior$posterior$Sigma[, , indices, drop = F]
  posterior$posterior$Theta0 = posterior$posterior$Theta0[, , indices, drop = F]
  posterior$posterior$shocks = posterior$posterior$shocks[, , indices, drop = F]
  posterior$posterior$ess    = ess
  
  return(posterior)
}

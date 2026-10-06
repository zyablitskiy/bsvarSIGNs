
#include <RcppArmadillo.h>

#include "utils.h"
#include "restrictions_narrative.h"
#include "restrictions_zero.h"
#include "compute.h"

using namespace Rcpp;
using namespace arma;


// Response bound restrictions
//
// Columns of response_bounds:
// 0: var_i, 1-based
// 1: shock_j, 1-based
// 2: horizon, 0-based
// 3: lower bound
// 4: upper bound
// 5: use_abs flag
//
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
bool match_response_bounds(
    const arma::mat& Q,
    const arma::cube& irf,
    const arma::mat& response_bounds
) {
  for (int r = 0; r < (int)response_bounds.n_rows; ++r) {
    
    int var_i    = (int)response_bounds(r, 0) - 1;
    int shock_j  = (int)response_bounds(r, 1) - 1;
    int horizon  = (int)response_bounds(r, 2);
    
    double lower  = response_bounds(r, 3);
    double upper  = response_bounds(r, 4);
    bool use_abs  = response_bounds(r, 5) > 0.5;
    
    if (var_i < 0 || var_i >= (int)irf.n_rows) return false;
    if (shock_j < 0 || shock_j >= (int)Q.n_cols) return false;
    if (horizon < 0 || horizon >= (int)irf.n_slices) return false;
    
    // Structural IRF at horizon h is irf.slice(h) * Q.
    // Element (var_i, shock_j) is computed by the dot product below.
    double val = arma::dot(irf.slice(horizon).row(var_i), Q.col(shock_j));
    
    if (use_abs) {
      val = std::abs(val);
    }
    
    if (val < lower || val > upper) {
      return false;
    }
  }
  
  return true;
}



// Short-run elasticity restrictions
// Columns of elasticity:
// 0: shock_j, 1-based
// 1: numerator variable, 1-based
// 2: denominator variable, 1-based
// 3: lower bound
// 4: upper bound
// 5: horizon, 0-based
// 6: use_abs flag
//
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
bool match_elasticity(
    const arma::mat& Q,
    const arma::cube& irf,
    const arma::mat& elasticity
) {
  for (int r = 0; r < (int)elasticity.n_rows; ++r) {
    
    int shock_j = (int)elasticity(r, 0) - 1;
    int num_i   = (int)elasticity(r, 1) - 1;
    int den_i   = (int)elasticity(r, 2) - 1;
    
    double lower  = elasticity(r, 3);
    double upper  = elasticity(r, 4);
    int horizon   = (int)elasticity(r, 5);
    bool use_abs  = elasticity(r, 6) > 0.5;
    
    if (shock_j < 0 || shock_j >= (int)Q.n_cols) return false;
    if (num_i   < 0 || num_i   >= (int)Q.n_rows) return false;
    if (den_i   < 0 || den_i   >= (int)Q.n_rows) return false;
    if (horizon < 0 || horizon >= (int)irf.n_slices) return false;
    
    // Structural impact/IRF response is irf.slice(horizon) * Q.
    // Element (i, j) is dot(irf.slice(horizon).row(i), Q.col(j)).
    double num = arma::dot(irf.slice(horizon).row(num_i), Q.col(shock_j));
    double den = arma::dot(irf.slice(horizon).row(den_i), Q.col(shock_j));
    
    // Elasticity is undefined if denominator is numerically zero.
    // You can adjust this tolerance depending on the scaling of your data.
    if (std::abs(den) < 1e-12) {
      Rcpp::stop("match_elasticity: denominator response is numerically zero "
                   "for shock %d, variable %d, horizon %d - check for a conflicting "
                   "zero restriction.", shock_j + 1, den_i + 1, horizon);
    }
    
    
    double ratio;
    if (use_abs) {
      ratio = std::abs(num) / std::abs(den);
    } else {
      ratio = num / den;
    }
    
    if (ratio < lower || ratio > upper) {
      return false;
    }
  }
  
  return true;
}


// If matches traditional sign restrictions on impulse response functions
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
bool match_sign_irf(
    const arma::mat& Q,
    const arma::cube& sign_irf,
    const arma::cube& irf
) {
  arma::uword N = Q.n_cols;
  
  // Loop over horizons
  for (arma::uword t = 0; t < sign_irf.n_slices; t++) {
    
    // Loop over columns (shocks)
    for (arma::uword j = 0; j < N; j++) {
      
      // Loop over rows (variables)
      for (arma::uword i = 0; i < N; i++) {
        
        double s = sign_irf(i, j, t);
        
        // Only do math if there is actually a sign restriction here!
        if (s != 0.0) {
          
          // Manually compute the (i, j) element of (irf.slice(t) * Q)
          double val = arma::dot(irf.slice(t).row(i), Q.col(j));
          // double val = 0.0;
          // for (arma::uword k = 0; k < N; k++) {
          //   val += irf(i, k, t) * Q(k, j);
          // }
          
          // Short-circuit instantly if the sign is violated
          if (s > 0.0 && val <= 0.0) return false;
          if (s < 0.0 && val >= 0.0) return false;
        }
      }
    }
  }
  
  return true;
}


// Sample rotation matrix Q and updates importance weight by reference
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::field<arma::mat> sample_Q(    
    const int&                    p,
    const arma::mat&              Y,
    const arma::mat&              X,
    arma::mat&                    B,
    arma::mat&                    h_invp,
    arma::mat&                    chol_Sigma,
    const Rcpp::List&             prior,
    const arma::cube&             sign_irf,
    const arma::cube&             sign_irf_cum,
    const arma::mat&              sign_narrative,
    const arma::mat&              sign_B,
    const arma::mat&              elasticity,
    const arma::mat&              elasticity_cum,
    const arma::mat&              response_bounds,
    const arma::mat&              response_bounds_cum,
    const arma::field<arma::mat>& Z,
    const int&                    max_tries,
    const arma::field<arma::mat>& W
) {
  
  const int N          = Y.n_cols;
  const int T          = Y.n_rows;
  
  int    n_tries       = 0;
  int h = 0;
  
  if (sign_narrative.n_rows > 0) {
    h = std::max(h, (int)sign_narrative.col(5).max());
  }
  
  if (sign_irf.n_slices > 0) {
    h = std::max(h, (int)sign_irf.n_slices - 1);
  }
  
  if (sign_irf_cum.n_slices > 0) {
    h = std::max(h, (int)sign_irf_cum.n_slices - 1);
  }
  
  if (elasticity.n_rows > 0) {
    h = std::max(h, (int)elasticity.col(5).max());
  }
  
  if (elasticity_cum.n_rows > 0) {
    h = std::max(h, (int)elasticity_cum.col(5).max());
  }
  
  if (response_bounds.n_rows > 0) {
    h = std::max(h, (int)response_bounds.col(2).max());
  }
  
  if (response_bounds_cum.n_rows > 0) {
    h = std::max(h, (int)response_bounds_cum.col(2).max());
  }
  
  
  bool has_narrative  = (sign_narrative.n_rows > 0) && (sign_narrative(0, 0) != 0);
  bool has_zero       = Z.n_elem > 0;
  bool has_sign_irf   = accu(abs(sign_irf)) > 0;
  bool has_sign_irf_cum    = accu(abs(sign_irf_cum)) > 0;
  bool has_sign_B     = accu(abs(sign_B)) > 0;
  bool has_elasticity = elasticity.n_rows > 0;
  bool has_elasticity_cum  = elasticity_cum.n_rows > 0;
  bool has_response_bounds = response_bounds.n_rows > 0;
  bool has_response_bounds_cum = response_bounds_cum.n_rows > 0;
  
  bool needs_cum_irf =
    has_sign_irf_cum ||
    has_elasticity_cum ||
    has_response_bounds_cum;
  
  mat    Q(N, N);
  
  cube   irf           = ir1_cpp(B, chol_Sigma, h, p);  // Cholesky structural irf
  cube irf_Q(irf.n_rows, irf.n_cols, irf.n_slices);
  
  // Cumulative Cholesky IRFs.
  cube irf_cum;
  if (needs_cum_irf) {
    irf_cum = cum_irf_cpp(irf);
  }
  
  
  field<mat> ZF;
  if (has_zero) {
    ZF = ZIRF(Z, irf.slice(0));
  }
  
  bool   success       = false;
  mat    shocks(N, T, fill::zeros);
  while (n_tries < max_tries && !success) {
    checkUserInterrupt();
    // Draw Q
    if (has_zero) {
      Q = rzeroQ(Z, ZF);
    } else {
      Q = rortho_cpp(N);
    }
    
    /// Horizon-specific IRF sign restrictions
    bool pass_irf = true;
    if (has_sign_irf) {
      pass_irf = match_sign_irf(Q, sign_irf, irf);
    }
    
    // Cumulative IRF sign restrictions.
    bool pass_irf_cum = true;
    if (has_sign_irf_cum) {
      pass_irf_cum = match_sign_irf(Q, sign_irf_cum, irf_cum);
    }
    
    // Traditional sign restrictions on contemporaneous matrix B
    bool pass_B = true;
    if (has_sign_B) {
      pass_B = match_sign(Q.t() * h_invp, sign_B);
    }
    
    // Horizon-specific elasticity restrictions
    bool pass_elasticity = true;
    if (has_elasticity) {
      pass_elasticity = match_elasticity(Q, irf, elasticity);
    }
    
    // Cumulative elasticity restrictions
    bool pass_elasticity_cum = true;
    if (has_elasticity_cum) {
      pass_elasticity_cum = match_elasticity(Q, irf_cum, elasticity_cum);
    }
    
    // Horizon-specific response bounds restrictions
    bool pass_response_bounds = true;
    if (has_response_bounds) {
      pass_response_bounds = match_response_bounds(Q, irf, response_bounds);
    }
    
    // Cumulative response bounds.
    bool pass_response_bounds_cum = true;
    if (has_response_bounds_cum) {
      pass_response_bounds_cum = match_response_bounds(Q, irf_cum, response_bounds_cum);
    }
    
    bool pass_traditional =
      pass_irf &&
      pass_irf_cum &&
      pass_B &&
      pass_elasticity &&
      pass_elasticity_cum &&
      pass_response_bounds &&
      pass_response_bounds_cum;
    
    if (pass_traditional) {
      shocks = Q.t() * h_invp * (Y - X * B).t();
      if (!has_narrative) {
        success = true;
      } else {
        for (uword t = 0; t < irf.n_slices; ++t) {
          irf_Q.slice(t) = irf.slice(t) * Q;
        }
        success = match_sign_narrative(shocks, sign_narrative, irf_Q);
      }
    }
    
    n_tries++;
  }
  
  double log_w_out = 0;
  if (!success) {
    log_w_out = -arma::datum::inf;
  } else {
    double log_w = 0;
    // Narrative restrictions: likelihood truncation correction
    if (has_narrative) {
      log_w += log_weight_narrative(T, sign_narrative, irf_Q);  
    }
    // Zero restrictions: volume-element correction
    if (has_zero) {
      log_w += log_weight_zero(Z, B, h_invp, Q, W);
    }
    log_w_out = log_w;
  }
  
  field<mat> result(3);
  result(0) = Q;
  result(1) = shocks;
  result(2) = log_w_out;
  
  return result;
}
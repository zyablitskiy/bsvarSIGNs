
#include <RcppArmadillo.h>

#include "utils.h"
#include "compute.h"

using namespace Rcpp;
using namespace arma;


// If matches narrative sign restrictions
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
bool match_sign_narrative(
    const arma::mat& Epsilon,
    const arma::mat& sign_narrative,
    const arma::cube& irf
) {
  int N = Epsilon.n_rows;
  
  for (int k = 0; k < (int)sign_narrative.n_rows; k++) {
    
    // sign_narrative columns:
    // 0: type
    // 1: sign
    // 2: var_i
    // 3: shock_j, unused for type C
    // 4: start
    // 5: horizon
    // 6 to 6+N-1: group membership dummies for type C
    
    int type        = (int)sign_narrative(k, 0);
    bool is_greater = sign_narrative(k, 1) == 1;
    int var_i       = (int)sign_narrative(k, 2);
    int shock_j     = (int)sign_narrative(k, 3) - 1;
    int t           = (int)sign_narrative(k, 4);
    int h           = (int)sign_narrative(k, 5);
    
    // ------------------------------------------------------------
    // Type 1: restriction on structural shock sign
    // ------------------------------------------------------------
    if (type == 1) {
      
      for (int col = t - 1; col < t + h; ++col) {
        if (is_greater) {
          if (Epsilon(shock_j, col) <= 0.0) return false;
        } else {
          if (Epsilon(shock_j, col) >= 0.0) return false;
        }
      }
      
      continue;
    }
    
    // For all historical-decomposition restrictions,
    // compute HD contributions for variable var_i.
    arma::mat hd = hd1_cpp(var_i, t, h, Epsilon, irf);
    
    // ------------------------------------------------------------
    // Type 4: group HD restriction
    //
    // sign = 1:
    //   sum_{j in group} |H_j| > sum_{j not in group} |H_j|
    //
    // sign = -1:
    //   sum_{j in group} |H_j| < sum_{j not in group} |H_j|
    // ------------------------------------------------------------
    if (type == 4) {
      
      if ((int)sign_narrative.n_cols < 6 + N) {
        return false;
      }
      
      arma::vec group_flag(N);
      int n_group = 0;
      
      for (int jj = 0; jj < N; ++jj) {
        group_flag(jj) = sign_narrative(k, 6 + jj);
        if (group_flag(jj) > 0.5) {
          n_group++;
        }
      }
      
      // Empty group is invalid.
      if (n_group == 0) {
        return false;
      }
      
      // If all shocks are in the group, the complement is empty.
      // For sign = 1, the restriction is trivially satisfied
      // except for the zero-probability equality case.
      // For sign = -1, it is impossible.
      if (n_group == N) {
        if (!is_greater) {
          return false;
        } else {
          continue;
        }
      }
      
      for (int col = 0; col <= h; ++col) {
        
        double sum_group = 0.0;
        double sum_other = 0.0;
        
        for (int jj = 0; jj < N; ++jj) {
          
          double val = std::abs(hd(jj, col));
          
          if (group_flag(jj) > 0.5) {
            sum_group += val;
          } else {
            sum_other += val;
          }
        }
        
        if (is_greater) {
          if (sum_group <= sum_other) return false;
        } else {
          if (sum_group >= sum_other) return false;
        }
      }
      
      continue;
    }
    
    // ------------------------------------------------------------
    // Existing Type 2 and Type 3 HD restrictions
    // ------------------------------------------------------------
    if (type != 2 && type != 3) {
      return false;
    }
    
    for (int col = 0; col <= h; ++col) {
      
      double val_j = std::abs(hd(shock_j, col));
      
      if (type == 2) {
        
        // Type A: compare against max/min of other individual shocks.
        
        if (is_greater) {
          
          double max_other = -1.0;
          
          for (int r = 0; r < (int)hd.n_rows; ++r) {
            if (r == shock_j) continue;
            
            double val_r = std::abs(hd(r, col));
            if (val_r > max_other) max_other = val_r;
          }
          
          if (val_j <= max_other) return false;
          
        } else {
          
          double min_other = arma::datum::inf;
          
          for (int r = 0; r < (int)hd.n_rows; ++r) {
            if (r == shock_j) continue;
            
            double val_r = std::abs(hd(r, col));
            if (val_r < min_other) min_other = val_r;
          }
          
          if (val_j >= min_other) return false;
        }
        
      } else if (type == 3) {
        
        // Type B: compare against sum of all other shocks.
        
        double sum_other = 0.0;
        
        for (int r = 0; r < (int)hd.n_rows; ++r) {
          if (r == shock_j) continue;
          sum_other += std::abs(hd(r, col));
        }
        
        if (is_greater) {
          if (val_j <= sum_other) return false;
        } else {
          if (val_j >= sum_other) return false;
        }
      }
    }
  }
  
  return true;
}


// approximate importance weight for narrative restrictions
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
double log_weight_narrative(
    const int&        T,
    arma::mat         sign_narrative,
    const arma::cube& irf
) {
  const int K = (int)sign_narrative.n_rows;
  const int N = (int)irf.n_rows;
  
  // ---- (1) are there any historical-decomposition restrictions? ----
  bool has_hd = false;
  for (int k = 0; k < K; k++) {
    if ((int)sign_narrative(k, 0) != 1) { has_hd = true; break; }
  }
  
  // ---- (2) shock-sign restrictions only: omega is analytic, = 2^{-n_cells} ----
  // The indicator does not depend on (B, Sigma, Q), so this weight is the SAME
  // for every draw and the resampling step is a no-op.
  if (!has_hd) {
    arma::mat cell(N, T, arma::fill::zeros);
    int n_cells = 0;
    for (int k = 0; k < K; k++) {
      int    j = (int)sign_narrative(k, 3) - 1;
      int    t = (int)sign_narrative(k, 4);
      int    h = (int)sign_narrative(k, 5);
      double s = (sign_narrative(k, 1) == 1) ? 1.0 : -1.0;
      for (int col = t - 1; col < t + h; ++col) {
        if (col < 0 || col >= T) return -arma::datum::inf;   // date out of sample
        if (cell(j, col) == 0.0)     { cell(j, col) = s; n_cells++; }
        else if (cell(j, col) != s)  return -arma::datum::inf; // contradictory signs
      }
    }
    return n_cells * std::log(2.0);
  }
  
  // ---- (3) HD restrictions present: Monte Carlo over the RESTRICTED periods only ----
  // ADRR eq. (16): the integral collapses to eps_v, the shocks that enter the
  // restrictions. Columns outside that set can stay at zero.
  arma::uvec flag(T, arma::fill::zeros);
  for (int k = 0; k < K; k++) {
    int t = (int)sign_narrative(k, 4);
    int h = (int)sign_narrative(k, 5);
    for (int col = t - 1; col < t + h; ++col) {
      if (col < 0 || col >= T) return -arma::datum::inf;
      flag(col) = 1;
    }
  }
  arma::uvec cols = arma::find(flag == 1);
  
  arma::mat Z_m(N, T, arma::fill::zeros);
  
  const int    M_block  = 10000;     // ADRR fn. 4: 1e3 often enough, 1e6 for many events
  const int    M_max    = 1000000;
  const double min_succ = 50.0;      // stop once 1/omega is stably estimated
  
  int    M         = 0;
  double n_success = 0.0;
  
  while (M < M_max && n_success < min_succ) {
    for (int m = 0; m < M_block; m++) {
      Z_m.cols(cols) = arma::randn<arma::mat>(N, cols.n_elem);
      if (match_sign_narrative(Z_m, sign_narrative, irf)) n_success++;
    }
    M += M_block;
  }
  
  if (n_success == 0.0) return -arma::datum::inf;
  
  return std::log((double)M / n_success);
}
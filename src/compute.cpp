
#include <RcppArmadillo.h>
#include "progress.hpp"
#include "bsvars.h"

using namespace Rcpp;
using namespace arma;

// cpp wrappers for bsvars functions

// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube bsvarSIGNs_structural_shocks (
    arma::cube&     posterior_B,    // (N, N, S)
    arma::cube&     posterior_A,    // (N, K, S)
    arma::cube&     posterior_Theta0, // (N, N, S)
    arma::mat&      Y,              // NxT dependent variables
    arma::mat&      X,               // KxT dependent variables
    const bool      standardise = false,
    arma::uvec      standardise_idx = 0
    
) {

  const int       N = Y.n_rows;
  const int       T = Y.n_cols;
  const int       S = posterior_B.n_slices;

  cube            structural_shocks(N, T, S);
  
  for (int s=0; s<S; s++) {
    // Raw shocks: u = B * epsilon = B * (Y - A*X), where B = Theta0^{-1} 
    structural_shocks.slice(s)    = posterior_B.slice(s) * (Y - posterior_A.slice(s) * X);
    
    // Standardize shocks if requested (to match standardized IRFs)
    if (standardise) {
      mat Theta0_s = posterior_Theta0.slice(s);
      if (standardise_idx.n_elem == 1 && standardise_idx(0) == 0) {
        structural_shocks.slice(s) = diagmat(abs(diagvec(Theta0_s))) * structural_shocks.slice(s);
      } else {
        mat std_mat = zeros(N, N);
        for (int i=0; i<N; i++) {
          std_mat(i, i) = abs(Theta0_s(standardise_idx(i), i));
        }
        structural_shocks.slice(s) = std_mat * structural_shocks.slice(s);
      }
    }
  } // END s loop

  return structural_shocks;
} // END bsvarSIGNs_structural_shocks



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube bsvarSIGNs_fitted_values (
    arma::cube&     posterior_A,        // NxKxS
    arma::cube&     posterior_B,        // NxNxS
    arma::cube&     posterior_sigma,    // NxTxS
    arma::mat&      X                   // KxT
) {
  cube fitted_values = bsvars::bsvars_fitted_values (posterior_A, posterior_B, posterior_sigma, X);
  return fitted_values;
} // END bsvarSIGNs_fitted_values




// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube ir1_cpp (
    const arma::mat& At,           // KxN
    const arma::mat& Theta0,      // NxN
    int              horizon,
    const int&       p
) {
  
  horizon++;
  int N         = Theta0.n_cols;
  
  cube irf      = zeros(N, N, horizon);
  irf.slice(0)  = Theta0;
  
  cube A(N, N, p);
  for(int j = 0; j < p; j++) {
    A.slice(j) = At.rows(j * N, (j + 1) * N - 1).t();
  }
  
  for(int t = 1; t < horizon; t++) {
    int max_lag = std::min(p, t);
    for(int j = 1; j <= max_lag; j++) {
      irf.slice(t) += A.slice(j - 1) * irf.slice(t - j);
    } // END j loop
  } // END t loop
  
  return irf;
} // END ir1_cpp


// Compute cumulative IRFs.
// If irf contains horizons 0, 1, ..., H,
// the returned cube contains:
// slice 0: IRF_0
// slice 1: IRF_0 + IRF_1
// slice h: sum_{l=0}^h IRF_l
//
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube cum_irf_cpp(const arma::cube& irf) {
  arma::cube cum_irf = irf;
  
  for (arma::uword t = 1; t < cum_irf.n_slices; ++t) {
    cum_irf.slice(t) += cum_irf.slice(t - 1);
  }
  
  return cum_irf;
}


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::field<arma::cube> bsvarSIGNs_ir ( // K=N*p+1
    arma::cube&   posterior_At,        // (K19, N3, S1000)
    arma::cube&   posterior_Theta0,   // (N3, N3, S1000)
    const int     horizon,
    const int     p,
    const bool    standardise = false,
    arma::uvec    standardise_idx = 0
) {
  
  const int       N = posterior_At.n_cols;
  const int       S = posterior_At.n_slices;
  
  cube            aux_irfs(N, N, horizon + 1);
  field<cube>     irfs(S);
  
  for (int s=0; s<S; s++) {
    mat irf_0           = posterior_Theta0.slice(s);
    if (standardise) {
      if (standardise_idx.n_elem == 1 && standardise_idx(0) == 0) {
        irf_0             = irf_0 * diagmat(pow(abs(diagvec(irf_0)), -1));
      } else {
        mat std_mat = zeros(N, N);
        for (int i=0; i<N; i++) {
          std_mat(i, i) = pow(abs(irf_0(standardise_idx(i), i)), -1);
        }
        irf_0 = irf_0 * std_mat;
      }
      
    }
    aux_irfs            = ir1_cpp( posterior_At.slice(s), irf_0, horizon, p );
    irfs(s)             = aux_irfs;
  } // END s loop
  
  return irfs;
} // END bsvarSIGNs_ir



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::field<arma::cube> bsvarSIGNs_hd ( // N - number of variables, S - number of draws, T - time length
    arma::field<arma::cube>&    posterior_irf_T,    // output of bsvars_irf with irfs at T horizons. NxNxT(horizon)xS
    arma::cube&                 structural_shocks,  // NxTxS output bsvars_structural_shocks
    const bool                  show_progress = true
) {
  
  const int       N = structural_shocks.n_rows;
  const int       T = structural_shocks.n_cols;
  const int       S = structural_shocks.n_slices;
  
  // Progress bar setup
  vec prog_rep_points = arma::round(arma::linspace(0, S, 50));
  if (show_progress) {
    Rcout << "**************************************************|" << endl;
    Rcout << "bsvars: Bayesian Structural Vector Autoregressions|" << endl;
    Rcout << "**************************************************|" << endl;
    Rcout << " Computing historical decomposition               |" << endl;
    Rcout << "**************************************************|" << endl;
    Rcout << " This might take a little while :)                " << endl;
    Rcout << "**************************************************|" << endl;
  }
  Progress p(50, show_progress);
  
  
  field<cube>     hds(S);
  mat             posterior_irf_t(N, N);
  
  for (int s=0; s<S; s++) {
    cube aux_hds(N, N, T, fill::zeros);
    
    // Increment progress bar
    if (any(prog_rep_points == s)) p.increment();
    // Check for user interrupts
    if (s % 200 == 0) checkUserInterrupt();
    
    for (int t=0; t<T; t++) {
      
      // cube        hds_at_t(N, N, t + 1);
      // for (int i=0; i<t; i++) {
      //   posterior_irf_t   = posterior_irf_T(s).slice(t - i - 1);
      //   hds_at_t.slice(i) = posterior_irf_t.each_col() % structural_shocks.slice(s).col(i);
      // } // END i loop
      for (int i=0; i<=t; i++) {
        posterior_irf_t   = posterior_irf_T(s).slice(i);
        aux_hds.slice(t) += posterior_irf_t.each_row() % structural_shocks.slice(s).col(t - i).t();
        // hds_at_t.slice(i) = posterior_irf_t.each_row() % structural_shocks.slice(s).col(t - i).t();
      } // END i loop
      // aux_hds.slice(t) = sum(hds_at_t, 2);
    } // END t loop
    
    hds(s)          = aux_hds;
  } // END s loop
  
  return hds;
} // END bsvarSIGNs_hd



// compute historical decomposition from t to t+h of the i-th variable
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat hd1_cpp(
    const int&        var_i,   // i-th variable
    const int&        t,       // start at period t
    const int&        h,       // number of horizons
    const arma::mat&  Epsilon, // structural shocks, NxT
    const arma::cube& irf
) {
  
  int ii = var_i - 1;
  int tt = t - 1;
  int N  = Epsilon.n_rows;
  
  mat hd(N, h + 1, arma::fill::zeros);
  
  // for each shock jj
  for(int jj = 0; jj < N; jj++) {
    // for each horizon hh
    for(int hh = 0; hh <= h; hh++) {
      // for each lag ll
      for(int ll = 0; ll <= hh; ll++) {
        hd(jj, hh) += irf(ii, jj, ll) * Epsilon(jj, tt + hh - ll);
      }
    }
  }
  return hd;
} // END hd1_cpp



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::field<arma::cube> bsvarSIGNs_fevd (
    arma::field<arma::cube>&    posterior_irf   // output of bsvars_irf
) {
  
  const int       N = posterior_irf(0).n_rows;
  const int       S = posterior_irf.n_rows;
  const int       horizon = posterior_irf(0).n_slices;
  
  field<cube>     fevds(S);
  
  for (int s=0; s<S; s++) {
    cube cum_fevds(N, N, horizon, fill::zeros);
    
    for (int h = 0; h < horizon; h++) {
      cum_fevds.slice(h) = square(posterior_irf(s).slice(h));
      if (h > 0) cum_fevds.slice(h) += cum_fevds.slice(h - 1);
    }
    
    for (int h = 0; h < horizon; h++) {
      // Normalize rows
      vec row_sums = sum(cum_fevds.slice(h), 1);
      for(int n = 0; n < N; n++) {
        if (row_sums(n) > 0) cum_fevds.slice(h).row(n) /= row_sums(n);
      }
    }
    
    cum_fevds *= 100;
    fevds(s) = cum_fevds;
    // Rcpp::Rcout << "Memory size aux_fevds: "
    //             << aux_fevds.n_elem * sizeof(double) / (1024*1024)
    //             << " MB" << std::endl;
  } // END s loop
  
  return fevds;
} // END bsvarSIGNs_fevd

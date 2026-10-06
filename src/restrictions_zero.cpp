
#include <RcppArmadillo.h>

#include "utils.h"
#include "compute.h"

using namespace Rcpp;
using namespace arma;

// Z_j * irf_0
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::field<arma::mat> ZIRF(
    const arma::field<arma::mat>& Z,
    const arma::mat& irf_0
) {
  arma::field<arma::mat> ZF(Z.n_elem);
  
  for (int j = 0; j < (int)Z.n_elem; j++) {
    if (Z(j).n_rows > 0) {
      ZF(j) = Z(j) * irf_0;
    } else {
      ZF(j) = arma::mat(0, irf_0.n_cols); // Skip math, return empty matrix
    }
  }
  
  return ZF;
}


// Zero restrictions
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::colvec zero_restrictions(
    const arma::field<arma::mat>& Z,
    const arma::colvec vec_structural
) {
  int N = Z(0).n_cols;
  
  int total_zeros = 0;
  for (int j = 0; j < (int)Z.n_elem; j++) {
    total_zeros += Z(j).n_rows;
  }
  
  arma::colvec z(total_zeros, arma::fill::none);
  if (total_zeros == 0) return z;
  
  arma::mat A0 = arma::reshape(vec_structural.rows(0, N * N - 1), N, N);
  
  // STABLE APPROACH:
  // Instead of inv(A0.t()), we solve A0.t() * X = I
  // This is mathematically identical but avoids squaring the condition number.
  arma::mat invA0t = arma::solve(A0.t(), arma::eye(N, N));
  
  int current_idx = 0;
  for (int j = 0; j < (int)Z.n_elem; j++) {
    int z_j = Z(j).n_rows;
    if (z_j > 0) {
      z.rows(current_idx, current_idx + z_j - 1) = Z(j) * invA0t.col(j);
      current_idx += z_j;
    }
  }
  
  return z;
}


// g ∘ f_h
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::colvec g_fh(
    const arma::field<arma::mat>& Z,
    const arma::mat& A0,
    const arma::mat& Aplus,
    const arma::field<arma::mat>& W
) {
  
  // FAST & STABLE: Solve B * A0 = Aplus 
  // Mathematically: B = Aplus * inv(A0)  =>  B * A0 = Aplus  =>  A0^T * B^T = Aplus^T
  // arma::solve solves X * Y = Z. So we solve A0^T * B^T = Aplus^T, then transpose the result.
  mat B = arma::solve(A0.t(), Aplus.t()).t();
  
  // S = A0 * A0^T is Symmetric Positive Definite (SPD)
  mat S = A0 * A0.t();
  
  // inv_sympd is highly optimized and stable specifically for SPD matrices.
  // It internally uses Cholesky decomposition, which is the mathematically correct way to invert SPD matrices.
  mat Sigma = arma::inv_sympd(S);
  mat irf_0 = arma::chol(Sigma, "lower");
  mat Q = irf_0.t() * A0;
  
  arma::field<arma::mat> ZF = ZIRF(Z, irf_0);
  
  colvec out = join_vert(vectorise(B), vectorise(Sigma));
  
  mat M_j, K_j;
  colvec w_j;
  int N = Q.n_cols;
  
  for (int j=0; j < (int)ZF.n_elem; j++) {
    if (j == 0) {
      M_j = ZF(j);
    } else {
      M_j = join_horiz(Q.cols(0, j-1), ZF(j).t()).t();
    }
    
    int s = W(j).n_rows;
    mat M_tilde_j = join_vert(M_j, W(j));
    
    mat K, R;
    qr(K, R, M_tilde_j.t());
    
    // Enforce strictly positive diagonals on R
    for(int i = N - s; i < N; i++) {
      if(R(i, i) < 0) {
        K.col(i) = -K.col(i);
      }
    }
    
    K_j = K.cols(N - s, N - 1);
    w_j = K_j.t() * Q.col(j);
    out = join_vert(out, w_j);
  }
  return out;
}


// g ∘ f_h with vectorized input
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::colvec g_fh_vec(
    const arma::field<arma::mat>& Z,
    const arma::colvec vec_structural,
    const arma::field<arma::mat>& W
) {
  int N = Z(0).n_cols;
  int M = vec_structural.n_rows;
  int K = (M - N*N)/N;
  mat A0 = reshape(vec_structural.rows(0, N*N-1), N, N);
  mat Aplus = reshape(vec_structural.rows(N*N, M-1), K, N);
  
  return g_fh(Z, A0, Aplus, W);
}


// log volume element
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
double log_volume_element(
    const arma::field<arma::mat>& Z,
    const arma::mat& A0,
    const arma::mat& Aplus,
    const arma::field<arma::mat>& W
) {
  colvec vec_structural = join_vert(vectorise(A0), vectorise(Aplus));
  
  mat Dz = Df([Z](const colvec& x) { return zero_restrictions(Z, x); },
              vec_structural);
  
  // Pass W into the Df lambda capture list so it stays constant during differentiation
  mat Dgf = Df([Z, W](const colvec& x) { return g_fh_vec(Z, x, W); },
               vec_structural);
  
  mat DN = Dgf * null(Dz);
  mat G  = arma::symmatu(DN.t() * DN);   // kill rounding asymmetry
  
  double log_det_G;
  if (!arma::log_det_sympd(log_det_G, G)) {
    // numerically singular: the draw is on a measure-zero set, give it weight 0
    return arma::datum::inf;
  }
  
  return 0.5 * log_det_G;
}


// importance weight for zero restrictions
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
double log_weight_zero(
    const arma::field<arma::mat>& Z,
    const arma::mat&              B,
    const arma::mat&              h_inv,
    const arma::mat&              Q,
    const arma::field<arma::mat>& W
) {
  int K     = B.n_rows;
  int N     = Q.n_cols;
  
  mat A0    = h_inv.t() * Q;
  mat Aplus = B * h_inv.t() * Q;
  
  double log_ve_f   = -(2*N+K+1) * log_det(A0).real();
  double log_ve_gfz = log_volume_element(Z, A0, Aplus, W);
  
  return log_ve_f - log_ve_gfz;
}


// draw Q conditional on zero restrictions
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat rzeroQ(
    const arma::field<arma::mat>& Z,
    const arma::field<arma::mat>& ZF
) {
  int N = Z.n_elem;
  
  // Pre-allocate Q to its final size to avoid memory thrashing
  arma::mat Q(N, N, arma::fill::none);
  
  for (int j = 0; j < (int)Z.n_elem; j++) {
    int z_j = Z(j).n_rows;
    
    // Draw and normalize isotropic normal vector
    arma::colvec x_j(N - j - z_j, arma::fill::randn);
    arma::colvec w_j = x_j / std::sqrt(arma::sum(arma::square(x_j)));
    
    arma::mat K_j;
    int m = j + z_j; // Number of rows in the constraint matrix
    
    if (m == 0) {
      // Fast path: first column and no zero restrictions
      K_j = arma::eye(N, N);
    } else {
      // Pre-allocate constraint matrix M_j perfectly
      arma::mat M_j(m, N, arma::fill::none);
      
      // Copy transposed columns of Q directly into M_j without join_horiz
      if (j > 0) {
        M_j.rows(0, j - 1) = Q.cols(0, j - 1).t();
      }
      // Copy zero restrictions
      if (z_j > 0) {
        M_j.rows(j, m - 1) = ZF(j);
      }
      
      // Fast Null Space via QR Decomposition (M_j^T = K * R)
      // The last (N - m) columns of K span the null space of M_j
      arma::mat K, R;
      arma::qr(K, R, M_j.t());
      K_j = K.cols(m, N - 1);
    }
    
    // Insert new vector directly into the pre-allocated Q
    Q.col(j) = K_j * w_j;
  }
  
  return Q;
}
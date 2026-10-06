
#include <functional>
#include <iostream>
#include <RcppArmadillo.h>

using namespace arma;


// QR decomposition, where the diagonal elements of R are positive
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat qr_sign_cpp(const arma::mat& A) {
  int N = A.n_rows;
  
  arma::mat Q(N, N), R(N, N);
  arma::qr_econ(Q, R, A);
  
  // Check and modify the diagonal elements of R
  for(arma::uword i = 0; i < R.n_cols; ++i) {
    if(R(i, i) < 0) {
      // R.row(i) *= -1;  // Change sign of the column
      Q.col(i) *= -1;  // Change sign of the corresponding row in Q
    }
  }
  
  return Q;
}


// Sample uniformly from the space of NxN orthogonal matrices
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat rortho_cpp(const int& N) {
  return qr_sign_cpp(arma::mat(N, N, fill::randn));
}


// If matches signs
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
bool match_sign(
    const arma::mat& A, 
    const arma::mat& sign
) {
  for(uword i = 0; i < sign.n_elem; ++i) {
    if (sign[i] > 0 && A[i] <= 0) return false;
    if (sign[i] < 0 && A[i] >= 0) return false;
  }
  return true;
}


// Numerical derivative of f: R^n -> R^m, at x (Rcpp::export is not possible)
// [[Rcpp:interface(cpp)]]
arma::mat Df(
    const std::function<arma::vec(const arma::vec&)>& f,
    const arma::vec& x,
    const double     h = 1e-5   // RELATIVE step; central differences give O(h^2)
)
{
  const int n = x.n_elem;
  const int m = f(x).n_elem;
  
  mat result(m, n);
  
  for (int i = 0; i < n; i++)
  {
    // scale the step to the magnitude of the coordinate; floor at 1 so that
    // near-zero entries of A0 / Aplus still get a usable perturbation
    const double h_i = h * std::max(1.0, std::abs(x(i)));
    
    vec x_p = x;  x_p(i) += h_i;
    vec x_m = x;  x_m(i) -= h_i;
    
    result.col(i) = (f(x_p) - f(x_m)) / (2.0 * h_i);
  }
  
  return result;
}

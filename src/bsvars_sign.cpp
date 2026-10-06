
#include <RcppArmadillo.h>
#include "progress.hpp"
#include "Rcpp/Rmath.h"
#include <bsvars.h>

#include "sample_hyper.h"
#include "sample_Q.h"
#include "sample_NIW.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::field<arma::mat> bsvar_sign_single_draw_cpp(
    const int &p,
    const arma::mat &Y,
    const arma::mat &X,
    const arma::cube &sign_irf,
    const arma::cube &sign_irf_cum,
    const arma::mat &sign_narrative,
    const arma::mat &sign_B,
    const arma::field<arma::mat> &Z,
    const arma::mat &elasticity,
    const arma::mat &elasticity_cum,
    const arma::mat &response_bounds,
    const arma::mat &response_bounds_cum,
    const Rcpp::List &prior,
    const arma::field<arma::mat>& W,
    const int &max_tries
) {
  const int T = Y.n_rows;
  const int N = Y.n_cols;
  const int K = X.n_cols;

  mat hypers = as<mat>(prior["hyper"]);
  int S_hyper  = hypers.n_cols - 1;
  int prior_nu = as<int>(prior["nu"]);
  int post_nu = prior_nu + T;
  int n_tries;

  double w, mu, delta, lambda, phi;

  vec hyper, psi;
  vec prior_v = as<mat>(prior["V"]).diag();

  mat B, Sigma, chol_Sigma, h_invp, Q, shocks;
  mat prior_V, prior_S, post_B, post_V, post_S;
  mat Ystar, Xstar, Yplus, Xplus;
  mat prior_B = as<mat>(prior["B"]);
  mat Ysoc = as<mat>(prior["Ysoc"]);
  mat Xsoc = as<mat>(prior["Xsoc"]);
  mat Ysur = as<mat>(prior["Ysur"]);
  mat Xsur = as<mat>(prior["Xsur"]);

  field<mat> result;

  hyper        = hypers.col(randi(distr_param(0, S_hyper)));
  mu           = hyper(0);
  delta        = hyper(1);
  lambda       = hyper(2);
  phi          = hyper(3);
  psi          = hyper.rows(4, N + 3);
  
  // update Minnesota prior
  vec v_all = prior_v % join_vert(lambda * lambda * repmat(1 / psi, p, 1),
                                  ones<vec>(K - N * p));
  uvec idx_dummy = as<uvec>(prior["idx_dummy"]);
  if (idx_dummy.n_elem > 0) {
    v_all.elem(idx_dummy).fill(1.0 / (phi * phi));   // Pandemic Priors
  }
  prior_V = diagmat(v_all);
  prior_S = diagmat(psi);

  // update dummy observation prior
  Ystar = join_vert(Ysoc / mu, Ysur / delta);
  Xstar = join_vert(Xsoc / mu, Xsur / delta);

  Yplus = join_vert(Ystar, Y);
  Xplus = join_vert(Xstar, X);

  // posterior parameters
  result = niw_cpp(Yplus, Xplus, prior_B, prior_V, prior_S, prior_nu);
  post_B = result(0);
  post_V = result(1);
  post_S = result(2);
  post_nu = as_scalar(result(3));
  mat post_V_chol_lower = result(4);

  double log_w = -arma::datum::inf;
  n_tries = 0;
  

  while (std::isinf(log_w) and (n_tries < max_tries or max_tries == 0)) {

    checkUserInterrupt();

    // sample reduced-form parameters
    Sigma = iwishrnd(post_S, post_nu);
    chol_Sigma = chol(Sigma, "lower");

    B = rmatnorm_cpp(post_B, post_V_chol_lower, chol_Sigma);

    h_invp = inv(trimatl(chol_Sigma)); // lower tri, h(Sigma) is upper tri

    result = sample_Q(p, Y, X, B, h_invp, chol_Sigma, prior,
                      sign_irf, sign_irf_cum, sign_narrative, sign_B, elasticity, elasticity_cum, response_bounds, response_bounds_cum, Z, 1, W);
    Q = result(0);
    shocks = result(1);
    double log_w_Q = as_scalar(result(2));

    if (!std::isinf(log_w_Q)) {
      log_w = log_w_Q;
    }
    n_tries++;
  }

  // w = std::exp(log_w);
  w = log_w; // Pass the log-weight directly, do NOT exponentiate here.

  arma::field<arma::mat> out(8);
  out(0) = arma::mat(1, 1, arma::fill::zeros); out(0)(0, 0) = w;
  out(1) = hyper;
  out(2) = B.t();
  out(3) = Q.t() * h_invp;
  out(4) = Q;
  out(5) = Sigma;
  out(6) = chol_Sigma * Q;
  out(7) = shocks;

  return out;
}

// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List bsvar_sign_par_cpp(
    const int &p,
    const arma::mat &Y,
    const arma::mat &X,
    const arma::cube &sign_irf,
    const arma::cube &sign_irf_cum,
    const arma::mat &sign_narrative,
    const arma::mat &sign_B,
    const arma::field<arma::mat> &Z,
    const arma::mat &elasticity,
    const arma::mat &elasticity_cum,
    const arma::mat &response_bounds,
    const arma::mat &response_bounds_cum,
    const Rcpp::List &prior,
    const arma::field<arma::mat>& W,
    const int &max_tries = 10000
) {
  arma::field<arma::mat> draw = bsvar_sign_single_draw_cpp(
    p, Y, X, sign_irf, sign_irf_cum, sign_narrative, sign_B, Z, elasticity, elasticity_cum, response_bounds, response_bounds_cum, prior, W, max_tries
  );

  return List::create(
      _["w"]      = as_scalar(draw(0)),
      _["hyper"]  = draw(1),
      _["A"]      = draw(2),
      _["B"]      = draw(3),
      _["Q"]      = draw(4),
      _["Sigma"]  = draw(5),
      _["Theta0"] = draw(6),
      _["shocks"] = draw(7)
  );
}

// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List bsvar_sign_cpp(
    const int&        S,
    const int&        p,
    const arma::mat&  Y,
    const arma::mat&  X,
    const arma::cube &sign_irf,
    const arma::cube &sign_irf_cum,
    const arma::mat &sign_narrative,
    const arma::mat &sign_B,
    const arma::field<arma::mat> &Z,
    const arma::mat &elasticity,
    const arma::mat &elasticity_cum,
    const arma::mat &response_bounds,
    const arma::mat &response_bounds_cum,
    const Rcpp::List& prior,
    const arma::field<arma::mat>& W,
    const bool        show_progress = true,
    const int&        max_tries = 10000
) {

  // Progress bar setup
  double num_threads = 1;
  vec prog_rep_points = arma::round(arma::linspace(0, S / num_threads, 50));
  if (show_progress) {
    Rcout << "**************************************************|" << endl;
    Rcout << " bsvarSIGNs: Bayesian Structural VAR with sign,   |" << endl;
    Rcout << "             zero and narrative restrictions      |" << endl;
    Rcout << "**************************************************|" << endl;
    Rcout << " Progress of simulation for " << S << " independent draws" << endl;
    Rcout << " Press Esc to interrupt the computations" << endl;
    Rcout << "**************************************************|" << endl;
  }
  Progress bar(50, show_progress);
  
  const int  T = Y.n_rows;
  const int  N = Y.n_cols;
  const int  K = X.n_cols;
  
  mat        hypers = as<mat>(prior["hyper"]);
  
  vec        posterior_w(S);
  mat        posterior_hyper(hypers.n_rows, S);
  cube       posterior_A(N, K, S);
  cube       posterior_B(N, N, S);
  cube       posterior_Q(N, N, S);
  cube       posterior_Sigma(N, N, S);
  cube       posterior_Theta0(N, N, S);
  cube       posterior_shocks(N, T, S);
  
  for (int s = 0; s < S; s++) {
    
    arma::field<arma::mat> draw = bsvar_sign_single_draw_cpp(
      p, Y, X, sign_irf, sign_irf_cum, sign_narrative, sign_B, Z, elasticity, elasticity_cum, response_bounds, response_bounds_cum, prior, W, max_tries
    );
    
    posterior_w(s)            = as_scalar(draw(0));
    posterior_hyper.col(s)    = draw(1);
    posterior_A.slice(s)      = draw(2);
    posterior_B.slice(s)      = draw(3);
    posterior_Q.slice(s)      = draw(4);
    posterior_Sigma.slice(s)  = draw(5);
    posterior_Theta0.slice(s) = draw(6);
    posterior_shocks.slice(s) = draw(7);
    
    // Increment progress bar
    if (any(prog_rep_points == s)) bar.increment();
    
  } // END s loop
  
  return List::create(
    _["posterior"]  = List::create(
      _["w"]        = posterior_w,
      _["hyper"]    = posterior_hyper,
      _["A"]        = posterior_A,
      _["B"]        = posterior_B,
      _["Q"]        = posterior_Q,
      _["Sigma"]    = posterior_Sigma,
      _["Theta0"]   = posterior_Theta0,
      _["shocks"]   = posterior_shocks
    )
  );
} // END bsvar_sign_cpp
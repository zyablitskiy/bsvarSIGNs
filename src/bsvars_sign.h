#ifndef _BSVARS_SIGN_H_
#define _BSVARS_SIGN_H_

#include <RcppArmadillo.h>

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
    const arma::field<arma::mat>& W
    const int &max_tries
);

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
);

Rcpp::List bsvar_sign_cpp(
    const int&        S,                  // number of draws from the posterior
    const int&        p,                  // number of lags
    const arma::mat&  Y,                  // TxN dependent variables
    const arma::mat&  X,                  // TxK dependent variables
    const arma::cube &sign_irf,
    const arma::cube &sign_irf_cum,
    const arma::mat &sign_narrative,
    const arma::mat &sign_B,
    const arma::field<arma::mat> &Z,
    const arma::mat &elasticity,
    const arma::mat &elasticity_cum,
    const arma::mat &response_bounds,
    const arma::mat &response_bounds_cum,
    const Rcpp::List& prior,              // a list of priors
    const arma::field<arma::mat>& W,
    const bool        show_progress = true,
    const int&        max_tries = 10000   // maximum tries for Q draw
);

#endif  // _BSVARS_SIGN_H_

#include <RcppArmadillo.h>
#include <cmath>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//[[Rcpp::export]]
List fast_wlm(const arma::mat& X, const arma::colvec& y, const arma::colvec& freqs) {
  int n = X.n_rows, k = X.n_cols;
  arma::mat w = arma::diagmat(freqs);
  arma::colvec coef = arma::solve(X.t()*w*X, X.t()*w*y);
  arma::colvec res  = y - X*coef;
  double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/(n - k);
  arma::colvec std_err = arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(X)*X)));
  return List::create(Named("coefficients") = coef, Named("sigma") = std_err);
}


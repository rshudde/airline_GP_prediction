#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat three_matrix_multiplication(const arma::mat& X, const arma::mat& Y, const arma::mat& Z)
{
  arma::mat XXt_return = X.t() * Y * Z;
  return(XXt_return);
}
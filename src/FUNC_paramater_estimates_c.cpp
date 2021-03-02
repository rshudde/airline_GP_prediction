#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

// [[Rcpp::export]]
arma::vec get_mu_c(const arma::mat y, const int n_datasets, const List g, const Rcpp::List V_mat,
                   const arma::uvec time_idx, const float sigma_2_mu, const float alpha_mu)
{
  // vector to be returned
  arma::vec mu = arma::ones<arma::vec>(n_datasets);
  
  for (int i = 0; i < n_datasets; i++)
  {
    NumericVector time (time_idx[i], 1);
    arma::mat V = V_mat[i]
    // float y_temp =y(i, time_idx[i]);
    // float sigma_2_mu_post_i = 1 / (time * V * time + 1 / sigma_2_mu);
    // float alpha_mu_post_i = sigma_2_mu_post_i * (alpha_mu / sigma_2_mu +
    //                                              (y_temp - g[i]).t() * arma::inv(V) * time );
    // alpha_mu_post_i = sqrt(alpha_mu_post_i);
    // mu[i] = R::rnorm(1, alpha_mu_post_i, sigma_2_mu_post_i);
  }
  
  return(mu);
}




// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/

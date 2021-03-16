// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include<cmath>

using namespace Rcpp;
using namespace arma;

// #include <RcppArmadillo.h>
// using namespace Rcpp;

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

typedef Rcpp::ListOf<Rcpp::NumericMatrix> MatList;
typedef Rcpp::ListOf<arma::vec> ListList;

// [[Rcpp::export]]
double numeric_to_arma( const Rcpp::NumericMatrix X ) { 
  return 0.0 ;
}

// // [[Rcpp :: export ()]
// arma::mat a3 (NumericMatrix x) {
//   arma::mat y = Rcpp::as<arma::mat>(x);
//   return(y);
// }

// [[Rcpp::export]]
arma::vec get_mu_c(const arma::mat& y, const int n_datasets, const List g, MatList V_mat,
                   const List time_idx, const float sigma_2_mu, const float alpha_mu)
{
  // vector to be returned - start with initializing a matrix of all 1
  arma::vec mu = arma::ones<arma::vec>(n_datasets);
  
  for (int i = 0; i < n_datasets; i++)
  {
    // NumericVector time (time_idx[i], 1);
    arma::vec t_temp = time_idx(i);
    int length = t_temp.n_elem;
    arma::vec time = arma::ones<arma::vec>(length);

    // // extract matrix we are inverting
    Rcpp::NumericMatrix V_temp = V_mat[i];
    arma::mat V = Rcpp::as<arma::mat>(wrap(V_temp));

    // // get sigma
    // arma::uvec y_temp = y(i, time_idx[i]);

    // this is getting 1V^(-1)1
    arma::mat temp1 = (time.t() * arma::inv(V));
    arma::mat tVt = temp1 * time;

    arma::mat sigma_2_mu_post_temp = 1 / (tVt + 1 / sigma_2_mu);
    float sigma_2_mu_post_i = sigma_2_mu_post_temp(0,0);

    // // get alpha
    // get y - g
    arma::uvec indices(as<arma::uvec>(wrap(time_idx(i))));
    // subtract 1 from indices
    for (int k = 0; k < indices.n_elem; k++)
    {
      indices[k] = indices[k] - 1;
    }
    
    arma::mat term_one = y.cols(indices); // get time_idx indices of the matrix
    arma::rowvec term_two = term_one.row(i); // get the specific row

    arma::uvec g_temp = g(i); // get g[i]

    arma::vec term_four = term_two.t() - g_temp; // y - g

    // // term_five is (y-g)^t %*% inv(V)
    arma::rowvec term_five = term_four.t() * arma::inv(V);
    arma::mat term_six = term_five * time;
    float term_seven = term_six(0,0);

    float alpha_mu_post_i = sigma_2_mu_post_i*(alpha_mu/sigma_2_mu + term_seven);
    mu[i] = sigma_2_mu_post_i * (arma::randn()) + alpha_mu_post_i;

    // std::cout << i << std::endl;
  }

  return(mu);
}

// [[Rcpp::export]]
arma::mat get_H_matrix_c(arma::vec w, arma::vec knots, int n_Knots)
{
  arma::mat H(w.n_elem, n_Knots, fill::zeros); // number of rows and columns 
  for (int i = 0; i < w.n_elem; i++)
  {
    arma::vec temp = 1 - abs((w[i] - knots)*(n_Knots - 1));
    for (int j = 0; j < temp.n_elem; j++)
    {
      temp[j] = std::fmax(temp[j], 0.0);
    }
    H.row(i) = temp.t();
  }
  
  return(H);
}

// [[Rcpp::export]]
Rcpp::List psi_alpha_c(arma::vec alpha, const arma::mat& y, const int n_datasets,
                       const List time_idx, MatList data, arma::vec mu, arma::vec xi,
                       MatList V_mat, arma::vec knots, const int n_Knots)
{
  // step to get the new beta
  float denomonator = 0;
  for (int i = 0; i < alpha.n_elem; i++)
  {
    denomonator += std::pow(alpha[i], 2);
  }

  arma::vec beta = arma::zeros<arma::vec>(alpha.n_elem);
  for (int i = 0; i <alpha.n_elem; i++)
  {
    beta[i] = alpha[i]/denomonator;
  }

  Rcpp::List w(n_datasets);
  Rcpp::List H_mat(n_datasets);
  Rcpp::List g(n_datasets);
  
  int negloglhood = 0;

  for (int i = 0; i < n_datasets; i++)
  {
    // get data
    // // extract matrix we are inverting
    Rcpp::NumericMatrix data_temp = data[i];
    arma::mat data_i = Rcpp::as<arma::mat>(wrap(data_temp));
    w[i] = (data_i.t() * beta + 1)/2;
    H_mat[i] = get_H_matrix_c(w[i], knots, n_Knots);

    Rcpp::NumericMatrix h_temp = H_mat[i];
    arma::mat h_i = Rcpp::as<arma::mat>(wrap(h_temp));
    g[i] = h_i * xi;
  }
  
  return Rcpp::List::create( Rcpp::Named("negloglhood") = negloglhood/2, 
                             Rcpp::Named("beta") = beta,
                             Rcpp::Named("w") = w,
                             Rcpp::Named("g") = g,
                             Rcpp::Named("H_mat") = H_mat);
}

// [[Rcpp::export]]
Rcpp::List get_alpha_c(arma::vec alpha_0, const arma::mat& y, const int n_datasets,
                       const List time_idx, MatList data, int n_covariates, arma::vec mu, 
                       arma::vec xi, MatList V_mat, arma::vec knots, const int n_Knots,
                       float c_prior = 1000.0)
{
  // step one 
  float theta = (2*M_PI)* arma::randu();
  arma::vec alpha_prior = c_prior * (arma::randn(n_covariates));
  arma::vec alpha_proposed = std::cos(theta) * alpha_0 + std::sin(theta) * alpha_prior;
  
  // step two 
  float theta_min = theta - 2 * M_PI;
  float theta_max = theta;
  
  // old negative loglikelihood
  Rcpp::List psi_out_old = psi_alpha_c(alpha_0, y, n_datasets, time_idx, data,mu, xi, V_mat, knots, n_Knots);
    
  // new negative loglikelihood
  float acceptance = -10.0;
  if (alpha_proposed[0] <= 0.0)
  {
      acceptance = 0.0;
  } 
  else
  {
    Rcpp::List psi_out_new = psi_alpha_c(alpha_proposed, y, n_datasets, time_idx, data, mu, xi, V_mat, knots, n_Knots);
    
    // calculate new acceptance value
    float old_loglike = psi_out_old["psi_out_old$negloglhood"];
    float new_loglike = psi_out_new["psi_out_old$negloglhood"];
    float value_new = exp(old_loglike - new_loglike);
    acceptance = std::fmin(1.0, value_new);
  }

  
  // step 3
  float zeta = randu();
  
  // continuation of step 3 - don't return until we get something we accept 
  while (acceptance <= zeta)
  {
  // step a
    if (theta < 0.0)
    {
      theta_min = theta;
    } else {
      theta_max = theta;
    }
    
    // step b 
    float a = theta_max - theta_min;
    float b = theta_min;
    theta = a*randu() + b;
      
    // step c
    alpha_proposed = std::cos(theta) * alpha_0 + std::sin(theta) * alpha_prior;
        
    // new negative loglikelihood
    // reject if alpha[1] is negative because we want beta[1] to be positive
    // new negative loglikelihood
    if (alpha_proposed[0] <= 0.0)
    {
      acceptance = 0.0;
    } 
    else
    {
      Rcpp::List psi_out_new = psi_alpha_c(alpha_proposed, y, n_datasets, time_idx, data, mu, xi, V_mat, knots, n_Knots);
     
      // calculate new acceptance value
      float old_loglike = psi_out_old["psi_out_old$negloglhood"];
      float new_loglike = psi_out_new["psi_out_old$negloglhood"];
      float value_new = exp(old_loglike - new_loglike);
      acceptance = std::fmin(1.0, value_new);
    }

  }
  
  // for return purposes
  Rcpp::List psi_out_new = psi_alpha_c(alpha_proposed, y, n_datasets, time_idx, data, mu, xi, V_mat, knots, n_Knots);

  return Rcpp::List::create( Rcpp::Named("negloglhood") = psi_out_new["negloglhood"], 
                             Rcpp::Named("beta") = psi_out_new["beta"],
                             Rcpp::Named("w") = psi_out_new["w"],
                             Rcpp::Named("g") = psi_out_new["g"],
                             Rcpp::Named("H_mat") = psi_out_new["H_mat"],
                             Rcpp::Named("alpha") = alpha_proposed);
}
  


  
    

      

        
      

        

// Rcpp::sourceCpp("src/FUNC_paramater_estimates_c.cpp")


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/

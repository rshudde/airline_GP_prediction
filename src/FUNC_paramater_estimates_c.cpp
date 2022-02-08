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
arma::mat inv_chol(arma::mat A){
  int n = A.n_rows;
  double nr = A(0,0);
  // Normalizing the matrix (if required):
  if(nr!=1){
    A = A/nr;
  }
  arma::rowvec av = A.row(0);
  arma::colvec avec = av.t();
  arma::colvec r = avec.tail(n-1);
  arma::colvec y = zeros(n);
  arma::mat R = eye(n,n);
  // Initial values:
  y(0) = - r(0);
  double a = - r(0);
  double b = 1;
  // Updating b:
  b = (1 - (a*a))*b;
  R(0,1) = y(0);
  R.col(1) = R.col(1)/pow(b,0.5);
  // Updating other columns:
  //int k = 0;
  for(int k = 0; k <= (n-3); k++){
    arma::colvec rsub = r.head(k+1);
    arma::colvec rrev = reverse(rsub);
    arma::colvec ysub = y.head(k+1);
    arma::colvec yrev = reverse(ysub);
    // Updating a:
    a = -(r(k+1) + sum(rrev.t()*ysub))/b;
    // Updating y:
    y.head(k+1) = ysub + a*yrev; y(k+1) = a;
    // Updating b:
    b = (1 - (a*a))*b;
    // Updating the columns:
    arma::colvec ysub2 = y.head(k+2);
    arma::colvec yrev2 = reverse(ysub2);
    R.submat(0,k+2,k+1,k+2) = yrev2;
    R.col(k+2) = R.col(k+2)/pow(b,0.5);
  }
  R = R/pow(nr,0.5);
  return R;
}

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
    for (arma::uword k = 0; k < indices.n_elem; k++)
    {
      indices[k] = indices[k] - 1;
    }
    
    arma::mat term_one = y.cols(indices); // get time_idx indices of the matrix
    arma::rowvec term_two = term_one.row(i); // get the specific row

    arma::vec g_temp = g(i); // get g[i]

    arma::vec term_four = term_two.t() - g_temp; // y - g

    // // term_five is (y-g)^t %*% inv(V)
    arma::rowvec term_five = term_four.t() * arma::inv(V);
    arma::mat term_six = term_five * time;
    float term_seven = term_six(0,0);

    float alpha_mu_post_i = sigma_2_mu_post_i*(alpha_mu/sigma_2_mu + term_seven);

    mu[i] = std::sqrt(sigma_2_mu_post_i) * arma::randn() + alpha_mu_post_i;
    
    // std::cout << i << std::endl;
  }

  return(mu);
}

// [[Rcpp::export]]
arma::mat get_H_matrix_c(arma::vec w, arma::vec knots, int n_Knots)
{
  arma::mat H(w.n_elem, n_Knots, fill::zeros); // number of rows and columns 
  for (arma::uword i = 0; i < w.n_elem; i++)
  {
    arma::vec temp = 1 - abs((w[i] - knots)*(n_Knots - 1));
    for (arma::uword j = 0; j < temp.n_elem; j++)
    {
      temp[j] = std::fmax(temp[j], 0.0);
    }
    H.row(i) = temp.t();
  }
  
  return(H);
}

// [[Rcpp::export]]
Rcpp::List psi_alpha_c(arma::vec alpha, const arma::mat& y, const int n_datasets,
                       const List time_idx, Rcpp::List data, arma::vec mu, arma::vec xi,
                       MatList V_mat, arma::vec knots, const int n_Knots)
{
  // step to get the new beta
  float denomonator = 0.0;
  for (arma::uword i = 0; i < alpha.n_elem; i++)
  {
    denomonator += std::pow(alpha[i], 2);
  }

  arma::vec beta = arma::zeros<arma::vec>(alpha.n_elem);
  for (arma::uword i = 0; i < alpha.n_elem; i++)
  {
    beta[i] = alpha[i]/std::sqrt(denomonator); 
  }
  
  Rcpp::List w(n_datasets);
  Rcpp::List H_mat(n_datasets);
  Rcpp::List g(n_datasets);
  
  float negloglhood = 0.0;

  for (int i = 0; i < n_datasets; i++)
  {
    // get data
    // // extract matrix we are inverting
    arma::mat data_i = data[i];
    arma::vec ones(data_i.n_rows);
    ones.fill(1.0);

    w[i] = ((data_i * beta) + ones)/2;
    
    H_mat[i] = get_H_matrix_c(w[i], knots, n_Knots);

    Rcpp::NumericMatrix h_temp = H_mat[i];
    arma::mat h_i = Rcpp::as<arma::mat>(wrap(h_temp));
    g[i] = h_i * xi;
    
    // getting negative loglikelihood
    Rcpp::NumericMatrix V_temp = V_mat[i];
    arma::mat V_i = Rcpp::as<arma::mat>(wrap(V_temp));
    
    // get time indices
    arma::uvec indices(as<arma::uvec>(wrap(time_idx(i))));
    // subtract 1 from indices
    for (arma::uword k = 0; k < indices.n_elem; k++)
    {
      indices[k] = indices[k] - 1;
    }
    
    arma::mat y_row = y.cols(indices); // get time_idx indices of the matrix
    arma::rowvec y_temp = y_row.row(i); // get the specific row
    
    // get y - mu[i]- g[[i]]
    arma::vec mu_temp = arma::ones<arma::vec>(y_temp.n_elem) * mu[i];
    mu_temp.fill(mu[i]);
    arma::vec g_temp = g(i);
    arma::rowvec term_two = y_temp - mu_temp.t() - g_temp.t();
    
    arma::mat temp_lik = term_two * arma::inv(V_i) * term_two.t();
    
    negloglhood += temp_lik(0,0);
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
  float theta = (2.0*M_PI) * arma::randu();
  arma::vec alpha_prior = c_prior * (arma::randn(n_covariates));
  arma::vec alpha_proposed = std::cos(theta) * alpha_0 + std::sin(theta) * alpha_prior;
  
  // step two 
  float theta_min = theta - 2.0 * M_PI;
  float theta_max = theta;
  
  // old negative loglikelihood
  Rcpp::List psi_out_old = psi_alpha_c(alpha_0, y, n_datasets, time_idx, data, mu, xi, V_mat, knots, n_Knots);
  Rcpp::List psi_out_new = psi_alpha_c(alpha_0, y, n_datasets, time_idx, data, mu, xi, V_mat, knots, n_Knots);
  
  // new negative loglikelihood
  float acceptance = 0.0;
  if (alpha_proposed[0] <= 0.0)
  {
      acceptance = 0.0;
  } 
  else
  {
    psi_out_new = psi_alpha_c(alpha_proposed, y, n_datasets, time_idx, data, mu, xi, V_mat, knots, n_Knots);
    
    // calculate new acceptance value
    float old_loglike = psi_out_old["negloglhood"];
    float new_loglike = psi_out_new["negloglhood"];
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
      psi_out_new = psi_alpha_c(alpha_proposed, y, n_datasets, time_idx, data, mu, xi, V_mat, knots, n_Knots);
     
      // calculate new acceptance value
      float old_loglike = psi_out_old["negloglhood"];
      float new_loglike = psi_out_new["negloglhood"];
      float value_new = exp(old_loglike - new_loglike);
      acceptance = std::fmin(1.0, value_new);
    }

  }
  
  // for return purposes
  // Rcpp::List psi_out_new = psi_alpha_c(alpha_proposed, y, n_datasets, time_idx, data, mu, xi, V_mat, knots, n_Knots);

  return Rcpp::List::create( Rcpp::Named("negloglhood") = psi_out_new["negloglhood"], 
                             Rcpp::Named("beta") = psi_out_new["beta"],
                             Rcpp::Named("w") = psi_out_new["w"],
                             Rcpp::Named("g") = psi_out_new["g"],
                             Rcpp::Named("H_mat") = psi_out_new["H_mat"],
                             Rcpp::Named("alpha") = alpha_proposed);
}
  
// [[Rcpp::export]]
Rcpp::List rfunc() {
  Rcpp::Environment invgamma("package:invgamma");
  Rcpp::Function dinvgamma = invgamma["dinvgamma"];
  Rcpp::Function rinvgamma = invgamma["rinvgamma"];
  Rcpp::NumericVector tmp = rinvgamma(5, 1);
  Rcpp::NumericVector a = dinvgamma(tmp, 1);
  return Rcpp::List::create(Rcpp::Named("tmp") = tmp,
                            Rcpp::Named("a") = a);
}


Rcpp::NumericVector rinvgamma(R_xlen_t n,
                              double shape,
                              double rate = 1.0) {
  return 1.0/Rcpp::rgamma(n, shape, rate);
}

// [[Rcpp::export]]
Rcpp::List get_sigma_2_c(const float a, const float b, const arma::mat& y, const int n_datasets,
                                   const int n_nonNA_y, const List time_idx, arma::vec mu, MatList M_mat,
                                   const List g)
{
  float rate_term = 0;
  
  for (int i = 0; i < n_datasets; i++)
  {
    
    // extract matrix we are inverting
    Rcpp::NumericMatrix M_temp = M_mat[i];
    arma::mat M = Rcpp::as<arma::mat>(wrap(M_temp));
    
    // create diagonal to add to M
    arma::vec t_temp = time_idx(i);
    int length = t_temp.n_elem;
    arma::mat I = arma::eye(length, length);

    arma::mat term_one = arma::inv(M + I);

    // // get term two
    arma::uvec indices(as<arma::uvec>(wrap(time_idx(i))));
    // subtract 1 from indices
    for (arma::uword k = 0; k < indices.n_elem; k++)
    {
      indices[k] = indices[k] - 1;
    }
    arma::mat y_row = y.cols(indices); // get time_idx indices of the matrix
    arma::rowvec y_temp = y_row.row(i); // get the specific row

    // get y - mu[i]- g[[i]]
    arma::vec mu_temp = arma::ones<arma::vec>(length) * mu[i];
    arma::vec g_temp = g(i);
    arma::rowvec term_two = y_temp - mu_temp.t() - g_temp.t();
    
    // arma::rowvec term_three = term_two * term_one;
    arma::mat term_three = term_two * term_one * term_two.t();
    rate_term += term_three(0,0);
    
  }
  // Rcpp::NumericVector sigma_2_temp = rinvgamma(1, a + n_nonNA_y/2, b + rate_term/2);
  Rcpp::NumericVector sigma_2_temp = Rcpp::rgamma(1, a + n_nonNA_y/2, 1/(b + rate_term/2));
  float sigma_2 = sigma_2_temp[0];
  
  // float sigma_2 = 0;
    
  return Rcpp::List::create( Rcpp::Named("sigma_2") = 1/sigma_2, 
                             Rcpp::Named("rate") = rate_term);
}
  
// [[Rcpp::export]] 
float get_matern_values_c(const float l_k, const float r_mj)
{
  float term_one = 1 + (sqrt(5) * r_mj) / l_k + (5 * pow(r_mj, 2)) / (3 * pow(l_k, 2));
  float exponential = exp(-(sqrt(5) * r_mj) / l_k);
  return(term_one * exponential);
}

// [[Rcpp::export]] 
arma::mat dist_c(const arma::vec input)
{
  arma::mat C(input.n_elem, input.n_elem, fill::zeros);
  
  for (arma::uword i = 0; i < input.n_elem; i++)
  {
    for (arma::uword j = 0; j < input.n_elem; j++)
    {
      C(i, j) = abs(input[i] - input[j]);
    }
  }
  
  return(C);
}

// [[Rcpp::export]]
arma::mat get_matern_c(const float l_k, const arma::vec time_points)
{
  // get the distance matrix
  arma::mat distance_matrix = dist_c(time_points);

  // get full matrix
  arma::mat M_i(time_points.n_elem, time_points.n_elem);
  
  for (arma::uword i = 0; i < M_i.n_rows; i++)
  {
    for (arma::uword j = 0; j < M_i.n_cols; j++)
    {
      M_i(i,j) = get_matern_values_c(l_k, distance_matrix(i,j));
    }
  }
  
  // return M matfix
  return(M_i);
}


// [[Rcpp::export]]
arma::mat get_K_i_c(const float sigma_2, const arma::mat M_i)
{
  arma::mat K_i = sigma_2 * M_i;
  
  return(K_i);
}


// [[Rcpp::export]]
arma::mat get_V_i_c(const float sigma_2, const arma::mat K_i)
{

  arma::mat A(K_i.n_rows, K_i.n_rows, fill::eye);
  arma::mat V_i = K_i + sigma_2 * A;
  return(V_i);
}



// [[Rcpp::export]]
arma::vec vector_differences_c(const arma::vec y, const float m_i, const arma::vec g_i)
{
  int Ti = y.n_elem;
  arma::vec fill = arma::vec(Ti);
  fill.fill(m_i);

  arma::vec to_return = y - fill - g_i.t();

  return(to_return);
}

// [[Rcpp::export]]
float lk_acceptance_c(const arma::mat& y, const arma::vec mu, const List g, const float sigma_2,
                    float lk_prime, float l_k, const List time)
{
  
  float to_return = 0.0;
  float ratio = 1.0;
  
  // indicator function part
  if (lk_prime < 0.1 || lk_prime > 1.0 || l_k < 0.1 ||  l_k > 1.0)
  {
    to_return = 0.0;
  }
  else
  {
    // calcualte first term outside of the product
    arma::rowvec y_noNA = y.row(0); // get the specific row
    arma::vec fill = arma::vec(y_noNA.n_elem);
    fill.fill(mu[0]);
    arma::vec g_fill = g[0];
    arma::rowvec term_one = y_noNA - fill.t() - g_fill.t();

    // get the two v terms necessary
    arma::vec time_temp = time[0];
    arma::mat M_temp = get_matern_c(l_k, time_temp);
    arma::mat M_prime = get_matern_c(lk_prime, time_temp);

    arma::mat K_temp = get_K_i_c(sigma_2, M_temp);
    arma::mat K_prime = get_K_i_c(sigma_2, M_prime);

    arma::mat V_temp = get_V_i_c(sigma_2, K_temp);
    arma::mat V_prime = get_V_i_c(sigma_2, K_prime);

    arma::mat term_two = arma::inv(V_prime) - arma::inv(V_temp);

    // do matrix multiplication
    arma::mat matrix_part1 = term_one * term_two;
    arma::mat matrix_part = matrix_part1.t() * term_one;

    // calculate ratio
    ratio = exp(-0.5 * matrix_part(0,0));
    
    // for loop goes here
    for (arma::uword i = 1; i < y.n_rows; i++)
    {
      arma::rowvec y_noNA = y.row(i); // get the specific row
      arma::vec fill = arma::vec(y_noNA.n_elem);
      fill.fill(mu[i]);
      arma::vec g_fill = g[i];
      arma::rowvec term_one = y_noNA - fill.t() - g_fill.t();

      // get the two v terms necessary
      arma::vec time_temp = time[i];
      arma::mat M_temp = get_matern_c(l_k, time_temp);
      arma::mat M_prime = get_matern_c(lk_prime, time_temp);

      arma::mat K_temp = get_K_i_c(sigma_2, M_temp);
      arma::mat K_prime = get_K_i_c(sigma_2, M_prime);

      arma::mat V_temp = get_V_i_c(sigma_2, K_temp);
      arma::mat V_prime = get_V_i_c(sigma_2, K_prime);

      arma::mat term_two = arma::inv(V_prime) - arma::inv(V_temp);

      // do matrix multiplication
      arma::mat matrix_part1 = term_one * term_two;
      arma::mat matrix_part = matrix_part1.t() * term_one;

      ratio = ratio * exp(-0.5 * matrix_part(0,0));
    }
    
    to_return = std::fmin(1.0, (lk_prime / l_k) * ratio);
  }
  
  return(to_return);
}

// [[Rcpp::export]]
float lb_acceptance_c(const arma::mat y, const float lb, const float lb_prime, const arma::vec xi, const arma::vec knots)
{
  float to_return = 0.0;
  if (lb_prime < 0.1 || lb_prime > 1.0 || lb < 0.1 ||  lb > 1.0)
  {
    to_return = 0.0;
  }
  else 
  {
    //  set term_one = xi 
    arma::vec term_one = xi;
    
    // # stuff we need to calcualte v_i
    arma::mat M_lb = get_matern_c(lb, knots);
    arma::mat M_lb_prime = get_matern_c(lb_prime, knots);
      
    // use tinv to invert M here
    arma::mat term_two = (M_lb.i()) - (M_lb_prime.i());
    
    
    // # matrix multiplication
    arma::mat matrix_part1 = term_one.t() * term_two;
    arma::mat matrix_part = matrix_part1 * term_one;
    float ratio = std::exp(-0.5 * matrix_part(0,0));
    
    to_return = std::fmin(1.0, (lb_prime / lb) * ratio);
  }
  
  return(to_return);
}

// [[Rcpp::export]]
float get_sigmaB_2_c(const float a, const float b, const arma::vec xi, const float lb,
                     const arma::vec knots, const int n_Knots)
{
  // calculating shape and rate of inv gamma
  arma::mat M = get_matern_c(lb, knots);
  arma::mat rate_term_temp = xi.t() * M.i() * xi;
  float rate_term = rate_term_temp(0, 0)/2;
  
  // do inverse gamma draw
  Rcpp::NumericVector sigmaB_2_temp = Rcpp::rgamma(1, a + n_Knots/2, 1/(b + rate_term));
  float sigmaB_2 = 1/sigmaB_2_temp[0];
  
  return(sigmaB_2);
}
       
// [[Rcpp::export]]
Rcpp::List psi_xi_c(const arma::vec xi, const arma::mat y, const int n_datasets, 
                    const List time_idx, const arma::vec mu, const List H_mat, 
                    const List V_mat)     
{
  List g(n_datasets);
  float negloglhood = 0.0;
  
  for (int i = 0; i < n_datasets; i++)
  {
    // // extract matrix we are inverting
    Rcpp::NumericMatrix H_temp = H_mat[i];
    arma::mat H = Rcpp::as<arma::mat>(wrap(H_temp));
    g(i) = H * xi;
    
    // get the loglikelihood
    Rcpp::NumericMatrix V_temp = V_mat[i];
    arma::mat V = Rcpp::as<arma::mat>(wrap(V_temp));

    // get the y term 
    arma::uvec indices(as<arma::uvec>(wrap(time_idx(i))));
    // subtract 1 from indices
    for (arma::uword k = 0; k < indices.n_elem; k++)
    {
      indices[k] = indices[k] - 1;
    }
    
    arma::mat y_row = y.cols(indices); // get time_idx indices of the matrix
    arma::rowvec y_temp = y_row.row(i); // get the specific row
    
    // get y - mu[i]- g[[i]]
    arma::vec mu_temp = arma::ones<arma::vec>(y_temp.n_elem) * mu[i];
    arma::vec g_temp = g(i);
    arma::rowvec term_two = y_temp - mu_temp.t() - g_temp.t();
    
    arma::mat temp_lik = term_two * arma::inv(V) * term_two.t();
    
    negloglhood += temp_lik(0,0)/2;
  }
  
  return Rcpp::List::create( Rcpp::Named("negloglhood") = negloglhood, 
                             Rcpp::Named("g") = g);
}
      
// [[Rcpp::export]]
Rcpp::List get_xi_backend_c(const arma::vec xi_0, const float sigmaB_2, const arma::mat y, const int n_datasets, 
                    const List time_idx, const arma::vec mu, const List H_mat, const List V_mat, 
                    const float lb, const arma::vec knots, const arma::vec xi_prior)
{
  // stuff before the while loop, assume wood and chan stuff is already passed in as xi_proposed
  
  // step one
  float theta = (2*M_PI)* arma::randu();
  arma::vec xi_proposed = std::cos(theta) * xi_0 + std::sin(theta) * xi_prior;
  
  // step two
  float theta_min = theta - 2 * M_PI;
  float theta_max = theta;
  
  // old negative loglikelihood
  Rcpp::List psi_out_old = psi_xi_c(xi_0, y, n_datasets, time_idx, mu, H_mat, V_mat);
  
  // new negative loglikelihood
  Rcpp::List psi_out_new = psi_xi_c(xi_proposed, y, n_datasets, time_idx, mu, H_mat, V_mat);
  
  // calculate new acceptance value

  // calculate new acceptance value
  float old_loglike = psi_out_old["negloglhood"];
  float new_loglike = psi_out_new["negloglhood"];
  float value_new = exp(old_loglike - new_loglike);
  float acceptance = std::fmin(1.0, value_new);
  
  // step 3
  float zeta = arma::randu();
  
  // continuation of step 3 - don't return until we get something we accept 
  while (acceptance <= zeta)
  {
    // step a
    if (theta < 0)
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
    xi_proposed = std::cos(theta) * xi_0 + std::sin(theta) * xi_prior;
    
    // new negative loglikelihood
    Rcpp::List psi_out_new = psi_xi_c(xi_proposed, y, n_datasets, time_idx, mu,H_mat, V_mat);
    
    // calculate new acceptance value
    new_loglike = psi_out_new["negloglhood"];
    value_new = exp(old_loglike - new_loglike);
    acceptance = std::fmin(1.0, value_new);  
  }
  
  psi_out_new = psi_xi_c(xi_proposed, y, n_datasets, time_idx, mu, H_mat, V_mat);
  
  return Rcpp::List::create(  Rcpp::Named("negloglhood") = psi_out_new["negloglhood"],
                              Rcpp::Named("xi") = xi_proposed,
                             Rcpp::Named("g") = psi_out_new["g"]);  
}

// [[Rcpp::export]]
Rcpp::List inner_gibbs_one(float sigma_2_post, Rcpp::List M_gibbs, int n_datasets)
{
  Rcpp::List K_gibbs(n_datasets);
  Rcpp::List V_gibbs(n_datasets);
  
  for (int i = 0; i < n_datasets; i++)
  {
    // extract matrix we are inverting
    Rcpp::NumericMatrix M_temp = M_gibbs[i];
    arma::mat M = Rcpp::as<arma::mat>(wrap(M_temp));
    
    arma:: mat K_temp = get_K_i_c(sigma_2_post, M);
    arma:: mat V_temp = get_V_i_c(sigma_2_post, K_temp);
    
    K_gibbs[i] = K_temp;
    V_gibbs[i] = V_temp;
  }
  return Rcpp::List::create(  Rcpp::Named("K_gibbs") = K_gibbs,
                              Rcpp::Named("V_gibbs") = V_gibbs);
    
}



// Rcpp::sourceCpp("src/FUNC_paramater_estimates_c.cpp")


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/


// [[Rcpp::export]]
float blah(float mean, float sd)
{
  float b = std::sqrt(sd)*arma::randn() + mean;
  return b;
}

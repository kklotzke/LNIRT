// [[Rcpp::depends(RcppArmadillo)]] 

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;




// MCMC functions used in LNIRT 

//'@export
// [[Rcpp::export]]
arma::mat Rcpp_DrawS_LNIRT(const arma::vec &alpha0, const arma::vec &beta0, const arma::vec &guess0, 
                           const arma::vec &theta0, const arma::mat &Y, const bool alt = false) {
  const int N = Y.n_rows;
  const int K = Y.n_cols;
  
 
  if(alt) {
    arma::mat a0(N, K);
    arma::mat b0(N, K);
    arma::mat g0(N, K);
    arma::mat t0(N, K);
    
    for (int i = 0; i < N; i++) {
      a0.row(i) = alpha0.t();
      b0.row(i) = beta0.t();
      g0.row(i) = guess0.t();
    }
    
    for (int i = 0; i < K; i++) {
      t0.col(i) = theta0;
    }
    
    NumericVector tmp = as<NumericVector>(wrap(a0 % t0 - b0));
    arma::mat eta(pnorm(tmp));
    eta.reshape(N, K);
    arma::mat probS = eta / (eta + g0 % (1 - eta));
    probS.reshape(N * K, 1);
    
    NumericVector randS(runif(N * K));
    NumericVector res = wrap(ifelse(randS > as<NumericVector>(wrap(probS)), 0, 1));

    arma::mat S(res);
    S.reshape(N, K);

    return(S % Y);
  }
  
  else {
    arma::mat S(N, K);
    
    for (int i = 0; i < N; i++) {
      //a0.row(i) = alpha0;
      //b0.row(i) = beta0;
      arma::vec tmp = (alpha0 * theta0(i) - beta0);
      arma::vec eta(pnorm(as<NumericVector>(wrap(tmp)), 0.0, 1.0, 1, 0)); // eta, mu, sd, lower.tail, log.p
      arma::vec probS = eta / (eta + guess0 % (1 - eta));
      NumericVector randS(runif(K));
      NumericVector res = wrap(ifelse(randS > as<NumericVector>(wrap(probS)), 0, 1));
      
      arma::vec resultS(res);
      S.row(i) = resultS.t();
      
    }
    
    return(S % Y);
  }
  
  // tmp = pnorm(a0 % t0 - b0) // element-wise multiplication and substraction
  
  //eta <- t(matrix(alpha0, ncol = N, nrow = K)) * matrix(theta0, ncol = K, nrow = N) - t(matrix(beta0, nrow = K, ncol = N))
  //eta <- matrix(pnorm(eta), ncol = K, nrow = N)
  //probS <- eta/(eta + t(matrix(guess0, nrow = K, ncol = N)) * (matrix(1, ncol = K, nrow = N) - eta))
  //S <- matrix(runif(N * K), ncol = K, nrow = N)
  //S <- matrix(ifelse(S > probS, 0, 1), ncol = K, nrow = N)
  //S <- S * Y

}


//'@export
// [[Rcpp::export]]
arma::fvec Rcpp_DrawC_LNIRT(const arma::mat &S, const arma::mat &Y) {
  
  const int N = Y.n_rows;
  const int K = Y.n_cols;

  arma::Col<int> Q1(K);
  arma::Col<int> Q2(K);
  arma::fvec guess(K);
  Q1.fill(20); // starting value 20
  Q2.fill(80); // starting value 80
  
  for (int i = 0; i < K; i++) {
    for (int j = 0; j < N; j++) {

      if(S(j, i) == 0) { // answer unknown
        Q2(i)++;
        if(Y(j, i) == 1) { // answer unknown and guess correctly
          Q1(i)++;
        }
      }
    }
    guess(i) = rbeta(1, Q1(i), Q2(i))(0);
  }
  
  return(guess);
}


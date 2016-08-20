// [[Rcpp::depends(RcppArmadillo)]] 

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;




// MCMC functions used in LNIRT 

//'@export
// [[Rcpp::export]]
arma::mat Rcpp_DrawS_LNIRT(const arma::vec &alpha0, const arma::vec &beta0, const arma::vec &guess0, 
                           const arma::vec &theta0, const arma::mat &Y) {
  const int N = Y.n_rows;
  const int K = Y.n_cols;
  arma::mat S(N, K);
  
  arma::mat randS(runif(N * K));
  randS.reshape(N, K);
  
  for (int i = 0; i < N; i++) {
    arma::vec tmp = (alpha0 * theta0(i) - beta0);
    arma::vec eta(pnorm(as<NumericVector>(wrap(tmp)), 0.0, 1.0, 1, 0)); // eta, mu, sd, lower.tail, log.p
    arma::vec probS = eta / (eta + guess0 % (1 - eta));
    arma::vec res(K);
    res.fill(1);
    for(int j = 0; j < K; j++) {
      if(randS(i, j) > probS(j))
        res(j) = 0;
    }
    
    S.row(i) = res.t();
  }
  
  return(S % Y);
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


// [[Rcpp::depends(RcppArmadillo)]] 

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

//'@export
// [[Rcpp::export]]
arma::fvec Rcpp_DrawC_LNIRT(const arma::mat &S, arma::mat &Y) {
  
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


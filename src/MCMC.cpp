// [[Rcpp::depends(RcppArmadillo)]] 

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;




// MCMC functions used in LNIRT 

//'@export
// [[Rcpp::export]]
List Rcpp_DrawS_LNIRT(const arma::vec &alpha0, const arma::vec &beta0, const arma::vec &guess0, 
                           const arma::vec &theta0, const arma::mat &Y) {
  const int N = Y.n_rows;
  const int K = Y.n_cols;
  arma::mat S(N, K);
  arma::mat etamat(N, K);
  
  arma::mat randS(runif(N * K));
  randS.reshape(N, K);
  
  for (int i = 0; i < N; i++) {
    arma::vec eta = (alpha0 * theta0(i) - beta0);
    arma::vec probEta(pnorm(as<NumericVector>(wrap(eta)), 0.0, 1.0, 1, 0)); // eta, mu, sd, lower.tail, log.p
    arma::vec probS = probEta / (probEta + guess0 % (1 - probEta));
    arma::vec res(K);
    res.fill(1);
    for (int j = 0; j < K; j++) {
      if(randS(i, j) > probS(j))
        res(j) = 0;
    }
    
    S.row(i) = res.t();
    etamat.row(i) = eta.t();
  }
  
  List ret; ret["S"] = S % Y; ret["eta"] = etamat;
  return(ret);
}


//'@export
// [[Rcpp::export]]
arma::mat Rcpp_DrawZ_LNIRT(const arma::vec &alpha0, const arma::vec &beta0, const arma::vec &theta0, 
                           const arma::mat &S, const arma::mat &D, const arma::mat &eta, const bool PNO = false) {
  const int N = S.n_rows;
  const int K = S.n_cols;
  
  arma::mat Z(N, K);
  arma::mat randU(runif(N * K));
  randU.reshape(N, K);

  for (int i = 0; i < N; i++) {
    arma::vec personEta(K);
    if (PNO)
      personEta = eta.row(i).t();
    else
      personEta = (alpha0 * theta0(i) - beta0);
    
    arma::vec BB(pnorm(as<NumericVector>(wrap(-personEta)), 0.0, 1.0, 1, 0)); // eta, mu, sd, lower.tail, log.p
    BB.elem(find(BB < 1e-05)).fill(1e-05);
    BB.elem(find(BB > (1 - 1e-05))).fill(1 - 1e-05);
    arma::vec tt = (BB % (1 - S.row(i).t()) + (1 - BB) % S.row(i).t()) % randU.row(i).t() + BB % S.row(i).t();
    
    for (int j = 0; j < K; j++) {
      if (D(i, j) == 1) {
        arma::vec tmp(qnorm(as<NumericVector>(wrap(tt(j))), 0.0, 1.0, 1, 0));
        Z(i, j) = tmp(0) + personEta(j);
      }
      else if (D(i, j) == 0) {
        arma::vec tmp(rnorm(1, 0.0, 1.0));
        Z(i, j) = tmp(0) + personEta(j);
      }
    }
  }
  return(Z);
}



//'@export
// [[Rcpp::export]]
arma::vec Rcpp_DrawTheta_LNIRT(const arma::vec &alpha0, const arma::vec &beta0, const arma::mat &Z,
                               const arma::vec &mu, const double sigma) {
  const int N = Z.n_rows;
  const int K = Z.n_cols;
  
  arma::mat b0(N, K);
  for (int i = 0; i < N; i++) {
    b0.row(i) = beta0.t();
  }
  
  arma::vec tmp = sum(pow(alpha0, 2));
  const double pvar = tmp(0) + 1 / sigma;
  const double sqrtpvar = sqrt(1/pvar);
  arma::vec thetahat = (Z + b0) * alpha0;
  arma::vec mutmp = (thetahat + mu / sigma) / pvar;
  arma::vec theta(N);
  for (int i = 0; i < N; i++) {
    arma::vec tmp = rnorm(1, mutmp(i), sqrtpvar);
    theta(i) = tmp(0);
  }
  
  return (theta);
}



//'@export
// [[Rcpp::export]]
arma::vec Rcpp_DrawZeta_LNIRT(const arma::mat &RT, const arma::vec &phi, const arma::vec &lambda,
                              const arma::vec &sigma2, const arma::vec &mu, const double sigmaz) {
  const int N = RT.n_rows;
  const int K = RT.n_cols;
  
  arma::mat l0(K, N);
  for (int i = 0; i < N; i++) {
    l0.col(i) = lambda;
  }
  
  arma::mat Z = l0 - RT.t();
  arma::mat sigma2inv = arma::diagmat(1 / sigma2);
  arma::mat ee(phi.t() * sigma2inv);
  arma::mat tmp = 1 / (ee * phi +  1 / sigmaz);
  
  const double vartheta = tmp(0, 0);
  const double sqrtvartheta = sqrt(vartheta); 
  arma::mat meantheta((ee * Z + (mu.t() / sigmaz)) * vartheta);
  arma::vec zeta(N);
  for (int i = 0; i < N; i++) {
    arma::vec tmp2 = rnorm(1, meantheta(i), sqrtvartheta);
    zeta(i) = tmp2(0);
  }
  
  return (zeta);
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


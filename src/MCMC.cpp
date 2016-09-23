// [[Rcpp::depends(RcppArmadillo)]] 

#include "RcppArmadillo.h"
#include "MCMC.h"
using namespace Rcpp;


// Common functions used in MCMC


//'@export
// [[Rcpp::export]]
arma::mat Rcpp_SimulateRT(const arma::mat &RT, const arma::vec &zeta, const arma::vec &lambda, 
                          const arma::vec &phi, const arma::vec &sigma2, const arma::mat DT) {
  const int N = RT.n_rows;
  const int K = RT.n_cols;
  arma::mat RT0 = RT;
  arma::mat DT0 = DT;
  
  arma::mat meanT = arma::repmat(lambda.t(), N, 1) - arma::repmat(phi.t(), N, 1) % arma::repmat(zeta, 1, K);
  meanT.reshape(N * K, 1);
  arma::mat sigmaL = arma::repmat(sqrt(sigma2.t()), N, 1);
  sigmaL.reshape(N * K, 1);
  RT0.reshape(N * K, 1);
  DT0.reshape(N * K, 1);
  
  for (int i = 0; i < N * K; i ++) {
   if (DT0(i, 0) == 0) {
     RT0(i, 0) = R::rnorm(meanT(i, 0), sigmaL(i, 0)); 
   }
  }
  RT0.reshape(N, K);
  
  return(RT0);
}


//'@export
// [[Rcpp::export]]
arma::mat Rcpp_SimulateY(const arma::mat &Y, const arma::vec &theta, const arma::vec &alpha0, 
                         const arma::vec &beta0, const arma::vec &guess0, const arma::mat &D) {
  const int N = Y.n_rows;
  const int K = Y.n_cols;
  arma::mat Y0 = Y;
  
  arma::mat G(N, K);
  G.fill(0);
  for (int kk = 0; kk < K; kk++) {
    G.col(kk) = arma::vec(rbinom(N, 1, guess0(kk)));
  }

  arma::mat rand(runif(N * K));
  rand.reshape(N, K);  
  arma::mat par = theta * alpha0.t() - arma::repmat(beta0.t(), N, 1);
  for (int i = 0; i < N; i++) {
    for (int kk = 0; kk < K; kk++) {
      double prob = R::pnorm(par(i, kk), 0.0, 1.0, 1, 0); // eta, mu, sd, lower.tail, log.p
      double tmp = 0;
      if (rand(i, kk) < prob) {
        tmp = 1;
      }
      
      if (D(i, kk) == 0) {
        if (G(i, kk) == 1) { 
          Y0(i, kk) = 1; // missing: guessed correctly
        }
        else {
          Y0(i, kk) = tmp; // missing: response generated
        }
      }
      
    }
  }
  
  return(Y0);
}

  
  
//'@export
// [[Rcpp::export]]
arma::vec Rcpp_DrawZeta(const arma::mat &RT, const arma::vec &phi, const arma::vec &lambda,
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
    zeta(i) = rnorm(1, meantheta(i), sqrtvartheta)(0);
  }
  
  return (zeta);
}



//'@export
// [[Rcpp::export]]
List Rcpp_SampleB(const arma::mat &Y, const arma::vec &X, const arma::mat &Sigma, const arma::vec &B0, 
                  const arma::mat &V0) {
  const int k = Y.n_cols;
  
  arma::mat tmp = 1/(X.t() * X);
  arma::mat Bvar = inv(inv(Sigma * tmp(0, 0)) + inv(V0));
  arma::vec Btilde = Bvar * (inv(Sigma) * (X.t() * Y).t() + inv(V0) * B0);
  arma::vec randN(rnorm(k, 0.0, 1.0));
  arma::vec B = Btilde + chol(Bvar) * randN;
  arma::mat pred = X * B.t();
  
  List ret; ret["B"] = B; ret["pred"] = pred;
  return(ret);
}



//'@export
// [[Rcpp::export]] 
arma::mat Rcpp_rwishart(const int nu, const arma::mat V) {
  const int m = V.n_rows;
  arma::vec df = (nu + nu - m + 1) - arma::regspace(nu - m + 1, nu);
  arma::mat RT;
  if (m > 1) {
    arma::vec tmp(m);
    for (int i = 0; i < m; i++) {
      tmp(i) = sqrt(rchisq(1, df(i))(0));
    }
    RT = arma::diagmat(tmp); 
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < m; j++) {
        if (i > j) {
          RT(i, j) = rnorm(1, 0.0, 1.0)(0);
        }
      }
    }
  }
  else {
    RT = arma::mat(1, 1);
    RT(0,0) = sqrt(rchisq(1, df(0))(0));
  }
  
  arma::mat U = chol(V);
  arma::mat C = arma::trimatu(RT.t() * U);
  arma::mat tmp(m, m);
  tmp.fill(0);
  tmp.diag() += 1;
  arma::mat CI = solve(C, tmp);
  arma::mat IW = CI * CI.t();
  
  return (IW);
}




// MCMC functions used in LNRT

//'@export
// [[Rcpp::export]] 
Rcpp::List Rcpp_DrawLambdaPhi_LNRT(const arma::mat &RT, const arma::vec &theta, const arma::vec &sigma2, 
                             const arma::mat muI, const arma::mat sigmaI, const arma::vec ingroup) {
  const int N = RT.n_rows;
  const int K = RT.n_cols;
  
  arma::mat invSigmaI = inv(sigmaI);
  arma::mat H(N, 2);
  H.col(0) = -theta % ingroup;
  H.col(1) = ingroup;
  
  arma::mat tmp1 = arma::diagmat(1 / sigma2);
  arma::mat tmp2 = H.t() * H;
  arma::vec ones(K);
  ones.fill(1);
  arma::mat diagones = arma::diagmat(ones);
  arma::mat varest = inv(kron(tmp1, tmp2) + kron(diagones, invSigmaI)); 
  
  tmp1 = H.t() * RT;
  tmp2 = arma::repmat(sigma2.t(), 2, 1);
  arma::mat tmp3 = arma::repmat((muI.row(0) * invSigmaI).t(), 1, K);
  arma::mat meanest = (tmp1 / tmp2 + tmp3).t();
  
  tmp1 = arma::repmat(meanest, 1, K);
  tmp2 = kron(diagones, arma::rowvec(2, arma::fill::ones));
  arma::rowvec meanest0 = sum((tmp1 * varest) % tmp2);

  tmp1 = arma::randn(1, 2 * K);
  arma::rowvec lambdaphi = arma::repmat(meanest0, 1, 1) + tmp1 * arma::chol(varest);
  arma::vec phi(K);
  arma::vec lambda(K);
  
  for (int i = 0; i <= (2 * K - 2); i = i + 2) {
    phi(i / 2) = lambdaphi(i);
    lambda(i / 2) = lambdaphi(i + 1);
    if (phi(i / 2) < 0.3)
      phi(i / 2) = 0.3;
  }
  
  List ret; ret["phi"] = phi; ret["lambda"] = lambda; 
  return(ret);
}



//'@export
// [[Rcpp::export]]
arma::vec Rcpp_SampleS_LNRT(const arma::mat &RT, const arma::vec &zeta, const arma::vec &lambda, 
                            const arma::vec &phi, const arma::vec &ingroup) {
  const int N = RT.n_rows;
  const int K = RT.n_cols;
  const int Nr = sum(ingroup);
  const int ss0 = 10;
  
  arma::mat Z0(N, K);
  arma::vec colsum(K);
  colsum.fill(0);
  for (int i = 0; i < N; i++) {
    Z0.row(i) = (RT.row(i) + phi.t() * zeta(i) - lambda.t()) * ingroup(i);
    Z0.row(i) %= Z0.row(i);
    
    for (int j = 0; j < K; j++) {
      colsum(j) += Z0(i, j); 
    }
  }
  
  arma::vec randChisq(rchisq(K, Nr));
  arma::vec sigma2 = (colsum + ss0) / randChisq;
  
  return (sigma2);
}



//'@export
// [[Rcpp::export]]
arma::vec Rcpp_DrawLambda_LNRT(const arma::mat &RT, const arma::vec &zeta, const arma::vec &sigma2, 
                               const arma::vec &mu, const double sigma) {
  const int N = RT.n_rows;
  const int K = RT.n_cols;
  
  arma::mat RT0(N, K);
  for (int i = 0; i < K; i++) {
    RT0.col(i) = RT.col(i) + zeta;
  }
  
  arma::mat XX(N, K);
  XX.fill(1);
  arma::mat tmp = (XX.t() * XX);
  arma::vec pvar = tmp.diag() / sigma2 + 1 / sigma;
  arma::vec sqrtpvar = sqrt(1 / pvar);
  tmp = (XX.t() * RT);
  arma::vec betahat = tmp.diag() / sigma2;
  arma::vec mu0 = (betahat + mu / sigma) / pvar;
  arma::vec beta(K);
  for (int i = 0; i < K; i++) {
    beta(i) = rnorm(1, mu0(i), sqrtpvar(i))(0);
  }
  
  return (beta);
}




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
                           const arma::mat &S, const arma::mat &D, const arma::mat &eta, const bool PNO) {
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
  arma::vec mu0 = (thetahat + mu / sigma) / pvar;
  arma::vec theta(N);
  for (int i = 0; i < N; i++) {
    theta(i) = rnorm(1, mu0(i), sqrtpvar)(0);
  }
  
  return (theta);
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




//'@export
// [[Rcpp::export]]
arma::vec Rcpp_DrawBeta_LNIRT(const arma::vec &theta, const arma::vec &alpha, const arma::mat &Z,
                              const arma::vec &mu, const double sigma) {
  const int N = Z.n_rows;
  const int K = Z.n_cols;
  
  arma::mat Z0(N, K);
  arma::mat a0(N, K);
  for (int i = 0; i < N; i++) {
    a0.row(i) = -alpha.t();
    Z0.row(i) = Z.row(i) - a0.row(i) * theta(i);
  }
  
  arma::mat tmp = (a0.t() * a0);
  arma::vec pvar = tmp.diag() + 1 / sigma;
  arma::vec sqrtpvar = sqrt(1 / pvar);
  tmp = (a0.t() * Z0);
  arma::vec betahat = tmp.diag();
  arma::vec mu0 = (betahat + mu / sigma) / pvar;
  arma::vec beta(K);
  for (int i = 0; i < K; i++) {
    beta(i) = rnorm(1, mu0(i), sqrtpvar(i))(0);
  }
  
  return (beta);
}



//'@export
// [[Rcpp::export]]
arma::vec Rcpp_DrawLambda_LNIRT(const arma::mat &RT, const arma::vec &phi, const arma::vec &zeta, const arma::vec &sigma2, 
                                const arma::vec &mu, const double sigmal) {
  const int N = RT.n_rows;
  const int K = RT.n_cols;
  
  arma::mat Z0(N, K);
  for (int i = 0; i < N; i++) {
    Z0.row(i) = RT.row(i) + phi.t() * zeta(i);
  }
  arma::vec X(N);
  X.fill(1);
  
  arma::mat sigma2diag = arma::diagmat(sigma2);
  arma::mat sigmaldiag = arma::mat(K, K);
  sigmaldiag.fill(0);
  sigmaldiag.diag() += sigmal;
  
  List sampleB = Rcpp_SampleB(Z0, X, sigma2diag, mu, sigmaldiag);
  arma::vec lambda = as<arma::vec>(sampleB["B"]);

  return (lambda);
}



//'@export
// [[Rcpp::export]]
arma::vec Rcpp_DrawAlpha_LNIRT(const arma::vec &theta, const arma::vec &beta, const arma::mat &Z,
                              const arma::vec &mu, const double sigma) {
  const int N = Z.n_rows;
  const int K = Z.n_cols;
  
  arma::mat b0(N, K);
  arma::mat XX(N, K);
  for (int i = 0; i < N; i++) {
    b0.row(i) = beta.t();
    XX.row(i) = theta(i) - b0.row(i);
  }
  
  arma::mat tmp = (XX.t() * XX);
  arma::vec pvar = tmp.diag() + 1 / sigma;
  arma::vec sqrtpvar = sqrt(1 / pvar);
  tmp = (XX.t() * Z);
  arma::vec alphahat = tmp.diag();
  arma::vec mu0 = (alphahat + mu / sigma) / pvar;
  arma::vec alpha(K);
  for (int i = 0; i < K; i++) {
    alpha(i) = rnorm(1, mu0(i), sqrtpvar(i))(0);
  }
  
  return (alpha);
}



//'@export
// [[Rcpp::export]]
arma::vec Rcpp_DrawPhi_LNIRT(const arma::mat &RT, const arma::vec &lambda, const arma::vec &zeta, const arma::vec &sigma2, 
                                const arma::vec &mu, const double sigmal) {
  const int N = RT.n_rows;
  const int K = RT.n_cols;
  
  arma::mat Z0(N, K);
  for (int i = 0; i < N; i++) {
    Z0.row(i) = -RT.row(i) + lambda.t();
  }
  arma::mat sigma2diag = arma::diagmat(sigma2);
  arma::mat sigmaldiag = arma::mat(K, K);
  sigmaldiag.fill(0);
  sigmaldiag.diag() += sigmal;
  
  List sampleB = Rcpp_SampleB(Z0, zeta, sigma2diag, mu, sigmaldiag);
  arma::vec phi = as<arma::vec>(sampleB["B"]);
  
  return (phi);
}



//'@export
// [[Rcpp::export]]
arma::vec Rcpp_SampleS2_LNIRT(const arma::mat &RT, const arma::vec &zeta, const arma::vec &lambda, 
                              const arma::vec &phi) {
  const int N = RT.n_rows;
  const int K = RT.n_cols;
  const int ss0 = 10;
  
  arma::mat Z0(N, K);
  arma::vec colsum(K);
  colsum.fill(0);
  for (int i = 0; i < N; i++) {
    Z0.row(i) = RT.row(i) + phi.t() * zeta(i) - lambda.t();
    Z0.row(i) %= Z0.row(i); // % Z0.row(i);
    
    for (int j = 0; j < K; j++) {
      colsum(j) += Z0(i, j); 
    }
  }
 
 arma::vec randChisq(rchisq(K, N));
 arma::vec sigma2 = (colsum + ss0) / randChisq;

  return (sigma2);
}







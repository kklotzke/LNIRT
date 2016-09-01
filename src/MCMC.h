#include "RcppArmadillo.h"

// Common functions used in MCMC


arma::vec Rcpp_DrawZeta(const arma::mat &RT, const arma::vec &phi, const arma::vec &lambda,
                        const arma::vec &sigma2, const arma::vec &mu, const double sigmaz);

Rcpp::List Rcpp_SampleB(const arma::mat &Y, const arma::vec &X, const arma::mat &Sigma, const arma::vec &mu, 
                        const arma::mat &V0);

arma::mat Rcpp_rwishart(const int nu, const arma::mat V);


// MCMC functions used in LNRT

Rcpp::List Rcpp_DrawLambdaPhi_LNRT(const arma::mat &RT, const arma::vec &theta, const arma::vec &sigma2, 
                             const arma::mat muI, const arma::mat sigmaI, const arma::vec ingroup);


arma::vec Rcpp_SampleS_LNRT(const arma::mat &RT, const arma::vec &zeta, const arma::vec &lambda, const arma::vec &phi, 
                            const arma::vec &ingroup);


arma::vec Rcpp_DrawLambda_LNRT(const arma::mat &RT, const arma::vec &zeta, const arma::vec &sigma2, 
                               const arma::vec &mu, const double sigma);


// MCMC functions used in LNIRT 

Rcpp::List Rcpp_DrawS_LNIRT(const arma::vec &alpha0, const arma::vec &beta0, const arma::vec &guess0, 
                            const arma::vec &theta0, const arma::mat &Y);

arma::mat Rcpp_DrawZ_LNIRT(const arma::vec &alpha0, const arma::vec &beta0, const arma::vec &theta0, 
                           const arma::mat &S, const arma::mat &D, const arma::mat &eta, const bool PNO = false);

arma::vec Rcpp_DrawTheta_LNIRT(const arma::vec &alpha0, const arma::vec &beta0, const arma::mat &Z,
                               const arma::vec &mu, const double sigma);

arma::fvec Rcpp_DrawC_LNIRT(const arma::mat &S, const arma::mat &Y);

arma::vec Rcpp_DrawBeta_LNIRT(const arma::vec &theta, const arma::vec &alpha, const arma::mat &Z,
                              const arma::vec &mu, const double sigma);

arma::vec Rcpp_DrawLambda_LNIRT(const arma::mat &RT, const arma::vec &phi, const arma::vec &zeta, const arma::vec &sigma2, 
                                const arma::vec &mu, const double sigmal);

arma::vec Rcpp_DrawAlpha_LNIRT(const arma::vec &theta, const arma::vec &beta, const arma::mat &Z,
                               const arma::vec &mu, const double sigma);

arma::vec Rcpp_DrawPhi_LNIRT(const arma::mat &RT, const arma::vec &lambda, const arma::vec &zeta, const arma::vec &sigma2, 
                             const arma::vec &mu, const double sigmal);

arma::vec Rcpp_SampleS2_LNIRT(const arma::mat &RT, const arma::vec &zeta, const arma::vec &lambda, const arma::vec &phi);


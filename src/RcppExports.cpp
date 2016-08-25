// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Rcpp_DrawS_LNIRT
List Rcpp_DrawS_LNIRT(const arma::vec& alpha0, const arma::vec& beta0, const arma::vec& guess0, const arma::vec& theta0, const arma::mat& Y);
RcppExport SEXP LNIRT_Rcpp_DrawS_LNIRT(SEXP alpha0SEXP, SEXP beta0SEXP, SEXP guess0SEXP, SEXP theta0SEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::vec& >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type guess0(guess0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta0(theta0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    __result = Rcpp::wrap(Rcpp_DrawS_LNIRT(alpha0, beta0, guess0, theta0, Y));
    return __result;
END_RCPP
}
// Rcpp_DrawZ_LNIRT
arma::mat Rcpp_DrawZ_LNIRT(const arma::vec& alpha0, const arma::vec& beta0, const arma::vec& theta0, const arma::mat& S, const arma::mat& D, const arma::mat& eta, const bool PNO);
RcppExport SEXP LNIRT_Rcpp_DrawZ_LNIRT(SEXP alpha0SEXP, SEXP beta0SEXP, SEXP theta0SEXP, SEXP SSEXP, SEXP DSEXP, SEXP etaSEXP, SEXP PNOSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::vec& >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta0(theta0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type D(DSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const bool >::type PNO(PNOSEXP);
    __result = Rcpp::wrap(Rcpp_DrawZ_LNIRT(alpha0, beta0, theta0, S, D, eta, PNO));
    return __result;
END_RCPP
}
// Rcpp_DrawTheta_LNIRT
arma::vec Rcpp_DrawTheta_LNIRT(const arma::vec& alpha0, const arma::vec& beta0, const arma::mat& Z, const arma::vec& mu, const double sigma);
RcppExport SEXP LNIRT_Rcpp_DrawTheta_LNIRT(SEXP alpha0SEXP, SEXP beta0SEXP, SEXP ZSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::vec& >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const double >::type sigma(sigmaSEXP);
    __result = Rcpp::wrap(Rcpp_DrawTheta_LNIRT(alpha0, beta0, Z, mu, sigma));
    return __result;
END_RCPP
}
// Rcpp_DrawZeta_LNIRT
arma::vec Rcpp_DrawZeta_LNIRT(const arma::mat& RT, const arma::vec& phi, const arma::vec& lambda, const arma::vec& sigma2, const arma::vec& mu, const double sigmaz);
RcppExport SEXP LNIRT_Rcpp_DrawZeta_LNIRT(SEXP RTSEXP, SEXP phiSEXP, SEXP lambdaSEXP, SEXP sigma2SEXP, SEXP muSEXP, SEXP sigmazSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type RT(RTSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const double >::type sigmaz(sigmazSEXP);
    __result = Rcpp::wrap(Rcpp_DrawZeta_LNIRT(RT, phi, lambda, sigma2, mu, sigmaz));
    return __result;
END_RCPP
}
// Rcpp_DrawC_LNIRT
arma::fvec Rcpp_DrawC_LNIRT(const arma::mat& S, const arma::mat& Y);
RcppExport SEXP LNIRT_Rcpp_DrawC_LNIRT(SEXP SSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    __result = Rcpp::wrap(Rcpp_DrawC_LNIRT(S, Y));
    return __result;
END_RCPP
}
// Rcpp_DrawBeta_LNIRT
arma::vec Rcpp_DrawBeta_LNIRT(const arma::vec& theta, const arma::vec& alpha, const arma::mat& Z, const arma::vec& mu, const double sigma);
RcppExport SEXP LNIRT_Rcpp_DrawBeta_LNIRT(SEXP thetaSEXP, SEXP alphaSEXP, SEXP ZSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const double >::type sigma(sigmaSEXP);
    __result = Rcpp::wrap(Rcpp_DrawBeta_LNIRT(theta, alpha, Z, mu, sigma));
    return __result;
END_RCPP
}
// Rcpp_DrawLambda_LNIRT
arma::vec Rcpp_DrawLambda_LNIRT(const arma::mat& RT, const arma::vec& phi, const arma::vec& zeta, const arma::vec& sigma2, const arma::vec& mu, const double sigmal);
RcppExport SEXP LNIRT_Rcpp_DrawLambda_LNIRT(SEXP RTSEXP, SEXP phiSEXP, SEXP zetaSEXP, SEXP sigma2SEXP, SEXP muSEXP, SEXP sigmalSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type RT(RTSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const double >::type sigmal(sigmalSEXP);
    __result = Rcpp::wrap(Rcpp_DrawLambda_LNIRT(RT, phi, zeta, sigma2, mu, sigmal));
    return __result;
END_RCPP
}
// Rcpp_DrawAlpha_LNIRT
arma::vec Rcpp_DrawAlpha_LNIRT(const arma::vec& theta, const arma::vec& beta, const arma::mat& Z, const arma::vec& mu, const double sigma);
RcppExport SEXP LNIRT_Rcpp_DrawAlpha_LNIRT(SEXP thetaSEXP, SEXP betaSEXP, SEXP ZSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const double >::type sigma(sigmaSEXP);
    __result = Rcpp::wrap(Rcpp_DrawAlpha_LNIRT(theta, beta, Z, mu, sigma));
    return __result;
END_RCPP
}
// Rcpp_DrawPhi_LNIRT
arma::vec Rcpp_DrawPhi_LNIRT(const arma::mat& RT, const arma::vec& lambda, const arma::vec& zeta, const arma::vec& sigma2, const arma::vec& mu, const double sigmal);
RcppExport SEXP LNIRT_Rcpp_DrawPhi_LNIRT(SEXP RTSEXP, SEXP lambdaSEXP, SEXP zetaSEXP, SEXP sigma2SEXP, SEXP muSEXP, SEXP sigmalSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type RT(RTSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const double >::type sigmal(sigmalSEXP);
    __result = Rcpp::wrap(Rcpp_DrawPhi_LNIRT(RT, lambda, zeta, sigma2, mu, sigmal));
    return __result;
END_RCPP
}
// Rcpp_SampleS2_LNIRT
arma::vec Rcpp_SampleS2_LNIRT(const arma::mat& RT, const arma::vec& zeta, const arma::vec& lambda, const arma::vec& phi);
RcppExport SEXP LNIRT_Rcpp_SampleS2_LNIRT(SEXP RTSEXP, SEXP zetaSEXP, SEXP lambdaSEXP, SEXP phiSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type RT(RTSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type phi(phiSEXP);
    __result = Rcpp::wrap(Rcpp_SampleS2_LNIRT(RT, zeta, lambda, phi));
    return __result;
END_RCPP
}
// Rcpp_SampleB_LNIRT
List Rcpp_SampleB_LNIRT(const arma::mat& Y, const arma::vec& X, const arma::mat& Sigma, const arma::vec& B0, const arma::mat& V0);
RcppExport SEXP LNIRT_Rcpp_SampleB_LNIRT(SEXP YSEXP, SEXP XSEXP, SEXP SigmaSEXP, SEXP B0SEXP, SEXP V0SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V0(V0SEXP);
    __result = Rcpp::wrap(Rcpp_SampleB_LNIRT(Y, X, Sigma, B0, V0));
    return __result;
END_RCPP
}
// Rcpp_rwishart_LNIRT
arma::mat Rcpp_rwishart_LNIRT(const int nu, const arma::mat V);
RcppExport SEXP LNIRT_Rcpp_rwishart_LNIRT(SEXP nuSEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const int >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type V(VSEXP);
    __result = Rcpp::wrap(Rcpp_rwishart_LNIRT(nu, V));
    return __result;
END_RCPP
}
// rcpp_hello
List rcpp_hello();
RcppExport SEXP LNIRT_rcpp_hello() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(rcpp_hello());
    return __result;
END_RCPP
}

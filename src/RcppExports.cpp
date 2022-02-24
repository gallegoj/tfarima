// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// acovtomaC
arma::colvec acovtomaC(const arma::colvec& g);
RcppExport SEXP _tfarima_acovtomaC(SEXP gSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type g(gSEXP);
    rcpp_result_gen = Rcpp::wrap(acovtomaC(g));
    return rcpp_result_gen;
END_RCPP
}
// decompHC
arma::mat decompHC(const arma::mat& T, const double mu);
RcppExport SEXP _tfarima_decompHC(SEXP TSEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type T(TSEXP);
    Rcpp::traits::input_parameter< const double >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(decompHC(T, mu));
    return rcpp_result_gen;
END_RCPP
}
// decompFC
arma::mat decompFC(const arma::mat& T, const double mu);
RcppExport SEXP _tfarima_decompFC(SEXP TSEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type T(TSEXP);
    Rcpp::traits::input_parameter< const double >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(decompFC(T, mu));
    return rcpp_result_gen;
END_RCPP
}
// deceffBC
arma::mat deceffBC(const arma::colvec& y, const bool& bc, const double& mu, const arma::colvec& phi, const arma::colvec& nabla, const arma::colvec& theta, double& sig2, const arma::mat& F, int type);
RcppExport SEXP _tfarima_deceffBC(SEXP ySEXP, SEXP bcSEXP, SEXP muSEXP, SEXP phiSEXP, SEXP nablaSEXP, SEXP thetaSEXP, SEXP sig2SEXP, SEXP FSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const bool& >::type bc(bcSEXP);
    Rcpp::traits::input_parameter< const double& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type nabla(nablaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double& >::type sig2(sig2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type F(FSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(deceffBC(y, bc, mu, phi, nabla, theta, sig2, F, type));
    return rcpp_result_gen;
END_RCPP
}
// seasadjC
arma::vec seasadjC(const arma::colvec& y, const bool& bc, const double& mu, const arma::colvec& phi, const arma::colvec& nabla, const arma::colvec& theta, double& sig2, const arma::cx_colvec& ariroots, int method);
RcppExport SEXP _tfarima_seasadjC(SEXP ySEXP, SEXP bcSEXP, SEXP muSEXP, SEXP phiSEXP, SEXP nablaSEXP, SEXP thetaSEXP, SEXP sig2SEXP, SEXP arirootsSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const bool& >::type bc(bcSEXP);
    Rcpp::traits::input_parameter< const double& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type nabla(nablaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double& >::type sig2(sig2SEXP);
    Rcpp::traits::input_parameter< const arma::cx_colvec& >::type ariroots(arirootsSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(seasadjC(y, bc, mu, phi, nabla, theta, sig2, ariroots, method));
    return rcpp_result_gen;
END_RCPP
}
// diffC
arma::colvec diffC(const arma::colvec& z, const arma::colvec& nabla, const bool& bc);
RcppExport SEXP _tfarima_diffC(SEXP zSEXP, SEXP nablaSEXP, SEXP bcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type nabla(nablaSEXP);
    Rcpp::traits::input_parameter< const bool& >::type bc(bcSEXP);
    rcpp_result_gen = Rcpp::wrap(diffC(z, nabla, bc));
    return rcpp_result_gen;
END_RCPP
}
// filterC
arma::colvec filterC(const arma::colvec& x, const arma::colvec& omega, const arma::colvec& delta, int b);
RcppExport SEXP _tfarima_filterC(SEXP xSEXP, SEXP omegaSEXP, SEXP deltaSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(filterC(x, omega, delta, b));
    return rcpp_result_gen;
END_RCPP
}
// forecastC
arma::mat forecastC(const arma::colvec& y, const bool bc, const double& mu, const arma::colvec& phi, const arma::colvec& nabla, const arma::colvec& theta, double sig2, int ori, const int hor);
RcppExport SEXP _tfarima_forecastC(SEXP ySEXP, SEXP bcSEXP, SEXP muSEXP, SEXP phiSEXP, SEXP nablaSEXP, SEXP thetaSEXP, SEXP sig2SEXP, SEXP oriSEXP, SEXP horSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const bool >::type bc(bcSEXP);
    Rcpp::traits::input_parameter< const double& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type nabla(nablaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type sig2(sig2SEXP);
    Rcpp::traits::input_parameter< int >::type ori(oriSEXP);
    Rcpp::traits::input_parameter< const int >::type hor(horSEXP);
    rcpp_result_gen = Rcpp::wrap(forecastC(y, bc, mu, phi, nabla, theta, sig2, ori, hor));
    return rcpp_result_gen;
END_RCPP
}
// backcastC
arma::colvec backcastC(const arma::colvec& y, const bool bc, const double& mu, const arma::colvec& phi, const arma::colvec& nabla, const arma::colvec& theta, double sig2, int ori, const int hor);
RcppExport SEXP _tfarima_backcastC(SEXP ySEXP, SEXP bcSEXP, SEXP muSEXP, SEXP phiSEXP, SEXP nablaSEXP, SEXP thetaSEXP, SEXP sig2SEXP, SEXP oriSEXP, SEXP horSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const bool >::type bc(bcSEXP);
    Rcpp::traits::input_parameter< const double& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type nabla(nablaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type sig2(sig2SEXP);
    Rcpp::traits::input_parameter< int >::type ori(oriSEXP);
    Rcpp::traits::input_parameter< const int >::type hor(horSEXP);
    rcpp_result_gen = Rcpp::wrap(backcastC(y, bc, mu, phi, nabla, theta, sig2, ori, hor));
    return rcpp_result_gen;
END_RCPP
}
// glsC
const arma::colvec glsC(const arma::colvec& y, const arma::mat& X, const arma::colvec& phi, const arma::colvec& theta);
RcppExport SEXP _tfarima_glsC(SEXP ySEXP, SEXP XSEXP, SEXP phiSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(glsC(y, X, phi, theta));
    return rcpp_result_gen;
END_RCPP
}
// ellarmaC
double ellarmaC(const arma::colvec& w, const arma::colvec& phi, const arma::colvec& theta);
RcppExport SEXP _tfarima_ellarmaC(SEXP wSEXP, SEXP phiSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(ellarmaC(w, phi, theta));
    return rcpp_result_gen;
END_RCPP
}
// gresC
arma::colvec gresC(const arma::colvec& w, const arma::colvec& phi, const arma::colvec& theta);
RcppExport SEXP _tfarima_gresC(SEXP wSEXP, SEXP phiSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(gresC(w, phi, theta));
    return rcpp_result_gen;
END_RCPP
}
// ssrC
double ssrC(const arma::colvec& w, const arma::colvec& phi, const arma::colvec& theta);
RcppExport SEXP _tfarima_ssrC(SEXP wSEXP, SEXP phiSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(ssrC(w, phi, theta));
    return rcpp_result_gen;
END_RCPP
}
// cllarmaC
double cllarmaC(const arma::colvec& w, const arma::colvec& phi, const arma::colvec& theta);
RcppExport SEXP _tfarima_cllarmaC(SEXP wSEXP, SEXP phiSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(cllarmaC(w, phi, theta));
    return rcpp_result_gen;
END_RCPP
}
// outliersC
arma::mat outliersC(const arma::colvec& z, bool bc, double mu, const arma::colvec& phi, const arma::colvec& nabla, const arma::colvec& theta, arma::ucolvec& timing, bool eres, double c);
RcppExport SEXP _tfarima_outliersC(SEXP zSEXP, SEXP bcSEXP, SEXP muSEXP, SEXP phiSEXP, SEXP nablaSEXP, SEXP thetaSEXP, SEXP timingSEXP, SEXP eresSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type z(zSEXP);
    Rcpp::traits::input_parameter< bool >::type bc(bcSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type nabla(nablaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::ucolvec& >::type timing(timingSEXP);
    Rcpp::traits::input_parameter< bool >::type eres(eresSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(outliersC(z, bc, mu, phi, nabla, theta, timing, eres, c));
    return rcpp_result_gen;
END_RCPP
}
// polyevalC
double polyevalC(const arma::colvec& pol, double z);
RcppExport SEXP _tfarima_polyevalC(SEXP polSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type pol(polSEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(polyevalC(pol, z));
    return rcpp_result_gen;
END_RCPP
}
// polyrootsC
arma::mat polyrootsC(const arma::colvec& pol);
RcppExport SEXP _tfarima_polyrootsC(SEXP polSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type pol(polSEXP);
    rcpp_result_gen = Rcpp::wrap(polyrootsC(pol));
    return rcpp_result_gen;
END_RCPP
}
// sortrootsC
arma::mat sortrootsC(const arma::cx_colvec& r);
RcppExport SEXP _tfarima_sortrootsC(SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cx_colvec& >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(sortrootsC(r));
    return rcpp_result_gen;
END_RCPP
}
// combinerootsC
arma::mat combinerootsC(arma::mat T);
RcppExport SEXP _tfarima_combinerootsC(SEXP TSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type T(TSEXP);
    rcpp_result_gen = Rcpp::wrap(combinerootsC(T));
    return rcpp_result_gen;
END_RCPP
}
// roots2polC
arma::mat roots2polC(arma::mat T);
RcppExport SEXP _tfarima_roots2polC(SEXP TSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type T(TSEXP);
    rcpp_result_gen = Rcpp::wrap(roots2polC(T));
    return rcpp_result_gen;
END_RCPP
}
// admregC
bool admregC(const arma::colvec& pol, bool ar);
RcppExport SEXP _tfarima_admregC(SEXP polSEXP, SEXP arSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type pol(polSEXP);
    Rcpp::traits::input_parameter< bool >::type ar(arSEXP);
    rcpp_result_gen = Rcpp::wrap(admregC(pol, ar));
    return rcpp_result_gen;
END_RCPP
}
// polymultC
arma::colvec polymultC(const arma::colvec& pol1, const arma::colvec& pol2);
RcppExport SEXP _tfarima_polymultC(SEXP pol1SEXP, SEXP pol2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type pol1(pol1SEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type pol2(pol2SEXP);
    rcpp_result_gen = Rcpp::wrap(polymultC(pol1, pol2));
    return rcpp_result_gen;
END_RCPP
}
// polydivC
arma::colvec polydivC(const arma::colvec& pol1, const arma::colvec& pol2, bool rem);
RcppExport SEXP _tfarima_polydivC(SEXP pol1SEXP, SEXP pol2SEXP, SEXP remSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type pol1(pol1SEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type pol2(pol2SEXP);
    Rcpp::traits::input_parameter< bool >::type rem(remSEXP);
    rcpp_result_gen = Rcpp::wrap(polydivC(pol1, pol2, rem));
    return rcpp_result_gen;
END_RCPP
}
// polygcdC
arma::colvec polygcdC(const arma::colvec& pol1, const arma::colvec& pol2);
RcppExport SEXP _tfarima_polygcdC(SEXP pol1SEXP, SEXP pol2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type pol1(pol1SEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type pol2(pol2SEXP);
    rcpp_result_gen = Rcpp::wrap(polygcdC(pol1, pol2));
    return rcpp_result_gen;
END_RCPP
}
// polyprsC
arma::colvec polyprsC(const arma::colvec& pol1, const arma::colvec& pol2);
RcppExport SEXP _tfarima_polyprsC(SEXP pol1SEXP, SEXP pol2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type pol1(pol1SEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type pol2(pol2SEXP);
    rcpp_result_gen = Rcpp::wrap(polyprsC(pol1, pol2));
    return rcpp_result_gen;
END_RCPP
}
// polyraiseC
arma::colvec polyraiseC(const arma::colvec& pol, int d);
RcppExport SEXP _tfarima_polyraiseC(SEXP polSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type pol(polSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(polyraiseC(pol, d));
    return rcpp_result_gen;
END_RCPP
}
// polyfactorsC
arma::mat polyfactorsC(const arma::colvec& pol);
RcppExport SEXP _tfarima_polyfactorsC(SEXP polSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type pol(polSEXP);
    rcpp_result_gen = Rcpp::wrap(polyfactorsC(pol));
    return rcpp_result_gen;
END_RCPP
}
// polyratioC
arma::colvec polyratioC(const arma::colvec& num, const arma::colvec& den, int d);
RcppExport SEXP _tfarima_polyratioC(SEXP numSEXP, SEXP denSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type num(numSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type den(denSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(polyratioC(num, den, d));
    return rcpp_result_gen;
END_RCPP
}
// condresC
arma::colvec condresC(const arma::colvec& w, const arma::colvec& phi, const arma::colvec& theta, const bool forward);
RcppExport SEXP _tfarima_condresC(SEXP wSEXP, SEXP phiSEXP, SEXP thetaSEXP, SEXP forwardSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const bool >::type forward(forwardSEXP);
    rcpp_result_gen = Rcpp::wrap(condresC(w, phi, theta, forward));
    return rcpp_result_gen;
END_RCPP
}
// inicondC
arma::colvec inicondC(const arma::colvec& w, const arma::colvec& phi, const arma::colvec& theta);
RcppExport SEXP _tfarima_inicondC(SEXP wSEXP, SEXP phiSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(inicondC(w, phi, theta));
    return rcpp_result_gen;
END_RCPP
}
// exactresC
arma::colvec exactresC(const arma::colvec& w, const arma::colvec& phi, const arma::colvec& theta);
RcppExport SEXP _tfarima_exactresC(SEXP wSEXP, SEXP phiSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(exactresC(w, phi, theta));
    return rcpp_result_gen;
END_RCPP
}
// cssrC
double cssrC(const arma::colvec& w, const arma::colvec& phi, const arma::colvec& theta);
RcppExport SEXP _tfarima_cssrC(SEXP wSEXP, SEXP phiSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(cssrC(w, phi, theta));
    return rcpp_result_gen;
END_RCPP
}
// simC
arma::colvec simC(arma::colvec a, const bool bc, const double mu, const arma::colvec& phi, const arma::colvec& nabla, const arma::colvec& theta, const arma::colvec& y0);
RcppExport SEXP _tfarima_simC(SEXP aSEXP, SEXP bcSEXP, SEXP muSEXP, SEXP phiSEXP, SEXP nablaSEXP, SEXP thetaSEXP, SEXP y0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type a(aSEXP);
    Rcpp::traits::input_parameter< const bool >::type bc(bcSEXP);
    Rcpp::traits::input_parameter< const double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type nabla(nablaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y0(y0SEXP);
    rcpp_result_gen = Rcpp::wrap(simC(a, bc, mu, phi, nabla, theta, y0));
    return rcpp_result_gen;
END_RCPP
}
// spectrumC
arma::mat spectrumC(const arma::colvec& phi, const arma::colvec& theta, double sigma2, int nfreq);
RcppExport SEXP _tfarima_spectrumC(SEXP phiSEXP, SEXP thetaSEXP, SEXP sigma2SEXP, SEXP nfreqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< int >::type nfreq(nfreqSEXP);
    rcpp_result_gen = Rcpp::wrap(spectrumC(phi, theta, sigma2, nfreq));
    return rcpp_result_gen;
END_RCPP
}
// pgramC
arma::mat pgramC(const arma::colvec& y, bool cpgram);
RcppExport SEXP _tfarima_pgramC(SEXP ySEXP, SEXP cpgramSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< bool >::type cpgram(cpgramSEXP);
    rcpp_result_gen = Rcpp::wrap(pgramC(y, cpgram));
    return rcpp_result_gen;
END_RCPP
}
// tacovC
arma::colvec tacovC(const arma::colvec& phi, const arma::colvec& theta, double sigma2, int nlags);
RcppExport SEXP _tfarima_tacovC(SEXP phiSEXP, SEXP thetaSEXP, SEXP sigma2SEXP, SEXP nlagsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< int >::type nlags(nlagsSEXP);
    rcpp_result_gen = Rcpp::wrap(tacovC(phi, theta, sigma2, nlags));
    return rcpp_result_gen;
END_RCPP
}
// pacorrC
const arma::colvec pacorrC(const arma::colvec& phi, const arma::colvec& theta, int nlags);
RcppExport SEXP _tfarima_pacorrC(SEXP phiSEXP, SEXP thetaSEXP, SEXP nlagsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type nlags(nlagsSEXP);
    rcpp_result_gen = Rcpp::wrap(pacorrC(phi, theta, nlags));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tfarima_acovtomaC", (DL_FUNC) &_tfarima_acovtomaC, 1},
    {"_tfarima_decompHC", (DL_FUNC) &_tfarima_decompHC, 2},
    {"_tfarima_decompFC", (DL_FUNC) &_tfarima_decompFC, 2},
    {"_tfarima_deceffBC", (DL_FUNC) &_tfarima_deceffBC, 9},
    {"_tfarima_seasadjC", (DL_FUNC) &_tfarima_seasadjC, 9},
    {"_tfarima_diffC", (DL_FUNC) &_tfarima_diffC, 3},
    {"_tfarima_filterC", (DL_FUNC) &_tfarima_filterC, 4},
    {"_tfarima_forecastC", (DL_FUNC) &_tfarima_forecastC, 9},
    {"_tfarima_backcastC", (DL_FUNC) &_tfarima_backcastC, 9},
    {"_tfarima_glsC", (DL_FUNC) &_tfarima_glsC, 4},
    {"_tfarima_ellarmaC", (DL_FUNC) &_tfarima_ellarmaC, 3},
    {"_tfarima_gresC", (DL_FUNC) &_tfarima_gresC, 3},
    {"_tfarima_ssrC", (DL_FUNC) &_tfarima_ssrC, 3},
    {"_tfarima_cllarmaC", (DL_FUNC) &_tfarima_cllarmaC, 3},
    {"_tfarima_outliersC", (DL_FUNC) &_tfarima_outliersC, 9},
    {"_tfarima_polyevalC", (DL_FUNC) &_tfarima_polyevalC, 2},
    {"_tfarima_polyrootsC", (DL_FUNC) &_tfarima_polyrootsC, 1},
    {"_tfarima_sortrootsC", (DL_FUNC) &_tfarima_sortrootsC, 1},
    {"_tfarima_combinerootsC", (DL_FUNC) &_tfarima_combinerootsC, 1},
    {"_tfarima_roots2polC", (DL_FUNC) &_tfarima_roots2polC, 1},
    {"_tfarima_admregC", (DL_FUNC) &_tfarima_admregC, 2},
    {"_tfarima_polymultC", (DL_FUNC) &_tfarima_polymultC, 2},
    {"_tfarima_polydivC", (DL_FUNC) &_tfarima_polydivC, 3},
    {"_tfarima_polygcdC", (DL_FUNC) &_tfarima_polygcdC, 2},
    {"_tfarima_polyprsC", (DL_FUNC) &_tfarima_polyprsC, 2},
    {"_tfarima_polyraiseC", (DL_FUNC) &_tfarima_polyraiseC, 2},
    {"_tfarima_polyfactorsC", (DL_FUNC) &_tfarima_polyfactorsC, 1},
    {"_tfarima_polyratioC", (DL_FUNC) &_tfarima_polyratioC, 3},
    {"_tfarima_condresC", (DL_FUNC) &_tfarima_condresC, 4},
    {"_tfarima_inicondC", (DL_FUNC) &_tfarima_inicondC, 3},
    {"_tfarima_exactresC", (DL_FUNC) &_tfarima_exactresC, 3},
    {"_tfarima_cssrC", (DL_FUNC) &_tfarima_cssrC, 3},
    {"_tfarima_simC", (DL_FUNC) &_tfarima_simC, 7},
    {"_tfarima_spectrumC", (DL_FUNC) &_tfarima_spectrumC, 4},
    {"_tfarima_pgramC", (DL_FUNC) &_tfarima_pgramC, 2},
    {"_tfarima_tacovC", (DL_FUNC) &_tfarima_tacovC, 4},
    {"_tfarima_pacorrC", (DL_FUNC) &_tfarima_pacorrC, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_tfarima(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

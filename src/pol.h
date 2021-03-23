#ifndef _ARIMA_POL_
#define _ARIMA_POL_

double polyevalC(const arma::colvec &pol, double z);
arma::mat polyrootsC(const arma::colvec &pol);
arma::mat sortrootsC(const arma::cx_colvec &roots);
arma::mat combinerootsC(arma::mat T);
arma::mat roots2polC(arma::mat T);
bool admregC(const arma::colvec &pol, bool ar);
arma::colvec polymultC(const arma::colvec &pol1, const arma::colvec &pol2);
arma::colvec polydivC(const arma::colvec &pol1, const arma::colvec &pol2, bool rem = false);
arma::colvec polygcdC(const arma::colvec &pol1, const arma::colvec &pol2);
arma::colvec polyprsC(const arma::colvec &pol1, const arma::colvec &pol2);
arma::colvec polyraiseC(const arma::colvec &pol, int d);
arma::colvec polyratioC(const arma::colvec &num, const arma::colvec &den, int d);
arma::mat polyfactorsC(const arma::colvec &pol);
bool simeqC(double x, double y, double tol = 1e-5);
bool ltC(double x, double y, double tol = 1e-5);

#endif //_ARIMA_POL_

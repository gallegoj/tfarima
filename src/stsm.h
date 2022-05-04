#ifndef __ARIMA_STSM__
#define __ARIMA_STSM__
double llrfC(const arma::colvec &w, const arma::colvec &d, const arma::mat &A,
             const arma::mat &Sv, double s2u, bool s2star);
#endif //__ARIMA_STSM__

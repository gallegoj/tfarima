#ifndef __ARIMA_TACOV__
#define __ARIMA_TACOV__

arma::colvec tacovC(const arma::colvec &phi, const arma::colvec &theta, double sigma2, int nlags);
const arma::colvec pacorrC(const arma::colvec &phi, const arma::colvec &theta, int nlags);

#endif //__ARIMA_TACOV__

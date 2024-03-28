#ifndef __ARIMA_LOGLIK__
#define __ARIMA_LOGLIK__

double ellarmaC(const arma::colvec &w, const arma::colvec &phi,const arma::colvec &theta);
arma::colvec gresC(const arma::colvec &w, const arma::colvec &phi, const arma::colvec &theta);
double ssrC(const arma::colvec &w, const arma::colvec &phi, const arma::colvec &theta);
double cllarmaC(const arma::colvec &w, const arma::colvec &phi,const arma::colvec &theta);
#endif //__ARIMA_LOGLIK__

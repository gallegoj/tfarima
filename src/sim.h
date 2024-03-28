#ifndef __ARIMA_SIM__
#define __ARIMA_SIM__

arma::colvec simC(arma::colvec a, const bool bc, const double mu,
                  const arma::colvec &phi, const arma::colvec &nabla,
                  const arma::colvec &theta, const arma::colvec &y0);


#endif //__ARIMA_SIM__

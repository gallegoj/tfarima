#ifndef __ARIMA_FILTER__
#define __ARIMA_FILTER__

arma::colvec filterC(const arma::colvec &x, const arma::colvec &omega,
                     const arma::colvec &delta, int b);

#endif //__ARIMA_FILTER__

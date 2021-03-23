#ifndef __ARIMA_OUTLIERS__
#define __ARIMA_OUTLIERS__

arma::mat outliersC(const arma::colvec &z, bool bc, double mu, const arma::colvec &phi,
                    const arma::colvec &nabla, const arma::colvec &theta, 
                    arma::ucolvec &timing, double c);  

#endif //__ARIMA_OUTLIERS__

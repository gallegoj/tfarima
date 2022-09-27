#ifndef __ARIMA_PRESAMPLE__
#define __ARIMA_PRESAMPLE__

arma::mat presampleCovC(const arma::colvec &phi, const arma::colvec &theta, 
                        bool fvf = true);

#endif //__ARIMA_ACOVTOMA__

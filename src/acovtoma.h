#ifndef __ARIMA_ACOVTOMA__
#define __ARIMA_ACOVTOMA__

arma::colvec acovtomaC(const arma::colvec &g, int &code, double tol = 1e-6, 
                       int max_iter = 100);
void checkmaC(arma::colvec &ma);
#endif //__ARIMA_ACOVTOMA__

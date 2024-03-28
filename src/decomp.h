#ifndef __ARIMA_DECOMP__
#define __ARIMA_DECOMP__

arma::mat decompFC(const arma::mat &T, const double mu);
arma::mat deceffBC(const arma::colvec &y, const bool &bc, const double &mu,
                     const arma::colvec &phi, const arma::colvec &nabla,
                     const arma::colvec &theta, double &sig2, const arma::mat &F, type = 1);

arma::vec seasadjC(const arma::colvec &y, const bool &bc, const double &mu,
                  const arma::colvec &phi, const arma::colvec &nabla,
                  const arma::colvec &theta, double &sig2,
                  const arma::cx_colvec &ariroots, int method = 3);

#endif //__ARIMA_DECOMP__


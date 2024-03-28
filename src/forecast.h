#ifndef __ARIMA_FORECAST__
#define __ARIMA_FORECAST__

arma::mat forecastC(const arma::colvec &y, const bool bc, const double &mu,
                    const arma::colvec &phi, const arma::colvec &nabla,
                    const arma::colvec &theta, double sig2,
                    int ori, const int hor);

arma::colvec backcastC(const arma::colvec &y, const bool bc, const double &mu,
                    const arma::colvec &phi, const arma::colvec &nabla,
                    const arma::colvec &theta, double sig2,
                    int ori, const int hor);

#endif //__ARIMA_FORECAST__


#ifndef __ARIMA_SPEC__
#define __ARIMA_SPEC__

arma::mat spectrumC(const arma::colvec &phi, const arma::colvec &theta, double sigma2, int nfreq);
arma::mat pgramC(const arma::colvec &y, bool cpgram = false);

#endif //__ARIMA_SPEC__

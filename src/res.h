#ifndef __MAUROOTS_CRES__
#define __MAUROOTS_CRES__

arma::colvec condresC(const arma::colvec &w, const arma::colvec &phi, const arma::colvec &theta, const bool forward = true);
arma::colvec inicondC(const arma::colvec &w, const arma::colvec &phi, const arma::colvec &theta);
arma::colvec exactresC(const arma::colvec &w, const arma::colvec &phi, const arma::colvec &theta);
double cssrC(const arma::colvec &w, const arma::colvec &phi, const arma::colvec &theta);
#endif //__MAUROOTS_CRES__

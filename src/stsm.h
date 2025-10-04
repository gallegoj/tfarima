#ifndef __ARIMA_STSM__
#define __ARIMA_STSM__
double llrfC(const arma::colvec &w, const arma::colvec &nabla,
             const arma::rowvec &b, const arma::mat &C,
             const arma::mat &S, arma::colvec &s2, bool cform);

double llucaC(arma::colvec &w, const arma::colvec &phi,
              const arma::mat &A, const arma::mat &S, arma::colvec &s2, 
              int res);

arma::colvec resrfC(const arma::colvec &w, const arma::colvec &nabla,
                    const arma::rowvec &b, const arma::mat &C,
                    const arma::mat &S, arma::colvec& s2, bool cform);


bool kf0C(const arma::colvec &y, const arma::colvec &b, const arma::mat &C,
          const arma::mat &S, const arma::colvec &x0, const arma::mat &P0,
          arma::colvec &v, arma::colvec &s2);
  
bool kfC(const arma::colvec &y, const arma::colvec &b, const arma::mat &C,
         const arma::mat &S, const arma::colvec &x0, const arma::mat &P0,
         arma::colvec &v, arma::colvec &s2, arma::mat &X, arma::mat &PX, 
         bool filtered);

bool ksC(const arma::colvec &y, const arma::colvec &b, const arma::mat &C,
         const arma::mat &S, arma::colvec &x0, arma::mat &P0,
         arma::mat &X, arma::mat &PX);

#endif //__ARIMA_STSM__

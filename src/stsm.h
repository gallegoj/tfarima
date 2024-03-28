#ifndef __ARIMA_STSM__
#define __ARIMA_STSM__
double llrfC(const arma::colvec &w, const arma::colvec &nabla,
             const arma::rowvec &b, const arma::mat &C,
             const arma::mat &S, arma::colvec &s2, bool cform);

arma::colvec resrfC(const arma::colvec &w, const arma::colvec &nabla,
                    const arma::rowvec &b, const arma::mat &C,
                    const arma::mat &S, arma::colvec& s2, bool cform);

bool kf0C(const arma::colvec &y, const arma::colvec &b, const arma::mat &C,
          const arma::mat &S, const arma::colvec &x1, const arma::mat &P1,
          arma::colvec &v, arma::colvec &s2);
  
bool kfC(const arma::colvec &y, const arma::colvec &b, const arma::mat &C,
         const arma::mat &S, arma::colvec &x1, arma::mat &P1,
         arma::colvec &v, arma::colvec &s2, arma::mat &X, arma::mat &PX, 
         bool cform, bool filtered, bool xn);

bool ksC(const arma::colvec &y, const arma::colvec &b, const arma::mat &C,
         const arma::mat &S, const arma::colvec &x1, const arma::mat &P1,
         arma::mat &X, arma::mat &PX, bool cform);

#endif //__ARIMA_STSM__

#include <Rcpp.h>
using namespace Rcpp;

//' Log of Multivarate Gamma Function - Gamma_p(a)
//' @reference https://en.wikipedia.org/wiki/Multivariate_gamma_function
// [[Rcpp::export]]
double lmvgamma(double a, int p){
  static const double pi = log(3.14159265); 
  double s=0;
  double x = pi*(p*(p-1.0))/2.0;
  for (int i=1; i<=p; i++){
    s += lgamma(a+(1.0-i)/2);
  }
  return(x+s);
}

//' Derivative of Log of Multivariate Gamma Function - Gamma_p(a)
//' @reference https://en.wikipedia.org/wiki/Multivariate_gamma_function
// [[Rcpp::export]]
double lmvgamma_deriv(double a, int p){
  double s=0;
  for (int i=1; i<=p; i++){
    s += R::digamma(a + 0.5*(1-i));
  }
  return s*lmvgamma(a,p);
}
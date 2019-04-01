#ifndef MONGREL_SPECIALFUNC_H
#define MONGREL_SPECIALFUNC_H

#include <Rcpp.h>
using namespace Rcpp;

//' Log of Multivarate Gamma Function
//' Gamma_p(a) - https://en.wikipedia.org/wiki/Multivariate_gamma_function
double lmvgamma(double a, int p);
//' Derivative of Log of Multivariate Gamma Function
//' https://en.wikipedia.org/wiki/Multivariate_gamma_function
//' Gamma_p(a)
double lmvgamma_deriv(double a, int p);

#endif
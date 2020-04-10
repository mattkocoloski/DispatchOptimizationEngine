// simple test function

// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sumCpp(NumericVector x) {
  
  double value = 0;
  
  for (int i = 0; i < x.length(); ++i) {
    value = value + x(i);
  }
  
  return(value);
  
}
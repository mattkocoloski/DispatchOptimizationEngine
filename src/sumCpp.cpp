// simple function to test whether package was installed and loaded correctly

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double sumCpp(NumericVector x) {
  
  double value = 0;
  
  for (int i = 0; i < x.length(); ++i) {
    value = value + x[i];
  }
  
  return(value);
  
}

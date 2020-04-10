// fills in on and off constraint matrices whenever gap between constraints falls below minimum gap time
// on gaps must be >= min downtime, off gaps >= min runtime

// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix fillConstraintMatrixGaps(NumericMatrix M,
                                       NumericVector MinGap) {
  
  // dimensions
  int num_hours = M.nrow();
  int num_units = M.ncol();
  
  // output matrix
  NumericMatrix Fill(num_hours, num_units);
  
  // matrices to track gap starts and stops
  NumericMatrix GapStart(num_hours, num_units);
  NumericMatrix GapStop(num_hours, num_units);
  NumericMatrix HoursSinceStart(num_hours, num_units);
  NumericMatrix HoursSinceStop(num_hours, num_units);
  NumericMatrix Violations(num_hours, num_units);
  
  for (int h = 0; h < num_units; ++h) {

    // starts
    for (int i = 0; i < num_hours; ++i) {

      if (i > 0) {

        if (M(i-1,h) == 1 && M(i,h) == 0) {
          GapStart(i,h) = 1;
        } else {
          GapStart(i,h) = 0;
        }

      } else {

        if (M(i,h) == 0) {
          GapStart(i,h) = 1;
        } else {
          GapStart(i,h) = 0;
        }

      }
    }

    // stops
    for (int i = 0; i < num_hours; ++i) {
      if (i > 0) {
        if (M(i-1,h) == 0 && M(i,h) == 1) {
          GapStop(i,h) = 1;
        } else {
          GapStop(i,h) = 0;
        }
      }
    }

    // update HoursSinceStop and HoursSinceStart
    for (int i = 0; i < num_hours; ++i) {

      if (i == 0) {
        HoursSinceStart(i,h) = 9999;
        HoursSinceStop(i,h)  = 9999;
      } else {
        HoursSinceStart(i,h) = HoursSinceStart(i-1,h) + 1;
        HoursSinceStop(i,h)  = HoursSinceStop(i-1,h) + 1;
      }

      if (GapStart(i,h) == 1) {
        HoursSinceStart(i,h) = 1;
      }

      if (GapStop(i,h) == 1) {
        HoursSinceStop(i,h) = 1;
      }
    }

    // violations
    for (int i = 0; i < num_hours; ++i) {
      
      if (GapStop(i,h) == 1 && HoursSinceStart(i,h) < MinGap(h)) {

        // flag violation
        Violations(i,h) = 1;

        // fill in preceding hours
        int first = i - HoursSinceStart(i,h) + 1;
        int last  = i - 1;
        for (int k = first; k <= last; ++k) {
          if (k >= 0) {
            Fill(k,h) = 1;
          }
        }
      }
    }
  }

  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //                                                                                              //
  // Return Results                                                                               //
  //                                                                                              //
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  // // combine into single list to return
  // List Out;
  // Out["GapStart"]         = GapStart;
  // Out["GapStop"]          = GapStop;
  // Out["HoursSinceStart"]  = HoursSinceStart;
  // Out["HoursSinceStop"]   = HoursSinceStop;
  // Out["Violations"]       = Violations;
  // Out["Fill"]             = Fill;
  // 
  // return(Out);
  
  // return matrix of values to fill
  return(Fill);
  
}
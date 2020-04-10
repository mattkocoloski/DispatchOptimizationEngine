// Rcpp engine for optimizing battery dispatch

// version smoothes load by assuming energy buy / sell prices = midpoints of load blocks

// version 2 includes initial and final states as inputs
// version 3 removes loop over sims
// promoted back to base

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List optimizeBatteryStorageDP(DataFrame BatteryData,
                              NumericVector LoadVector,
                              NumericVector ChargeStates,
                              DataFrame TransitionList,
                              int initial_state,
                              int final_state) {
  
  // optimization tuning parameters
  // int max_iterations_hard                = OptParameters["MaxIterationsHard"];
  
  // individual battery parameters
  double charge_efficiency               = BatteryData["ChargeEfficiency"];
  double discharge_efficiency            = BatteryData["DischargeEfficiency"];
  // double max_charge                      = BatteryData["MaxChargeMW"];
  // double max_discharge                   = BatteryData["MaxDischargeMW"];
  
  // dimensions
  int num_states                         = ChargeStates.length();
  int num_hours                          = LoadVector.length();
  int first                              = 0;
  int last                               = num_hours - 1;
  // int nPerHour                           = 1;
  
  // create solution matrices
  NumericVector Path(num_hours);            // unitless index
  NumericVector BatteryEnergy(num_hours);   // MWh
  NumericVector Charge(num_hours);          // MW seen by grid
  NumericVector Discharge(num_hours);       // MW seen by grid
  NumericVector NewLoad(num_hours);         // MW

  // create state transition matrices
  NumericMatrix TransitionValue(num_states, num_states);
  NumericMatrix YMatrix(num_states, num_states);
  
  for (int j = 0; j < num_states; ++j) {
    for (int k = 0; k < num_states; ++k) {
      
      TransitionValue(j,k) = -999999999;
      YMatrix(j,k)         = -999999999;
      
    }
  }
  
  
  //////// Create Transition Vectors ////////
  
  
  NumericVector TransitionStartIndexR        = TransitionList["StartIndex"];
  NumericVector TransitionEndIndexR          = TransitionList["EndIndex"];
  NumericVector TransitionLossFactor         = TransitionList["LossFactor"];
  NumericVector TransitionBatteryEnergyDelta = TransitionList["BatteryEnergyDelta"];
  NumericVector TransitionLoadDelta          = TransitionList["LoadDelta"];
  
  NumericVector TransitionStart              = TransitionStartIndexR - 1;
  NumericVector TransitionEnd                = TransitionEndIndexR - 1;
  
  int num_transitions                        = TransitionStart.length();
  
  // create matrix for future best state for each current state
  NumericMatrix FutureBestState(num_hours, num_states);
  
  // create dynamic program variables across states (update each hour)
  NumericVector FutureBestValue(num_states);
  NumericVector FutureBestIter(num_states);
  NumericVector Y(num_states);
  
  // initialize future best value
  for (int j = 0; j < num_states; ++j) {
    FutureBestValue(j)  = -999999999;
  }
  
  // set value for desired final state
  FutureBestValue(final_state) = 999999999;
  
  
  //////// Run Optimization Model ////////
  
  
  // perform backward pass for optimization model
  for (int i = last; i >= first; --i) {
    
    // reset future best value for iteration
    for (int j = 0; j < num_states; ++j) {
      FutureBestIter(j) = -999999999;
    }
    
    // calculate transition values and update future best values
    for (int m = 0; m < num_transitions; ++m) {
      
      int j = TransitionStart(m);
      int k = TransitionEnd(m);
      
      // calculate value in transitioning from j to k
      if (i < last) {
        
        double old_load = LoadVector(i+1);
        double new_load = old_load + TransitionLoadDelta(m);
        
        TransitionValue(j,k) = -TransitionLoadDelta(m) * (new_load + old_load) / 2;
        
      } else {
        
        TransitionValue(j,k) = 0;
        
      }
      
      // add future best value for destination state
      YMatrix(j,k) = TransitionValue(j,k) + FutureBestValue(k);
      
      // update best transition
      if (YMatrix(j,k) > FutureBestIter(j)) {
        FutureBestIter(j)    = YMatrix(j,k);
        FutureBestState(i,j) = k;
      }
      
    }
    
    // update future value for each current state
    for (int j = 0; j < num_states; ++j) {
      FutureBestValue(j) = FutureBestIter(j);
    }
    
  }
  
  // add transition cost from initial state to first period state
  for (int j = 0; j < num_states; ++j) {
    FutureBestValue(j) = FutureBestValue(j) + TransitionValue(initial_state, j);
  }
  
  // select optimal state for first time period
  Path(first) = which_max(FutureBestValue);
  
  // forward pass to determine optimal path
  for (int i = (first + 1); i <= last; ++i) {
    Path(i) = FutureBestState(i-1, Path(i-1));
  }
  
  
  //////// Update Solution Matrices ////////
  
  
  for (int i = first; i <= last; ++i) {
    
    // battery energy level
    BatteryEnergy(i) = ChargeStates(Path(i));
    
    // charge and discharge
    if (i == 0) {
      
      // charging
      if (Path(i) > initial_state) {
        Charge(i) = (ChargeStates[Path(i)] - ChargeStates[initial_state]) / charge_efficiency;
      } else {
        Charge(i) = 0;
      }
      
      // discharging
      if (Path(i) < initial_state) {
        Discharge(i) = (ChargeStates[initial_state] - ChargeStates[Path(i)]) * discharge_efficiency;
      }
      
    } else {
      
      // charging
      if (Path(i) > Path(i-1)) {
        Charge(i) = (ChargeStates[Path(i)] - ChargeStates[Path(i-1)]) / charge_efficiency;
      } else {
        Charge(i) = 0;
      }
      
      // discharging
      if (Path(i) < Path(i-1)) {
        Discharge(i) = (ChargeStates[Path(i-1)] - ChargeStates[Path(i)]) * discharge_efficiency;
      } else {
        Discharge(i) = 0;
      }
    }
    
    // change in load, buy, and sell values
    NewLoad(i) = LoadVector(i) + Charge(i) - Discharge(i);
    
  }
  
  // create list to return
  List Out;
  
  Out["Path"]                = Path;
  Out["BatteryEnergy"]       = BatteryEnergy;
  Out["Charge"]              = Charge;
  Out["Discharge"]           = Discharge;
  Out["NewLoad"]             = NewLoad;
  
  return(Out);
  
}

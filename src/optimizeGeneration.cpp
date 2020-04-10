// dispatch optimization engine

// version New used pre-optimized unit commitments and simply optimizes
// generation through greedy search

// create list of all transitions for on units (1-2, 2-3, 3-4) (including buys/sells)
// sort by lowest marginal cost

// version 2 includes power purchases
// promoted back to base

// version 2 includes regulation constraint
// promoted back to base

// to do:
// check for ramp violations

// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List optimizeGeneration(NumericMatrix CommitMatrix,
                        IntegerVector AlwaysOn,
                        DataFrame PlantParameters,
                        DataFrame UnitList,
                        NumericMatrix StateGenerationMatrix,
                        NumericMatrix StateHeatRateMatrix,
                        NumericVector PowerPriceVector,
                        NumericVector LoadVector,
                        DataFrame HourIndexList) {
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //                                                                                              //
  // Process Inputs                                                                               //
  //                                                                                              //
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  // create month index variable
  NumericVector MonthIndexR          = HourIndexList["MonthIndex"];
  NumericVector MonthIndex           = MonthIndexR - 1;
  
  // dimensions
  int num_hours                      = HourIndexList.nrow();
  int num_months                     = max(MonthIndex) + 1;
  int num_units                      = UnitList.nrow();
  int max_num_states                 = StateGenerationMatrix.ncol();
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //                                                                                              //
  // Reshape Inputs Across Units                                                                  //
  //                                                                                              //
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  // create individual vectors from PlantParameters data frame
  NumericVector PPUnitIndexR            = PlantParameters["UnitIndex"];
  NumericVector PPMonthIndexR           = PlantParameters["MonthIndex"];

  NumericVector PPUnitIndex             = PPUnitIndexR - 1;
  NumericVector PPMonthIndex            = PPMonthIndexR - 1;

  NumericVector PPNumStatesUnit         = PlantParameters["NumStates"];
  
  NumericVector PPMinGeneration         = PlantParameters["MinGeneration"];
  NumericVector PPMaxGeneration         = PlantParameters["MaxGeneration"];
  NumericVector PPFuelPrice             = PlantParameters["FuelPrice"];
  NumericVector PPFuelPriceAdder        = PlantParameters["FuelPriceAdder"];
  NumericVector PPStartFuelPrice        = PlantParameters["StartFuelPrice"];
  NumericVector PPStartFuelPriceAdder   = PlantParameters["StartFuelPriceAdder"];
  NumericVector PPVOM                   = PlantParameters["VOM"];
  
  NumericVector PPRampUpRate            = PlantParameters["RampUpRate"];
  NumericVector PPRampDownRate          = PlantParameters["RampDownRate"];
  
  NumericVector PPMaxRegulationUp       = PlantParameters["AnServiceCapMax1"];
  NumericVector PPMaxRegulationDown     = PlantParameters["AnServiceCapMax2"];
  
  NumericVector PPPowerPriceFactor      = PlantParameters["PowerPriceFactor"];
  
  int num_rows = PPUnitIndex.length();

  // create initial empty variables
  NumericMatrix NumStatesUnit(num_months, num_units);
  
  NumericMatrix MinGeneration(num_months, num_units);
  NumericMatrix MaxGeneration(num_months, num_units);
  NumericMatrix FuelPrice(num_months, num_units);
  NumericMatrix FuelPriceAdder(num_months, num_units);
  NumericMatrix StartFuelPrice(num_months, num_units);
  NumericMatrix StartFuelPriceAdder(num_months, num_units);
  NumericMatrix VOM(num_months, num_units);
  NumericMatrix RampUpRate(num_months, num_units);
  NumericMatrix RampDownRate(num_months, num_units);
  NumericMatrix MaxRegulationUp(num_months, num_units);
  NumericMatrix MaxRegulationDown(num_months, num_units);
  NumericMatrix PowerPriceFactor(num_months, num_units);
  
  arma::cube StateGenerationMonth     = arma::zeros<arma::cube>(num_months, max_num_states, num_units);
  arma::cube StateHeatRateMonth       = arma::zeros<arma::cube>(num_months, max_num_states, num_units);
  arma::cube StateRegulationUpMonth   = arma::zeros<arma::cube>(num_months, max_num_states, num_units);
  arma::cube StateRegulationDownMonth = arma::zeros<arma::cube>(num_months, max_num_states, num_units);
  
  // fill in parameters values by unit
  for (int k = 0; k < num_rows; ++k) {

    int i = PPMonthIndex(k);
    int j = PPUnitIndex(k);
    
    NumStatesUnit(i,j)        = PPNumStatesUnit(k);
    
    MinGeneration(i,j)        = PPMinGeneration(k);
    MaxGeneration(i,j)        = PPMaxGeneration(k);
    FuelPrice(i,j)            = PPFuelPrice(k);
    FuelPriceAdder(i,j)       = PPFuelPriceAdder(k);
    StartFuelPrice(i,j)       = PPStartFuelPrice(k);
    StartFuelPriceAdder(i,j)  = PPStartFuelPriceAdder(k);
    VOM(i,j)                  = PPVOM(k);
    RampUpRate(i,j)           = PPRampUpRate(k);
    RampDownRate(i,j)         = PPRampDownRate(k);
    
    MaxRegulationUp(i,j)      = PPMaxRegulationUp(k);
    MaxRegulationDown(i,j)    = PPMaxRegulationDown(k);
    PowerPriceFactor(i,j)     = PPPowerPriceFactor(k);

    // state generation and state heat rate
    int num_states = NumStatesUnit(i,j);
    for (int m = 0; m < num_states; ++m) {
      
      StateGenerationMonth(i,m,j)   = StateGenerationMatrix(k,m);
      StateHeatRateMonth(i,m,j)     = StateHeatRateMatrix(k,m);
      
    }
    
    // state regulation up and down
    for (int m = 1; m < num_states; ++m) {
      
      double a1 = MaxRegulationUp(i,j);
      double a2 = MaxGeneration(i,j);
      double a3 = StateGenerationMonth(i,m,j);
      StateRegulationUpMonth(i,m,j)      = std::min(a1, a2 - a3);
      
      double b1 = MaxRegulationDown(i,j);
      double b2 = MinGeneration(i,j);
      double b3 = StateGenerationMonth(i,m,j);
      
      StateRegulationDownMonth(i,m,j)    = std::min(b1, b3 - b2);
      
    }

  }
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //                                                                                              //
  // Calculate Cost by State                                                                      //
  //                                                                                              //
  //////////////////////////////////////////////////////////////////////////////////////////////////


  // solution matrices (by hour, unit, and sim)
  arma::cube StateGeneration                = arma::zeros<arma::cube>(num_hours, max_num_states, num_units);
  arma::cube StateHeatRate                  = arma::zeros<arma::cube>(num_hours, max_num_states, num_units);
  arma::cube StateCost                      = arma::zeros<arma::cube>(num_hours, max_num_states, num_units);
  arma::cube MarginalStateGeneration        = arma::zeros<arma::cube>(num_hours, max_num_states, num_units);
  arma::cube MarginalStateCost              = arma::zeros<arma::cube>(num_hours, max_num_states, num_units);
  arma::cube IncrementalCost                = arma::zeros<arma::cube>(num_hours, max_num_states, num_units);
  arma::cube DecrementalCost                = arma::zeros<arma::cube>(num_hours, max_num_states, num_units);
  
  arma::cube StateRegulationUp              = arma::zeros<arma::cube>(num_hours, max_num_states, num_units);
  arma::cube MarginalStateRegulationUp      = arma::zeros<arma::cube>(num_hours, max_num_states, num_units);
  arma::cube IncrementalCostRegulationUp    = arma::zeros<arma::cube>(num_hours, max_num_states, num_units);
  arma::cube DecrementalCostRegulationUp    = arma::zeros<arma::cube>(num_hours, max_num_states, num_units);

  arma::cube StateRegulationDown            = arma::zeros<arma::cube>(num_hours, max_num_states, num_units);
  arma::cube MarginalStateRegulationDown    = arma::zeros<arma::cube>(num_hours, max_num_states, num_units);
  arma::cube IncrementalCostRegulationDown  = arma::zeros<arma::cube>(num_hours, max_num_states, num_units);
  // arma::cube DecrementalCostRegulationDown  = arma::zeros<arma::cube>(num_hours, max_num_states, num_units);
  
  for (int i = 0; i < num_hours; ++i) {
    
    int month_i = MonthIndex(i);

    for (int k = 0; k < num_units; ++k) {

      int num_states = NumStatesUnit(month_i, k);

      if ((CommitMatrix(i,k) == 1) || (AlwaysOn(k) == 1)) {
        
        // calculate cost, generation, and regulation
        for (int j = 0; j < num_states; ++j) {
          
          StateGeneration(i,j,k)  = StateGenerationMonth(month_i,j,k);
          StateHeatRate(i,j,k)    = StateHeatRateMonth(month_i,j,k);
          StateCost(i,j,k)        = ((FuelPrice(month_i,k) + FuelPriceAdder(month_i,k)) * 
            StateHeatRateMonth(month_i,j,k) + VOM(month_i,k)) * StateGeneration(i,j,k);
          
          // add price for purchases
          StateCost(i,j,k)        = StateCost(i,j,k) + 
            PowerPriceFactor(month_i,k) * PowerPriceVector(i) * StateGeneration(i,j,k);
          
          StateRegulationUp(i,j,k)    = StateRegulationUpMonth(month_i,j,k);
          StateRegulationDown(i,j,k)  = StateRegulationDownMonth(month_i,j,k);
          
        }
        
        for (int j = 0; j < (num_states - 1); ++j) {
          
          MarginalStateGeneration(i,j,k)      = StateGeneration(i,j+1,k) - StateGeneration(i,j,k);
          MarginalStateCost(i,j,k)            = StateCost(i,j+1,k) - StateCost(i,j,k);
          MarginalStateRegulationUp(i,j,k)    = StateRegulationUp(i,j+1,k) - StateRegulationUp(i,j,k);
          MarginalStateRegulationDown(i,j,k)  = StateRegulationDown(i,j+1,k) - StateRegulationDown(i,j,k);
          
          // calculate marginal cost of generation ($/MWh generation)
          if (MarginalStateGeneration(i,j,k) > 0) {
            IncrementalCost(i,j,k)        = MarginalStateCost(i,j,k) / MarginalStateGeneration(i,j,k);
            DecrementalCost(i,j+1,k)      = IncrementalCost(i,j,k);
          }
          

          // calculate marginal cost of regulation down ($/MWh regulation down)
          if (MarginalStateRegulationDown(i,j,k) > 0) {
            IncrementalCostRegulationDown(i,j,k)    = MarginalStateCost(i,j,k) / MarginalStateRegulationDown(i,j,k);
          } else {
            IncrementalCostRegulationDown(i,j,k)   = 9999;
          }
        }
        
        // incremental cost from max state
        IncrementalCost(i, num_states-1, k) = 9999;
        IncrementalCostRegulationDown(i, num_states-1, k) = 9999;
        
        // decremental cost from off state
        DecrementalCost(i,0,k) = 9999;
        
      } else {
        
        for (int j = 0; j < num_states; ++j) {
          IncrementalCost(i,j,k) = 9999;
          DecrementalCost(i,j,k) = 9999;
          
          IncrementalCostRegulationDown(i,j,k) = 9999;
        }
      }
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  //                                                                                              //
  // Determine Unit Generation by Hour                                                            //
  //                                                                                              //
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  // results by unit
  NumericMatrix StateMatrix(num_hours, num_units);
  NumericMatrix Generation(num_hours, num_units);
  NumericMatrix HeatRate(num_hours, num_units);
  NumericMatrix RegulationUp(num_hours, num_units);
  NumericMatrix RegulationDown(num_hours, num_units);
  
  // results for entire system
  NumericVector TotalGeneration(num_hours);
  NumericVector TotalRegulationUp(num_hours);
  NumericVector TotalRegulationDown(num_hours);
  NumericVector BuyFlag(num_hours);
  NumericVector SellFlag(num_hours);
  
  NumericVector BuyMW(num_hours);
  NumericVector SellMW(num_hours);
  
  NumericVector MarginalCost(num_hours);
  
  NumericMatrix IncrementalMatrix(num_hours, num_units);
  NumericMatrix DecrementalMatrix(num_hours, num_units);
  
  NumericMatrix TrimGeneration(num_hours, num_units);
  NumericMatrix UpGeneration(num_hours, num_units);
  
  NumericMatrix IncrementalMatrixRegulationDown(num_hours, num_units);
  
  
  for (int i = 0; i < num_hours; ++i) {
    
    // set initial generation levels
    for (int k = 0; k < num_units; ++k) {
      
      StateMatrix(i,k)        = CommitMatrix(i,k);
      Generation(i,k)         = StateGeneration(i, StateMatrix(i,k), k);
      HeatRate(i,k)           = StateHeatRate(i, StateMatrix(i,k), k);
      RegulationUp(i,k)       = StateRegulationUp(i, StateMatrix(i,k), k);
      RegulationDown(i,k)     = StateRegulationDown(i, StateMatrix(i,k), k);
      
      IncrementalMatrix(i,k)  = IncrementalCost(i, StateMatrix(i,k), k);
      DecrementalMatrix(i,k)  = DecrementalCost(i, StateMatrix(i,k), k);
      
      TotalGeneration(i)      = TotalGeneration(i) + Generation(i,k);
      TotalRegulationUp(i)    = TotalRegulationUp(i) + RegulationUp(i,k);
      TotalRegulationDown(i)  = TotalRegulationDown(i) + RegulationDown(i,k);
      
      IncrementalMatrixRegulationDown(i,k) = IncrementalCostRegulationDown(i, StateMatrix(i,k), k);
      
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                              //
    // Meet Regulation Down Constraint                                                              //
    //                                                                                              //
    //////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    // adjust for additional regulation down requirement
    NumericVector IncrementalVectorRegulation = IncrementalMatrixRegulationDown(i,_);
    
    int iteration = 0;
    
    while ((TotalRegulationDown(i) < 200) && iteration < 10) {
      
      iteration       = iteration + 1;
      double min_cost = min(IncrementalVectorRegulation);
      
      // check for upper bound
      if (min_cost == 9999) {break;}
      
      // proposed parameters to change
      int k                   = which_min(IncrementalVectorRegulation);
      int current_state       = StateMatrix(i,k);
      int new_state           = current_state + 1;
      
      double current_gen      = StateGeneration(i,current_state,k);
      double new_gen          = StateGeneration(i,new_state,k);
      
      double current_reg_up   = StateRegulationUp(i, current_state, k);
      double new_reg_up       = StateRegulationUp(i, new_state, k);
      double current_reg_down = StateRegulationDown(i, current_state, k);
      double new_reg_down     = StateRegulationDown(i, new_state, k);
      
      // compare to regulation up constraint
      if (TotalRegulationUp(i) >= 200) {
        
        if ((TotalRegulationUp(i) + new_reg_up - current_reg_up) < 200) {
          
          // set incremental cost values to rule out unit
          IncrementalVectorRegulation(k)        = 9999;
          IncrementalMatrixRegulationDown(i,k)  = 9999;
          
          // skip to next unit
          continue;
          
        }
      }
      
      // update generation
      StateMatrix(i,k)        = new_state;
      Generation(i,k)         = new_gen;
      HeatRate(i,k)           = StateHeatRate(i, new_state, k);
      TotalGeneration(i)      = TotalGeneration(i) + new_gen - current_gen;
      UpGeneration(i,k)       = UpGeneration(i,k) + new_gen - current_gen;
      
      // update regulation
      RegulationUp(i,k)       = new_reg_up;
      RegulationDown(i,k)     = new_reg_down;
      TotalRegulationUp(i)    = TotalRegulationUp(i) + new_reg_up - current_reg_up;
      TotalRegulationDown(i)  = TotalRegulationDown(i) + new_reg_down - current_reg_down;
      
      // update incremental and decremental costs
      IncrementalMatrix(i,k)  = IncrementalCost(i,new_state,k);
      DecrementalMatrix(i,k)  = DecrementalCost(i,new_state,k);
      
      IncrementalVectorRegulation(k)        = IncrementalCostRegulationDown(i,new_state,k);
      IncrementalMatrixRegulationDown(i,k)  = IncrementalCostRegulationDown(i,new_state,k);
      
      // set marginal cost
      MarginalCost(i) = IncrementalCost(i,new_state,k);
      
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                              //
    // Meet Load Constraint                                                                         //
    //                                                                                              //
    //////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    // initial value for incremental costs across units for hour
    NumericVector IncrementalVector = IncrementalMatrix(i,_);
    NumericVector DecrementalVector = IncrementalMatrix(i,_);
    
    while (TotalGeneration(i) < LoadVector(i)) {

      double min_cost         = min(IncrementalVector);
      
      // check for upper bound
      if (min_cost == 9999) {break;}
      
      // proposed parameters to change
      int k                   = which_min(IncrementalVector);
      int current_state       = StateMatrix(i,k);
      int new_state           = current_state + 1;

      double current_gen      = StateGeneration(i,current_state,k);
      double new_gen          = StateGeneration(i,new_state,k);

      double current_reg_up   = StateRegulationUp(i, current_state, k);
      double new_reg_up       = StateRegulationUp(i, new_state, k);
      double current_reg_down = StateRegulationDown(i, current_state, k);
      double new_reg_down     = StateRegulationDown(i, new_state, k);
      
      // compare to regulation up constraint
      if (TotalRegulationUp(i) >= 200) {
        
        if ((TotalRegulationUp(i) + new_reg_up - current_reg_up) < 200) {
          
          // set incremental cost values to rule out unit
          IncrementalVector(k)    = 9999;
          IncrementalMatrix(i,k)  = 9999;
          
          // skip to next unit
          continue;
          
        }
      }
      
      // update generation
      StateMatrix(i,k)        = new_state;
      Generation(i,k)         = new_gen;
      HeatRate(i,k)           = StateHeatRate(i, new_state, k);
      TotalGeneration(i)      = TotalGeneration(i) + new_gen - current_gen;
      UpGeneration(i,k)       = UpGeneration(i,k) + new_gen - current_gen;
      
      // update regulation
      RegulationUp(i,k)       = new_reg_up;
      RegulationDown(i,k)     = new_reg_down;
      TotalRegulationUp(i)    = TotalRegulationUp(i) + new_reg_up - current_reg_up;
      TotalRegulationDown(i)  = TotalRegulationDown(i) + new_reg_down - current_reg_down;
      
      // update incremental and decremental costs
      IncrementalVector(k)    = IncrementalCost(i,new_state,k);
      IncrementalMatrix(i,k)  = IncrementalCost(i,new_state,k);
      DecrementalVector(k)    = DecrementalCost(i,new_state,k);
      DecrementalMatrix(i,k)  = DecrementalCost(i,new_state,k);
      
      IncrementalMatrixRegulationDown(i,k) = IncrementalCostRegulationDown(i,new_state,k);
      
      // trim generation to exactly match load
      if (TotalGeneration(i) > LoadVector(i)) {
        
        double over         = TotalGeneration(i) - LoadVector(i);
        
        TrimGeneration(i,k) = TrimGeneration(i,k) + over;
        Generation(i,k)     = Generation(i,k) - over;
        TotalGeneration(i)  = TotalGeneration(i) - over;
          
      }
      
      // set marginal cost
      MarginalCost(i) = min_cost;

    }
    
  }
  
  

  // combine into single list to return
  List Out;
  
  // solution variables
  Out["StateMatrix"]              = StateMatrix;
  Out["Generation"]               = Generation;
  Out["HeatRate"]                 = HeatRate;
  Out["RegulationUp"]             = RegulationUp;
  Out["RegulationDown"]           = RegulationDown;
  Out["MarginalCost"]             = MarginalCost;
  
  // development variables
  Out["MarginalStateGeneration"]  = MarginalStateGeneration;
  Out["MarginalStateCost"]        = MarginalStateCost;
  
  Out["IncrementalCost"]          = IncrementalCost;
  Out["IncrementalMatrix"]        = IncrementalMatrix;
  Out["DecrementalCost"]          = DecrementalCost;
  Out["DecrementalMatrix"]        = DecrementalMatrix;
  
  // Out["IncrementalCostRegulationUp"]    = IncrementalCostRegulationUp;
  
  
  Out["TotalGeneration"]          = TotalGeneration;
  Out["UpGeneration"]             = UpGeneration;
  Out["TrimGeneration"]           = TrimGeneration;
  Out["BuyFlag"]                  = BuyFlag;
  Out["SellFlag"]                 = SellFlag;
  Out["BuyMW"]                    = BuyMW;
  Out["SellMW"]                   = SellMW;
  
  Out["StateRegulationUp"]        = StateRegulationUp;
  Out["StateRegulationDown"]      = StateRegulationDown;

  Out["IncrementalCostRegulationDown"]    = IncrementalCostRegulationDown;
  Out["IncrementalMatrixRegulationDown"]  = IncrementalMatrixRegulationDown;
  
  Out["TotalRegulationUp"]        = TotalRegulationUp;
  Out["TotalRegulationDown"]      = TotalRegulationDown;
  
  return(Out);
  
}
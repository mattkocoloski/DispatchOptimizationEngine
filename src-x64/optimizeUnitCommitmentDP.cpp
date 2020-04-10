// unit commitment optimization engine

// optimizes on/off decisions across multiple units over time, ensuring that range of generation
// converges to load

// simplifies generation by assuming unit is either at max or min generation

// version 2 removes iteration over DTL shadow prices
// promoted back to base

// version 2 separates must run, outage, and stack on/off constraints
// promoted back to base

// version 2 removes min runtime penalty if following hour is outage
// promoted back to base

// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List optimizeUnitCommitmentDP(List OptParameters,
                              DataFrame PlantParameters,
                              DataFrame UnitList,
                              NumericVector PowerPriceVector,
                              NumericVector LoadVector,
                              NumericMatrix MustRunMatrix,
                              NumericMatrix OutageMatrix,
                              NumericVector StackOnMatrix,
                              NumericVector StackOffMatrix,
                              DataFrame HourIndexList,
                              IntegerVector InitialState) {
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //                                                                                              //
  // Process Inputs                                                                               //
  //                                                                                              //
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  // OptParameters values
  String initial_start_temp           = OptParameters["InitialStartTemp"];

  int max_iterations_hard             = OptParameters["MaxIterationsHard"];

  int min_downtime_penalty            = OptParameters["MinDowntimePenalty"];
  int min_runtime_penalty             = OptParameters["MinRuntimePenalty"];

  double min_downtime_decay_rate      = OptParameters["MinDowntimeDecayRate"];
  double min_runtime_decay_rate       = OptParameters["MinRuntimeDecayRate"];

  int must_run_penalty                = OptParameters["MustRunPenalty"];
  int outage_penalty                  = OptParameters["OutagePenalty"];
  
  int stack_on_penalty                = OptParameters["StackOnPenalty"];
  // int stack_off_penalty               = OptParameters["StackOffPenalty"];
  
  // create year, month, and day index variables
  NumericVector YearIndexR            = HourIndexList["YearIndex"];
  NumericVector MonthIndexR           = HourIndexList["MonthIndex"];
  NumericVector DayIndexR             = HourIndexList["DayIndex"];
  
  // create unit-specific dummy variables
  IntegerVector IncludeUnit           = UnitList["IncludeUnit"];
  
  // adjust for zero-index C++
  NumericVector YearIndex             = YearIndexR - 1;
  NumericVector MonthIndex            = MonthIndexR - 1;
  NumericVector DayIndex              = DayIndexR - 1;

  // dimensions
  int num_hours                       = HourIndexList.nrow();
  int num_months                      = max(MonthIndex) + 1;
  int num_units                       = IncludeUnit.length();
  int num_states                      = 2;
  
  
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

  NumericVector PPMinGeneration         = PlantParameters["MinGeneration"];
  NumericVector PPMinGenHeatRate        = PlantParameters["MinGenHeatRate"];
  NumericVector PPMaxGeneration         = PlantParameters["MaxGeneration"];
  NumericVector PPMaxGenHeatRate        = PlantParameters["MaxGenHeatRate"];

  NumericVector PPFuelPrice             = PlantParameters["FuelPrice"];
  NumericVector PPFuelPriceAdder        = PlantParameters["FuelPriceAdder"];
  NumericVector PPStartFuelPrice        = PlantParameters["StartFuelPrice"];
  NumericVector PPStartFuelPriceAdder   = PlantParameters["StartFuelPriceAdder"];

  NumericVector PPHotStartCostFixed     = PlantParameters["HotStartCostFixed"];
  NumericVector PPWarmStartCostFixed    = PlantParameters["WarmStartCostFixed"];
  NumericVector PPColdStartCostFixed    = PlantParameters["ColdStartCostFixed"];

  NumericVector PPHotStartFuel          = PlantParameters["HotStartFuel"];
  NumericVector PPWarmStartFuel         = PlantParameters["WarmStartFuel"];
  NumericVector PPColdStartFuel         = PlantParameters["ColdStartFuel"];

  NumericVector PPHotWarmStartSplit     = PlantParameters["HotWarmStartSplit"];
  NumericVector PPWarmColdStartSplit    = PlantParameters["WarmColdStartSplit"];

  NumericVector PPVOM                   = PlantParameters["VOM"];
  NumericVector PPMinRuntime            = PlantParameters["MinRuntime"];
  NumericVector PPMinDowntime           = PlantParameters["MinDowntime"];

  NumericVector PPPowerPriceFactor      = PlantParameters["PowerPriceFactor"];

  int num_rows = PPUnitIndex.length();

  // create initial empty variables
  NumericMatrix NumStatesUnit(num_months, num_units);
  NumericMatrix MinGeneration(num_months, num_units);
  NumericMatrix MinGenHeatRate(num_months, num_units);
  NumericMatrix MaxGeneration(num_months, num_units);
  NumericMatrix MaxGenHeatRate(num_months, num_units);

  NumericMatrix FuelPrice(num_months, num_units);
  NumericMatrix FuelPriceAdder(num_months, num_units);
  NumericMatrix StartFuelPrice(num_months, num_units);
  NumericMatrix StartFuelPriceAdder(num_months, num_units);

  NumericMatrix HotStartCostFixed(num_months, num_units);
  NumericMatrix WarmStartCostFixed(num_months, num_units);
  NumericMatrix ColdStartCostFixed(num_months, num_units);

  NumericMatrix HotStartFuel(num_months, num_units);
  NumericMatrix WarmStartFuel(num_months, num_units);
  NumericMatrix ColdStartFuel(num_months, num_units);

  NumericMatrix HotStartTime(num_months, num_units);
  NumericMatrix WarmStartTime(num_months, num_units);
  NumericMatrix ColdStartTime(num_months, num_units);

  NumericMatrix HotWarmStartSplit(num_months, num_units);
  NumericMatrix WarmColdStartSplit(num_months, num_units);

  NumericMatrix VOM(num_months, num_units);

  NumericMatrix MinRuntime(num_months, num_units);
  NumericMatrix MinDowntime(num_months, num_units);

  NumericMatrix PowerPriceFactor(num_months, num_units);
  
  // fill in parameters values by unit
  for (int k = 0; k < num_rows; ++k) {

    int i = PPMonthIndex(k);
    int j = PPUnitIndex(k);
    
    MinGeneration(i,j)       = PPMinGeneration(k);
    MinGenHeatRate(i,j)      = PPMinGenHeatRate(k);
    MaxGeneration(i,j)       = PPMaxGeneration(k);
    MaxGenHeatRate(i,j)      = PPMaxGenHeatRate(k);

    FuelPrice(i,j)           = PPFuelPrice(k);
    FuelPriceAdder(i,j)      = PPFuelPriceAdder(k);
    StartFuelPrice(i,j)      = PPStartFuelPrice(k);
    StartFuelPriceAdder(i,j) = PPStartFuelPriceAdder(k);

    HotStartCostFixed(i,j)   = PPHotStartCostFixed(k);
    WarmStartCostFixed(i,j)  = PPWarmStartCostFixed(k);
    ColdStartCostFixed(i,j)  = PPColdStartCostFixed(k);

    HotStartFuel(i,j)        = PPHotStartFuel(k);
    WarmStartFuel(i,j)       = PPWarmStartFuel(k);
    ColdStartFuel(i,j)       = PPColdStartFuel(k);

    // HotStartTime(i,j)        = PPHotStartTime(k);
    // WarmStartTime(i,j)       = PPWarmStartTime(k);
    // ColdStartTime(i,j)       = PPColdStartTime(k);

    HotWarmStartSplit(i,j)   = PPHotWarmStartSplit(k);
    WarmColdStartSplit(i,j)  = PPWarmColdStartSplit(k);

    VOM(i,j)                 = PPVOM(k);

    MinRuntime(i,j)          = PPMinRuntime(k);
    MinDowntime(i,j)         = PPMinDowntime(k);

    PowerPriceFactor(i,j)    = PPPowerPriceFactor(k);

  }
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //                                                                                              //
  // Create Solution Matrices and Lists                                                           //
  //                                                                                              //
  //////////////////////////////////////////////////////////////////////////////////////////////////


  // solution matrices (by hour, unit, and sim)
  NumericMatrix OnGeneration(num_hours, num_units);
  NumericMatrix OnHeatRate(num_hours, num_units);
  
  NumericMatrix On(num_hours, num_units);
  NumericMatrix PathGeneration(num_hours, num_units);
  NumericMatrix PathHeatRate(num_hours, num_units);
  NumericMatrix Start(num_hours, num_units);
  NumericMatrix Stop(num_hours, num_units);
  NumericMatrix HoursSinceStart(num_hours, num_units);
  NumericMatrix HoursSinceStop(num_hours, num_units);
  NumericMatrix StartTime(num_hours, num_units);
  NumericMatrix StartFuelMMBtu(num_hours, num_units);
  NumericMatrix StartCostFuel(num_hours, num_units);
  NumericMatrix StartCostFixed(num_hours, num_units);

  // matrices for tracking total generation for dispatch to load
  NumericVector TotalPathGeneration(num_hours);
  NumericVector MinOnGeneration(num_hours);
  NumericVector MaxOnGeneration(num_hours);
  NumericVector GenerationInRange(num_hours);
  
  // load buffers for market purchases and sales
  NumericVector BuyBuffer(num_hours);
  NumericVector SellBuffer(num_hours);
  NumericVector DumpBuffer(num_hours);
  
  // tracking number of iterations
  NumericVector NumIterationsUnit(num_units);

  // variables for tracking periods to optimize by unit
  NumericMatrix StartHour(num_hours, num_units);
  NumericMatrix EndHour(num_hours, num_units);
  NumericVector NumPeriods(num_units);
  
  for (int h = 0; h < num_units; ++h) {
    StartHour(0,h) = 0;
    EndHour(0,h)   = num_hours - 1;
    NumPeriods(h)  = 1;
  }
  
  // variables for tracking shadow prices across units
  NumericMatrix SPStartMinDowntime(num_hours, num_units);
  NumericMatrix SPStopMinRuntime(num_hours, num_units);
  
  // dynamic programming variables
  arma::cube FutureBestStateCube  = arma::zeros<arma::cube>(num_hours, num_states, num_units);
  arma::cube FutureBestValueCube  = arma::zeros<arma::cube>(num_hours, num_states, num_units);
  NumericMatrix TransitionCost(num_states, num_states);
  NumericVector PathSum(num_units);

  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //                                                                                              //
  // Calculate State Margins, Start Costs, and Must Run / Locked States                           //
  //                                                                                              //
  //////////////////////////////////////////////////////////////////////////////////////////////////


  //// Calculate State Margins Across Hours, States, and Units (Excluding Shadow Prices) ////
  
  arma::cube StateMargin = arma::zeros<arma::cube>(num_hours, num_states, num_units);
  
  // cost and margin assuming that unit is on and chooses
  NumericMatrix MinGenFuelCostPerMWh(num_hours, num_units);
  NumericMatrix MaxGenFuelCostPerMWh(num_hours, num_units);
  
  NumericMatrix MinGenCost(num_hours, num_units);
  NumericMatrix MaxGenCost(num_hours, num_units);
  
  NumericMatrix MinGenMargin(num_hours, num_units);
  NumericMatrix MaxGenMargin(num_hours, num_units);
  
  // states that are locked by constraints
  NumericMatrix StateIsLocked(num_hours, num_units);
  
  // calculate StateMargin
  for (int k = 0; k < num_units; ++k) {
    
    // check for units that are excluded from optimization
    if (IncludeUnit(k) == 0) {continue;}
    
    for (int i = 0; i < num_hours; ++i) {
      
      int month_i = MonthIndex(i);
      
      // calculate fuel cost
      MinGenFuelCostPerMWh(i,k) = (FuelPrice(month_i,k) + FuelPriceAdder(month_i,k)) * MinGenHeatRate(month_i,k);
      MaxGenFuelCostPerMWh(i,k) = (FuelPrice(month_i,k) + FuelPriceAdder(month_i,k)) * MaxGenHeatRate(month_i,k);

      // calculate total cost
      MinGenCost(i,k) = (MinGenFuelCostPerMWh(i,k) + VOM(month_i,k) + 
        PowerPriceFactor(month_i,k) * PowerPriceVector(i)) * MinGeneration(month_i,k);
      MaxGenCost(i,k) = (MaxGenFuelCostPerMWh(i,k) + VOM(month_i,k) + 
        PowerPriceFactor(month_i,k) * PowerPriceVector(i)) * MaxGeneration(month_i,k);
      
      // calculate running margin
      // MinGenMargin(i,k) = PowerPriceVector(i) * MinGeneration(month_i,k) - MinGenCost(i,k);
      // MaxGenMargin(i,k) = PowerPriceVector(i) * MaxGeneration(month_i,k) - MaxGenCost(i,k);
      MinGenMargin(i,k) = -MinGenCost(i,k);
      MaxGenMargin(i,k) = -MaxGenCost(i,k);
      

      // set OnGeneration, OnHeatRate, and StateMargin to higher of min and max values
      if (MinGenMargin(i,k) > MaxGenMargin(i,k)) {

        OnGeneration(i,k)   = MinGeneration(month_i,k);
        OnHeatRate(i,k)     = MinGenHeatRate(month_i,k);
        StateMargin(i,1,k)  = MinGenMargin(i,k);

      } else {

        OnGeneration(i,k)   = MaxGeneration(month_i,k);
        OnHeatRate(i,k)     = MaxGenHeatRate(month_i,k);
        StateMargin(i,1,k)  = MaxGenMargin(i,k);

      }
      
      // adjust for outage hours
      if (OutageMatrix(i,k) == 1) {
        
        MinGenCost(i,k)     = outage_penalty;
        MaxGenCost(i,k)     = outage_penalty;
        
        MinGenMargin(i,k)   = -outage_penalty;
        MaxGenMargin(i,k)   = -outage_penalty;
        
        StateMargin(i,1,k)  = -outage_penalty;
        StateIsLocked(i,k)  = 1;
        On(i,k)             = 0;
        
      }
      
      // // adjust for stack off hours
      // if (StackOffMatrix(i,k) == 1) {
      //   
      //   MinGenCost(i,k)     = stack_off_penalty;
      //   MaxGenCost(i,k)     = stack_off_penalty;
      //   
      //   MinGenMargin(i,k)   = -stack_off_penalty;
      //   MaxGenMargin(i,k)   = -stack_off_penalty;
      //   
      //   StateMargin(i,1,k)  = -stack_off_penalty;
      //   On(i,k)             = 0;
      //   
      // }
      
      // adjust for must run hours
      if (MustRunMatrix(i,k) == 1) {
        
        StateMargin(i,0,k)  = -must_run_penalty;
        StateIsLocked(i,k)  = 1;
        On(i,k)             = 1;
        
        OnGeneration(i,k)   = MaxGeneration(month_i,k);
        OnHeatRate(i,k)     = MaxGenHeatRate(month_i,k);
        PathGeneration(i,k) = OnGeneration(i,k);
        
      }
      
      // adjust for stack on hours
      if (StackOnMatrix(i,k) == 1) {
        
        StateMargin(i,0,k)  = -stack_on_penalty;
        On(i,k)             = 1;
        
        OnGeneration(i,k)   = MaxGeneration(month_i,k);
        OnHeatRate(i,k)     = MaxGenHeatRate(month_i,k);
        PathGeneration(i,k) = OnGeneration(i,k);
        
      }
    }
  }
  
  
  //// Calculate Start and Stop Variables ////
  
  
  // starts
  for (int i = 0; i < num_hours; ++i) {
    
    for (int h = 0; h < num_units; ++h) {
      
      if (i > 0) {
        
        if (On(i-1,h) == 0 && On(i,h) == 1) {
          Start(i,h) = 1;
        } else {
          Start(i,h) = 0;
        }
        
      } else {
        
        if (On(i,h) == 1) {
          Start(i,h) = 1;
        } else {
          Start(i,h) = 0;
        }
        
      }
    }
  }
  
  // stops
  for (int i = 0; i < num_hours; ++i) {
    for (int h = 0; h < num_units; ++h) {
      if (i > 0) {
        if (On(i-1,h) == 1 && On(i,h) == 0) {
          Stop(i,h) = 1;
        } else {
          Stop(i,h) = 0;
        }
      }
    }
  }
  
  
  //// Set Initial StartCost and StartTime Values Across Hours and Units ////
  
  
  if (initial_start_temp == "Hot") {
    
    for (int i = 0; i < num_hours; ++i) {
      
      int month_i = MonthIndex(i);
      
      for (int k = 0; k < num_units; ++k) {
        
        StartTime(i,k)      = HotStartTime(month_i,k);
        StartFuelMMBtu(i,k) = HotStartFuel(month_i,k);
        StartCostFuel(i,k)  = HotStartFuel(month_i,k) * (StartFuelPrice(month_i,k) + StartFuelPriceAdder(month_i,k));
        StartCostFixed(i,k) = HotStartCostFixed(month_i,k);
        
      }
    }
    
  } else if (initial_start_temp == "Warm") {
    
    for (int i = 0; i < num_hours; ++i) {
      
      int month_i = MonthIndex(i);
      
      for (int k = 0; k < num_units; ++k) {
        
        StartTime(i,k)      = WarmStartTime(month_i,k);
        StartFuelMMBtu(i,k) = WarmStartFuel(month_i,k);
        StartCostFuel(i,k)  = WarmStartFuel(month_i,k) * (StartFuelPrice(month_i,k) + StartFuelPriceAdder(month_i,k));
        StartCostFixed(i,k) = WarmStartCostFixed(month_i,k);
        
      }
    }
    
  } else if (initial_start_temp == "Cold") {
    
    for (int i = 0; i < num_hours; ++i) {
      
      int month_i = MonthIndex(i);
      
      for (int k = 0; k < num_units; ++k) {
        
        StartTime(i,k)      = ColdStartTime(month_i,k);
        StartFuelMMBtu(i,k) = ColdStartFuel(month_i,k);
        StartCostFuel(i,k)  = ColdStartFuel(month_i,k) * (StartFuelPrice(month_i,k) + StartFuelPriceAdder(month_i,k));
        StartCostFixed(i,k) = ColdStartCostFixed(month_i,k);
        
      }
    }
  }
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //                                                                                              //
  // Determine Periods to Optimize                                                                //
  //                                                                                              //
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  // // only reoptimize over unit-hours that are not locked
  // for (int k = 0; k < num_units; ++k) {
  //   
  //   // check for included unit
  //   if (IncludeUnit(k) == 0) {continue;}
  //   
  //   // reset number of periods
  //   NumPeriods(k) = 0;
  //   int is_locked = 1;
  //   
  //   for (int i = 0; i < num_hours; ++i) {
  //     
  //     if (StateIsLocked(i,k) == 0) {
  //       
  //       // if previous hour is locked, create new period and set start hour
  //       if (is_locked == 1) {
  //         NumPeriods(k) = NumPeriods(k) + 1;
  //         StartHour(NumPeriods(k) - 1, k) = i;
  //         is_locked = 0;
  //       }
  //       
  //       // update end hour
  //       EndHour(NumPeriods(k) - 1, k) = i;
  //       
  //     } else {
  //       
  //       // set variable to locked
  //       is_locked = 1;
  //       
  //     }
  //   }
  // }
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //                                                                                              //
  // Begin Iterating Over Units                                                                   //
  //                                                                                              //
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  for (int h = 0; h < num_units; ++h) {
    
    // variables to track progress and request reiteration for unit
    int unit_iteration             = 0;
    int same_solution_count        = 0;
    int max_same_solution_count    = 10;
    
    // int new_solution_reiterate     = 1;
    int min_runtime_reiterate      = 0;
    int min_downtime_reiterate     = 0;
    
    // check for units that are excluded from optimization
    if (IncludeUnit(h) == 0) {continue;}
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                              //
    // Begin Tuning Start Costs and Unit Shadow Prices                                              //
    //                                                                                              //
    //////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    // only reiterate for constraint violations
    while (unit_iteration         == 0 ||
           min_runtime_reiterate  == 1 ||
           min_downtime_reiterate == 1) {
      
      // check for hard conditions to exit
      if (unit_iteration >= max_iterations_hard) {break;}
      if (same_solution_count >= max_same_solution_count) {break;}
      
      // reset reiteration switches
      min_runtime_reiterate      = 0;
      min_downtime_reiterate     = 0;
      
      // create dynamic program variables across states (update each hour)
      NumericVector FutureBestValue(num_states);
      NumericVector FutureBestIter(num_states);
      NumericMatrix FutureBestState(num_hours, num_states);
      
      
      //////////////////////////////////////////////////////////////////////////////////////////////////
      //                                                                                              //
      // Begin Looping Over Periods                                                                   //
      //                                                                                              //
      //////////////////////////////////////////////////////////////////////////////////////////////////
      
      
      int num_periods   = NumPeriods(h);
      int initial_state = 0;
      
      for (int p = 0; p < num_periods; ++p) {
        
        int first = StartHour(p,h);
        int last  = EndHour(p,h);
        
        // update initial state
        if (first == 0) {
          initial_state = InitialState(h);
        } else {
          initial_state = On(first-1, h);
        }
        
        // update future best value for destination state
        if (last < (num_hours - 1)) {
          if (On(last+1,h) == 1) {
            
            FutureBestValue(0) = -999999999;
            FutureBestValue(1) = 999999999;
            
          } else {
            
            FutureBestValue(0) = 999999999;
            FutureBestValue(1) = -999999999;
            
          }
        }
        
        
        //////// Run Optimization Over Hours in Time Period //////
        
        
        // perform backward pass for optimization model
        for (int i = last; i >= first; --i) {
          
          // update startup costs to include shadow prices
          TransitionCost(0,1) = StartCostFuel(i,h) + StartCostFixed(i,h) + SPStartMinDowntime(i,h);

          // update shutdown costs to include shadow prices
          TransitionCost(1,0) = SPStopMinRuntime(i,h);

          // run optimization
          for (int j = 0; j < num_states; ++j) {
            
            // calculate costs for all future states
            NumericVector YY(num_states);
            
            for (int k = 0; k < num_states; ++k) {
              YY(k) = StateMargin(i,j,h) - TransitionCost(j,k) + FutureBestValue(k);
            }
            
            // select best future state for each current state
            FutureBestIter(j)          = max(YY);
            FutureBestState(i,j)       = which_max(YY);
            FutureBestStateCube(i,j,h) = which_max(YY);
            
          }
          
          // update future value for each current state
          for (int j = 0; j < num_states; ++j) {
            FutureBestValue(j)          = FutureBestIter(j);
            FutureBestValueCube(i,j,h)  = FutureBestIter(j);
          }
          
        }
        
        // add transition cost from initial state to first period state
        if (first > 0) {
          for (int j = 0; j < num_states; ++j) {
            FutureBestValue(j) = FutureBestValue(j) - TransitionCost(initial_state, j);
            FutureBestValueCube(first,j,h) = FutureBestValueCube(first,j,h) - TransitionCost(initial_state, j);
          }
        }
        
        // select optimal state for first time period
        On(first,h) = which_max(FutureBestValue);
        
        // forward pass to determine optimal path
        for (int i = (first + 1); i <= last; ++i) {
          On(i,h) = FutureBestState(i-1, On(i-1,h));
        }
        
      } // close out loop over time periods
      
      
      //////// Update Generation, Start, Stop, and Run Hour Variables ////////
      
      
      // path generation
      for (int i = 0; i < num_hours; ++i) {
        PathGeneration(i,h) = On(i,h) * OnGeneration(i,h);
      }
      
      // starts
      for (int i = 0; i < num_hours; ++i) {
        
        if (i > 0) {
          
          if (On(i-1,h) == 0 && On(i,h) == 1) {
            Start(i,h) = 1;
          } else {
            Start(i,h) = 0;
          }
          
        } else {
          
          if (On(i,h) == 1) {
            Start(i,h) = 1;
          } else {
            Start(i,h) = 0;
          }
          
        }
      }
      
      // stops
      for (int i = 0; i < num_hours; ++i) {
        if (i > 0) {
          if (On(i-1,h) == 1 && On(i,h) == 0) {
            Stop(i,h) = 1;
          } else {
            Stop(i,h) = 0;
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
        
        if (Start(i,h) == 1) {
          HoursSinceStart(i,h) = 1;
        }
        
        if (Stop(i,h) == 1) {
          HoursSinceStop(i,h) = 1;
        }
      }
      
      
      //////// Recalculate Start Cost and Time by Hour ////////
      
      
      // update start cost total and start time
      for (int i = 0; i < num_hours; ++i) {
        
        int month_i = MonthIndex(i);
        
        if (HoursSinceStop(i,h) < HotWarmStartSplit(month_i,h)) {
          
          StartTime(i,h)      = HotStartTime(month_i,h);
          StartFuelMMBtu(i,h) = HotStartFuel(month_i,h);
          StartCostFuel(i,h)  = HotStartFuel(month_i,h) * (StartFuelPrice(month_i,h) + StartFuelPriceAdder(month_i,h));
          StartCostFixed(i,h) = HotStartCostFixed(month_i,h);
          
        } else if (HoursSinceStop(i,h) < WarmColdStartSplit(month_i,h)) {
          
          StartTime(i,h)      = WarmStartTime(month_i,h);
          StartFuelMMBtu(i,h) = WarmStartFuel(month_i,h);
          StartCostFuel(i,h)  = WarmStartFuel(month_i,h) * (StartFuelPrice(month_i,h) + StartFuelPriceAdder(month_i,h));
          StartCostFixed(i,h) = WarmStartCostFixed(month_i,h);
          
        } else {
          
          StartTime(i,h)      = ColdStartTime(month_i,h);
          StartFuelMMBtu(i,h) = ColdStartFuel(month_i,h);
          StartCostFuel(i,h)  = ColdStartFuel(month_i,h) * (StartFuelPrice(month_i,h) + StartFuelPriceAdder(month_i,h));
          StartCostFixed(i,h) = ColdStartCostFixed(month_i,h);
          
        }
        
      }
      
      
      //////// Update Minimum Runtime/Downtime Shadow Prices ////////
      
      
      // SPStartMinDowntime: penalty for starting within MinDowntime of stop
      for (int i = 0; i < num_hours; ++i) {
        
        int month_i = MonthIndex(i);
        
        // // check to see if unit is constrained to be on in next hour
        // if (i < (num_hours - 1)) {
        //   if ((MustRunMatrix(i+1,h) == 1) || (StackOnMatrix(i+1,h) == 1)) {
        //     continue;
        //   }
        // }
        
        // decay previous penalties
        SPStartMinDowntime(i,h) = SPStartMinDowntime(i,h) * (1 - min_downtime_decay_rate);
        
        // update penalties
        if (HoursSinceStop(i,h) <= (MinDowntime(month_i,h) + StartTime(i,h) - 1)) {
          
          // update value
          // SPStartMinDowntime(i,h) = min_downtime_penalty;
          SPStartMinDowntime(i,h) = min_downtime_penalty * (1 - MustRunMatrix(i,h)) * (1 - StackOnMatrix(i,h));
          
          // check for constraint violations
          if (Start(i,h) == 1) {
            min_downtime_reiterate = 1;
          }
        }
      }
      
      // SPStopMinRuntime: penalty for stopping within MinRuntime of start
      for (int i = 0; i < num_hours; ++i) {
        
        int month_i = MonthIndex(i);
        
        // check to see if unit is constrained to be off in following hour
        if (i < (num_hours - 1)) {
          if (OutageMatrix(i+1,h) == 1) {
            continue;
          }
        }
        
        // decay previous penalties
        SPStopMinRuntime(i,h) = SPStopMinRuntime(i,h) * (1 - min_runtime_decay_rate);
        
        // update penalties
        if (HoursSinceStart(i,h) < MinRuntime(month_i,h)) {
          
          // update value
          SPStopMinRuntime(i,h) = min_runtime_penalty;
          
          // check for constraint violations
          if (Stop(i,h) == 1) {
            min_runtime_reiterate = 1;
          }
        }
      }
      
      
      //////// Check for Solution Convergence ////////
      
      
      // calculate new sum of all states
      int new_sum = 0;
      for (int i = 0; i < num_hours; ++i) {
        new_sum = new_sum + On(i,h);
      }
      
      // check for identical solution to previous iteration
      if (new_sum == PathSum(h)) {
        same_solution_count    = same_solution_count + 1;
      } else {
        same_solution_count    = 0;
        // new_solution_reiterate = 1; // currently does not trigger reiteration
      }
      
      // update path_sum
      PathSum(h) = new_sum;
      
      // increment iteration
      unit_iteration = unit_iteration + 1;
      
      // store number of iterations
      NumIterationsUnit(h) = unit_iteration;
      
    } // close loop tuning start costs and unit constraints
    
  } // close loop iterating over units
  
  
  //// Calculate Total Generation Across Units ////
  
  
  for (int i = 0; i < num_hours; ++i) {
    
    // reset total generation value
    TotalPathGeneration(i) = 0;
    MinOnGeneration(i)     = 0;
    MaxOnGeneration(i)     = 0;
    int month_i            = MonthIndex(i);
    
    // calculate total generation across units
    for (int k = 0; k < num_units; ++k) {
      TotalPathGeneration(i) = TotalPathGeneration(i) + PathGeneration(i,k);
      MinOnGeneration(i)     = MinOnGeneration(i) + On(i,k) * MinGeneration(month_i,k);
      MaxOnGeneration(i)     = MaxOnGeneration(i) + On(i,k) * MaxGeneration(month_i,k);
    }
    
    if ((LoadVector(i) >= MinOnGeneration(i)) && (LoadVector(i) <= MaxOnGeneration(i))) {
      GenerationInRange(i) = 1;
    } else {
      GenerationInRange(i) = 0;
    }
    
  }
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //                                                                                              //
  // Create Output List                                                                           //
  //                                                                                              //
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  // necessary solution parameters
  List Solution;
  Solution["On"]                = On;
  Solution["Generation"]        = PathGeneration;
  Solution["HeatRate"]          = PathHeatRate;
  Solution["Start"]             = Start;
  Solution["Stop"]              = Stop;
  Solution["StartFuelMMBtu"]    = StartFuelMMBtu;
  Solution["StartCostFuel"]     = StartCostFuel;
  Solution["StartCostFixed"]    = StartCostFixed;
  
  Solution["HoursSinceStart"]   = HoursSinceStart;
  Solution["HoursSinceStop"]    = HoursSinceStop;
  
  Solution["TotalPathGeneration"] = TotalPathGeneration;
  Solution["OnGeneration"]        = OnGeneration;
  Solution["MinOnGeneration"]     = MinOnGeneration;
  Solution["MaxOnGeneration"]     = MaxOnGeneration;
  Solution["GenerationInRange"]   = GenerationInRange;
  
  // // add solution diagnostics
  // if (output_solution_diagnostics) {
  //   
  //   Solution["Path"]               = Path;
  //   Solution["Stop"]               = Stop;
  //   Solution["HeatRate"]           = PathHeatRate;
  //   Solution["HoursSinceStart"]    = HoursSinceStart;
  //   Solution["HoursSinceStop"]     = HoursSinceStop;
  //   Solution["StartTime"]          = StartTime;
  //   
  //   Solution["RunHoursYear"]       = RunHoursYearList;
  //   Solution["RunHoursMonth"]      = RunHoursMonthList;
  //   Solution["RunHoursDay"]        = RunHoursDayList;
  //   Solution["StartsYear"]         = StartsYearList;
  //   Solution["StartsMonth"]        = StartsMonthList;
  //   Solution["StartsDay"]          = StartsDayList;
  //   
  // }
  // 
  // // update shadow price list
  // List ShadowPrices;
  // 
  // if (output_sp_min_runtime) {
  //   
  //   ShadowPrices["MinDowntime"]    = SPStartMinDowntimeList;
  //   ShadowPrices["MinRuntime"]     = SPStopMinRuntimeList;
  //   
  // }
  // 
  // // update dispatch to load list
  // List DTLParameters;
  // 
  // if (output_sp_dispatch_to_load) {
  //   
  //   DTLParameters["DispatchToLoad"]       = SPDispatchToLoadList;
  //   DTLParameters["DTLStateLock"]         = SPDTLStateLockList;
  //   DTLParameters["DTLGenerationInRange"] = DTLGenerationInRange;
  //   DTLParameters["DTLStateIsLocked"]     = DTLStateIsLocked;
  //   DTLParameters["DTLTriedSearch"]       = DTLTriedSearch;
  //   DTLParameters["DTLSmoothGenDelta"]    = DTLSmoothGenDelta;
  //   
  // }
  
  // combine into single list to return
  List Out;
  
  Out["Solution"]          = Solution;
  
  Out["StartHour"]    = StartHour;
  Out["EndHour"]      = EndHour;
  Out["NumPeriods"]   = NumPeriods;
  Out["MinGenMargin"] = MinGenMargin;
  Out["MaxGenMargin"] = MaxGenMargin;
  Out["StateMargin"]  = StateMargin;
  
  // Out["NumIterations"] = dtl_iteration;
  Out["BuyBuffer"]     = BuyBuffer;
  Out["SellBuffer"]    = SellBuffer;
  Out["DumpBuffer"]    = DumpBuffer;
  
  Out["SPStopMinRuntime"] = SPStopMinRuntime;
  Out["SPStartMinDowntime"] = SPStartMinDowntime;
  Out["FutureBestValueCube"] = FutureBestValueCube;
  
  // Out["SPDispatchToLoad"] = SPDispatchToLoad;
  
  // Out["ShadowPrices"]      = ShadowPrices;
  // Out["DispatchToLoad"]    = DTLParameters;
  // 
  // Out["NumIterationsUnit"] = NumIterationsUnit;
  // Out["NumIterationsDTL"]  = NumIterationsDTL;
  // 
  // Out["FutureBestState"]   = FutureBestStateList;
  // Out["FutureBestValue"]   = FutureBestValueList;
  
  return(Out);
  
}

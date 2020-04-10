// estimates marginal cost for dispatch to load using stack and cut algorithm

// version 2 initializes cubes with zeros
// promoted back to base

// version 2 recognizes must run hours
// promoted back to base

// version 2 includes market purchases
// promoted back to base

// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List stackCut(DataFrame PlantParameters,
              NumericMatrix LoadMatrix,
              NumericMatrix OutageMatrix,
              NumericMatrix MustRunMatrix,
              NumericVector PowerPrice,
              DataFrame HourIndexList,
              IntegerVector IncludeUnit) {
  
  
  
  //////// Create Parameters Consistent Across Runs ///////
  
  
  
  // create year, month, and day index variables
  NumericVector MonthIndexR          = HourIndexList["MonthIndex"];
  
  // adjust for zero-index C++
  NumericVector MonthIndex           = MonthIndexR - 1;
  
  // dimensions
  int num_hours                      = HourIndexList.nrow();
  int num_sims                       = LoadMatrix.ncol();
  int num_months                     = max(MonthIndex) + 1;
  int num_units                      = IncludeUnit.length();
  
  
  //////// Reorganize Plant Parameters into Matrices ////////
  
  
  // create individual vectors of plant parameter data
  NumericVector PPUnitIndexR        = PlantParameters["UnitIndex"];
  NumericVector PPMonthIndexR       = PlantParameters["MonthIndex"];
  
  NumericVector PPUnitIndex         = PPUnitIndexR - 1;
  NumericVector PPMonthIndex        = PPMonthIndexR - 1;
  
  NumericVector PPFuelPrice         = PlantParameters["FuelPrice"];
  NumericVector PPFuelPriceAdder    = PlantParameters["FuelPriceAdder"];
  NumericVector PPMinGeneration     = PlantParameters["MinCapacity"];
  NumericVector PPMaxGeneration     = PlantParameters["MaxGeneration"];
  NumericVector PPMaxGenHeatRate    = PlantParameters["MaxGenHeatRate"];
  NumericVector PPVOM               = PlantParameters["VOM"];
  NumericVector PPPowerPriceFactor  = PlantParameters["PowerPriceFactor"];
  
  // create initial empty values
  // NumericMatrix MustRun(num_hours, num_units);
  NumericMatrix FuelPrice(num_months, num_units);
  NumericMatrix FuelPriceAdder(num_months, num_units);
  NumericMatrix MinGeneration(num_months, num_units);
  NumericMatrix MaxGeneration(num_months, num_units);
  NumericMatrix MaxGenHeatRate(num_months, num_units);
  NumericMatrix VOM(num_months, num_units);
  NumericMatrix PowerPriceFactor(num_months, num_units);
  
  int num_rows = PPUnitIndex.length();
  
  // fill in parameters values by unit
  for (int k = 0; k < num_rows; ++k) {
    
    int i = PPMonthIndex(k);
    int j = PPUnitIndex(k);
    
    FuelPrice(i,j)        = PPFuelPrice(k);
    FuelPriceAdder(i,j)   = PPFuelPriceAdder(k);
    MinGeneration(i,j)    = PPMinGeneration(k);
    MaxGeneration(i,j)    = PPMaxGeneration(k);
    MaxGenHeatRate(i,j)   = PPMaxGenHeatRate(k);
    VOM(i,j)              = PPVOM(k);
    PowerPriceFactor(i,j) = PPPowerPriceFactor(k);
    
  }
  
  // calculate hourly maximum generation including outages
  NumericMatrix MinGenerationHour(num_hours, num_units);
  NumericMatrix MaxGenerationHour(num_hours, num_units);

  for (int i = 0; i < num_hours; ++i) {

    int month_i = MonthIndex(i);

    for (int k = 0; k < num_units; ++k) {
      MinGenerationHour(i,k) = MinGeneration(month_i, k) * (1 - OutageMatrix(i,k));
      MaxGenerationHour(i,k) = MaxGeneration(month_i, k) * (1 - OutageMatrix(i,k));
    }
  }
  
  // solution matrices
  NumericMatrix MarginalUnitCost(num_hours, num_sims);
  NumericMatrix TotalCapacity(num_hours, num_sims);
  arma::cube StackCutCommit     = arma::zeros<arma::cube>(num_hours, num_units, num_sims);
  arma::cube StackCumulativeMW  = arma::zeros<arma::cube>(num_hours, num_units, num_sims);

  // dispatch cost to save during development
  arma::cube DispatchCostSave   = arma::zeros<arma::cube>(num_hours, num_units, num_sims);

  // iterate over price/load simulation
  for (int n = 0; n < num_sims; ++n) {

    // iterate over hours
    for (int i = 0; i < num_hours; ++i) {

      int month_i = MonthIndex(i);

      // calculate dispatch cost by unit
      NumericVector FuelCostPerMWh(num_units);
      NumericVector EnvironmentalCostPerMWh(num_units);
      NumericVector AverageStartCost(num_units);
      NumericVector DispatchCost(num_units);
      NumericVector Generation(num_units);

      for (int k = 0; k < num_units; ++k) {

        if (IncludeUnit(k) == 1) {

          FuelCostPerMWh(k) = (FuelPrice(month_i,k) + FuelPriceAdder(month_i,k)) * MaxGenHeatRate(month_i,k);

          // calculate dispatch cost and generation
          DispatchCost(k)         = FuelCostPerMWh(k) + VOM(month_i,k) + PowerPriceFactor(month_i,k) * PowerPrice(i);
          Generation(k)           = MaxGenerationHour(i,k);
          DispatchCostSave(i,k,n) = DispatchCost(k);

          // check for must run constraint
          if (MustRunMatrix(i,k) == 1) {
            
            DispatchCost(k)         = -999;
            // Generation(k)           = MinGenerationHour(i,k);
            Generation(k)           = MaxGenerationHour(i,k);
            DispatchCostSave(i,k,n) = DispatchCost(k);
            
          }
          
        } else {

          Generation(k)           = 0;
          DispatchCost(k)         = 9999;
          DispatchCostSave(i,k,n) = 9999;

        }

      }

      // sort generation by dispatch cost
      IntegerVector idx = seq_along(Generation) - 1;
      std::sort(idx.begin(), idx.end(), [&](int i, int j){return DispatchCost[i] < DispatchCost[j];});
      NumericVector SortedGeneration   = Generation[idx];
      NumericVector SortedDispatchCost = DispatchCost.sort();

      // select marginal unit
      double capacity = 0;
      int below_load = 1;
      for (int k = 0; k < num_units; ++k) {

        if (IncludeUnit(idx(k)) == 1) {

          // commit unit
          StackCutCommit(i, idx(k), n)    = below_load;
          StackCumulativeMW(i, idx(k), n) = capacity;

          // update cumulative capacity
          capacity = capacity + SortedGeneration(k);

          // check to see if load is met
          if (capacity > LoadMatrix(i,n) && below_load == 1) {
            MarginalUnitCost(i,n) = DispatchCost(k);
            below_load            = 0;
          }
        }
      }
      
      // save total generation capacity
      TotalCapacity(i,n) = capacity;
      
      // check for unmet load
      if (MarginalUnitCost(i,n) == 0) {
        MarginalUnitCost(i,n) = 9999;
      }
    }
  }
  
  // create output list to return
  List Out;
  
  Out["DispatchCost"]      = DispatchCostSave;
  Out["MarginalUnitCost"]  = MarginalUnitCost;
  Out["StackCutCommit"]    = StackCutCommit;
  Out["StackCumulativeMW"] = StackCumulativeMW;
  Out["TotalCapacity"]     = TotalCapacity;
  
  return(Out);
  
}
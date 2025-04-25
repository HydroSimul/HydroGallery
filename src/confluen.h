// Defines a header file containing function for CONFLUEN/
#ifndef CONFLUEN
#define CONFLUEN

#include "00utilis.h"

NumericVector riverout_LinearResorvoir(
    NumericVector RIVER_water_m3,
    NumericVector RIVER_inflow_m3,
    NumericVector RIVER_velocity_km,
    NumericVector RIVER_length_km
);
NumericVector riverlakout_LinearResorvoir(
    NumericVector Riverlak_water_m3,
    NumericVector Riverlak_inflow_m3,
    NumericVector Riverlak_capacity_m3,
    NumericVector param_Riverlak_lin_storeFactor
);
NumericVector reservoireleas_Hanasaki(
    NumericVector Reservoi_water_m3,
    NumericVector Reservoi_inflow_m3,
    NumericVector Reservoi_demand_m3,
    NumericVector Reservoi_capacity_m3,
    NumericVector Reservoi_meanInflow_m3,
    NumericVector Reservoi_meanDemand_m3,
    NumericVector Reservoi_releaseCoefficient_1,
    LogicalVector Reservoi_isIrrigate_01
);
#endif
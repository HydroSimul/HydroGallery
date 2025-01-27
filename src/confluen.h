// Defines a header file containing function for CONFLUEN/
#ifndef CONFLUEN
#define CONFLUEN

#include "00utilis.h"

NumericVector river_LinearResorvoir(
    NumericVector RIVER_water_m3,
    NumericVector RIVER_inflow_m3,
    NumericVector RIVER_velocity_km,
    NumericVector RIVER_length_km
);
NumericVector riverlak_LinearResorvoir(
    NumericVector Riverlak_water_m3,
    NumericVector Riverlak_inflow_m3,
    NumericVector Riverlak_capacity_m3,
    NumericVector param_Riverlak_lin_storeFactor
);
#endif
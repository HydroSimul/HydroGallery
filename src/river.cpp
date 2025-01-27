#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]

//' **river outflow**
//' @name river
//' @inheritParams all_vari
//' @description
//'
//' The concept of river estimates the waterbody outflow for waternet concentation
//' @return outflow (m3)
//' @export
// [[Rcpp::export]]
NumericVector river_LinearResorvoir(
   NumericVector RIVER_water_m3,
   NumericVector RIVER_inflow_m3,
   NumericVector RIVER_velocity_km,
   NumericVector RIVER_length_km
)
{
 
 NumericVector RIVER_paramK_TS = RIVER_length_km / RIVER_velocity_km; // pmax(RIVER_length_km / RIVER_velocity_km, 1.);
 
 // return RIVER_water_m3 * (1 / (RIVER_paramK_TS + 0.5)) + RIVER_inflow_m3 * (0.5 / (RIVER_paramK_TS + 0.5));
 return RIVER_water_m3 * (1 - exp(-1. / RIVER_paramK_TS)) + RIVER_inflow_m3 * (1 - RIVER_paramK_TS * (1 - exp(-1. / RIVER_paramK_TS)));
}

//' @rdname river
//' @param param_Riverlak_lin_storeFactor <uknow> parameter for [riverlak_LinearResorvoir()],
//' @export
// [[Rcpp::export]]
NumericVector riverlak_LinearResorvoir(
   NumericVector Riverlak_water_m3,
   NumericVector Riverlak_inflow_m3,
   NumericVector Riverlak_capacity_m3,
   NumericVector param_Riverlak_lin_storeFactor
)
{
 
 
 
 NumericVector Riverlak_outflow_m3 = Riverlak_water_m3 * (1 - exp(-1. / param_Riverlak_lin_storeFactor)) + Riverlak_inflow_m3 * (1 - param_Riverlak_lin_storeFactor * (1 - exp(-1. / param_Riverlak_lin_storeFactor)));
 
 
 NumericVector Riverlak_water_New = pmin(Riverlak_water_m3 + Riverlak_inflow_m3 - Riverlak_outflow_m3, Riverlak_capacity_m3);
 Riverlak_water_New = pmax(Riverlak_water_New, 0);
 Riverlak_outflow_m3 = Riverlak_water_m3 + Riverlak_inflow_m3 - Riverlak_water_New;
 
 return (Riverlak_outflow_m3);
}




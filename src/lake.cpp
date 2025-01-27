#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]

//' **lake outflow**
//' @name lakeout
//' @description
//'
//' The concept of lake estimates the waterbody outflow for waternet concentation
//' @inheritParams all_vari
//' @param param_Lakeout_acp_storeFactor <uknow> parameter for [lake_AcceptPow()],
//' @param param_Lakeout_acp_gamma <uknow> parameter for [lake_AcceptPow()],
//' @return outflow (m3)
//' @export
// [[Rcpp::export]]
NumericVector lakeout_AcceptPow(
   NumericVector Lake_water_m3,
   NumericVector Lake_capacity_m3,
   NumericVector param_Lakeout_acp_storeFactor,
   NumericVector param_Lakeout_acp_gamma
)
{
 
 
 NumericVector Lake_outflow_m3 = (1 / param_Lakeout_acp_storeFactor) * Lake_water_m3 * vecpow(Lake_water_m3 / Lake_capacity_m3, param_Lakeout_acp_gamma);
 
 Lake_outflow_m3 = pmin(Lake_outflow_m3, Lake_water_m3);
 return (Lake_outflow_m3);
}







#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]

//' **potential evapotranspiration**
//' @name evatransPotential
//' @description 
//' 
//' The concept 
//'  
//' @references
//' \insertAllCited{}
//' @inheritParams all_vari
//' @export
// [[Rcpp::export]]
NumericMatrix landLeafAreaIndex_WaterGAP3(NumericMatrix ATMOS_temperature_Cel,
                                         NumericMatrix ATMOS_precipitation_mm,
                                         NumericVector CELL_latitude_deg,
                                         IntegerVector LAND_growUpDay_d,
                                         NumericVector LAND_leafAreaIndexMin_,
                                         NumericVector LAND_leafAreaIndexMax_,
                                         IntegerVector Time_dayOfYear_d) {
 int n_Days = ATMOS_temperature_Cel.nrow();
 int n_Grids = ATMOS_temperature_Cel.ncol();
 NumericMatrix LAND_leafAreaIndex_(n_Days, n_Grids);
 NumericVector LAND_leafAreaRatio_(n_Days);
 NumericVector range_LAI = LAND_leafAreaIndexMax_ - LAND_leafAreaIndexMin_;
 NumericVector num_Growup(366, 1.0);
 NumericVector num_Droop(366, 0.0);
 
 // Initialize growth and droop vectors
 for (int i = 0; i < 30; i++) {
   num_Growup[i] = (i + 1) / 30.0;
   num_Droop[i] = (30 - i) / 30.0;
 }
 
 for (int g = 0; g < n_Grids; ++g) {
   IntegerVector tag_Day;
   if (CELL_latitude_deg[g] >= 0) {
     tag_Day = find_locations(Time_dayOfYear_d, 1);
   } else {
     tag_Day = find_locations(Time_dayOfYear_d, 183);
   }
   
   int n_Year = tag_Day.size();
   tag_Day.push_back(n_Days + 1); // Add a placeholder for end of data
   
   for (int y = 0; y < n_Year; ++y) {
     int start_day = tag_Day[y];
     int next_year_start = tag_Day[y + 1];
     
     double cumsum_Perc_Temp = 0.0;
     for (int d = start_day; d < next_year_start && d < n_Days; ++d) {
       cumsum_Perc_Temp += ATMOS_precipitation_mm(d, g);
       
       // Temperature check over LAND_growUpDay_d[g] days
       int start_idx = std::max(0, d - LAND_growUpDay_d[g]);
       double cum_Temp_Temp = min(NumericVector(ATMOS_temperature_Cel(_, g).begin() + start_idx,
                                                ATMOS_temperature_Cel(_, g).begin() + d + 1));
       
       if (cumsum_Perc_Temp > 40 && cum_Temp_Temp > 8) {
         for (int i = d; i < next_year_start && i < n_Days; ++i) {
           LAND_leafAreaRatio_(i) = num_Growup[std::min(i - d, 365)];
         }
         break;
       }
     }
     
     for (int d = start_day + 182; d < next_year_start && d < n_Days; ++d) {
       // Temperature check for drooping
       int start_idx = std::max(0, d - LAND_growUpDay_d[g]);
       double cum_Temp_Temp = max(NumericVector(ATMOS_temperature_Cel(_, g).begin() + start_idx,
                                                ATMOS_temperature_Cel(_, g).begin() + d + 1));
       
       if (cum_Temp_Temp < 8) {
         for (int i = d; i < next_year_start && i < n_Days; ++i) {
           LAND_leafAreaRatio_(i) = num_Droop[std::min(i - d, 365)];
         }
         break;
       }
     }
   }
   
   
   LAND_leafAreaIndex_(_,g) = LAND_leafAreaIndexMin_(g) + LAND_leafAreaRatio_ * range_LAI(g);
   
 }
 
 return LAND_leafAreaIndex_;
}

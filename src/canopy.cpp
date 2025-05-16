#include "utils.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

//' **potential evapotranspiration**
//' @name evatransPotential
//' @description 
//' The concept 
//'  
//' @references
//' \insertAllCited{}
//' @inheritParams all_vari
//' @export
// [[Rcpp::export]]
arma::mat landLeafAreaIndex_WaterGAP3(const arma::mat& ATMOS_temperature_Cel,
                                      const arma::mat& ATMOS_precipitation_mm,
                                      const arma::vec& CELL_latitude_deg,
                                      const arma::uvec& LAND_growUpDay_d,
                                      const arma::vec& LAND_leafAreaIndexMin_,
                                      const arma::vec& LAND_leafAreaIndexMax_,
                                      const arma::uvec& Time_dayOfYear_d) {
  
  int n_Days = ATMOS_temperature_Cel.n_cols; // (grid, day)
  int n_Grids = ATMOS_temperature_Cel.n_rows; // number of grids
  
  arma::mat LAND_leafAreaIndex_(n_Grids, n_Days, arma::fill::zeros); // (grid, day)
  arma::vec LAND_leafAreaRatio_(n_Days, arma::fill::zeros);
  arma::vec range_LAI = LAND_leafAreaIndexMax_ - LAND_leafAreaIndexMin_;
  
  arma::vec num_Growup(366, arma::fill::ones);
  arma::vec num_Droop(366, arma::fill::zeros);
  
  // Initialize growth and droop vectors
  for (int i = 0; i < 30; ++i) {
    num_Growup(i) = (i + 1.0) / 30.0;
    num_Droop(i) = (30.0 - i) / 30.0;
  }
  
  for (int g = 0; g < n_Grids; ++g) {
    arma::uvec tag_Day;
    
    if (CELL_latitude_deg(g) >= 0) {
      tag_Day = arma::find(Time_dayOfYear_d == 1); // Northern Hemisphere
    } else {
      tag_Day = arma::find(Time_dayOfYear_d == 183); // Southern Hemisphere
    }
    
    int n_Year = tag_Day.n_elem;
    tag_Day.resize(n_Year + 1); // Add a placeholder for end of data
    tag_Day(n_Year) = n_Days;   // Set last element to total days
    
    for (int y = 0; y < n_Year; ++y) {
      int start_day = tag_Day(y);
      int next_year_start = tag_Day(y + 1);
      
      double cumsum_Perc_Temp = 0.0;
      
      for (int d = start_day; d < next_year_start && d < n_Days; ++d) {
        cumsum_Perc_Temp += ATMOS_precipitation_mm(g, d);
        
        // Check temperature over LAND_growUpDay_d(g) days
        int start_idx = std::max(0, d - static_cast<int>(LAND_growUpDay_d(g)));
        double cum_Temp_Temp = arma::min(ATMOS_temperature_Cel.row(g).subvec(start_idx, d));
        
        if (cumsum_Perc_Temp > 40 && cum_Temp_Temp > 8) {
          for (int i = d; i < next_year_start && i < n_Days; ++i) {
            LAND_leafAreaRatio_(i) = num_Growup(std::min(i - d, 365));
          }
          break;
        }
      }
      
      for (int d = start_day + 182; d < next_year_start && d < n_Days; ++d) {
        // Check temperature for drooping
        int start_idx = std::max(0, d - static_cast<int>(LAND_growUpDay_d(g)));
        double cum_Temp_Temp = arma::max(ATMOS_temperature_Cel.row(g).subvec(start_idx, d));
        
        if (cum_Temp_Temp < 8) {
          for (int i = d; i < next_year_start && i < n_Days; ++i) {
            LAND_leafAreaRatio_(i) = num_Droop(std::min(i - d, 365));
          }
          break;
        }
      }
    }
    
    // Assign final leaf area index
    LAND_leafAreaIndex_.row(g) = (LAND_leafAreaIndexMin_(g) +  range_LAI(g) * LAND_leafAreaRatio_).t();
    
  }
  
  return LAND_leafAreaIndex_;
}
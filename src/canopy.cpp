#include <RcppArmadillo.h>
#include <algorithm>  // for std::clamp
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

//' **potential evapotranspiration**
//' @name evatransPotential
//' @description Leaf Area Index model from WaterGAP3
//' @inheritParams all_vari
//' @export
// [[Rcpp::export]]
arma::mat landLeafAreaIndex_WaterGAP3(const arma::mat& ATMOS_temperature_Cel,
                                      const arma::mat& ATMOS_precipitation_mm,
                                      const arma::vec& CELL_latitude_deg,
                                      const arma::ivec& LAND_growUpDay_d,
                                      const arma::vec& LAND_leafAreaIndexMin_,
                                      const arma::vec& LAND_leafAreaIndexMax_,
                                      const arma::ivec& Time_dayOfYear_d) {
  
  int n_Days = ATMOS_temperature_Cel.n_rows;
  int n_Grids = ATMOS_temperature_Cel.n_cols;
  
  arma::mat LAND_leafAreaIndex_(n_Days, n_Grids, arma::fill::zeros);
  arma::vec LAND_leafAreaRatio_(n_Days, arma::fill::zeros);
  arma::vec range_LAI = LAND_leafAreaIndexMax_ - LAND_leafAreaIndexMin_;
  
  arma::vec num_Growup(366, arma::fill::ones);
  arma::vec num_Droop(366, arma::fill::zeros);
  
  for (int i = 0; i < 30; ++i) {
    num_Growup(i) = (i + 1.0) / 30.0;
    num_Droop(i) = (30.0 - i) / 30.0;
  }
  
  for (int g = 0; g < n_Grids; ++g) {
    arma::uvec tag_Day = (CELL_latitude_deg(g) >= 0) ? 
      arma::find(Time_dayOfYear_d == 1) : 
      arma::find(Time_dayOfYear_d == 183);
    
    std::vector<unsigned int> tag_vec = arma::conv_to<std::vector<unsigned int>>::from(tag_Day);
    tag_vec.push_back(n_Days);  // not n_Days + 1, avoids overflow
    
    int n_Year = tag_vec.size() - 1;
    
    for (int y = 0; y < n_Year; ++y) {
      int start_day = tag_vec[y];
      int next_year_start = tag_vec[y + 1];
      
      double cumsum_Prec = 0.0;
      for (int d = start_day; d < next_year_start && d < n_Days; ++d) {
        cumsum_Prec += ATMOS_precipitation_mm(d, g);
        
        int start_idx = std::max(0, d - LAND_growUpDay_d(g));
        double min_Temp = ATMOS_temperature_Cel(arma::span(start_idx, d), g).min();
        
        if (cumsum_Prec > 40 && min_Temp > 8) {
          for (int i = d; i < next_year_start && i < n_Days; ++i) {
            LAND_leafAreaRatio_(i) = num_Growup(std::clamp(i - d, 0, 365));
          }
          break;
        }
      }
      
      for (int d = start_day + 182; d < next_year_start && d < n_Days; ++d) {
        int start_idx = std::max(0, d - LAND_growUpDay_d(g));
        double max_Temp = ATMOS_temperature_Cel(arma::span(start_idx, d), g).max();
        
        if (max_Temp < 8) {
          for (int i = d; i < next_year_start && i < n_Days; ++i) {
            LAND_leafAreaRatio_(i) = num_Droop(std::clamp(i - d, 0, 365));
          }
          break;
        }
      }
    }
    
    LAND_leafAreaIndex_.col(g) = LAND_leafAreaIndexMin_(g) + LAND_leafAreaRatio_ * range_LAI(g);
  }
  
  return LAND_leafAreaIndex_;
}

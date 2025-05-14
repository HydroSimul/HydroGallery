#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::interfaces(r, cpp)]]

//' **river outflow**
//' @name river
//' @inheritParams all_vari
//' @description
//' The concept of river estimates the waterbody outflow for waternet concentration
//' @return outflow (m3)
//' @export
// [[Rcpp::export]]
arma::vec riverout_LinearResorvoir(
    const arma::vec& RIVER_water_m3,
    const arma::vec& RIVER_inflow_m3,
    const arma::vec& RIVER_velocity_km,
    const arma::vec& RIVER_length_km)
{
    // Calculate time parameter
    arma::vec RIVER_paramK_TS = RIVER_length_km / RIVER_velocity_km;
    
    // Calculate exponential terms
    arma::vec exp_term = arma::exp(-1.0 / RIVER_paramK_TS);
    
    // Compute outflow using linear reservoir equation
    return RIVER_water_m3 % (1.0 - exp_term) + 
           RIVER_inflow_m3 % (1.0 - RIVER_paramK_TS % (1.0 - exp_term));
}

//' @rdname river
//' @param param_Riverlak_lin_storeFactor <unknown> parameter for [riverlak_LinearResorvoir()]
//' @export
// [[Rcpp::export]]
arma::vec riverlakout_LinearResorvoir(
    const arma::vec& Riverlak_water_m3,
    const arma::vec& Riverlak_inflow_m3,
    const arma::vec& Riverlak_capacity_m3,
    const arma::vec& param_Riverlak_lin_storeFactor)
{
    // Calculate exponential terms
    arma::vec exp_term = arma::exp(-1.0 / param_Riverlak_lin_storeFactor);
    
    // Compute initial outflow estimate
    arma::vec Riverlak_outflow_m3 = Riverlak_water_m3 % (1.0 - exp_term) + 
                                   Riverlak_inflow_m3 % (1.0 - param_Riverlak_lin_storeFactor % (1.0 - exp_term));
    
    // Calculate new water volume with constraints
    arma::vec Riverlak_water_New = arma::min(Riverlak_water_m3 + Riverlak_inflow_m3 - Riverlak_outflow_m3, 
                                           Riverlak_capacity_m3);
    Riverlak_water_New.transform([](double val) { return std::max(val, 0.0); });
    
    // Adjust outflow to maintain mass balance
    Riverlak_outflow_m3 = Riverlak_water_m3 + Riverlak_inflow_m3 - Riverlak_water_New;
    
    return Riverlak_outflow_m3;
}
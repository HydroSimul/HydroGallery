#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::interfaces(r, cpp)]]

//' **interception** water from land goes into the soil.
//' @name intercep
//' @inheritParams all_vari
//' @description 
//' \loadmathjax
//' 
//' In hydrological modeling, interception refers to the process by which water from precipitation is temporarily retained on the surfaces of vegetation, such as leaves and branches, before being returned to the atmosphere through evaporation or drip.
//' 
//' Under the concept of the conceptual HM, the interception will simply be calculated with the maximal interception of the land.
//' And the interception water will also not go to the land, but will be evaporated.
//' The maximal Interception of the canopy is maybe difficult to estimate 
//' but the process is really simple and there are also not so many methods to describe it. 
//' 
//' @details
//' # **_Full** : 
//' 
//' consider only the radiation and temperature as the main factors. 
//' \mjsdeqn{F_{itcp} = C_{icpt} - W_{icpt}}
//' where
//'   - \mjseqn{F_{icp}} is `intercept_water_mm`
//'   - \mjseqn{C_{icpt}} is `LAND_intercepCapaciy_mm`
//'   - \mjseqn{W_{icpt}} is `LAND_intercepWater_mm`
//' @return intercept_water_mm (mm/m2) intercepted water in this timestep
//' @export
// [[Rcpp::export]]
arma::vec intercep_Full(
    const arma::vec& ATMOS_precipitation_mm,
    const arma::vec& LAND_interceptWater_mm,
    const arma::vec& LAND_interceptCapacity_mm
)
{
  arma::vec water_diff_mm = LAND_interceptCapacity_mm - LAND_interceptWater_mm;
  return arma::min(water_diff_mm, ATMOS_precipitation_mm);
}

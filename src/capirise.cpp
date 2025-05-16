#include "utils.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]


//' **capilarise**
//' @name capirise
//' @inheritParams all_vari
//' @description
//' 
//' \loadmathjax
//' 
//' In hydrological modeling, capillary rise refers to the process by which water is drawn upward from groundwater (table) through the soil due to the force of capillary action.
//' In conceptual watershed models, the capillary rise term often refers to a process that moves water from lower to higher soil water stores, 
//' which may also implicitly include lateral groundwater flow processes in a sloping domain.
//' 
//' It can be calculated by the water in the ground layer \mjseqn{W_{grnd}}, which can also be treated as part of \mjseqn{W_{grnd}}. 
//' There are not many methods to describe this process. Most HMs ignore this process, 
//' perhaps because it is not significant in most situations, or because the process of percolation can deal with this process at the same time.
//' 
//' So we can give the function from:
//' 
//' \mjsdeqn{F_{capi} = f_{capirise}(D_{grnd}, D_{soil})}
//' 
//' to:
//' 
//' \mjsdeqn{F_{capi} = f_{capirise}(W_{grnd}, W_{soil}, C_{soil}, ...)}
//' \mjsdeqn{F_{capi} \leq W_{grnd}}
//' \mjsdeqn{F_{capi} \leq C_{soil} - W_{soil}}
//' 
//' where
//' - \mjseqn{F_{capi}} is `GROUND_capirise_mm`
//' - \mjseqn{W_{grnd}} is `GROUND_water_mm`
//' - \mjseqn{W_{soil}} is `water_SOIL_mm`
//' - \mjseqn{C_{soil}} is `capacity_SOIL_mm`
//' 
//' The output density distribution from 4 methods:
//'
//' @references
//' \insertAllCited{}
//' @return GROUND_capirise_mm (mm/m2/TS) capillary rise
//' 
//' @details
//' # **_HBV** \insertCite{HBV_Lindstrom_1997}{HydroGallery}: 
//' 
//' \mjsdeqn{F_{capi} = M_{capi} \left( 1 - \frac{W_{soil}}{C_{soil}} \right)}
//' where
//'   - \mjseqn{M_{capi}} is `SOIL_potentialCapirise_mm`
//'   
//' @export
// [[Rcpp::export]]
arma::vec capirise_HBV(
   const arma::vec& GROUND_water_mm, 
   const arma::vec& SOIL_water_mm,
   const arma::vec& SOIL_capacity_mm, 
   const arma::vec& SOIL_potentialCapirise_mm)
{
  arma::vec SOIL_diff_mm = SOIL_capacity_mm - SOIL_water_mm;
  SOIL_diff_mm.elem(arma::find(SOIL_diff_mm < 0)).zeros();

  arma::vec limit_mm = arma::min(SOIL_diff_mm, GROUND_water_mm);

  return arma::min(SOIL_potentialCapirise_mm % (SOIL_diff_mm / SOIL_capacity_mm), limit_mm);
}

//' @rdname capirise
//' @details
//' # **_HBVfix** \insertCite{HBV_Lindstrom_1997}{HydroGallery}: 
//' 
//' \mjsdeqn{F_{capi} = M_{capi} \left( 1 - \frac{W_{soil}}{k_{fc}C_{soil}} \right), \quad W_{soil} < k_{fc}C_{soil}}
//' where
//'   - \mjseqn{k_{fc}} is `SOIL_fieldCapacityPerc_1`
//' @export
// [[Rcpp::export]]
arma::vec capirise_HBVfix(
   const arma::vec& GROUND_water_mm, 
   const arma::vec& SOIL_water_mm,
   const arma::vec& SOIL_capacity_mm, 
   const arma::vec& SOIL_fieldCapacityPerc_1,
   const arma::vec& SOIL_potentialCapirise_mm)
{
  arma::vec SOIL_diff_mm = SOIL_capacity_mm % (1 - SOIL_fieldCapacityPerc_1) - SOIL_water_mm;
  SOIL_diff_mm.elem(arma::find(SOIL_diff_mm < 0)).zeros();

  arma::vec limit_mm = arma::min(SOIL_diff_mm, GROUND_water_mm);

  return arma::min(SOIL_potentialCapirise_mm % (SOIL_diff_mm / SOIL_capacity_mm), limit_mm);
}

//' @rdname capirise
//' @details
//' # **_AcceptRatio**: 
//' 
//' \mjsdeqn{F_{capi} = k \left( W_{soil} - k_{fc}C_{soil} \right), \quad W_{soil} < k_{fc}C_{soil}}
//' where
//'   - \mjseqn{k} is `param_CAPIRISE_acr_k`
//'   - \mjseqn{k_{fc}} is `SOIL_fieldCapacityPerc_1`
//' @param param_CAPIRISE_acr_k <0.01, 1> coefficient parameter [capirise_AcceptRatio()]
//' @export
// [[Rcpp::export]]
arma::vec capirise_AcceptRatio(
   const arma::vec& GROUND_water_mm, 
   const arma::vec& SOIL_water_mm,
   const arma::vec& SOIL_capacity_mm, 
   const arma::vec& SOIL_fieldCapacityPerc_1,
   const arma::vec& param_CAPIRISE_acr_k)
{
  arma::vec SOIL_diff_mm = SOIL_capacity_mm % (1 - SOIL_fieldCapacityPerc_1) - SOIL_water_mm;
  SOIL_diff_mm.elem(arma::find(SOIL_diff_mm < 0)).zeros();

  arma::vec limit_mm = arma::min(SOIL_diff_mm, GROUND_water_mm);

  return arma::min(SOIL_diff_mm % param_CAPIRISE_acr_k, limit_mm);
}

//' @rdname capirise
//' @details
//' # **_AcceptPow**: 
//' 
//' \mjsdeqn{F_{capi} = k \left( W_{soil} - k_{fc}C_{soil} \right)^\gamma, \quad W_{soil} < k_{fc}C_{soil}}
//' where
//'   - \mjseqn{k} is `param_CAPIRISE_acp_k`
//'   - \mjseqn{\gamma} is `param_CAPIRISE_acp_gamma`
//'   
//' @param param_CAPIRISE_acp_k <0.01, 1> coefficient parameter for [capirise_AcceptPow()]
//' @param param_CAPIRISE_acp_gamma <0.01, 1> exponential parameter for [capirise_AcceptPow()]
//' @export
// [[Rcpp::export]]
arma::vec capirise_AcceptPow(
   const arma::vec& GROUND_water_mm, 
   const arma::vec& SOIL_water_mm,
   const arma::vec& SOIL_capacity_mm,
   const arma::vec& SOIL_fieldCapacityPerc_1,
   const arma::vec& param_CAPIRISE_acp_k,
   const arma::vec& param_CAPIRISE_acp_gamma)
{
  arma::vec SOIL_diff_mm = SOIL_capacity_mm % (1 - SOIL_fieldCapacityPerc_1) - SOIL_water_mm;
  SOIL_diff_mm.elem(arma::find(SOIL_diff_mm < 0)).zeros();

  arma::vec capirise_mm = (param_CAPIRISE_acp_k % pow(SOIL_diff_mm / (SOIL_capacity_mm % (1 - SOIL_fieldCapacityPerc_1)), param_CAPIRISE_acp_gamma)) % SOIL_diff_mm;
  capirise_mm.elem(arma::find(capirise_mm < 0)).zeros();

  arma::vec limit_mm = arma::min(SOIL_diff_mm, GROUND_water_mm);
  return arma::min(capirise_mm, limit_mm);
}

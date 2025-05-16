#include "utils.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

//' **actuall evapotranspiration**
//' @name evatransActual
//' @inheritParams all_vari
//' @description
//' \loadmathjax
//' Actual ET, or actual evapotranspiration, is a measure of the amount of water that is lost from the land surface through evaporation and transpiration by plants.
//' 
//' Under the concept of the conceptual HM, the actual ET is always calculated by the potential ET \mjseqn{E_p}, which evaluates the meteorological and landuse (vegetation) situations. 
//' The second point to consider is the water availability of the land or soil.
//' 
//' So we can give the function from:
//' 
//' \mjsdeqn{E_a = f_{evatransActual}(D_{atms}, D_{lssg})}
//' 
//' 
//' to:
//' 
//' \mjsdeqn{E_a = f_{evatransActual}(E_p, W_{lssg}, ...) = k^* E_p}
//' 
//' where
//' - \mjseqn{E_a} is `LAND_evatrans_mm` or `soil_evatrans_mm`
//' - \mjseqn{E_p} is `ATMOS_potentialEvatrans_mm`
//' - \mjseqn{k^*} is estimated ratio.
//' 
//' Then the different `evatransActual` methods will estimate the ratio \mjseqn{k^*}.
//' 
//' The output density distribution from 7 methods:
//'
//' @references
//' \insertAllCited{}
//' @return actually ET in (mm/m2/TS)
//' - evaporation in interception (landLy), `LAND_evatrans_mm`
//' - transpiration in root
//' - evaporation in soil (soilLy), `soil_evatrans_mm`
//' @details
//' # **_SupplyRatio**: 
//'
//' 
//' The water content (the ratio to the maximal capacity) 
//' is considered as th main factors for the ratio \mjseqn{k^*}.
//' \mjsdeqn{k^* = k  \frac{W}{C}}
//' where
//'   - \mjseqn{W} is water volume in (mm/m2/TS), `water_mm`, `LAND_interceptWater_mm`, `soil_water_mm`
//'   - \mjseqn{C} is water capacity in (mm/m2), `capacity_mm`, `LAND_interceptCapacity_mm`, `soil_capacity_mm`
//'   - \mjseqn{k} is `param_EVATRANS_sur_k`
//' @param param_EVATRANS_sur_k <0.1, 1> parameter for [evatransActual_SupplyRatio()], ratio of potential ET, that is estimated as actually ET  
//' @export
// [[Rcpp::export]]
arma::vec evatransActual_SupplyRatio(
    const arma::vec& ATMOS_potentialEvatrans_mm,
    const arma::vec& water_mm,
    const arma::vec& capacity_mm,
    const arma::vec& param_EVATRANS_sur_k
)
{
  arma::vec AET = ATMOS_potentialEvatrans_mm % ((water_mm / capacity_mm) % param_EVATRANS_sur_k);
  return arma::min(AET, water_mm);
}

//' @rdname evatransActual
//' @details
//' # **_SupplyPow**: 
//'
//' 
//' The water content (the ratio to the maximal capacity) 
//' is considered as th main factors for the ratio \mjseqn{k^*}.
//' \mjsdeqn{k^* = k  \left(\frac{W}{C}\right)^\gamma}
//' where
//'   - \mjseqn{k} is `param_EVATRANS_sup_k`
//'   - \mjseqn{\gamma} is `param_EVATRANS_sup_gamma`
//' @param param_EVATRANS_sup_k <0.1, 1> parameter for [evatransActual_SupplyPow()], ratio of this method
//' @param param_EVATRANS_sup_gamma <1, 5> parameter for [evatransActual_SupplyPow()], exponent of this method
//' @export
// [[Rcpp::export]]
arma::vec evatransActual_SupplyPow(
    const arma::vec& ATMOS_potentialEvatrans_mm,
    const arma::vec& water_mm,
    const arma::vec& capacity_mm,
    const arma::vec& param_EVATRANS_sup_k,
    const arma::vec& param_EVATRANS_sup_gamma
)
{
  arma::vec AET = ATMOS_potentialEvatrans_mm % (param_EVATRANS_sup_k % arma::pow(water_mm / capacity_mm, param_EVATRANS_sup_gamma));
  return arma::min(AET, water_mm);
}

//' @rdname evatransActual
//' @details
//' # **_VIC** \insertCite{VIC_Wood_1992}{HydroGallery}: 
//'
//' 
//' This method is similar with [evatransActual_SupplyPow()], estimate the water content in the storage.
//' \mjsdeqn{k^* = 1-\left(1-\frac{W}{C}\right)^{\gamma}}
//' where
//'   - \mjseqn{\gamma} is `param_EVATRANS_vic_gamma`
//' @param param_EVATRANS_vic_gamma <0.2, 5> parameter for [evatransActual_VIC()]
//' @export
// [[Rcpp::export]]
arma::vec evatransActual_VIC(
    const arma::vec& ATMOS_potentialEvatrans_mm,
    const arma::vec& water_mm,
    const arma::vec& capacity_mm,
    const arma::vec& param_EVATRANS_vic_gamma
)
{
  arma::vec AET = ATMOS_potentialEvatrans_mm % (1 - arma::pow(1 - water_mm / capacity_mm, param_EVATRANS_vic_gamma));
  return arma::min(AET, water_mm);
}

//' @rdname evatransActual
//' @details
//' # **_GR4J** \insertCite{GR4J_Perrin_2003}{HydroGallery}: 
//'
//' 
//' It is a little different than other method, it estimate not the ratio \mjseqn{f},
//' rather dieectly a equation with potential ET and water content.
//' And it need **no parameter**.
//' \mjsdeqn{E_a = \frac{W\left(2-\frac{W}{C}\right)\tanh \left(\frac{E_p}{C}\right)}{1 + \left(1-\frac{W}{C}\right)\tanh \left(\frac{E_p}{C}\right)}}
//' @export
// [[Rcpp::export]]
arma::vec evatransActual_GR4J(
    const arma::vec& ATMOS_potentialEvatrans_mm,
    const arma::vec& water_mm,
    const arma::vec& capacity_mm
)
{
  arma::vec AET = water_mm % (2 - water_mm / capacity_mm) % arma::tanh(ATMOS_potentialEvatrans_mm / capacity_mm) /
                  (1 + (1 - water_mm / capacity_mm) % arma::tanh(ATMOS_potentialEvatrans_mm / capacity_mm));
  return arma::min(AET, water_mm);
}

//' @rdname evatransActual
//' @details
//' # **_UBC** \insertCite{UBC_Quick_1977}{HydroGallery}: 
//'
//' 
//' It estimates the water content in the storage. 
//' (This is a little different than original, the parameter `P0AGEN` is replaced by \mjseqn{\frac{C}{\gamma}}.)
//' \mjsdeqn{k^* = 10^{\gamma \frac{W-C}{C}}}
//' where
//'   - \mjseqn{\gamma} is `param_EVATRANS_ubc_gamma`
//' @param param_EVATRANS_ubc_gamma <0.5, 2> parameter for [evatransActual_UBC()]
//' @export
// [[Rcpp::export]]
arma::vec evatransActual_UBC(
    const arma::vec& ATMOS_potentialEvatrans_mm,
    const arma::vec& water_mm,
    const arma::vec& capacity_mm,
    const arma::vec& param_EVATRANS_ubc_gamma
) {
  arma::vec diff_mm = capacity_mm - water_mm;
  arma::vec k_ = exp(-diff_mm / (param_EVATRANS_ubc_gamma(0) * capacity_mm) * std::log(10.0));
  arma::vec AET = ATMOS_potentialEvatrans_mm % k_;
  return arma::min(AET, water_mm);
}



//' @rdname evatransActual
//' 
//' @details
//' # **_LiangLand** \insertCite{VIC2_Liang_1994}{HydroGallery}: 
//'
//' 
//' It is also a similar method like [evatransActual_SupplyPow()], 
//' but it will estimate the supply ability agian, whwn the water is still not enough.
//' \mjsdeqn{E_a^* = \left(\frac{W}{C}\right)^\gamma E_p}
//' \mjsdeqn{E_a = \min \left(1, \frac{W}{E_a^*}\right) E_a^*}
//' where
//'   - \mjseqn{E_l^*} is the first estimated actuall ET
//'   - \mjseqn{E_l} is actuall ET from land, `LAND_evatrans_mm`
//'   - \mjseqn{\gamma} is `param_EVATRANS_lia_gamma`
//' @param param_EVATRANS_lia_gamma <0.4, 1> parameter for [evatransActual_LiangLand()]
//' @export
// [[Rcpp::export]]
arma::vec evatransActual_LiangLand(
    const arma::vec& ATMOS_potentialEvatrans_mm,
    const arma::vec& water_mm,
    const arma::vec& capacity_mm,
    const arma::vec& param_EVATRANS_lia_gamma
)
{
  arma::vec AET = ATMOS_potentialEvatrans_mm % (arma::pow(water_mm / capacity_mm, param_EVATRANS_lia_gamma));
  AET = arma::min(water_mm / AET, 1) % AET;
  return arma::min(AET, water_mm);
}

//' @rdname evatransActual
//' @details
//' # **_LiangSoil** \insertCite{VIC2_Liang_1994}{HydroGallery}: 
//'
//' 
//' It estimates the water content in the storage. 
//' (This is a little different than original, the parameter `P0AGEN` is replaced by \mjseqn{\frac{C}{\gamma}}.)
//' \mjsdeqn{k^* = \int_{0}^{A_{s}} {\rm d} A + \int_{A_{s}}^{1} \frac{i_{0}}{i_{m} [1-(1-A)^{1 / B} ]} {\rm d} A }
//' where
//'   - \mjseqn{B} is `param_EVATRANS_lia_B`
//'   - \mjseqn{A} is fraction of area
//' 
//' ![](liang_evatransSoil.png)
//' @param param_EVATRANS_lia_B <0.01, 3> parameter for [evatransActual_LiangSoil()]
//' @export
// [[Rcpp::export]]
arma::vec evatransActual_LiangSoil(
    const arma::vec& ATMOS_potentialEvatrans_mm,
    const arma::vec& water_mm,
    const arma::vec& capacity_mm,
    const arma::vec& param_EVATRANS_lia_B
)
{
  arma::vec i_m = capacity_mm * (param_EVATRANS_lia_B + 1);
  arma::vec i_0 = i_m % (1 - arma::pow(1 - water_mm / capacity_mm, 1 / (param_EVATRANS_lia_B + 1)));
  arma::vec A_s = 1 - arma::pow((1 - i_0 / i_m), param_EVATRANS_lia_B);
  arma::vec A_s_1 = 1 - A_s;
  
  arma::vec k_ = A_s + i_0 / i_m % A_s_1 % (1 + param_EVATRANS_lia_B / (1 + param_EVATRANS_lia_B) % arma::pow(A_s_1, 1 / param_EVATRANS_lia_B) +
                                             param_EVATRANS_lia_B / (2 + param_EVATRANS_lia_B) % arma::pow(A_s_1, 2 / param_EVATRANS_lia_B) +
                                             param_EVATRANS_lia_B / (3 + param_EVATRANS_lia_B) % arma::pow(A_s_1, 3 / param_EVATRANS_lia_B));
  
  arma::vec AET = ATMOS_potentialEvatrans_mm % k_;
  return arma::min(AET, water_mm);
}

//' @rdname evatransActual
//' @param param_EVATRANS_wat_petmax <10, 20> parameter for [evatransActual_WaterGAP()], 10 for arid area, 20 for humid area
//' @export
// [[Rcpp::export]]
arma::vec evatransActual_WaterGAP3(
    const arma::vec& ATMOS_potentialEvatrans_mm,
    const arma::vec& water_mm,
    const arma::vec& capacity_mm,
    const arma::vec& param_EVATRANS_wat_petmax
)
{
  arma::vec AET = arma::min(water_mm / capacity_mm % param_EVATRANS_wat_petmax, ATMOS_potentialEvatrans_mm);
  return arma::min(AET, water_mm);
}

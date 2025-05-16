#include "utils.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

//' caculate **snowfall**
//' @name atmosSnow
//' @inheritParams all_vari
//' @description
//' \loadmathjax
//' 
//' Under the concept of the conceptual HM, the amount of snowfall is always calculated by the temperature \mjseqn{T} and the precipitation \mjseqn{P} availability. 
//' The proportion of snowfall is always determined by the air temperature.
//' 
//' So we can give the function from:
//' 
//' \mjsdeqn{P_s = f_{atmosSnow}(D_{atms})}
//' 
//' to:
//' 
//' \mjsdeqn{P_s = f_{atmosSnow}(P, T) = k^*P}
//' \mjsdeqn{0 \leq k^* \leq 1}
//' where
//'   - \mjseqn{P} is `ATMOS_precpitation_mm`
//'   - \mjseqn{T} is `ATMOS_teperature_Cel`
//' - \mjseqn{k^*} is estimated portion
//' 
//' Then the different `atmosSnow` methods will estimate the portion \mjseqn{k^*}.
//' 
//' The output density distribution from 2 methods:
//' 
//' @references
//' \insertAllCited{}
//' @return ATMOS_snow_mm (mm/m2/TS) snowfall volume
//' @details
//' # **_ThresholdT**: 
//' 
//' Only a temperature is as the threshold defined, so estimate the portion \mjseqn{k^*} as: 
//' \mjsdeqn{k^{*}=1, \quad T \leq T_s}
//' where
//'   - \mjseqn{T_s} is `param_ATMOS_thr_Ts`
//' 
//' @param param_ATMOS_thr_Ts <-1, 3> (Cel) threshold air temperature that snow, parameter for [atmosSnow_ThresholdT()]
//' @export
// [[Rcpp::export]]
arma::vec atmosSnow_ThresholdT(
   const arma::vec& ATMOS_precipitation_mm, 
   const arma::vec& ATMOS_temperature_Cel, 
   const arma::vec& param_ATMOS_thr_Ts)
{
 arma::vec result = ATMOS_precipitation_mm;
 result.elem(arma::find(ATMOS_temperature_Cel > param_ATMOS_thr_Ts)).zeros();
 return result;
}

//' @rdname atmosSnow
//' @details
//' # **_UBC** \insertCite{UBC_Quick_1977}{HydroGallery}: 
//' 
//' estimate the portion \mjseqn{k^*}{} as:
//' \mjsdeqn{k^* = 1- \frac{T}{T_0}}
//' \mjsdeqn{k^* \geq 0}
//' where
//'   - \mjseqn{T_0} is `param_ATMOS_ubc_A0FORM`
//' 
//' @param param_ATMOS_ubc_A0FORM <0.01, 3> (Cel) threshold air temperature that snow, it can not equal or small than 0, parameter for [atmosSnow_UBC()]
//' @export
// [[Rcpp::export]]
arma::vec atmosSnow_UBC(
   const arma::vec& ATMOS_precipitation_mm, 
   const arma::vec& ATMOS_temperature_Cel, 
   const arma::vec& param_ATMOS_ubc_A0FORM)
{
 arma::vec ATMOS_snow_mm = (1 - ATMOS_temperature_Cel / param_ATMOS_ubc_A0FORM) % ATMOS_precipitation_mm;
 
 ATMOS_snow_mm.elem(arma::find(ATMOS_temperature_Cel <= 0)) = 
   ATMOS_precipitation_mm.elem(arma::find(ATMOS_temperature_Cel <= 0));
 
 ATMOS_snow_mm.elem(arma::find(ATMOS_temperature_Cel > param_ATMOS_ubc_A0FORM)).zeros();
 
 return ATMOS_snow_mm;
}

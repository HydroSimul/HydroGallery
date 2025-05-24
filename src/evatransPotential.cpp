#include "utils.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

//' **potential evapotranspiration**
//' @name evatransPotential
//' @description 
//' The concept of potential evapotranspiration (ET) estimates the ability of water lost from the soil and vegetation in an area due to evaporation and transpiration. 
//' It assumes that there is always enough water in the ET area to meet the demand for evapotranspiration.
//' However, the characteristics of the ET area, such as whether it is covered with vegetation or bare soil, can affect the amount of evapotranspiration that occurs. 
//' In order to accurately estimate potential ET, we need to consider these characteristics. 
//' 
//' But we may not always have access to the necessary information or effective methods to do this.
//' In these cases, we can use a simplified method known as **reference ET**. 
//' This method defines the ET area using certain fixed characteristics, such as those provided by the [evatransPotential_FAO56()] function. 
//' In this situation, we need to provide factors to account for the differences between the actual ET area and the reference ET area.
//' 
//' @references
//' \insertAllCited{}
//' @inheritParams all_vari
//' @details
//' - **_TurcWendling** \insertCite{ET_TurcWendling_1991}{HydroGallery}: consider only the radiation and temperature as the main factors. 
//' \mjsdeqn{E_p = \frac{(100 R_s + 3.875 t_h k)\cdot(T + 22)}{150 (T + 123)}}
//' where
//'   - \mjseqn{E_p} is potential ET, `ATMOS_potentialEvatrans_mm`
//'   - \mjseqn{R_s} is solar radiation, `ATMOS_solarRadiat_MJ`
//'   - \mjseqn{t_h} is time step in hour, `time_step_h`
//'   - \mjseqn{T} is average air temperature, `ATMOS_temperature_Cel`
//'   - \mjseqn{k} is `param_EVATRANS_tur_k`
//' @param param_EVATRANS_tur_k <0.6, 1> parameter for [evatransPotential_TurcWendling()], higher value when closer to the sea
//' @return potential evapotranspiration (mm/m2)
//' @export
// [[Rcpp::export]]
arma::vec evatransPotential_TurcWendling(
    const arma::vec& ATMOS_temperature_Cel, 
    const arma::vec& ATMOS_solarRadiat_MJ, 
    const arma::vec& param_EVATRANS_tur_k 
)
{
  arma::vec PET = (ATMOS_solarRadiat_MJ * 100 + 3.875 * 24 * param_EVATRANS_tur_k) % 
    (ATMOS_temperature_Cel + 22) / 150 / (ATMOS_temperature_Cel + 123);
  PET.transform([](double val) { return std::max(val, 0.0); });
  return PET;
}

//' @rdname evatransPotential
//' @details
//' - **_Linacre** \insertCite{ET_Linacre_1977}{HydroGallery}: consider only the temperature as the main factors. 
//' \mjsdeqn{E_p = \frac{\frac{100(0.75 - \alpha)(T + 0.006 z)}{100 - \phi} + 15(T - T_d)}{80 - T}}
//' \mjsdeqn{T_d = T - 20 (1-H_R)}
//' where
//'   - \mjseqn{\alpha} is albedo, `LAND_albedo_1`
//'   - \mjseqn{z} is elevation, `LAND_elevation_m`
//'   - \mjseqn{T_d} is dewpoint temperature,
//'   - \mjseqn{H_R} is relative humidity, `ATMOS_relativeHumidity_1`
//' @export
// [[Rcpp::export]]
arma::vec evatransPotential_Linacre(
    const arma::vec& ATMOS_temperature_Cel,
    const arma::vec& ATMOS_relativeHumidity_1,
    const arma::vec& LAND_latitude_Degree,
    const arma::vec& LAND_elevation_m,
    const arma::vec& LAND_albedo_1
)
{// Compute the full expression
  arma::vec PET = ((0.75 - LAND_albedo_1) * 100 % (ATMOS_temperature_Cel + 0.006 * LAND_elevation_m) / 
    (100 - LAND_latitude_Degree) + 3 * 100 * (1 - ATMOS_relativeHumidity_1)) / 
    (80 - ATMOS_temperature_Cel);
  
  // Efficiently clamp negatives to 0 (no temporary copies)
  PET.transform([](double val) { return std::max(val, 0.0); });
  return PET;
}

//' @rdname evatransPotential
//' @details
//' - **_FAO56** \insertCite{ET_FAO56_1998}{HydroGallery}: consider not only radiation and temperature but also other variables like wind speed
//' as the main factors. 
//' \mjsdeqn{E_p =\frac{0.408 \Delta\left(R_n - G\right)+\gamma \frac{900}{T+273} {u}_{2}\left({e}_{{s}}-{e}_{{a}}\right)}{\Delta+\gamma\left(1+0.34 {u}_{2}\right)}}
//' where
//'   - \mjseqn{\Delta} is slope vapour pressure curve (kPa °C-1)
//'   - \mjseqn{R_n} is net radiation, `ATMOS_netRadiat_MJ`
//'   - \mjseqn{G} is soil heat flux density
//'   - \mjseqn{u_2} is wind speed at 2 m height, `ATMOS_windSpeed2m_m_s`
//'   - \mjseqn{e_s} is saturation vapour pressure, `ATMOS_saturatVaporPress_hPa`
//'   - \mjseqn{e_a} is actual vapour pressure, `ATMOS_vaporPress_hPa`
//'   - \mjseqn{\gamma} is psychrometric constant
//' @export
// [[Rcpp::export]]
arma::vec evatransPotential_FAO56(
    const arma::vec& ATMOS_temperature_Cel, 
    const arma::vec& ATMOS_vaporPress_hPa, 
    const arma::vec& ATMOS_saturatVaporPress_hPa, 
    const arma::vec& ATMOS_netRadiat_MJ,
    const arma::vec& ATMOS_windSpeed2m_m_s,
    const arma::vec& LAND_elevation_m
)
{
  arma::vec Delta_, e_s, e_a, R_n, gamma_, u_2, ET_o;

  R_n = ATMOS_netRadiat_MJ;
  u_2 = ATMOS_windSpeed2m_m_s;
  e_a = ATMOS_vaporPress_hPa;
  e_s = ATMOS_saturatVaporPress_hPa;

  Delta_ = 4098 * (0.6108 * exp(17.27 * ATMOS_temperature_Cel / (ATMOS_temperature_Cel + 237.3))) / ((ATMOS_temperature_Cel + 237.3) % (ATMOS_temperature_Cel + 237.3));

  gamma_ = 0.665e-3 * 101.3 * pow(((293 - 0.0065 * LAND_elevation_m) / 293), 5.26);

  ET_o = (0.408 * Delta_ % (R_n - 0.) + gamma_ * 90 % u_2 % (e_s - e_a) / (ATMOS_temperature_Cel + 273)) / (Delta_ + gamma_ % (1 + 0.34 * u_2));
  ET_o.transform([](double val) { return std::max(val, 0.0); });
  return ET_o;
}

//' @rdname evatransPotential
//' @param param_EVATRANS_prt_alpha <1, 2> parameter for [evatransPotential_PriestleyTaylor()], higher value when closer to the tropical
//' @export
// [[Rcpp::export]]
arma::vec evatransPotential_PriestleyTaylor(
    const arma::vec& ATMOS_temperature_Cel, 
    const arma::vec& ATMOS_netRadiat_MJ,
    const arma::vec& LAND_elevation_m,
    const arma::vec& param_EVATRANS_prt_alpha
)
{
  arma::vec Delta_, P_, gamma_, ET_o, lat_heat;

  lat_heat = 2.501 - 0.002361 * ATMOS_temperature_Cel;
  lat_heat = arma::clamp(lat_heat, 2.835, arma::datum::inf);

  Delta_ = 4098 * (0.6108 * exp(17.27 * ATMOS_temperature_Cel / (ATMOS_temperature_Cel + 237.3))) / ((ATMOS_temperature_Cel + 237.3) % (ATMOS_temperature_Cel + 237.3));

  P_ = 101.3 * pow(((293 - 0.0065 * LAND_elevation_m) / 293), 5.26);
  gamma_ = 0.665e-3 * P_;

  ET_o = param_EVATRANS_prt_alpha % Delta_ % (ATMOS_netRadiat_MJ - 0.) / (Delta_ + gamma_);
  ET_o /= lat_heat;  // Element-wise division (if lat_heat is scalar or same-sized)
  ET_o.transform([](double val) { return std::max(val, 0.0); });
  return ET_o;
}

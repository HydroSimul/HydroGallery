#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::interfaces(r, cpp)]]

#include "meteo.h"

//' **lake outflow**
//' @name lakeout
//' @description
//' The concept of lake estimates the waterbody outflow for waternet concentration
//' @inheritParams all_vari
//' @param param_Lakeout_sup_storeFactor <unknown> parameter for [lakeout_SupplyPow()],
//' @param param_Lakeout_sup_gamma <unknown> parameter for [lakeout_SupplyPow()],
//' @return outflow (m3)
//' @export
// [[Rcpp::export]]
arma::vec lakeout_SupplyPow(
   const arma::vec& Lake_water_m3,
   const arma::vec& Lake_capacity_m3,
   const arma::vec& param_Lakeout_sup_storeFactor,
   const arma::vec& param_Lakeout_sup_gamma
)
{
  arma::vec Lake_outflow_m3 = (1.0 / param_Lakeout_sup_storeFactor) % Lake_water_m3 % arma::pow(Lake_water_m3 / Lake_capacity_m3, param_Lakeout_sup_gamma);
  
  // Return the minimum of the computed outflow and the current water volume
  return arma::min(Lake_outflow_m3, Lake_water_m3);
}

//' **lake evaporation**
//' @name lakevap
//' @description
//' The concept of lake estimates the waterbody outflow for waternet concentration
//' @inheritParams all_vari
//' @return evaporation of lake area (mm / day)
//' @export
// [[Rcpp::export]]
arma::vec lakeevap_Zhao(const arma::vec& ATMOS_solarRadiat_MJ,
                        const arma::vec& ATMOS_temperature_Cel,
                        const arma::vec& ATMOS_vaporPress_kPa,
                        const arma::vec& ATMOS_windSpeed2m_m_s,
                        const arma::vec& LAND_latitude_Degree,
                        const arma::vec& LAND_elevation_m,
                        arma::vec& Lake_temperature_Cel,
                        arma::vec& Lake_depth_m,
                        const arma::vec& Lake_area_km2,
                        const arma::vec& Lake_fetchLength_m,
                        const arma::vec& Time_dayOfYear)
{
  // Constants
  const double const_waterDensity = 1000.0;
  const double const_waterHeatCapacity = 0.0042;
  const double const_stefanBoltzmann = 4.9e-9;
  const double const_tempAbs = 273.15;
  const double const_albedo = 0.1;
  const double const_waterEmissivity = 0.97;

  // Ensure minimum values for wind speed and vapor pressure
  arma::vec ATMOS_windSpeed2m_m_s_modified = arma::max(ATMOS_windSpeed2m_m_s, 0.01);
  arma::vec ATMOS_vaporPress_kPa_modified = arma::max(ATMOS_vaporPress_kPa, 0.0001);

  // Vapor pressure calculations
  arma::vec num_SaturatVaporPress = meteo_saturatVaporPress_kPa(ATMOS_temperature_Cel);
  arma::vec ATMOS_vaporPress_kPa_mod = arma::min(ATMOS_vaporPress_kPa_modified, num_SaturatVaporPress * 0.99);

  // Atmospheric emissivity
  arma::vec ATMOS_atmosEmissivity_ = meteo_atmosEmissivity_UNKNOW(
    Time_dayOfYear,
    ATMOS_temperature_Cel,
    ATMOS_vaporPress_kPa_mod,
    ATMOS_solarRadiat_MJ,
    LAND_latitude_Degree,
    LAND_elevation_m);

  // Net radiation balance
  arma::vec num_Net_Radiat = (1.0 - const_albedo) * ATMOS_solarRadiat_MJ  +
    (ATMOS_atmosEmissivity_ - const_waterEmissivity) % (const_stefanBoltzmann *
    arma::pow(ATMOS_temperature_Cel + const_tempAbs, 4.0));

  // Wet-bulb temperature
  arma::vec num_WetBulbTemperature = meteo_wetBulbTemperature(ATMOS_vaporPress_kPa_mod,
                                                              ATMOS_temperature_Cel);

  // Latent heat of vaporization and psychrometric constant
  arma::vec num_Lambda_Air = 2.501 - ATMOS_temperature_Cel * 2.361e-3;
  arma::vec num_Gamma_TEMP = 101.3 * arma::pow((const_tempAbs + ATMOS_temperature_Cel -
    0.0065 * LAND_elevation_m) / (const_tempAbs + ATMOS_temperature_Cel), 5.26);
  arma::vec num_Gamma = 0.00163 * num_Gamma_TEMP / num_Lambda_Air;

  // Slope of the saturation vapor pressure curve
  arma::vec num_Delta_Tair = meteo_saturatDelta(ATMOS_temperature_Cel);
  arma::vec num_Delta_TwetBulb = meteo_saturatDelta(num_WetBulbTemperature);

  // Wind function
  arma::vec num_Factor_Wind = (2.33 + 1.65 * ATMOS_windSpeed2m_m_s_modified) % arma::pow(Lake_fetchLength_m, -0.1) % num_Lambda_Air;

  // Equilibrium temperature
  arma::vec num_T_Equilibrium = ((0.46 * ATMOS_atmosEmissivity_ + num_Factor_Wind % (num_Delta_Tair + num_Gamma)) % ATMOS_temperature_Cel +
    (1.0 - const_albedo) * ATMOS_solarRadiat_MJ -
    28.38 * (const_waterEmissivity - ATMOS_atmosEmissivity_) - num_Factor_Wind % (num_SaturatVaporPress - ATMOS_vaporPress_kPa_mod)) /
    (0.46 * const_waterEmissivity + num_Factor_Wind % (num_Delta_Tair + num_Gamma));

  // Lake Depth limit
  arma::vec num_LakeDepthLimit = 4.6 * arma::pow(Lake_area_km2, 0.205);
  Lake_depth_m = arma::min(Lake_depth_m, num_LakeDepthLimit);

  // Time constant and water temperature
  arma::vec num_Lake_Heat_LagTime = (const_waterDensity * const_waterHeatCapacity * Lake_depth_m) /
    (4.0 * const_stefanBoltzmann * arma::pow(num_WetBulbTemperature + const_tempAbs, 3.0) +
      num_Factor_Wind % (num_Delta_TwetBulb + num_Gamma));

  arma::vec Lake_newTemperature_Cel = num_T_Equilibrium + (Lake_temperature_Cel - num_T_Equilibrium) % arma::exp(-1.0 / num_Lake_Heat_LagTime);
  Lake_newTemperature_Cel = arma::max(Lake_newTemperature_Cel, 0.0);

  // Heat storage change
  arma::vec num_HeatChange_Lake = const_waterDensity * const_waterHeatCapacity *
    Lake_depth_m % (Lake_newTemperature_Cel - Lake_temperature_Cel);

  // Update temperature
  Lake_temperature_Cel = Lake_newTemperature_Cel;

  // Latent heat flux and evaporation
  arma::vec num_Latent_Heat = (num_Delta_Tair % (num_Net_Radiat - num_HeatChange_Lake) + num_Gamma % num_Factor_Wind % (num_SaturatVaporPress - ATMOS_vaporPress_kPa_mod)) /
    (num_Delta_Tair + num_Gamma);

  return arma::max(num_Latent_Heat / num_Lambda_Air, 0.0);
}

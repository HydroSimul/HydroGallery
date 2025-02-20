#include "00utilis.h"
#include "meteo.h"
// [[Rcpp::interfaces(r, cpp)]]

//' **lake outflow**
//' @name lakeout
//' @description
//'
//' The concept of lake estimates the waterbody outflow for waternet concentation
//' @inheritParams all_vari
//' @param param_Lakeout_sup_storeFactor <uknow> parameter for [lakeout_SupplyPow()],
//' @param param_Lakeout_sup_gamma <uknow> parameter for [lakeout_SupplyPow()],
//' @return outflow (m3)
//' @export
// [[Rcpp::export]]
NumericVector lakeout_SupplyPow(
   NumericVector Lake_water_m3,
   NumericVector Lake_capacity_m3,
   NumericVector param_Lakeout_sup_storeFactor,
   NumericVector param_Lakeout_sup_gamma
)
{
 
 
 NumericVector Lake_outflow_m3 = (1 / param_Lakeout_sup_storeFactor) * Lake_water_m3 * vecpow(Lake_water_m3 / Lake_capacity_m3, param_Lakeout_sup_gamma);
 
 Lake_outflow_m3 = pmin(Lake_outflow_m3, Lake_water_m3);
 return (Lake_outflow_m3);
}

//' **lake evaporation**
//' @name lakevap
//' @description
//'
//' The concept of lake estimates the waterbody outflow for waternet concentation
//' @inheritParams all_vari
//' @return evaporation o flake area (mm / day)
//' @export
// [[Rcpp::export]]
NumericVector lakeevap_Zhao(NumericVector ATMOS_solarRadiat_MJ,
                            NumericVector ATMOS_temperature_Cel,
                            NumericVector ATMOS_vaporPress_kPa,
                            NumericVector ATMOS_windSpeed2m_m_s,
                            NumericVector LAND_latitude_Degree,
                            NumericVector LAND_elevation_m,
                            NumericVector& Lake_temperature_Cel,
                            NumericVector Lake_depth_m,
                            NumericVector Lake_area_km2,
                            NumericVector Lake_fetchLength_m,
                            NumericVector Time_dayOfYear) {
  
  // Constants
  const double const_waterDensity = 1000.0;
  const double const_waterHeatCapacity = 0.0042;
  // const double const_airHeatCapacity = 1.013;
  const double const_stefanBoltzmann = 4.9e-9;
  const double const_tempAbs = 273.15;
  const double const_albedo = 0.1;
  const double const_waterEmissivity = 0.97;
  
  // Ensure minimum values for wind speed and vapor pressure
  ATMOS_windSpeed2m_m_s = pmax(ATMOS_windSpeed2m_m_s, 0.01);
  ATMOS_vaporPress_kPa = pmax(ATMOS_vaporPress_kPa, 0.0001);
  
  // Vapor pressure calculations
  NumericVector num_SaturatVaporPress = meteo_saturatVaporPress_kPa(ATMOS_temperature_Cel);
  NumericVector ATMOS_vaporPress_kPa_mod = pmin(ATMOS_vaporPress_kPa, num_SaturatVaporPress * 0.99);
  
  // Atmospheric emissivity
  NumericVector ATMOS_atmosEmissivity_ = meteo_atmosEmissivity_UNKNOW(
    Time_dayOfYear,
    ATMOS_temperature_Cel,
    ATMOS_vaporPress_kPa_mod,
    ATMOS_solarRadiat_MJ,
    LAND_latitude_Degree,
    LAND_elevation_m);
  
  // Net radiation balance
  NumericVector num_Net_Radiat = ATMOS_solarRadiat_MJ * (1.0 - const_albedo) +
    (ATMOS_atmosEmissivity_ - const_waterEmissivity) * (const_stefanBoltzmann *
    pow(ATMOS_temperature_Cel + const_tempAbs, 4.0));
  
  // Wet-bulb temperature
  NumericVector num_WetBulbTemperature = meteo_wetBulbTemperature(ATMOS_vaporPress_kPa_mod,
                                                                  ATMOS_temperature_Cel);
  
  // Latent heat of vaporization and psychrometric constant
  NumericVector num_Lambda_Air = 2.501 - ATMOS_temperature_Cel * 2.361e-3;
  NumericVector num_Gamma_TEMP = 101.3 * pow((const_tempAbs + ATMOS_temperature_Cel -
    0.0065 * LAND_elevation_m) / (const_tempAbs + ATMOS_temperature_Cel), 5.26);
  NumericVector num_Gamma = 0.00163 * num_Gamma_TEMP / num_Lambda_Air;
  
  // Slope of the saturation vapor pressure curve
  NumericVector num_Delta_Tair = meteo_saturatDelta(ATMOS_temperature_Cel);
  NumericVector num_Delta_TwetBulb = meteo_saturatDelta(num_WetBulbTemperature);
  
  // Wind function
  NumericVector num_Factor_Wind = (2.33 + 1.65 * ATMOS_windSpeed2m_m_s) * pow(Lake_fetchLength_m, -0.1) * num_Lambda_Air;
  
  // Equilibrium temperature
  NumericVector num_T_Equilibrium = ((0.46 * ATMOS_atmosEmissivity_ + num_Factor_Wind * (num_Delta_Tair + num_Gamma)) * ATMOS_temperature_Cel +
    (1.0 - const_albedo) * ATMOS_solarRadiat_MJ -
    28.38 * (const_waterEmissivity - ATMOS_atmosEmissivity_) - num_Factor_Wind * (num_SaturatVaporPress - ATMOS_vaporPress_kPa_mod)) /
      (0.46 * const_waterEmissivity + num_Factor_Wind * (num_Delta_Tair + num_Gamma));
  
  
  // Lake Depth
  NumericVector num_LakeDepthLimit = 4.6*pow(Lake_area_km2, 0.205);
  Lake_depth_m =  ifelse (Lake_depth_m < num_LakeDepthLimit, Lake_depth_m, num_LakeDepthLimit);
  
  // Time constant and water temperature
  NumericVector num_Lake_Heat_LagTime = (const_waterDensity * const_waterHeatCapacity * Lake_depth_m) /
    (4.0 * const_stefanBoltzmann * pow(num_WetBulbTemperature + const_tempAbs, 3.0) +
      num_Factor_Wind * (num_Delta_TwetBulb + num_Gamma));
  
  NumericVector Lake_newTemperature_Cel = num_T_Equilibrium + (Lake_temperature_Cel - num_T_Equilibrium) * exp(-1. / num_Lake_Heat_LagTime); //30//
  Lake_newTemperature_Cel = pmax(Lake_newTemperature_Cel, 0.0);
  
  // Heat storage change
  NumericVector num_HeatChange_Lake = const_waterDensity * const_waterHeatCapacity *
    Lake_depth_m * (Lake_newTemperature_Cel - Lake_temperature_Cel); //30//
  
  // Lake_temperature_Cel = Lake_newTemperature_Cel;
  std::copy(Lake_newTemperature_Cel.begin(), Lake_newTemperature_Cel.end(), Lake_temperature_Cel.begin());
  // Latent heat flux and evaporation
  NumericVector num_Latent_Heat = (num_Delta_Tair * (num_Net_Radiat - num_HeatChange_Lake) + num_Gamma * num_Factor_Wind * (num_SaturatVaporPress - ATMOS_vaporPress_kPa_mod)) /
    (num_Delta_Tair + num_Gamma);
  
  return pmax(num_Latent_Heat / num_Lambda_Air, 0.0);
  // return ATMOS_atmosEmissivity_;
}



#ifndef METEO
#define METEO

#include <RcppArmadillo.h>

arma::vec meteo_solarRadiatClearSky_FAO56(
    const arma::vec& Time_dayOfYear_,
    const arma::vec& LAND_latitude_Degree,
    const arma::vec& LAND_elevation_m);

arma::vec meteo_saturatVaporPress(const arma::vec& ATMOS_temperature_Cel);

arma::vec meteo_saturatVaporPress_kPa(const arma::vec& ATMOS_temperature_Cel);

arma::vec meteo_vaporPress(const arma::vec& ATMOS_temperature_Cel, const arma::vec& ATMOS_relativeHumidity_1);

arma::vec meteo_nettoRadiat_FAO56(
    const arma::vec& Time_dayOfYear_,
    const arma::vec& ATMOS_temperature_Cel,
    const arma::vec& ATMOS_temperatureMax_Cel,
    const arma::vec& ATMOS_temperatureMin_Cel,
    const arma::vec& ATMOS_relativeHumidity_1,
    const arma::vec& ATMOS_solarRadiat_MJ,
    const arma::vec& LAND_latitude_Degree,
    const arma::vec& LAND_elevation_m);

arma::vec meteo_atmosEmissivity_FAO56(
    const arma::vec& Time_dayOfYear_,
    const arma::vec& ATMOS_temperature_Cel,
    const arma::vec& ATMOS_relativeHumidity_1,
    const arma::vec& ATMOS_solarRadiat_MJ,
    const arma::vec& LAND_latitude_Degree,
    const arma::vec& LAND_elevation_m);

arma::vec meteo_cloudFactor_UNKNOW(
    const arma::vec& ATMOS_solarRadiat_MJ,
    const arma::vec& Time_dayOfYear_,
    const arma::vec& LAND_latitude_Degree,
    const arma::vec& LAND_elevation_m);

arma::vec meteo_atmosEmissivity_UNKNOW(
    const arma::vec& Time_dayOfYear_,
    const arma::vec& ATMOS_temperature_Cel,
    const arma::vec& ATMOS_vaporPress_kPa,
    const arma::vec& ATMOS_solarRadiat_MJ,
    const arma::vec& LAND_latitude_Degree,
    const arma::vec& LAND_elevation_m);

arma::vec meteo_nettoRadiat_FAO56Simplify(
    const arma::vec& Time_dayOfYear_,
    const arma::vec& ATMOS_temperature_Cel,
    const arma::vec& ATMOS_relativeHumidity_1,
    const arma::vec& ATMOS_solarRadiat_MJ,
    const arma::vec& LAND_latitude_Degree,
    const arma::vec& LAND_elevation_m);

arma::vec meteo_windSpeed2m(const arma::vec& ATMOS_windSpeed_m_s, const arma::vec& ATMOS_windMeasureHeight_m);

arma::vec meteo_airDensity(const arma::vec& ATMOS_temperature_Cel, const arma::vec& LAND_elevation_m);

arma::vec meteo_saturatDelta(const arma::vec& ATMOS_temperature_Cel);

arma::vec meteo_wetBulbTemperature(
    const arma::vec& ATMOS_vaporPress_kPa, 
    const arma::vec& ATMOS_temperature_Cel);

#endif

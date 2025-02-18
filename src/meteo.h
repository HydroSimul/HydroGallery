// Defines a header file containing function of METEO/
#ifndef METEO
#define METEO

#include "00utilis.h"

NumericVector meteo_solarRadiatClearSky_FAO56(
    NumericVector Time_dayOfYear_,
    NumericVector LAND_latitude_Degree,
    NumericVector LAND_elevation_m);
NumericVector meteo_saturatVaporPress(NumericVector ATMOS_temperature_Cel);
NumericVector meteo_saturatVaporPress_kPa(NumericVector ATMOS_temperature_Cel);
NumericVector meteo_vaporPress(NumericVector ATMOS_temperature_Cel, NumericVector ATMOS_relativeHumidity_1);
NumericVector meteo_nettoRadiat_FAO56(
    NumericVector Time_dayOfYear_,
    NumericVector ATMOS_temperature_Cel,
    NumericVector ATMOS_temperatureMax_Cel,
    NumericVector ATMOS_temperatureMin_Cel,
    NumericVector ATMOS_relativeHumidity_1,
    NumericVector ATMOS_solarRadiat_MJ,
    NumericVector LAND_latitude_Degree,
    NumericVector LAND_elevation_m);
NumericVector meteo_atmosEmissivity_FAO56(
    NumericVector Time_dayOfYear_,
    NumericVector ATMOS_temperature_Cel,
    NumericVector ATMOS_relativeHumidity_1,
    NumericVector ATMOS_solarRadiat_MJ,
    NumericVector LAND_latitude_Degree,
    NumericVector LAND_elevation_m);
NumericVector meteo_cloudFactor_UNKNOW(NumericVector ATMOS_solarRadiat_MJ,
                                       NumericVector Time_dayOfYear_,
                                       NumericVector LAND_latitude_Degree,
                                       NumericVector LAND_elevation_m);
NumericVector meteo_atmosEmissivity_UNKNOW(
    NumericVector Time_dayOfYear_,
    NumericVector ATMOS_temperature_Cel,
    NumericVector ATMOS_vaporPress_kPa,
    NumericVector ATMOS_solarRadiat_MJ,
    NumericVector LAND_latitude_Degree,
    NumericVector LAND_elevation_m);
NumericVector meteo_nettoRadiat_FAO56Simplify(
    NumericVector Time_dayOfYear_,
    NumericVector ATMOS_temperature_Cel,
    NumericVector ATMOS_relativeHumidity_1,
    NumericVector ATMOS_solarRadiat_MJ,
    NumericVector LAND_latitude_Degree,
    NumericVector LAND_elevation_m);
NumericVector meteo_windSpeed2m(NumericVector ATMOS_windSpeed_m_s, NumericVector ATMOS_windMeasureHeight_m);
NumericVector meteo_airDensity(NumericVector ATMOS_temperature_Cel,
                               NumericVector LAND_elevation_m);
NumericVector meteo_saturatDelta(NumericVector ATMOS_temperature_Cel);
NumericVector meteo_wetBulbTemperature(NumericVector ATMOS_vaporPress_kPa, 
                                       NumericVector ATMOS_temperature_Cel);



#endif
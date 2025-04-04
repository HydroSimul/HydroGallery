#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]

//' **meteological variables**
//' some functions to calculate the meteological variables
//' @name meteo
//' @export
// [[Rcpp::export]]
NumericVector meteo_extraterreSolarRadiat_FAO56(
   NumericVector Time_dayOfYear_,
   NumericVector LAND_latitude_Degree) {
 
 NumericVector phi_ = M_PI / 180 * LAND_latitude_Degree;
 NumericVector d_r = 1 + 0.033 * cos(2 * M_PI / 365 * Time_dayOfYear_);
 NumericVector delta_ = 0.409 * sin(2 * M_PI / 365 * Time_dayOfYear_ - 1.39);
 NumericVector omega_s = acos(-tan(phi_) * tan(delta_));
 NumericVector R_a = 37.58603 * d_r * (omega_s * sin(phi_) * sin(delta_) + cos(phi_) * cos(delta_) * sin(omega_s));
 
 return pmax(R_a, 0);
}

//' @rdname meteo
//' @inheritParams all_vari
//' @return meteological variables
//' @export
// [[Rcpp::export]]
NumericVector meteo_solarRadiatClearSky_FAO56(
   NumericVector Time_dayOfYear_,
   NumericVector LAND_latitude_Degree,
   NumericVector LAND_elevation_m) {
 
 NumericVector R_a = meteo_extraterreSolarRadiat_FAO56(
    Time_dayOfYear_,
    LAND_latitude_Degree); //21
 NumericVector R_so = (0.75 + 2e-5 * LAND_elevation_m) * R_a; // eq37
 
 return R_so;
}


//' @rdname meteo
//' @export
// [[Rcpp::export]]
NumericVector meteo_saturatVaporPress(NumericVector ATMOS_temperature_Cel) {
 return 6.1078 * exp(17.27 * ATMOS_temperature_Cel / (ATMOS_temperature_Cel + 237.3));
}

//' @rdname meteo
//' @export
// [[Rcpp::export]]
NumericVector meteo_saturatVaporPress_kPa(NumericVector ATMOS_temperature_Cel) {
 return .61078 * exp(17.27 * ATMOS_temperature_Cel / (ATMOS_temperature_Cel + 237.3));
}

//' @rdname meteo
//' @export
// [[Rcpp::export]]
NumericVector meteo_vaporPress(NumericVector ATMOS_temperature_Cel, NumericVector ATMOS_relativeHumidity_1) {
 return 6.1078 * exp(17.27 * ATMOS_temperature_Cel / (ATMOS_temperature_Cel + 237.3)) * ATMOS_relativeHumidity_1;
}


//' @rdname meteo
//' @export
// [[Rcpp::export]]
NumericVector meteo_nettoRadiat_FAO56(
    NumericVector Time_dayOfYear_,
    NumericVector ATMOS_temperature_Cel,
    NumericVector ATMOS_temperatureMax_Cel,
    NumericVector ATMOS_temperatureMin_Cel,
    NumericVector ATMOS_relativeHumidity_1,
    NumericVector ATMOS_solarRadiat_MJ,
    NumericVector LAND_latitude_Degree,
    NumericVector LAND_elevation_m) {
  
  const double alpha_ = 0.23;
  const double sigma_ = 4.903e-09;
  NumericVector e_s = 0.3054 * exp(17.27 * ATMOS_temperatureMax_Cel / (ATMOS_temperatureMax_Cel + 237.3)) +
    0.3054 * exp(17.27 * ATMOS_temperatureMin_Cel / (ATMOS_temperatureMin_Cel + 237.3));
  NumericVector e_a = e_s * ATMOS_relativeHumidity_1;
  NumericVector phi_ = M_PI / 180 * LAND_latitude_Degree;
  NumericVector d_r = 1 + 0.033 * cos(2 * M_PI / 365 * Time_dayOfYear_);
  NumericVector delta_ = 0.409 * sin(2 * M_PI / 365 * Time_dayOfYear_ - 1.39);
  NumericVector omega_s = acos(-tan(phi_) * tan(delta_));
  NumericVector R_a = 37.58603 * d_r * (omega_s * sin(phi_) * sin(delta_) + cos(phi_) * cos(delta_) * sin(omega_s));
  NumericVector R_so = (0.75 + 2e-5 * LAND_elevation_m) * R_a;
  NumericVector R_ns = (1 - alpha_) * ATMOS_solarRadiat_MJ;
  NumericVector R_nl = sigma_ * (((ATMOS_temperatureMax_Cel + 273.16) * (ATMOS_temperatureMax_Cel + 273.16) *
    (ATMOS_temperatureMax_Cel + 273.16) * (ATMOS_temperatureMax_Cel + 273.16) +
    (ATMOS_temperatureMin_Cel + 273.16) * (ATMOS_temperatureMin_Cel + 273.16) *
    (ATMOS_temperatureMin_Cel + 273.16) * (ATMOS_temperatureMin_Cel + 273.16)) / 2) *
    (0.34 - 0.14 * sqrt(e_a)) * (1.35 * ATMOS_solarRadiat_MJ / R_so - 0.35);
  return R_ns - R_nl;
}

//' @rdname meteo
//' @export
// [[Rcpp::export]]
NumericVector meteo_atmosEmissivity_FAO56(
   NumericVector Time_dayOfYear_,
   NumericVector ATMOS_temperature_Cel,
   NumericVector ATMOS_relativeHumidity_1,
   NumericVector ATMOS_solarRadiat_MJ,
   NumericVector LAND_latitude_Degree,
   NumericVector LAND_elevation_m) {
 
 NumericVector e_a = meteo_vaporPress(ATMOS_temperature_Cel, ATMOS_relativeHumidity_1);
 NumericVector R_so = meteo_solarRadiatClearSky_FAO56(Time_dayOfYear_, LAND_latitude_Degree, LAND_elevation_m);
 NumericVector epsilon_a = (0.34 - 0.14 * sqrt(e_a)) * (1.35 * ATMOS_solarRadiat_MJ / R_so - 0.35); // eq39
 
 return epsilon_a;
}


//' @rdname meteo
//' @export
// [[Rcpp::export]]
NumericVector meteo_cloudFactor_UNKNOW(NumericVector ATMOS_solarRadiat_MJ,
                                      NumericVector Time_dayOfYear_,
                                      NumericVector LAND_latitude_Degree,
                                      NumericVector LAND_elevation_m) {
 // Convert latitude to radians
 NumericVector lat_r = LAND_latitude_Degree * M_PI / 180.0;
 
 // Calculate solar declination
 NumericVector delta = 0.409 * sin(2.0 * M_PI * Time_dayOfYear_ / 365.0 - 1.39);
 
 // Calculate sunset hour angle with validity check
 NumericVector omega_input = -tan(lat_r) * tan(delta);
 LogicalVector valid_omega = abs(omega_input) < 1.0;
 
 // Calculate clear sky radiation
 NumericVector Kso = meteo_solarRadiatClearSky_FAO56(Time_dayOfYear_, LAND_latitude_Degree, LAND_elevation_m);
 
 // Calculate relative shortwave radiation (bounded between 0 and 1)
 NumericVector Kr = pmax(pmin(ATMOS_solarRadiat_MJ / Kso, 1.0), 0.0);
 
 // Calculate cloud factor
 return ifelse(valid_omega, 1.0 - Kr, 0.0);
}

//' @rdname meteo
//' @export
// [[Rcpp::export]]
NumericVector meteo_atmosEmissivity_UNKNOW(
   NumericVector Time_dayOfYear_,
   NumericVector ATMOS_temperature_Cel,
   NumericVector ATMOS_vaporPress_kPa,
   NumericVector ATMOS_solarRadiat_MJ,
   NumericVector LAND_latitude_Degree,
   NumericVector LAND_elevation_m) {
 
 NumericVector num_CloudFactor = meteo_cloudFactor_UNKNOW(ATMOS_solarRadiat_MJ, Time_dayOfYear_,
                                                   LAND_latitude_Degree, LAND_elevation_m);
 
 NumericVector epsilon_a = 1.08 * (1.0 - exp(-vecpow(ATMOS_vaporPress_kPa * 10.0,
                                                (ATMOS_temperature_Cel + 273.15) / 2016.0))) *
                                                  (1.0 + 0.22 * vec_const_pow(num_CloudFactor, 2.75));
 
 return pmax(pmin(epsilon_a, 1.0), 0.0);
}

//' @rdname meteo
//' @export
// [[Rcpp::export]]
NumericVector meteo_atmosEmissivity_Idso(
   NumericVector ATMOS_temperature_Cel) {
 

 NumericVector epsilon_a = 0.261 * exp(-0.000777 * ATMOS_temperature_Cel * ATMOS_temperature_Cel) -0.02;
 
 return pmax(pmin(epsilon_a, 1.0), 0.0);
}

//' @rdname meteo
//' @export
// [[Rcpp::export]]
NumericVector meteo_nettoRadiat_WaterGAP3(
   NumericVector ATMOS_temperature_Cel,
   NumericVector ATMOS_solarRadiat_MJ,
   NumericVector ATMOS_solarRadiatClearSky_MJ,
   NumericVector LAND_albedo_1) {
  
  const double sigma_ = 4.903e-09;
  
  NumericVector R_ns = (1 - LAND_albedo_1) * ATMOS_solarRadiat_MJ;
  NumericVector epsilon_a = meteo_atmosEmissivity_Idso(ATMOS_temperature_Cel);
  NumericVector factor_Cloud = ATMOS_solarRadiat_MJ / ATMOS_solarRadiatClearSky_MJ;
  factor_Cloud = pmin(factor_Cloud, 1);  // Ensure max value is 1
  NumericVector ATMOS_temperature_T = ATMOS_temperature_Cel + 273.16;
  NumericVector factor_Cloud_Rnl = pmax(1.35 * factor_Cloud - 0.35, 0);
  NumericVector R_nl = sigma_ * ATMOS_temperature_T * ATMOS_temperature_T * ATMOS_temperature_T * ATMOS_temperature_T *
    epsilon_a * factor_Cloud_Rnl; // 
  
  return pmax(R_ns - R_nl, 0);
  
}

//' @rdname meteo
//' @export
// [[Rcpp::export]]
NumericVector meteo_nettoRadiat_FAO56Simplify(
   NumericVector Time_dayOfYear_,
   NumericVector ATMOS_temperature_Cel,
   NumericVector ATMOS_relativeHumidity_1,
   NumericVector ATMOS_solarRadiat_MJ,
   NumericVector LAND_latitude_Degree,
   NumericVector LAND_elevation_m) {
 
 const double alpha_ = 0.23;
 const double sigma_ = 4.903e-09;
 
 NumericVector e_a = meteo_vaporPress(ATMOS_temperature_Cel, ATMOS_relativeHumidity_1);
 NumericVector R_so = meteo_solarRadiatClearSky_FAO56(Time_dayOfYear_, LAND_latitude_Degree, LAND_elevation_m);
 NumericVector R_ns = (1 - alpha_) * ATMOS_solarRadiat_MJ;
 NumericVector factor_Cloud = ifelse(R_so > 0, ATMOS_solarRadiat_MJ / R_so, 1);
 factor_Cloud = pmin(factor_Cloud, 1);  // Ensure max value is 1
 NumericVector factor_Cloud_Rnl = pmax(1.35 * factor_Cloud - 0.35, 0);
 NumericVector R_nl = sigma_ * vec_const_pow(ATMOS_temperature_Cel + 273.16, 4.0) *
   (0.34 - 0.14 * sqrt(e_a)) * factor_Cloud_Rnl; // eq39
 
 return pmax(R_ns - R_nl, 0);
}





//' @rdname meteo
//' @export
// [[Rcpp::export]]
NumericVector meteo_windSpeed2m(NumericVector ATMOS_windSpeed_m_s, NumericVector ATMOS_windMeasureHeight_m) {
  return ATMOS_windSpeed_m_s * 4.87 / (log(67.8 * ATMOS_windMeasureHeight_m - 5.42));
}


//' @rdname meteo
//' @export
// [[Rcpp::export]]
NumericVector meteo_airDensity(NumericVector ATMOS_temperature_Cel,
                               NumericVector LAND_elevation_m) {
  // Constants
  const double const_P0 = 101.325;           // sea level standard atmospheric pressure (kPa)
  const double const_gravAccel = 9.80665;    // earth-surface gravitational acceleration (m/s2)
  const double const_gasConstant = 8.31447;  // ideal gas constant (J/(mol K))
  const double const_molarMassDryAir = 0.0289644;  // molar mass of dry air (kg/mol)
  const double const_tempLapseRate = 0.0065; // temperature lapse rate (K/m)
  const double const_maxDensity = 1.225;     // maximum allowed air density (kg/m3)
  
  // Convert temperature to Kelvin and calculate T0
  NumericVector ta = ATMOS_temperature_Cel + 273.15;
  NumericVector t0 = ta + const_tempLapseRate * LAND_elevation_m;
  
  // Calculate the base for power operation
  NumericVector base = 1.0 - (const_tempLapseRate * LAND_elevation_m) / t0;
  
  // Create vector of constant exponent
  NumericVector exponent(base.size(),
                         (const_gravAccel * const_molarMassDryAir) / (const_gasConstant * const_tempLapseRate));
  
  // Calculate pressure (in Pa) using vectorized power operation
  NumericVector p = const_P0 * vecpow(base, exponent) * 1000.0;
  
  // Calculate air density
  NumericVector airds = p * const_molarMassDryAir / (const_gasConstant * ta);
  
  // Apply maximum density constraint using ifelse
  return ifelse(airds > const_maxDensity, const_maxDensity, airds);
}





//' @rdname meteo
//' @export
// [[Rcpp::export]]
NumericVector meteo_saturatDelta(NumericVector ATMOS_temperature_Cel) {
  //// Delta,
  NumericVector Delta = 4098 * (0.6108 * exp(17.27 * ATMOS_temperature_Cel / (ATMOS_temperature_Cel + 237.3))) / ((ATMOS_temperature_Cel + 237.3) * (ATMOS_temperature_Cel + 237.3)); // Eq.13
  return(Delta);
}


// [[Rcpp::export]]
NumericVector meteo_wetBulbTemperature(NumericVector ATMOS_vaporPress_kPa, 
                                       NumericVector ATMOS_temperature_Cel) {
  NumericVector t_d = (116.9 + 237.3 * log(ATMOS_vaporPress_kPa)) / (16.78 - log(ATMOS_vaporPress_kPa));
  NumericVector twb = (0.00066 * 100.0 * ATMOS_temperature_Cel +
    4098.0 * ATMOS_vaporPress_kPa / pow(t_d + 237.3, 2.0) * t_d) /
      (0.00066 * 100.0 + 4098.0 * ATMOS_vaporPress_kPa / pow(t_d + 237.3, 2.0));
  
  return pmin(twb, ATMOS_temperature_Cel);
}



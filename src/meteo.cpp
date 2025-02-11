#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]

//' **meteological variables**
//' some functions to calculate the meteological variables
//' @name meteo
//' @inheritParams all_vari
//' @return meteological variables
//' @export
// [[Rcpp::export]]
NumericVector meteo_nettoRadiat(
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
NumericVector meteo_nettoRadiatSimplify(
    NumericVector Time_dayOfYear_,
    NumericVector ATMOS_temperature_Cel,
    NumericVector ATMOS_relativeHumidity_1,
    NumericVector ATMOS_solarRadiat_MJ,
    NumericVector LAND_latitude_Degree,
    NumericVector LAND_elevation_m) {
  
  const double alpha_ = 0.23;
  const double sigma_ = 4.903e-09;
  NumericVector e_s = 0.6108 * exp(17.27 * ATMOS_temperature_Cel / (ATMOS_temperature_Cel + 237.3));
  NumericVector e_a = e_s * ATMOS_relativeHumidity_1;
  NumericVector phi_ = M_PI / 180 * LAND_latitude_Degree;
  NumericVector d_r = 1 + 0.033 * cos(2 * M_PI / 365 * Time_dayOfYear_);
  NumericVector delta_ = 0.409 * sin(2 * M_PI / 365 * Time_dayOfYear_ - 1.39);
  NumericVector omega_s = acos(-tan(phi_) * tan(delta_));
  NumericVector R_a = 37.58603 * d_r * (omega_s * sin(phi_) * sin(delta_) + cos(phi_) * cos(delta_) * sin(omega_s));
  NumericVector R_so = (0.75 + 2e-5 * LAND_elevation_m) * R_a;
  NumericVector R_ns = (1 - alpha_) * ATMOS_solarRadiat_MJ;
  NumericVector num_Temp_T = (ATMOS_temperature_Cel + 273.16);
  NumericVector R_nl = sigma_ * (num_Temp_T * num_Temp_T * num_Temp_T * num_Temp_T) *
    (0.34 - 0.14 * sqrt(e_a)) * (1.35 * ATMOS_solarRadiat_MJ / R_so - 0.35);
  return R_ns - R_nl;
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
NumericVector meteo_vaporPress(NumericVector ATMOS_temperature_Cel, NumericVector ATMOS_relativeHumidity_1) {
  return 6.1078 * exp(17.27 * ATMOS_temperature_Cel / (ATMOS_temperature_Cel + 237.3)) * ATMOS_relativeHumidity_1;
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
NumericVector meteo_cloudFactor(NumericVector ATMOS_solarRadiat_MJ,
                                NumericVector Time_dayOfYear,
                                NumericVector LAND_latitude_Degree,
                                NumericVector LAND_elevation_m) {
  // Constants
  const double const_solarConstant = 4.92;    // MJ/m2/h
  const double const_elevCoef = 2e-5;         // elevation coefficient
  
  // Convert latitude to radians
  NumericVector lat_r = LAND_latitude_Degree * M_PI / 180.0;
  
  // Calculate solar declination
  NumericVector delta = 0.409 * sin(2.0 * M_PI * Time_dayOfYear / 366.0 - 1.39);
  
  // Calculate omega input
  NumericVector omega_input = -tan(lat_r) * tan(delta);
  
  
  // Create condition vector for polar night check
  LogicalVector valid_omega = abs(omega_input) < 1.0;
  
  // Calculate sunset hour angle (use acos(1) for invalid cases)
  NumericVector omega = acos(omega_input);
  // Calculate intermediate variables
  NumericVector dr = 1.0 + 0.033 * cos(2.0 * M_PI * Time_dayOfYear / 366.0);
  
  // Calculate extraterrestrial radiation
  NumericVector Ket = 24.0 / M_PI * const_solarConstant * dr *
    (omega * sin(lat_r) * sin(delta) +
    cos(lat_r) * cos(delta) * sin(omega));
  
  // Calculate clear sky radiation
  NumericVector Kso = (0.75 + const_elevCoef * LAND_elevation_m) * Ket;
  
  // Calculate relative shortwave radiation (handle division by zero)
  NumericVector Kr = ATMOS_solarRadiat_MJ / Kso;
  
  // Bound Kr between 0 and 1
  Kr = pmax(pmin(Kr, 1.0), 0.0);
  
  // Calculate cloud factor
  return ifelse(valid_omega, 1.0 - Kr, 0.0);
}


//' @rdname meteo
//' @export
// [[Rcpp::export]]
NumericVector meteo_saturatDelta(NumericVector ATMOS_temperature_Cel) {
  //// Delta,
  NumericVector Delta = 4098 * (0.6108 * exp(17.27 * ATMOS_temperature_Cel / (ATMOS_temperature_Cel + 237.3))) / ((ATMOS_temperature_Cel + 237.3) * (ATMOS_temperature_Cel + 237.3)); // Eq.13
  return(Delta);
}


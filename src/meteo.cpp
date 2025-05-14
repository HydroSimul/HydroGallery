#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::interfaces(r, cpp)]]

//' **meteological variables**
//' some functions to calculate the meteological variables
//' @inheritParams all_vari
//' @name meteo
//' @export
// [[Rcpp::export]]
arma::vec meteo_extraterreSolarRadiat_FAO56(
   const arma::vec& Time_dayOfYear_,
   const arma::vec& LAND_latitude_Degree) {
  
  arma::vec phi_ = M_PI / 180 * LAND_latitude_Degree; //eq22
  arma::vec d_r = 1 + 0.033 * cos(2 * M_PI / 365 * Time_dayOfYear_); // eq23
  arma::vec delta_ = 0.409 * sin(2 * M_PI / 365 * Time_dayOfYear_ - 1.39); // 24
  arma::vec omega_s_TEMP = arma::clamp(-tan(phi_) % tan(delta_), -1, 1);
  arma::vec omega_s = acos(omega_s_TEMP); //25
  arma::vec R_a = 37.58603 * d_r % (omega_s % sin(phi_) % sin(delta_) + cos(phi_) % cos(delta_) % sin(omega_s)); //21
  
  R_a.elem(arma::find(R_a < 0)).zeros();  // Sets negative values to 0
  return R_a;
}

//' @rdname meteo
//' @return meteological variables
//' @export
// [[Rcpp::export]]
arma::vec meteo_solarRadiatClearSky_FAO56(
   const arma::vec& Time_dayOfYear_,
   const arma::vec& LAND_latitude_Degree,
   const arma::vec& LAND_elevation_m) {
 
  arma::vec R_a = meteo_extraterreSolarRadiat_FAO56(Time_dayOfYear_, LAND_latitude_Degree); //21
  arma::vec R_so = (0.75 + 2e-5 * LAND_elevation_m) % R_a; // eq37
 
  return R_so;
}


//' @rdname meteo
//' @export
// [[Rcpp::export]]
arma::vec meteo_saturatVaporPress(const arma::vec& ATMOS_temperature_Cel) {
  return 6.1078 * exp(17.27 * ATMOS_temperature_Cel / (ATMOS_temperature_Cel + 237.3));
}

//' @rdname meteo
//' @export
// [[Rcpp::export]]
arma::vec meteo_saturatVaporPress_kPa(const arma::vec& ATMOS_temperature_Cel) {
  return .61078 * exp(17.27 * ATMOS_temperature_Cel / (ATMOS_temperature_Cel + 237.3));
}

//' @rdname meteo
//' @export
// [[Rcpp::export]]
arma::vec meteo_vaporPress(const arma::vec& ATMOS_temperature_Cel, const arma::vec& ATMOS_relativeHumidity_1) {
  return 6.1078 * exp(17.27 * ATMOS_temperature_Cel / (ATMOS_temperature_Cel + 237.3)) % ATMOS_relativeHumidity_1;
}


//' @rdname meteo
//' @export
// [[Rcpp::export]]
arma::vec meteo_nettoRadiat_FAO56(
    const arma::vec& Time_dayOfYear_,
    const arma::vec& ATMOS_temperature_Cel,
    const arma::vec& ATMOS_temperatureMax_Cel,
    const arma::vec& ATMOS_temperatureMin_Cel,
    const arma::vec& ATMOS_relativeHumidity_1,
    const arma::vec& ATMOS_solarRadiat_MJ,
    const arma::vec& LAND_latitude_Degree,
    const arma::vec& LAND_elevation_m) {
  
  const double alpha_ = 0.23;
  const double sigma_ = 4.903e-09;
  arma::vec e_s = 0.3054 * exp(17.27 * ATMOS_temperatureMax_Cel / (ATMOS_temperatureMax_Cel + 237.3)) +
    0.3054 * exp(17.27 * ATMOS_temperatureMin_Cel / (ATMOS_temperatureMin_Cel + 237.3));
  arma::vec e_a = e_s % ATMOS_relativeHumidity_1;
  arma::vec phi_ = M_PI / 180 * LAND_latitude_Degree;
  arma::vec d_r = 1 + 0.033 * cos(2 * M_PI / 365 * Time_dayOfYear_);
  arma::vec delta_ = 0.409 * sin(2 * M_PI / 365 * Time_dayOfYear_ - 1.39);
  arma::vec omega_s = acos(-tan(phi_) % tan(delta_));
  arma::vec R_a = 37.58603 * d_r % (omega_s % sin(phi_) % sin(delta_) + cos(phi_) % cos(delta_) % sin(omega_s));
  arma::vec R_so = (0.75 + 2e-5 * LAND_elevation_m) % R_a;
  arma::vec R_ns = (1 - alpha_) * ATMOS_solarRadiat_MJ;
  arma::vec R_nl = sigma_ * (((ATMOS_temperatureMax_Cel + 273.16) % (ATMOS_temperatureMax_Cel + 273.16) %
    (ATMOS_temperatureMax_Cel + 273.16) % (ATMOS_temperatureMax_Cel + 273.16) +
    (ATMOS_temperatureMin_Cel + 273.16) % (ATMOS_temperatureMin_Cel + 273.16) %
    (ATMOS_temperatureMin_Cel + 273.16) % (ATMOS_temperatureMin_Cel + 273.16)) / 2) *
    (0.34 - 0.14 * sqrt(e_a)) % (1.35 * ATMOS_solarRadiat_MJ / R_so - 0.35);
  return R_ns - R_nl;
}

//' @rdname meteo
//' @export
// [[Rcpp::export]]
arma::vec meteo_atmosEmissivity_FAO56(
   const arma::vec& Time_dayOfYear_,
   const arma::vec& ATMOS_temperature_Cel,
   const arma::vec& ATMOS_relativeHumidity_1,
   const arma::vec& ATMOS_solarRadiat_MJ,
   const arma::vec& LAND_latitude_Degree,
   const arma::vec& LAND_elevation_m) {
 
  arma::vec e_a = meteo_vaporPress(ATMOS_temperature_Cel, ATMOS_relativeHumidity_1);
  arma::vec R_so = meteo_solarRadiatClearSky_FAO56(Time_dayOfYear_, LAND_latitude_Degree, LAND_elevation_m);
  arma::vec epsilon_a = (0.34 - 0.14 * sqrt(e_a)) % (1.35 * ATMOS_solarRadiat_MJ / R_so - 0.35); // eq39
 
  return epsilon_a;
}



//' @rdname meteo
//' @export
// [[Rcpp::export]]
arma::vec meteo_cloudFactor_UNKNOW(
    const arma::vec& ATMOS_solarRadiat_MJ,
    const arma::vec& Time_dayOfYear_,
    const arma::vec& LAND_latitude_Degree,
    const arma::vec& LAND_elevation_m) {
    // Convert latitude to radians
    arma::vec lat_r = LAND_latitude_Degree * M_PI / 180.0;

    // Calculate solar declination
    arma::vec delta = 0.409 * arma::sin(2.0 * M_PI * Time_dayOfYear_ / 365.0 - 1.39);

    // Calculate sunset hour angle with validity check
    arma::vec omega_input = -arma::tan(lat_r) % arma::tan(delta);
    arma::uvec valid_omega = arma::abs(omega_input) < 1.0;

    // Calculate clear sky radiation
    arma::vec Kso = meteo_solarRadiatClearSky_FAO56(Time_dayOfYear_, LAND_latitude_Degree, LAND_elevation_m);

    // Calculate relative shortwave radiation (bounded between 0 and 1)
    arma::vec Kr = arma::clamp(ATMOS_solarRadiat_MJ / Kso, 0.0, 1.0);

    // Calculate cloud factor
    arma::vec cloudFactor = arma::zeros<arma::vec>(ATMOS_solarRadiat_MJ.n_elem);
    cloudFactor.elem(valid_omega) = 1.0 - Kr.elem(valid_omega);

    return cloudFactor;
}

//' @rdname meteo
//' @export
// [[Rcpp::export]]
arma::vec meteo_atmosEmissivity_UNKNOW(
    const arma::vec& Time_dayOfYear_,
    const arma::vec& ATMOS_temperature_Cel,
    const arma::vec& ATMOS_vaporPress_kPa,
    const arma::vec& ATMOS_solarRadiat_MJ,
    const arma::vec& LAND_latitude_Degree,
    const arma::vec& LAND_elevation_m) {

    arma::vec num_CloudFactor = meteo_cloudFactor_UNKNOW(ATMOS_solarRadiat_MJ, Time_dayOfYear_,
                                                          LAND_latitude_Degree, LAND_elevation_m);

    arma::vec epsilon_a = 1.08 * (1.0 - arma::exp(-arma::pow(ATMOS_vaporPress_kPa * 10.0,
                                                           (ATMOS_temperature_Cel + 273.15) / 2016.0))) *
                          (1.0 + 0.22 * arma::pow(num_CloudFactor, 2.75));

    return arma::clamp(epsilon_a, 0.0, 1.0);
}

//' @rdname meteo
//' @export
// [[Rcpp::export]]
arma::vec meteo_atmosEmissivity_Idso(const arma::vec& ATMOS_temperature_Cel) {
    arma::vec epsilon_a = 0.261 * arma::exp(-0.000777 * arma::pow(ATMOS_temperature_Cel, 2)) - 0.02;
    return arma::clamp(epsilon_a, 0.0, 1.0);
}

//' @rdname meteo
//' @export
// [[Rcpp::export]]
arma::vec meteo_nettoRadiat_WaterGAP3(
   const arma::vec& ATMOS_temperature_Cel,
   const arma::vec& ATMOS_solarRadiat_MJ,
   const arma::vec& ATMOS_solarRadiatClearSky_MJ,
   const arma::vec& LAND_albedo_1) {
 
 const double sigma_ = 4.903e-09;
 
 arma::vec R_ns = (1 - LAND_albedo_1) % ATMOS_solarRadiat_MJ;
 arma::vec epsilon_a = meteo_atmosEmissivity_Idso(ATMOS_temperature_Cel);
 arma::vec factor_Cloud = ATMOS_solarRadiat_MJ / ATMOS_solarRadiatClearSky_MJ;
 factor_Cloud.transform([](double val) { return std::min(val, 1.0); });

 arma::vec ATMOS_temperature_T = ATMOS_temperature_Cel + 273.16;
 arma::vec factor_Cloud_Rnl = 1.35 * factor_Cloud - 0.35;
 factor_Cloud_Rnl.transform([](double val) { return std::max(val, 0.0); });
 
 
 arma::vec R_nl = sigma_ * arma::pow(ATMOS_temperature_T, 4) % epsilon_a % factor_Cloud_Rnl;
 R_nl = arma::clamp(R_nl, 0.0, arma::datum::inf); // R_nl > 0 if condition
 
 arma::vec diff = R_ns - R_nl;
 diff.transform([](double val) { return std::max(val, 0.0); });
 return diff;
}

//' @rdname meteo
//' @export
// [[Rcpp::export]]
arma::vec meteo_nettoRadiat_FAO56Simplify(
   const arma::vec& Time_dayOfYear_,
   const arma::vec& ATMOS_temperature_Cel,
   const arma::vec& ATMOS_relativeHumidity_1,
   const arma::vec& ATMOS_solarRadiat_MJ,
   const arma::vec& LAND_latitude_Degree,
   const arma::vec& LAND_elevation_m) {
 
 const double alpha_ = 0.23;
 const double sigma_ = 4.903e-09;
 
 arma::vec e_a = meteo_vaporPress(ATMOS_temperature_Cel, ATMOS_relativeHumidity_1);
 arma::vec R_so = meteo_solarRadiatClearSky_FAO56(Time_dayOfYear_, LAND_latitude_Degree, LAND_elevation_m);
 arma::vec R_ns = (1 - alpha_) * ATMOS_solarRadiat_MJ;
 arma::vec factor_Cloud = ATMOS_solarRadiat_MJ / R_so;                      // Division
 factor_Cloud.elem(arma::find(R_so <= 0)).fill(1.0);                       // Set factor_Cloud to 1 where R_so <= 0
 factor_Cloud.transform([](double val) { return std::min(val, 1.0); });    // Clamp max to 1
 
 arma::vec factor_Cloud_Rnl = 1.35 * factor_Cloud - 0.35;
 factor_Cloud_Rnl.transform([](double val) { return std::max(val, 0.0); });
 arma::vec R_nl = sigma_ * arma::pow(ATMOS_temperature_Cel + 273.16, 4.0) %
   (0.34 - 0.14 * arma::sqrt(e_a)) % factor_Cloud_Rnl;  // eq39
 
 arma::vec diff = R_ns - R_nl;
 diff.transform([](double val) { return std::max(val, 0.0); });
 return diff;
}

//' @rdname meteo
//' @export
// [[Rcpp::export]]
arma::vec meteo_windSpeed2m(const arma::vec& ATMOS_windSpeed_m_s, const arma::vec& ATMOS_windMeasureHeight_m) {
    return ATMOS_windSpeed_m_s * 4.87 / (arma::log(67.8 * ATMOS_windMeasureHeight_m - 5.42));
}

//' @rdname meteo
//' @export
// [[Rcpp::export]]
arma::vec meteo_airDensity(const arma::vec& ATMOS_temperature_Cel,
                           const arma::vec& LAND_elevation_m) {
    // Constants
    const double const_P0 = 101.325;           // sea level standard atmospheric pressure (kPa)
    const double const_gravAccel = 9.80665;    // earth-surface gravitational acceleration (m/s2)
    const double const_gasConstant = 8.31447;  // ideal gas constant (J/(mol K))
    const double const_molarMassDryAir = 0.0289644;  // molar mass of dry air (kg/mol)
    const double const_tempLapseRate = 0.0065; // temperature lapse rate (K/m)
    const double const_maxDensity = 1.225;     // maximum allowed air density (kg/m3)

    // Convert temperature to Kelvin and calculate T0
    arma::vec ta = ATMOS_temperature_Cel + 273.15;
    arma::vec t0 = ta + const_tempLapseRate * LAND_elevation_m;

    // Calculate the base for power operation
    arma::vec base = 1.0 - (const_tempLapseRate * LAND_elevation_m) / t0;

    // Create vector of constant exponent
    arma::vec exponent = arma::ones<arma::vec>(base.size()) *
                         (const_gravAccel * const_molarMassDryAir) / (const_gasConstant * const_tempLapseRate);

    // Calculate pressure (in Pa) using vectorized power operation
    arma::vec p = const_P0 * arma::pow(base, exponent) * 1000.0;

    // Calculate air density
    arma::vec airds = p * const_molarMassDryAir / (const_gasConstant * ta);

    // Apply maximum density constraint using ifelse
    return arma::clamp(airds, 0.0, const_maxDensity);
}

//' @rdname meteo
//' @export
// [[Rcpp::export]]
arma::vec meteo_saturatDelta(const arma::vec& ATMOS_temperature_Cel) {
    // Delta
    arma::vec Delta = 4098 * (0.6108 * arma::exp(17.27 * ATMOS_temperature_Cel / (ATMOS_temperature_Cel + 237.3))) /
                      arma::pow((ATMOS_temperature_Cel + 237.3), 2.0); // Eq.13
    return Delta;
}

//' @rdname meteo
//' @export
// [[Rcpp::export]]
arma::vec meteo_wetBulbTemperature(
    const arma::vec& ATMOS_vaporPress_kPa, 
    const arma::vec& ATMOS_temperature_Cel) {
    arma::vec t_d = (116.9 + 237.3 * arma::log(ATMOS_vaporPress_kPa)) /
                    (16.78 - arma::log(ATMOS_vaporPress_kPa));
    arma::vec twb = (0.00066 * 100.0 * ATMOS_temperature_Cel +
                     4098.0 * ATMOS_vaporPress_kPa / arma::pow(t_d + 237.3, 2.0) * t_d) /
                    (0.00066 * 100.0 + 4098.0 * ATMOS_vaporPress_kPa / arma::pow(t_d + 237.3, 2.0));

    return arma::min(twb, ATMOS_temperature_Cel);
}

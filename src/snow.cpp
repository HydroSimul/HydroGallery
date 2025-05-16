#include "utils.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

//' **snow**
//' @name snowMelt
//' @inheritParams all_vari
//' @description
//' \loadmathjax
//' 
//' Under the concept of the conceptual HM, the melt of snowpack is always calculated by 
//' the energy availability (the state-variable temperature \mjseqn{T} or flux-variable (nett-) radiation \mjseqn{Rn}) 
//' and the solid water (snow or ice) availability \mjseqn{W_{snow}} of the snowpack. 
//' 
//' Some more complex processes, such as refrozen and residual water, will be ignored. 
//' To simplify the model, the layer snowLy will store only the solid water and will melt it as much as possible when the energy is sufficient.
//' 
//' So we can give the function from:
//' 
//' \mjsdeqn{F_{melt} = f_{snowMelt}(D_{atms}, D_{snow})}
//' 
//' 
//' to:
//' \mjsdeqn{F_{melt}  = f_{snowMelt}(T, ...)}
//' \mjsdeqn{F_{melt} \leq W_{snow} }
//' 
//' where
//'   - \mjseqn{F_{melt}} is `SNOW_melt_mm`
//'   - \mjseqn{W_{snow}} is `SNOW_ice_mm`
//'   - \mjseqn{T} is average temperature
//' 
//' Then the different `snowMelt` methods will estimate the maximal snow melt \mjseqn{M_{max}}.
//' 
//' The output density distribution from 2 methods:
//'
//' @references
//' \insertAllCited{}
//' @return SNOW_melt_mm (mm/m2) melted snow
//' 
//' @details
//' # **_Kustas** \insertCite{SNOW_kustas_1994}{HydroGallery}: 
//' 
//'
//' 
//' \mjsdeqn{F_{melt}  = m_T T + m_E R_n}
//' but due to the temperature is one energy-state-variable, 
//' in order to adjust to subday scale we need to add a new time interval \mjseqn{t_h} from 1 to 24 hour
//' \mjsdeqn{F_{melt}  = m_T T t_h + m_E R_n}
//' where
//'   - \mjseqn{m_T} is `param_SNOWMELT_kus_fT`
//'   - \mjseqn{m_E} is `param_SNOWMELT_kus_fE`
//'   - \mjseqn{R_n} is daily net radiation
//' 
//' @param param_SNOWMELT_kus_fE <0.0005, 0.003> (mm/m2/MJ) snow melt temperature parameter for [snowMelt_Factor()]
//' @param param_SNOWMELT_kus_fT <0.05, 1> (mm/m2/h/Cel) potential melt volum per Cel per hour parameter for [snowMelt_Factor()]
//' @export
// [[Rcpp::export]]
arma::vec snowMelt_Kustas(
    const arma::vec& SNOW_ice_mm,
    const arma::vec& ATMOS_temperature_Cel,
    const arma::vec& ATMOS_netRadiat_MJ,
    const arma::vec& param_SNOWMELT_kus_fE,
    const arma::vec& param_SNOWMELT_kus_fT
)
{
  // Apply the formula for snow melt based on temperature and net radiation
  arma::vec SNOW_melt_mm = arma::clamp(ATMOS_temperature_Cel, 0.0, arma::datum::inf) % param_SNOWMELT_kus_fT * 24 + param_SNOWMELT_kus_fE % ATMOS_netRadiat_MJ;
  
  // Ensure that the snow melt does not exceed the available snow
  return arma::min(SNOW_melt_mm, SNOW_ice_mm);
}

//' @rdname snowMelt
//' @details
//' # **_Factor** \insertCite{phyHydro_dingman_2014}{HydroGallery}: 
//'
//' 
//' \mjsdeqn{F_{melt}  = m_T (T-T_b), \quad T > T_b}
//' where
//'   - \mjseqn{m_T} is `param_SNOWMELT_fac_f`
//'   - \mjseqn{T_b} is `param_SNOWMELT_fac_Tmelt`
//' 
//' @param param_SNOWMELT_fac_Tmelt <0, 3> (Cel) snow melt temperature parameter for [snowMelt_Factor()]
//' @param param_SNOWMELT_fac_f <0.5, 2> (mm/m2/day/Cel) potential melt volum per Cel per hour parameter for [snowMelt_Factor()]
//' @export
// [[Rcpp::export]]
arma::vec snowMelt_Factor(
    const arma::vec& SNOW_ice_mm,
    const arma::vec& ATMOS_temperature_Cel,
    const arma::vec& param_SNOWMELT_fac_f,
    const arma::vec& param_SNOWMELT_fac_Tmelt
)
{
  arma::vec diff_T = ATMOS_temperature_Cel - param_SNOWMELT_fac_Tmelt;
  
  // Apply condition that only temperatures above T_b are considered for melt
  diff_T = arma::clamp(diff_T, 0.0, arma::datum::inf);
  
  // Calculate snow melt based on temperature difference
  arma::vec SNOW_melt_mm = param_SNOWMELT_fac_f % diff_T;
  
  // Ensure that the snow melt does not exceed the available snow
  return arma::min(SNOW_melt_mm, SNOW_ice_mm);
}

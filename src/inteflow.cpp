#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]



//' **interflow**
//' @name inteflow
//' @inheritParams all_vari
//' @description
//' \loadmathjax
//' 
//' In hydrological modeling, interflow refers to the movement of water that is transported horizontally through the soil or aquifer.
//' Like [baseflow], the impact of other RUs (response units) on the route to the river will be ignored.
//' 
//' It can be calculated by the water in the soil layer \mjseqn{W_{soil}},
//' which can also be tread as the part of the \mjseqn{W_{soil}}.
//' 
//' So we can give the function from:
//' 
//' \mjsdeqn{F_{itfl} = f_{inteflow}(D_{grnd}, D_{soil})}
//' 
//' 
//' to:
//' 
//' \mjsdeqn{F_{itfl} = f_{inteflow}(W_{soil}, C_{soil}, ...) = k^* W_{soil}}
//' \mjsdeqn{F_{itfl} \leq W_{soil}}
//' 
//' 
//' where
//' - \mjseqn{F_{itfl}} is `SOIL_inteflow_mm`
//' - \mjseqn{W_{soil}} is `water_SOIL_mm`
//' - \mjseqn{C_{soil}} is `capacity_SOIL_mm`
//' - \mjseqn{k^*} is estimated ratio
//' 
//' The output density distribution from 8 methods:
//'
//' @references
//' \insertAllCited{}
//' @return inteflow_mm (mm/m2)
//' @details
//' # **_GR4Jfix** \insertCite{GR4J_Perrin_2003}{HydroGallery}: 
//'
//' 
//' \mjsdeqn{k^* = 1 - \left[ 1 + \left(k \frac{W_{soil}}{C_{soil}} \right)^\gamma \right]^{-1/\gamma}}
//' where
//'   - \mjseqn{k} is `param_INTEFLOW_grf_k`
//'   - \mjseqn{\gamma} is `param_baseflow_grf_gamma`
//' @param param_INTEFLOW_grf_k <0.01, 1> coefficient parameter for [inteflow_GR4Jfix()]
//' @param param_INTEFLOW_grf_gamma <2, 7> exponential parameter for [baseflow_GR4Jfix()]
//' @export
// [[Rcpp::export]]
NumericVector inteflow_GR4Jfix(
    NumericVector SOIL_water_mm,
    NumericVector SOIL_capacity_mm,
    NumericVector param_INTEFLOW_grf_k,
    NumericVector param_INTEFLOW_grf_gamma
) 
{
  return SOIL_water_mm * (1 - vecpow((1 + vecpow(param_INTEFLOW_grf_k * SOIL_water_mm / SOIL_capacity_mm, param_INTEFLOW_grf_gamma)), -1.0 / param_INTEFLOW_grf_gamma));
}


//' @rdname inteflow
//' @details
//' # **_MaxPow**: 
//'
//' 
//' \mjsdeqn{F_{itfl} = M_{itfl} \left(\frac{W_{soil}}{C_{soil}} \right)^\gamma}
//' where
//'   - \mjseqn{M_{itfl}} is `SOIL_potentialInteflow_mm`
//'   - \mjseqn{\gamma} is `param_INTEFLOW_map_gamma`
//' @param param_INTEFLOW_map_gamma <0.1, 5> exponential parameter for [inteflow_MaxPow()]
//' @export
// [[Rcpp::export]]
NumericVector inteflow_MaxPow(
    NumericVector SOIL_water_mm,
    NumericVector SOIL_capacity_mm,
    NumericVector SOIL_potentialInteflow_mm,
    NumericVector param_INTEFLOW_map_gamma
)
{
  NumericVector inteflow_;
  
  inteflow_ = SOIL_potentialInteflow_mm * vecpow(SOIL_water_mm / SOIL_capacity_mm, param_INTEFLOW_map_gamma);
  
  return ifelse(inteflow_ > SOIL_water_mm, SOIL_water_mm, inteflow_) ;
}

//' @rdname inteflow
//' @details
//' # **_ThreshPow** 
//'
//' 
//' based on the `_MaxPow` and add the one threshold \mjseqn{\phi_b}: 
//' \mjsdeqn{F_{itfl} = 0, \quad \frac{W_{soil}}{C_{soil}} < \phi_b}
//' \mjsdeqn{F_{itfl} = M_{itfl} \left(\frac{\frac{W_{soil}}{C_{soil}} - \phi_b}{1-\phi_b} \right)^\gamma, \quad \frac{W_{soil}}{C_{soil}} \geq \phi_b}
//' where
//'   - \mjseqn{\phi_b} is `param_INTEFLOW_thp_thresh`
//'   - \mjseqn{\gamma} is `param_INTEFLOW_thp_gamma`
//' @param param_INTEFLOW_thp_thresh <0.1, 0.9> coefficient parameter for [inteflow_ThreshPow()]
//' @param param_INTEFLOW_thp_gamma <0.1, 5> exponential parameter for [inteflow_ThreshPow()]
//' @export
// [[Rcpp::export]]
NumericVector inteflow_ThreshPow(
    NumericVector SOIL_water_mm,
    NumericVector SOIL_capacity_mm,
    NumericVector SOIL_potentialInteflow_mm,
    NumericVector param_INTEFLOW_thp_thresh,
    NumericVector param_INTEFLOW_thp_gamma
)
{
  NumericVector inteflow_, inteflow_temp;
  inteflow_temp = (SOIL_water_mm / SOIL_capacity_mm - param_INTEFLOW_thp_thresh);
  inteflow_temp = ifelse(inteflow_temp < 0, 0, inteflow_temp);
  inteflow_ = SOIL_potentialInteflow_mm * vecpow(inteflow_temp / (1 - param_INTEFLOW_thp_thresh), param_INTEFLOW_thp_gamma);
  inteflow_ = ifelse(inteflow_ > SOIL_potentialInteflow_mm, SOIL_potentialInteflow_mm, inteflow_);
  return ifelse(inteflow_ > SOIL_water_mm, SOIL_water_mm, inteflow_) ;
}


//' @rdname inteflow
//' @details
//' # **_Arno** \insertCite{baseflow_Arno_1991,VIC2_Liang_1994}{HydroGallery} 
//'
//' 
//' has also in two cases divided by a threshold water content \mjseqn{\phi_b}:
//' (*This method is actually not the original method, but an analogy with `inteflow_Arno`*) 
//' \mjsdeqn{F_{itfl} = k M_{itfl} \frac{W_{soil}}{C_{soil}}, \quad \frac{W_{soil}}{C_{soil}} < \phi_b}
//' \mjsdeqn{F_{itfl} = k M_{itfl} \frac{W_{soil}}{C_{soil}} + (1-k) M_{itfl} \left(\frac{W_{soil} - W_s}{C_{soil} - W_s} \right)^2, \quad \frac{W_{soil}}{C_{soil}} \geq \phi_b}
//' \mjsdeqn{W_s = k C_{soil}}
//' where
//'   - \mjseqn{\phi_b} is `param_INTEFLOW_arn_thresh`
//'   - \mjseqn{k} is `param_INTEFLOW_arn_k`
//' @param param_INTEFLOW_arn_thresh <0.1, 0.9> coefficient parameter for [inteflow_ThreshPow()]
//' @param param_INTEFLOW_arn_k <0.1, 1> exponential parameter for [inteflow_ThreshPow()]
//' @export
// [[Rcpp::export]]
NumericVector inteflow_Arno(
    NumericVector SOIL_water_mm,
    NumericVector SOIL_capacity_mm,
    NumericVector SOIL_potentialInteflow_mm,
    NumericVector param_INTEFLOW_arn_thresh,
    NumericVector param_INTEFLOW_arn_k
)
{
  NumericVector inteflow_, inteflow_1, inteflow_2, Ws_Wc;
  Ws_Wc = SOIL_capacity_mm * param_INTEFLOW_arn_thresh;
  inteflow_1 = param_INTEFLOW_arn_k * SOIL_potentialInteflow_mm / (SOIL_capacity_mm) * SOIL_water_mm;
  inteflow_2 = param_INTEFLOW_arn_k * SOIL_potentialInteflow_mm / (SOIL_capacity_mm) * SOIL_water_mm + SOIL_potentialInteflow_mm * (1 - param_INTEFLOW_arn_k) * pow((SOIL_water_mm - Ws_Wc) / (SOIL_capacity_mm - Ws_Wc),2);
  inteflow_ = ifelse(SOIL_water_mm < Ws_Wc, inteflow_1, inteflow_2);
  inteflow_ = ifelse(SOIL_potentialInteflow_mm > Ws_Wc, SOIL_water_mm, inteflow_);
  inteflow_ = ifelse(inteflow_ > SOIL_potentialInteflow_mm, SOIL_potentialInteflow_mm, inteflow_);
  return ifelse(inteflow_ > SOIL_water_mm, SOIL_water_mm, inteflow_) ;
}


//' @rdname inteflow
//' @details
//' # **_BevenWood** \insertCite{percola_BevenWood_1983,TOPMODEL_Beven_1995}{HydroGallery}: 
//'
//' 
//' \mjsdeqn{k =  \frac{W_{soil}}{C_{soil} - W_{soil}} \quad {\rm and} \quad k \leq 1}
//' \mjsdeqn{F_{itfl} = k M_{itfl}}
//' where
//'   - \mjseqn{k_{fc}} is `SOIL_fieldCapacityPerc_1`
//'   - \mjseqn{\gamma} is `param_INTEFLOW_sup_gamma`
//' @export
// [[Rcpp::export]]
NumericVector inteflow_BevenWood(
    NumericVector SOIL_water_mm,
    NumericVector SOIL_capacity_mm,
    NumericVector SOIL_fieldCapacityPerc_1,
    NumericVector SOIL_potentialInteflow_mm
)
{
  NumericVector SOIL_inteflow_mm, SOIL_inteflowAvilibale_mm, SOIL_diff_mm, k_;
  SOIL_inteflowAvilibale_mm = SOIL_water_mm - SOIL_capacity_mm * (1-SOIL_fieldCapacityPerc_1);
  SOIL_inteflowAvilibale_mm = ifelse(SOIL_inteflowAvilibale_mm < 0, 0, SOIL_inteflowAvilibale_mm);
  SOIL_diff_mm = SOIL_capacity_mm - SOIL_water_mm;
  SOIL_diff_mm = ifelse(SOIL_diff_mm < SOIL_water_mm, SOIL_water_mm, SOIL_diff_mm);
  k_ = SOIL_water_mm / SOIL_diff_mm;
  SOIL_inteflow_mm = k_ * SOIL_potentialInteflow_mm;
  SOIL_inteflow_mm = ifelse(SOIL_water_mm > SOIL_inteflowAvilibale_mm, SOIL_inteflow_mm, 0.0);
  return ifelse(SOIL_inteflow_mm > SOIL_inteflowAvilibale_mm, SOIL_inteflowAvilibale_mm, SOIL_inteflow_mm) ;
}

//' @rdname inteflow
//' @details
//' # **_SupplyPow0**: 
//'
//' 
//' \mjsdeqn{F_{base} = k(W_{grnd})^\gamma}
//' where
//'   - \mjseqn{k} is `param_INTEFLOW_sup_k`
//'   - \mjseqn{\gamma} is `param_INTEFLOW_sup_gamma`
//' @param param_INTEFLOW_sp0_k <0.01, 1> coefficient parameter for [inteflow_SupplyPow0()]
//' @param param_INTEFLOW_sp0_gamma <0, 1> exponential parameter for [inteflow_SupplyPow0()]
//' @export
// [[Rcpp::export]]
NumericVector inteflow_SupplyPow0(
    NumericVector SOIL_water_mm,
    NumericVector param_INTEFLOW_sp0_k,
    NumericVector param_INTEFLOW_sp0_gamma
)
{
  NumericVector inteflow_;
  
  inteflow_ = param_INTEFLOW_sp0_k * vecpow(ceil(SOIL_water_mm), param_INTEFLOW_sp0_gamma);
  
  return ifelse(inteflow_ > SOIL_water_mm, SOIL_water_mm, inteflow_) ;
}


//' @rdname inteflow
//' @details
//' # **_SupplyPow**: 
//'
//' 
//' \mjsdeqn{k^* = k \left(\frac{W_{soil}}{C_{soil}} \right)^\gamma}
//' where
//'   - \mjseqn{k} is `param_INTEFLOW_sup_k`
//'   - \mjseqn{\gamma} is `param_INTEFLOW_sup_gamma`
//' @param param_INTEFLOW_sup_k <0.01, 1> coefficient parameter for [inteflow_SupplyPow()]
//' @param param_INTEFLOW_sup_gamma <0, 7> parameter for [inteflow_SupplyPow()]
//' @export
// [[Rcpp::export]]
NumericVector inteflow_SupplyPow(
    NumericVector SOIL_water_mm,
    NumericVector SOIL_capacity_mm,
    NumericVector param_INTEFLOW_sup_k,
    NumericVector param_INTEFLOW_sup_gamma
)
{
  NumericVector SOIL_inteflow_mm, k_;
  
  k_ = param_INTEFLOW_sup_k * vecpow((SOIL_water_mm / SOIL_capacity_mm), param_INTEFLOW_sup_gamma);
  SOIL_inteflow_mm = k_ * SOIL_water_mm;
  return ifelse(SOIL_inteflow_mm > SOIL_water_mm, SOIL_water_mm, SOIL_inteflow_mm) ;
}

//' @rdname inteflow
//' @details
//' # **_SupplyRatio**: 
//'
//' 
//' \mjsdeqn{k^* = k}
//' where
//'   - \mjseqn{k} is `param_INTEFLOW_sur_k`
//' @param param_INTEFLOW_sur_k <0.01, 1> coefficient parameter for [inteflow_SupplyRatio()]
//' @export
// [[Rcpp::export]]
NumericVector inteflow_SupplyRatio(
    NumericVector SOIL_water_mm,
    NumericVector param_INTEFLOW_sur_k
)
{
  
  return param_INTEFLOW_sur_k * SOIL_water_mm;
}

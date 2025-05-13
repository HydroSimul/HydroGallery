#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::interfaces(r, cpp)]]

//' **percolation**
//' @name percola
//' @inheritParams all_vari
//' @description
//' \loadmathjax
//' 
//' In hydrological modeling, percolation refers to the process by which water from the soil moves downward through the pores and cracks in the soil or rock.
//' This process is physically driven by a moisture gradient, but this is often simplified in conceptual percolation models \insertCite{Raven_Manual_35}{HydroGallery}.
//' It can be calculated by the water in the soil layer \mjseqn{W_{soil}},
//' it can also be tread as the part of the \mjseqn{W_{soil}}.
//' 
//' So we can give the function from:
//' 
//' \mjsdeqn{F_{prcl} = f_{percola}(D_{grnd}, D_{soil})}
//' 
//' 
//' to:
//' 
//' \mjsdeqn{F_{prcl} = f_{percola}(W_{soil}, C_{soil}, W_{grnd}, C_{grnd}, ...) = k^* W_{soil}}
//' \mjsdeqn{F_{prcl} \leq W_{soil}}
//' \mjsdeqn{F_{prcl} \leq C_{grnd} - W_{grnd}}
//' 
//' 
//' where
//' - \mjseqn{F_{prcl}} is `SOIL_percola_mm`
//' - \mjseqn{W_{soil}} is `water_SOIL_mm`
//' - \mjseqn{C_{soil}} is `capacity_SOIL_mm`
//' - \mjseqn{W_{grnd}} is `ground_water_mm`
//' - \mjseqn{C_{grnd}} is `capacity_water_mm`
//' - \mjseqn{k^*} is estimated ratio
//' 
//' The output density distribution from 8 methods:
//'
//' @references
//' \insertAllCited{}
//' @return percola_mm (mm/m2)
//' @details
//' # **_GR4J** \insertCite{GR4J_Perrin_2003}{HydroGallery}: 
//'
//' 
//' \mjsdeqn{k^* = 1 - \left[ 1 + \left(\frac{4}{9} \frac{W_{soil}}{C_{soil}} \right)^4 \right]^{-1/4}}
//' where
//'   - \mjseqn{k^*} is estimated ratio
//' @export
// [[Rcpp::export]]
arma::vec percola_GR4J(
    const arma::vec& SOIL_water_mm,
    const arma::vec& SOIL_capacity_mm
) {
    return SOIL_water_mm % (1 - arma::pow((1 + arma::pow(4.0 / 9.0 * SOIL_water_mm / SOIL_capacity_mm, 4)), -0.25));
}

//' @rdname percola
//' @details
//' # **_GR4Jfix** \insertCite{GR4J_Perrin_2003}{HydroGallery}: 
//'
//' 
//' \mjsdeqn{k^* = 1 - \left[ 1 + \left(k \frac{W_{soil}}{C_{soil}} \right)^4 \right]^{-1/4}}
//' where
//'   - \mjseqn{k} is `param_PERCOLA_grf_k`
//' @param param_PERCOLA_grf_k <0.01, 1> coefficient parameter for [percola_GR4Jfix()]
//' @export
// [[Rcpp::export]]
arma::vec percola_GR4Jfix(
    const arma::vec& SOIL_water_mm,
    const arma::vec& SOIL_capacity_mm,
    const arma::vec& param_PERCOLA_grf_k
) {
    return SOIL_water_mm % (1 - arma::pow((1 + arma::pow(param_PERCOLA_grf_k % SOIL_water_mm / SOIL_capacity_mm, 4)), -0.25));
}

//' @rdname percola
//' @details
//' # **_MaxPow**: 
//'
//' 
//' \mjsdeqn{F_{prcl} = M_{prcl} \left(\frac{W_{soil}}{C_{soil}} \right)^\gamma}
//' where
//'   - \mjseqn{M_{prcl}} is `SOIL_potentialPercola_mm`
//'   - \mjseqn{\gamma} is `param_baseflow_map_gamma`
//' @param param_PERCOLA_map_gamma <0.1, 5> exponential parameter for [percola_MaxPow()]
//' @export
// [[Rcpp::export]]
arma::vec percola_MaxPow(
    const arma::vec& SOIL_water_mm,
    const arma::vec& SOIL_capacity_mm,
    const arma::vec& SOIL_potentialPercola_mm,
    const arma::vec& param_PERCOLA_map_gamma
) {
    arma::vec percola_ = SOIL_potentialPercola_mm % arma::pow(SOIL_water_mm / SOIL_capacity_mm, param_PERCOLA_map_gamma);
    return arma::min(percola_, SOIL_water_mm);
}

//' @rdname percola
//' @details
//' # **_ThreshPow** 
//'
//' 
//' based on the `_MaxPow` and add the one threshold \mjseqn{\phi_b}: 
//' \mjsdeqn{F_{prcl} = 0, \quad \frac{W_{soil}}{C_{soil}} < \phi_b}
//' \mjsdeqn{F_{prcl} = M_{prcl} \left(\frac{\frac{W_{soil}}{C_{soil}} - \phi_b}{1-\phi_b} \right)^\gamma, \quad \frac{W_{soil}}{C_{soil}} \geq \phi_b}
//' where
//'   - \mjseqn{\phi_b} is `param_PERCOLA_thp_thresh`
//'   - \mjseqn{\gamma} is `param_PERCOLA_thp_gamma`
//' @param param_PERCOLA_thp_thresh <0.1, 0.9> coefficient parameter for [percola_ThreshPow()]
//' @param param_PERCOLA_thp_gamma <0.1, 5> exponential parameter for [percola_ThreshPow()]
//' @export
// [[Rcpp::export]]
arma::vec percola_ThreshPow(
    const arma::vec& SOIL_water_mm,
    const arma::vec& SOIL_capacity_mm,
    const arma::vec& SOIL_potentialPercola_mm,
    const arma::vec& param_PERCOLA_thp_thresh,
    const arma::vec& param_PERCOLA_thp_gamma
) {
    arma::vec percola_temp = arma::max(arma::zeros<arma::vec>(SOIL_water_mm.n_elem), SOIL_water_mm / SOIL_capacity_mm - param_PERCOLA_thp_thresh);
    arma::vec percola_ = SOIL_potentialPercola_mm % arma::pow(percola_temp / (1 - param_PERCOLA_thp_thresh), param_PERCOLA_thp_gamma);
    percola_ = arma::min(percola_, SOIL_potentialPercola_mm);
    return arma::min(percola_, SOIL_water_mm);
}

//' @rdname percola
//' @details
//' # **_Arno** \insertCite{baseflow_Arno_1991,VIC2_Liang_1994}{HydroGallery} 
//'
//' 
//' has also in two cases divided by a threshold water content \mjseqn{\phi_b}:
//' (*This method is actually not the original method, but an analogy with `baseflow_Arno`*) 
//' \mjsdeqn{F_{prcl} = k M_{prcl} \frac{W_{soil}}{C_{soil}}, \quad \frac{W_{soil}}{C_{soil}} < \phi_b}
//' \mjsdeqn{F_{prcl} = k M_{prcl} \frac{W_{soil}}{C_{soil}} + (1-k) M_{prcl} \left(\frac{W_{soil} - W_s}{C_{soil} - W_s} \right)^2, \quad \frac{W_{soil}}{C_{soil}} \geq \phi_b}
//' \mjsdeqn{W_s = k C_{soil}}
//' where
//'   - \mjseqn{\phi_b} is `param_PERCOLA_arn_thresh`
//'   - \mjseqn{k} is `param_PERCOLA_arn_k`
//' @param param_PERCOLA_arn_thresh <0.1, 0.9> coefficient parameter for [percola_ThreshPow()]
//' @param param_PERCOLA_arn_k <0.1, 1> exponential parameter for [percola_ThreshPow()]
//' @export
// [[Rcpp::export]]
arma::vec percola_Arno(
    const arma::vec& SOIL_water_mm,
    const arma::vec& SOIL_capacity_mm,
    const arma::vec& SOIL_potentialPercola_mm,
    const arma::vec& param_PERCOLA_arn_thresh,
    const arma::vec& param_PERCOLA_arn_k)
{
  // Calculate threshold water content
  arma::vec Ws_Wc = SOIL_capacity_mm % param_PERCOLA_arn_thresh;
  
  // Calculate percolation components
  arma::vec percola_1 = param_PERCOLA_arn_k % SOIL_potentialPercola_mm / 
                       SOIL_capacity_mm % SOIL_water_mm;
  
  arma::vec percola_2 = percola_1 + SOIL_potentialPercola_mm % 
                       (1.0 - param_PERCOLA_arn_k) % 
                       arma::pow((SOIL_water_mm - Ws_Wc) / 
                       (SOIL_capacity_mm - Ws_Wc), 2);
  
  // Combine components based on conditions
  arma::vec percola_ = arma::zeros<arma::vec>(SOIL_water_mm.n_elem);
  
  // First condition: SOIL_water_mm < Ws_Wc
  arma::uvec cond1 = arma::find(SOIL_water_mm < Ws_Wc);
  percola_.elem(cond1) = percola_1.elem(cond1);
  
  // Second condition: SOIL_water_mm >= Ws_Wc
  arma::uvec cond2 = arma::find(SOIL_water_mm >= Ws_Wc);
  percola_.elem(cond2) = percola_2.elem(cond2);
  
  // Additional constraints
  percola_.elem(arma::find(SOIL_potentialPercola_mm <= Ws_Wc)) = 
    arma::min(arma::min(percola_, SOIL_potentialPercola_mm), SOIL_water_mm);
  
  return percola_;
}

//' @rdname percola
//' @details
//' # **_BevenWood** \insertCite{percola_BevenWood_1983,TOPMODEL_Beven_1995}{HydroGallery}: 
//'
//' 
//' \mjsdeqn{k =  \frac{W_{soil}}{C_{soil} - W_{soil}} \quad {\rm and} \quad k \leq 1}
//' \mjsdeqn{F_{prcl} = k M_{prcl}}
//' where
//'   - \mjseqn{k_{fc}} is `SOIL_fieldCapacityPerc_1`
//'   - \mjseqn{\gamma} is `param_PERCOLA_sup_gamma`
//' @export
// [[Rcpp::export]]
arma::vec percola_BevenWood(
    const arma::vec& SOIL_water_mm,
    const arma::vec& SOIL_capacity_mm,
    const arma::vec& SOIL_fieldCapacityPerc_1,
    const arma::vec& SOIL_potentialPercola_mm
)
{
    arma::vec SOIL_percola_mm, SOIL_percolaAvilibale_mm, SOIL_diff_mm, k_;
    SOIL_percolaAvilibale_mm = SOIL_water_mm - SOIL_capacity_mm % (1 - SOIL_fieldCapacityPerc_1);
    SOIL_percolaAvilibale_mm = arma::clamp(SOIL_percolaAvilibale_mm, 0.0, arma::datum::inf);
    SOIL_diff_mm = SOIL_capacity_mm - SOIL_water_mm;
    SOIL_diff_mm = arma::max(SOIL_diff_mm, SOIL_water_mm);
    k_ = SOIL_water_mm / SOIL_diff_mm;
    SOIL_percola_mm = k_ % SOIL_potentialPercola_mm;
    SOIL_percola_mm = arma::min(SOIL_percola_mm, SOIL_percolaAvilibale_mm);
    return SOIL_percola_mm;
}

//' @rdname percola
//' @details
//' # **_SupplyPow**: 
//'
//' \mjsdeqn{k^* = k \left(\frac{W_{soil}}{C_{soil}} \right)^\gamma}
//' where
//'   - \mjseqn{k} is `param_PERCOLA_sup_k`
//'   - \mjseqn{\gamma} is `param_PERCOLA_sup_gamma`
//' @param param_PERCOLA_sup_k <0.01, 1> coefficient parameter for [percola_SupplyPow()]
//' @param param_PERCOLA_sup_gamma <0, 7> parameter for [percola_SupplyPow()]
//' @export
// [[Rcpp::export]]
arma::vec percola_SupplyPow(
    const arma::vec& SOIL_water_mm,
    const arma::vec& SOIL_capacity_mm,
    const arma::vec& param_PERCOLA_sup_k,
    const arma::vec& param_PERCOLA_sup_gamma
)
{
    arma::vec k_, SOIL_percola_mm;
    k_ = param_PERCOLA_sup_k % arma::pow(SOIL_water_mm / SOIL_capacity_mm, param_PERCOLA_sup_gamma);
    SOIL_percola_mm = k_ % SOIL_water_mm;
    return arma::min(SOIL_percola_mm, SOIL_water_mm);
}

//' @rdname percola
//' @details
//' # **_SupplyRatio**: 
//'
//' \mjsdeqn{k^* = k}
//' where
//'   - \mjseqn{k} is `param_PERCOLA_sur_k`
//' @param param_PERCOLA_sur_k <0.01, 1> coefficient parameter for [percola_SupplyRatio()]
//' @export
// [[Rcpp::export]]
arma::vec percola_SupplyRatio(
    const arma::vec& SOIL_water_mm,
    const arma::vec& param_PERCOLA_sur_k
)
{
    return param_PERCOLA_sur_k % SOIL_water_mm;
}

//' @rdname percola
//' @param param_PERCOLA_wat_01 <0, 1> 0: percolation, 1: non percolation [percola_WaterGAP3()]
//' @param param_PERCOLA_wat_thresh <10, 15> coefficient parameter for [percola_WaterGAP3()]
//' @param param_PERCOLA_wat_k <0.1, 1> exponential parameter for [percola_WaterGAP3()]
//' @export
// [[Rcpp::export]]
arma::vec percola_WaterGAP3(
    const arma::vec& LAND_water_mm,
    const arma::vec& SOIL_potentialPercola_mm,
    const arma::uvec& param_PERCOLA_wat_01,
    const arma::vec& param_PERCOLA_wat_thresh,
    const arma::vec& param_PERCOLA_wat_k)
{
    // Calculate potential percolation for all cells
    arma::vec potential_percola = arma::min(SOIL_potentialPercola_mm,
                                          param_PERCOLA_wat_k % LAND_water_mm);
    
    // Logical mask for arid regions
    arma::uvec arid_mask = arma::find(param_PERCOLA_wat_01 == 1 && 
                                     LAND_water_mm < param_PERCOLA_wat_thresh);
    
    // Set percolation to 0 for arid regions
    potential_percola.elem(arid_mask).zeros();
    
    // Apply final constraint
    return arma::min(potential_percola, LAND_water_mm);
}
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::interfaces(r, cpp)]]

//' **lateral flux**
//' @name lateral
//' @inheritParams all_vari
//' @description
//' \loadmathjax
//' 
//' In hydrological modeling, lateral flow refers to the process by which water flows horizontally through the soil or aquifer, rather than vertically.
//' It is typically represented by a loss term in the water balance equation, so it also named as groundwater exchange (e.g. GR4J \insertCite{GR4J_Perrin_2003}{HydroGallery}).
//' The flux of lateral exchange is always calculated (only) by the water in the ground layer \mjseqn{W_{grnd}}. 
//' Unlike other fluxes, the lateral exchange can be positive or negative, 
//' with positive indicating a supply from other regions and negative indicating distribution to other regions.
//' 
//' This process is so flexible that we must carefully use it, 
//' because it can easily destroy the waster balance in the research catchment.
//' 
//' \mjsdeqn{F_{ltrl} = f_{lateral}(D_{grnd})}
//' 
//' 
//' to:
//' 
//' \mjsdeqn{F_{ltrl} = f_{lateral}(W_{grnd}, C_{grnd}, ...)}
//' 
//' 
//' where
//' - \mjseqn{W_{grnd}} is `GROUND_water_mm`
//' - \mjseqn{C_{grnd}} is `GROUND_capacity_mm`, but not all the methods need the \mjseqn{C_{grnd}}
//' 
//' The output density distribution from 6 methods:
//'
//' @references
//' \insertAllCited{}
//' @return lateral_mm (mm/m2)
//' 
//' 
//' 
//' @details
//' # **_SupplyPow**: 
//'
//' 
//' \mjsdeqn{F_{ltrl} = k \left( \frac{W_{grnd}}{C_{grnd}} \right)^\gamma  W_{grnd}}
//' where
//'   - \mjseqn{k} is `param_LATERAL_sup_k`
//'   - \mjseqn{\gamma} is `param_LATERAL_sup_gamma`
//' @param param_LATERAL_sup_k <-1, 1> coefficient parameter for [lateral_SupplyPow()]
//' @param param_LATERAL_sup_gamma <0.01, 5> parameters for [lateral_SupplyPow()]
//' @export
// [[Rcpp::export]]
arma::vec lateral_SupplyPow(
    const arma::vec& GROUND_water_mm,
    const arma::vec& GROUND_capacity_mm,
    const arma::vec& param_LATERAL_sup_k,
    const arma::vec& param_LATERAL_sup_gamma
)
{
  arma::vec GROUND_lateral_mm, k_;
  arma::vec GROUND_diff_mm = (GROUND_capacity_mm - GROUND_water_mm);
  
  k_ = param_LATERAL_sup_k % arma::pow((GROUND_water_mm / GROUND_capacity_mm), param_LATERAL_sup_gamma);
  GROUND_lateral_mm = k_ % GROUND_water_mm;
  GROUND_lateral_mm = arma::min(GROUND_lateral_mm, GROUND_diff_mm);
  return arma::min(GROUND_lateral_mm, -GROUND_water_mm);
}

//' @rdname lateral
//' @details
//' # **_SupplyRatio**: 
//'
//' 
//' \mjsdeqn{F_{ltrl} = k * W_{grnd}}
//' where
//'   - \mjseqn{k} is `param_LATERAL_sur_k`
//' @param param_LATERAL_sur_k <-2, 1> coefficient parameter for [lateral_SupplyRatio()]
//' @export
// [[Rcpp::export]]
arma::vec lateral_SupplyRatio(
    const arma::vec& GROUND_water_mm,
    const arma::vec& param_LATERAL_sur_k
)
{
  arma::vec GROUND_lateral_mm = param_LATERAL_sur_k % GROUND_water_mm;
  return arma::min(GROUND_lateral_mm, -GROUND_water_mm);
}

//' @rdname lateral
//' @details
//' # **_GR4J** \insertCite{GR4J_Perrin_2003}{HydroGallery}: 
//'
//' 
//' \mjsdeqn{F_{ltrl} = M_{ltrl} \left( \frac{W_{grnd}}{C_{grnd}} \right)^{7/2}  }
//' where
//'   - \mjseqn{M_{ltrl}} is `GROUND_potentialLateral_mm`
//' @export
// [[Rcpp::export]]
arma::vec lateral_GR4J(
    const arma::vec& GROUND_water_mm,
    const arma::vec& GROUND_capacity_mm,
    const arma::vec& GROUND_potentialLateral_mm
) 
{
  arma::vec GROUND_lateral_mm;
  arma::vec GROUND_diff_mm = (GROUND_capacity_mm - GROUND_water_mm);
  GROUND_lateral_mm = GROUND_potentialLateral_mm % arma::pow((GROUND_water_mm / GROUND_capacity_mm), 3.5);
  GROUND_lateral_mm = arma::min(GROUND_lateral_mm, GROUND_diff_mm);
  return arma::min(GROUND_lateral_mm, -GROUND_water_mm);
}

//' @rdname lateral
//' @details
//' # **_GR4Jfix** \insertCite{GR4J_Perrin_2003}{HydroGallery} 
//'
//' 
//' based on `_GR4J` use a new parameter to replace the numer 4: 
//' \mjsdeqn{F_{ltrl} = M_{ltrl} \left( \frac{W_{grnd}}{C_{grnd}} \right)^\gamma  }
//' where
//'   - \mjseqn{\gamma} is `param_LATERAL_grf_gamma`
//' @param param_LATERAL_grf_gamma <0.01, 5> parameter for [lateral_GR4Jfix()]
//' @export
// [[Rcpp::export]]
arma::vec lateral_GR4Jfix(
    const arma::vec& GROUND_water_mm,
    const arma::vec& GROUND_capacity_mm,
    const arma::vec& GROUND_potentialLateral_mm,
    const arma::vec& param_LATERAL_grf_gamma
) 
{
  arma::vec GROUND_lateral_mm;
  arma::vec GROUND_diff_mm = (GROUND_capacity_mm - GROUND_water_mm);
  GROUND_lateral_mm = GROUND_potentialLateral_mm % arma::pow((GROUND_water_mm / GROUND_capacity_mm), param_LATERAL_grf_gamma);
  GROUND_lateral_mm = arma::min(GROUND_lateral_mm, GROUND_diff_mm);
  return arma::min(GROUND_lateral_mm, -GROUND_water_mm);
}

//' @rdname lateral
//' @details
//' # **_ThreshPow** 
//'
//' 
//' based on the `_GR4Jfix` and add the one threshold \mjseqn{\phi_b}: 
//' \mjsdeqn{F_{ltrl} = 0, \quad \frac{W_{grnd}}{C_{grnd}} < \phi_b}
//' \mjsdeqn{F_{ltrl} = M_{ltrl} \left(\frac{\frac{W_{grnd}}{C_{grnd}} - \phi_b}{1-\phi_b} \right)^\gamma, \quad \frac{W_{grnd}}{C_{grnd}} \geq \phi_b}
//' where
//'   - \mjseqn{\phi_b} is `param_LATERAL_thp_thresh`
//'   - \mjseqn{\gamma} is `param_LATERAL_thp_gamma`
//' @param param_LATERAL_thp_thresh <0.1, 0.9> coefficient parameter for [lateral_ThreshPow()]
//' @param param_LATERAL_thp_gamma <0.1, 5> exponential parameter for [lateral_ThreshPow()]
//' @export
// [[Rcpp::export]]
arma::vec lateral_ThreshPow(
    const arma::vec& GROUND_water_mm,
    const arma::vec& GROUND_capacity_mm,
    const arma::vec& GROUND_potentialLateral_mm,
    const arma::vec& param_LATERAL_thp_thresh,
    const arma::vec& param_LATERAL_thp_gamma
)
{
  arma::vec GROUND_lateral_mm, lateral_temp;
  arma::vec GROUND_diff_mm = (GROUND_capacity_mm - GROUND_water_mm);
  lateral_temp = (GROUND_water_mm / GROUND_capacity_mm - param_LATERAL_thp_thresh);
  lateral_temp = arma::clamp(lateral_temp, 0.0, arma::datum::inf);
  
  GROUND_lateral_mm = GROUND_potentialLateral_mm % arma::pow(lateral_temp / (1 - param_LATERAL_thp_thresh), param_LATERAL_thp_gamma);
  GROUND_lateral_mm = arma::min(GROUND_lateral_mm, GROUND_diff_mm);
  return arma::min(GROUND_lateral_mm, -GROUND_water_mm);
}

//' @rdname lateral
//' @details
//' # **_Arno** \insertCite{baseflow_Arno_1991,VIC2_Liang_1994}{HydroGallery} 
//'
//' 
//' has also in two cases divided by a threshold water content \mjseqn{\phi_b}: 
//' \mjsdeqn{F_{ltrl} = k M_{ltrl} \frac{W_{grnd}}{C_{grnd}}, \quad \frac{W_{grnd}}{C_{grnd}} < \phi_b}
//' \mjsdeqn{F_{ltrl} = k M_{ltrl} \frac{W_{grnd}}{C_{grnd}} + (1-k) M_{ltrl} \left(\frac{W_{grnd} - W_s}{C_{grnd} - W_s} \right)^2, \quad \frac{W_{grnd}}{C_{grnd}} \geq \phi_b}
//' \mjsdeqn{W_s = k C_{grnd}}
//' where
//'   - \mjseqn{\phi_b} is `param_LATERAL_arn_thresh`
//'   - \mjseqn{k} is `param_LATERAL_arn_k`
//' @param param_LATERAL_arn_thresh <0.1, 0.9> coefficient parameter for [lateral_ThreshPow()]
//' @param param_LATERAL_arn_k <0.1, 1> exponential parameter for [lateral_ThreshPow()]
//' @export
// [[Rcpp::export]]
arma::vec lateral_Arno(
    const arma::vec& GROUND_water_mm,
    const arma::vec& GROUND_capacity_mm,
    const arma::vec& GROUND_potentialLateral_mm,
    const arma::vec& param_LATERAL_arn_thresh,
    const arma::vec& param_LATERAL_arn_k
)
{
  arma::vec GROUND_lateral_mm, lateral_1, lateral_2, Ws_Wc;
  arma::vec GROUND_diff_mm = GROUND_capacity_mm - GROUND_water_mm;

  // Threshold water level
  Ws_Wc = GROUND_capacity_mm % param_LATERAL_arn_thresh;

  // Case 1: Below threshold
  lateral_1 = param_LATERAL_arn_k % GROUND_potentialLateral_mm / GROUND_capacity_mm % GROUND_water_mm;

  // Case 2: Above threshold
  arma::vec ratio = (GROUND_water_mm - Ws_Wc) / (GROUND_capacity_mm - Ws_Wc);
  lateral_2 = lateral_1 + GROUND_potentialLateral_mm % (1.0 - param_LATERAL_arn_k) % arma::square(ratio);

  // Combine cases
  GROUND_lateral_mm = arma::vec(GROUND_water_mm.n_elem);
  GROUND_lateral_mm.elem(arma::find(GROUND_water_mm < Ws_Wc)) = lateral_1.elem(arma::find(GROUND_water_mm < Ws_Wc));
  GROUND_lateral_mm.elem(arma::find(GROUND_water_mm >= Ws_Wc)) = lateral_2.elem(arma::find(GROUND_water_mm >= Ws_Wc));

  // Cap to water content if potential is greater than available
  arma::uvec idx = arma::find(GROUND_potentialLateral_mm > Ws_Wc);
  GROUND_lateral_mm.elem(idx) = GROUND_water_mm.elem(idx);

  // Bounds enforcement for positive potential lateral flow
  arma::uvec idx_pos = arma::find((GROUND_lateral_mm > GROUND_potentialLateral_mm) % (GROUND_potentialLateral_mm > 0));
  GROUND_lateral_mm.elem(idx_pos) = GROUND_potentialLateral_mm.elem(idx_pos);

  // Bounds enforcement for negative potential lateral flow
  arma::uvec idx_neg = arma::find((GROUND_lateral_mm < GROUND_potentialLateral_mm) % (GROUND_potentialLateral_mm < 0));
  GROUND_lateral_mm.elem(idx_neg) = GROUND_potentialLateral_mm.elem(idx_neg);

  // Ensure it doesn't exceed drainable water
  // Ensure lateral flow does not exceed available water or go below -GROUND_water_mm
  GROUND_lateral_mm = arma::min(GROUND_lateral_mm, GROUND_diff_mm);
  GROUND_lateral_mm = arma::max(GROUND_lateral_mm, -GROUND_water_mm);
  
  return GROUND_lateral_mm;
}

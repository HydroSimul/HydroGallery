#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::interfaces(r, cpp)]]

//' **infiltration**
//' @name infilt
//' @inheritParams all_vari
//' @description
//' \loadmathjax
//' In hydrological modeling, infiltration refers to the process by which water from precipitation snowmelt or irrigation enters the soil \insertCite{Handbook_Hydrology_1993}{HydroGallery}. 
//' 
//' Under the concept of the conceptual HM, the flux of infiltration is always calculated by the amount of water on the land \mjseqn{W_{land}}, 
//' which can be precipitation, precipitation after interception, or precipitation with snowmelt, among others. 
//' The second point to consider is the water acceptability of the soil layer (\mjseqn{C_{soil} - W_{soil}}).
//' 
//' So we can give the function from:
//' 
//' \mjsdeqn{F_{iflt} = f_{infilt}(D_{land}, D_{soil})}
//' 
//' 
//' to:
//' 
//' \mjsdeqn{F_{iflt} = f_{infilt}(W_{land}, W_{soil}, C_{soil}, ...)}
//' 
//' some methods will tread the infiltartion as the part of th pounded water so there is also:
//' 
//' \mjsdeqn{F_{iflt} = k^* W_{land}}
//' 
//' 
//' where
//' - \mjseqn{F_{iflt}} is `infilt_mm`
//' - \mjseqn{W_{land}} is `LAND_water_mm`
//' - \mjseqn{W_{soil}} is `SOIL_water_mm`
//' - \mjseqn{C_{soil}} is `SOIL_capacity_mm`
//' - \mjseqn{k^*} is estimated ratio.
//' 
//' The output density distribution from 9 methods:
//'
//' @references
//' \insertAllCited{}
//' @return flux of infiltration from land surface to soil layer
//' 
//' @details
//' # **_GR4J** \insertCite{GR4J_Perrin_2003}{HydroGallery}: 
//'
//' 
//' \mjsdeqn{F_{iflt}=\frac{C_{soil}\left(1-\left(\frac{W_{soil}}{C_{soil}}\right)^{2}\right) \tanh \left(\frac{W_{land}}{C_{soil}}\right)}{1+\frac{W_{soil}}{C_{soil}} \tanh \left(\frac{W_{land}}{C_{soil}}\right)}}
//' @export
// [[Rcpp::export]]
arma::vec infilt_GR4J(
    const arma::vec& LAND_water_mm,
    const arma::vec& SOIL_water_mm,
    const arma::vec& SOIL_capacity_mm
) 
{
  arma::vec SOIL_diff_mm = SOIL_capacity_mm - SOIL_water_mm;
  arma::vec limit_mm = arma::min(arma::join_cols(SOIL_diff_mm, LAND_water_mm), 1);
  
  arma::vec tanh_pn_x1 = arma::tanh(LAND_water_mm / SOIL_capacity_mm);
  arma::vec s_x1 = SOIL_water_mm / SOIL_capacity_mm;
  arma::vec infilt_water_mm = SOIL_capacity_mm % (1 - square(s_x1)) % tanh_pn_x1 / (1 + s_x1 % tanh_pn_x1); 
  
  return arma::min(infilt_water_mm, limit_mm);
}

//' @rdname infilt
//' @param param_INFILT_ubc_P0AGEN <0.1, 4> coefficient parameter for [infilt_UBC()]
//' @details
//' # **_UBC** \insertCite{UBC_Quick_1977}{HydroGallery}: 
//'
//' 
//' estimate the ratio \mjseqn{k^*} as:
//' \mjsdeqn{k^* = p_{imper} 10^{\frac{W_{soil}-C_{soil}}{p_{AGEN}}}}
//' where
//'   - \mjseqn{p_{imper}} is `LAND_impermeableFrac_1`
//'   - \mjseqn{p_{AGEN}} is `param_INFILT_ubc_P0AGEN`
//' @export
// [[Rcpp::export]]
arma::vec infilt_UBC(
    const arma::vec& LAND_water_mm, 
    const arma::vec& LAND_impermeableFrac_1, 
    const arma::vec& SOIL_water_mm,
    const arma::vec& SOIL_capacity_mm, 
    const arma::vec& param_INFILT_ubc_P0AGEN
)
{
  arma::vec SOIL_diff_mm = SOIL_capacity_mm - SOIL_water_mm;
  arma::vec limit_mm = arma::min(arma::join_cols(SOIL_diff_mm, LAND_water_mm), 1);
  
  arma::vec k_ = 1 - LAND_impermeableFrac_1 % arma::exp(-SOIL_diff_mm / (SOIL_capacity_mm % param_INFILT_ubc_P0AGEN) * std::log(10.0));
  arma::vec infilt_water_mm = LAND_water_mm % k_;
  
  return arma::min(infilt_water_mm, limit_mm);
}

//' @rdname infilt
//' @details
//' # **_SupplyRatio**: 
//'
//' 
//' is a very simple method, which estimate only the pounded water:
//' \mjsdeqn{k^* = k}
//' where
//'   - \mjseqn{k} is `param_INFILT_sur_k`
//' @param param_INFILT_sur_k <0.01, 1> coefficient parameter for [infilt_SupplyRatio()]
//' @return infilt_mm (mm/m2) 
//' @export
// [[Rcpp::export]]
arma::vec infilt_SupplyRatio(
    const arma::vec& LAND_water_mm,
    const arma::vec& SOIL_water_mm,
    const arma::vec& SOIL_capacity_mm, 
    const arma::vec& param_INFILT_sur_k
)
{
  arma::vec SOIL_diff_mm = SOIL_capacity_mm - SOIL_water_mm;
  arma::vec limit_mm = arma::min(arma::join_cols(SOIL_diff_mm, LAND_water_mm), 1);
  
  arma::vec infilt_water_mm = param_INFILT_sur_k % LAND_water_mm;
  
  return arma::min(infilt_water_mm, limit_mm);
}

//' @rdname infilt
//' @details
//' # **_AcceptRatio**: 
//'
//' 
//' \mjsdeqn{F_{iflt} = k (C_{soil} - W_{soil})}
//' where
//'   - \mjseqn{k} is `param_INFILT_acr_k`
//' @param param_INFILT_acr_k <0.01, 1> coefficient parameter for [infilt_AcceptRatio()]
//' @export
// [[Rcpp::export]]
arma::vec infilt_AcceptRatio(
    const arma::vec& LAND_water_mm, 
    const arma::vec& SOIL_water_mm,
    const arma::vec& SOIL_capacity_mm, 
    const arma::vec& param_INFILT_acr_k
)
{
  arma::vec SOIL_diff_mm = SOIL_capacity_mm - SOIL_water_mm;
  arma::vec limit_mm = arma::min(arma::join_cols(SOIL_diff_mm, LAND_water_mm), 1);
  
  arma::vec infilt_water_mm = SOIL_diff_mm % param_INFILT_acr_k;
  
  return arma::min(infilt_water_mm, limit_mm);
}

//' @rdname infilt
//' @details
//' # **_SupplyPow**: 
//'
//' 
//' is a very simple method, which estimate only the pounded water:
//' \mjsdeqn{F_{iflt} = kW_{land}^{\gamma}}
//' where
//'   - \mjseqn{k} is `param_INFILT_sup_k`
//'   - \mjseqn{\gamma} is `param_INFILT_sup_gamma`
//' @param param_INFILT_sup_k <0.01, 1> coefficient parameter for [infilt_SupplyPow()]
//' @param param_INFILT_sup_gamma <0, 1> parameters for [infilt_SupplyPow()]
//' @export
// [[Rcpp::export]]
arma::vec infilt_SupplyPow(
    const arma::vec& LAND_water_mm,
    const arma::vec& SOIL_water_mm,
    const arma::vec& SOIL_capacity_mm, 
    const arma::vec& param_INFILT_sup_k,
    const arma::vec& param_INFILT_sup_gamma
)
{
  arma::vec SOIL_diff_mm = SOIL_capacity_mm - SOIL_water_mm;
  arma::vec limit_mm = arma::min(arma::join_cols(SOIL_diff_mm, LAND_water_mm), 1);
  
  arma::vec infilt_water_mm = param_INFILT_sup_k % pow(arma::ceil(LAND_water_mm), param_INFILT_sup_gamma);
  infilt_water_mm = arma::min(infilt_water_mm, LAND_water_mm);
  
  return arma::min(infilt_water_mm, limit_mm);
}

//' @rdname infilt
//' @details
//' # **_AcceptPow**: 
//'
//' 
//' \mjsdeqn{F_{iflt} = k \left(\frac{C_{soil} - W_{soil}}{C_{soil}} \right)^{\gamma}}
//' where
//'   - \mjseqn{k} is `param_INFILT_acp_k`
//'   - \mjseqn{\gamma} is `param_INFILT_acp_gamma`
//' @param param_INFILT_acp_k <0.01, 1> coefficient parameter for [infilt_AcceptPow()]
//' @param param_INFILT_acp_gamma <0.001, 5> parameters for [infilt_AcceptPow()]
//' @export
// [[Rcpp::export]]
arma::vec infilt_AcceptPow(
    const arma::vec& LAND_water_mm, 
    const arma::vec& SOIL_water_mm,
    const arma::vec& SOIL_capacity_mm,
    const arma::vec& param_INFILT_acp_k,
    const arma::vec& param_INFILT_acp_gamma
)
{
  arma::vec SOIL_diff_mm = SOIL_capacity_mm - SOIL_water_mm;
  arma::vec limit_mm = arma::min(arma::join_cols(SOIL_diff_mm, LAND_water_mm), 1);
  
  arma::vec k_ = param_INFILT_acp_k % pow(SOIL_diff_mm / SOIL_capacity_mm, param_INFILT_acp_gamma);
  arma::vec infilt_water_mm = k_ % SOIL_diff_mm;
  
  return arma::min(infilt_water_mm, LAND_water_mm);
}

//' @rdname infilt
//' @details
//' # **_HBV** \insertCite{HBV_Lindstrom_1997}{HydroGallery}: 
//'
//' 
//' estimate the ratio \mjseqn{k^*} as:
//' \mjsdeqn{k^* = 1-\left(\frac{W_{soil}}{C_{soil}}\right)^{\beta}}
//' where
//'   - \mjseqn{\beta} is `param_INFILT_hbv_beta`
//' @param param_INFILT_hbv_beta <0.001, 5> parameters for [infilt_HBV()]
//' @export
// [[Rcpp::export]]
arma::vec infilt_HBV(
    const arma::vec& LAND_water_mm, 
    const arma::vec& SOIL_water_mm,
    const arma::vec& SOIL_capacity_mm, 
    const arma::vec& param_INFILT_hbv_beta 
)
{
  arma::vec SOIL_diff_mm = SOIL_capacity_mm - SOIL_water_mm;
  arma::vec limit_mm = arma::min(arma::join_cols(SOIL_diff_mm, LAND_water_mm), 1);
  
  arma::vec k_ = (1 - pow(SOIL_water_mm / SOIL_capacity_mm, param_INFILT_hbv_beta));
  arma::vec infilt_water_mm = LAND_water_mm % k_;
  
  return arma::min(infilt_water_mm, limit_mm);
}

//' @rdname infilt
//' @details
//' # **_XAJ** \insertCite{XAJ_Zhao_1992}{HydroGallery}:
//'
//' 
//' \mjsdeqn{F_{iflt} = MM  \frac{\left( \frac{MM - AU}{MM} \right)^{B+1} - \left( \frac{MM - AU - W_{land}}{MM} \right)^{B+1}}{B+1}}
//' \mjsdeqn{AU = MM - \left( \frac{(1 - W_{soil})(B+1)}{MM} \right)^{1 / B - 1}  }
//' \mjsdeqn{MM = C_{soil}(B+1)  }
//' where
//'   - \mjseqn{B} is `param_INFILT_xaj_B`
//' 
//' ![](xaj_infilt.png)
//' 
//' @param param_INFILT_xaj_B <0.01, 3> parameters for [infilt_XAJ()]
//' @export
// [[Rcpp::export]]
arma::vec infilt_XAJ(
    const arma::vec& LAND_water_mm, 
    const arma::vec& SOIL_water_mm,
    const arma::vec& SOIL_capacity_mm, 
    const arma::vec& param_INFILT_xaj_B
)
{
  arma::vec SOIL_diff_mm = SOIL_capacity_mm - SOIL_water_mm;
  arma::vec limit_mm = arma::min(arma::join_cols(SOIL_diff_mm, LAND_water_mm), 1);
  
  arma::vec MM_ = SOIL_capacity_mm * (param_INFILT_xaj_B + 1);
  arma::vec B_p_1 = (param_INFILT_xaj_B + 1);
  arma::vec B_B_1 = param_INFILT_xaj_B / B_p_1;
  arma::vec B_1 = 1 / param_INFILT_xaj_B;
  
  arma::vec AU_ = MM_ % (1 - pow(1 - SOIL_water_mm * B_p_1 / MM_, B_1));
  
  arma::vec AU_L_MM = (MM_ - AU_ - LAND_water_mm) / MM_;
  AU_L_MM = arma::max(AU_L_MM, arma::zeros<arma::vec>(AU_L_MM.n_elem));
  
  arma::vec MM_AU = (MM_ - AU_) / MM_;
  
  arma::vec infilt_water_mm = - MM_ % (pow(AU_L_MM, B_p_1) - pow(MM_AU, B_p_1)) / B_p_1;
  
  return arma::min(infilt_water_mm, limit_mm);
}

//' @rdname infilt
//' @details
//' # **_VIC** \insertCite{VIC_Wood_1992}{HydroGallery}:
//'
//' 
//' \mjsdeqn{F_{infilt} = \int_{i_{0}}^{i_{0}+P} A(i) {\rm d} i}
//' \mjsdeqn{i = C_{soil}(B+1) \left[ 1 - (1-A)^{1/B} \right]}
//' where
//'   - \mjseqn{B} is `param_INFILT_vic_B`
//' @param param_INFILT_vic_B <0.01, 3> parameters for [infilt_VIC()]
//' @export
// [[Rcpp::export]]
arma::vec infilt_VIC(
    const arma::vec& LAND_water_mm,
    const arma::vec& SOIL_water_mm,
    const arma::vec& SOIL_capacity_mm,
    const arma::vec& param_INFILT_vic_B)
{
    // Calculate intermediate variables
    arma::vec B_p_1 = param_INFILT_vic_B + 1;
    arma::vec B_1 = 1 / B_p_1;
    arma::vec i_m = SOIL_capacity_mm % B_p_1;
    
    // Calculate i_0
    arma::vec soil_ratio = 1 - SOIL_water_mm / SOIL_capacity_mm;
    arma::vec i_0 = i_m % (1 - arma::pow(soil_ratio, B_1));
    
    // Calculate infiltration
    arma::vec SOIL_diff_mm = SOIL_capacity_mm - SOIL_water_mm;
    
    // Create condition vector (uvec for find)
    arma::uvec condition = (i_0 + LAND_water_mm) > i_m;
    
    // Initialize result vector
    arma::vec infilt_water_mm = arma::zeros<arma::vec>(SOIL_diff_mm.n_elem);
    
    // Apply conditions
    infilt_water_mm.elem(arma::find(condition)) = SOIL_diff_mm.elem(arma::find(condition));
    
    arma::uvec else_cond = arma::find(condition == false);
    infilt_water_mm.elem(else_cond) = SOIL_diff_mm.elem(else_cond) - 
        (SOIL_capacity_mm.elem(else_cond) % 
         arma::pow(1 - (i_0.elem(else_cond) + LAND_water_mm.elem(else_cond)) / 
                  i_m.elem(else_cond), B_p_1.elem(else_cond)));
    
    // Calculate limit
    arma::vec limit_mm = arma::min(SOIL_diff_mm, LAND_water_mm);
    
    // Return final result
    return arma::min(infilt_water_mm, limit_mm);
}

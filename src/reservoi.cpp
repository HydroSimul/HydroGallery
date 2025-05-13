#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::interfaces(r, cpp)]]

//' **Reservoir release**
//' @name reservoi
//' @inheritParams all_vari
//' @description
//' The concept of river estimates the waterbody outflow for waternet concentration
//' @return outflow (m3)
//' @export
// [[Rcpp::export]]
arma::vec reservoireleas_Hanasaki(
    arma::vec Reservoi_water_m3,
    const arma::vec& Reservoi_inflow_m3,
    const arma::vec& Reservoi_demand_m3,
    const arma::vec& Reservoi_capacity_m3,
    const arma::vec& Reservoi_meanInflow_m3,
    const arma::vec& Reservoi_meanDemand_m3,
    const arma::vec& Reservoi_releaseCoefficient_1,
    const arma::uvec& Reservoi_isIrrigate_01)
{
    // Update reservoir water with inflow
    Reservoi_water_m3 += Reservoi_inflow_m3;
    
    // Initialize provisional release (eq-4)
    arma::vec Reservoi_releaseProvis_m3 = Reservoi_meanInflow_m3;
    
    // Apply irrigation conditions (eq-5)
    arma::uvec irrigate_idx = arma::find(Reservoi_isIrrigate_01 == 1);
    arma::uvec special_case_idx = arma::find(Reservoi_isIrrigate_01 == 1 && 
                                           Reservoi_meanDemand_m3 < 0.5 * Reservoi_meanInflow_m3);
    
    // Regular irrigation case
    Reservoi_releaseProvis_m3.elem(irrigate_idx) = 
        0.5 * Reservoi_meanInflow_m3.elem(irrigate_idx) % 
        (1.0 + Reservoi_demand_m3.elem(irrigate_idx) / 
         Reservoi_meanDemand_m3.elem(irrigate_idx));
    
    // Special irrigation case
    Reservoi_releaseProvis_m3.elem(special_case_idx) = 
        Reservoi_meanInflow_m3.elem(special_case_idx) + 
        Reservoi_demand_m3.elem(special_case_idx) - 
        Reservoi_meanDemand_m3.elem(special_case_idx);
    
    // Calculate capacity inflow ratio (eq-7)
    arma::vec Reservoi_ratioCapacityInflow_1 = Reservoi_capacity_m3 / Reservoi_meanInflow_m3;
    arma::vec temp_inflowRatio_ = 4.0 * Reservoi_ratioCapacityInflow_1 % Reservoi_ratioCapacityInflow_1;
    
    // Determine final release (eq-7)
    arma::uvec large_ratio_idx = arma::find(Reservoi_ratioCapacityInflow_1 > 0.5);
    arma::vec release = Reservoi_releaseCoefficient_1 % Reservoi_releaseProvis_m3;
    
    arma::uvec small_ratio_idx = arma::find(Reservoi_ratioCapacityInflow_1 <= 0.5);
    release.elem(small_ratio_idx) = 
        temp_inflowRatio_.elem(small_ratio_idx) % 
        Reservoi_releaseCoefficient_1.elem(small_ratio_idx) % 
        Reservoi_releaseProvis_m3.elem(small_ratio_idx) + 
        (1.0 - temp_inflowRatio_.elem(small_ratio_idx)) % 
        Reservoi_inflow_m3.elem(small_ratio_idx);
    
    return release;
}

//' @rdname reservoi
//' @param param_Reservoi_han_alpha <0,1> 0.85 parameter for [reservoireleasCoefficent_Hanasaki()],
//' @return new Reservoi_releaseCoefficient_1
//' @export
// [[Rcpp::export]]
arma::vec reservoiReleasCoefficent_Hanasaki(
    const arma::vec& Reservoi_water_m3,
    const arma::vec& Reservoi_capacity_m3,
    const arma::vec& Reservoi_releaseCoefficient_1,
    const arma::uvec& Reservoi_isOperateStart_01,
    const arma::vec& param_Reservoi_han_alpha)
{
    // Calculate new coefficient (eq-3)
    arma::vec Reservoi_releaseCoefficient_1_NEW = 
        Reservoi_water_m3 / (param_Reservoi_han_alpha % Reservoi_capacity_m3);
    
    // Apply minimum threshold
    Reservoi_releaseCoefficient_1_NEW = arma::clamp(Reservoi_releaseCoefficient_1_NEW, 0.1, arma::datum::inf);
    
    // Update coefficients where operation starts (eq-2)
    arma::uvec operate_idx = arma::find(Reservoi_isOperateStart_01 == 1);
    arma::vec result = Reservoi_releaseCoefficient_1;
    result.elem(operate_idx) = Reservoi_releaseCoefficient_1_NEW.elem(operate_idx);
    
    return result;
}
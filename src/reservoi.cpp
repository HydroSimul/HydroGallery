#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]

//' **Reservoir releas**
//' @name reservoi
//' @inheritParams all_vari
//' @description
//'
//' The concept of river estimates the waterbody outflow for waternet concentation
//' @return outflow (m3)
//' @export
// [[Rcpp::export]]
NumericVector reservoireleas_Hanasaki(
    NumericVector Reservoi_water_m3,
    NumericVector Reservoi_inflow_m3,
    NumericVector Reservoi_demand_m3,
    NumericVector Reservoi_capacity_m3,
    NumericVector Reservoi_meanInflow_m3,
    NumericVector Reservoi_meanDemand_m3,
    NumericVector Reservoi_releaseCoefficient_1,
    LogicalVector Reservoi_isIrrigate_01
)
{
  Reservoi_water_m3 += Reservoi_inflow_m3;


  NumericVector Reservoi_releaseProvis_m3 = Reservoi_meanInflow_m3; // eq-4
  Reservoi_releaseProvis_m3 = ifelse(Reservoi_isIrrigate_01,
                                      0.5 * Reservoi_meanInflow_m3 * (1 + Reservoi_demand_m3 / Reservoi_meanDemand_m3),
                                      Reservoi_releaseProvis_m3);  // eq-5
  Reservoi_releaseProvis_m3 = ifelse(Reservoi_isIrrigate_01 & (Reservoi_meanDemand_m3 < 0.5 * Reservoi_meanInflow_m3),
                                      Reservoi_meanInflow_m3 + Reservoi_demand_m3 - Reservoi_meanDemand_m3,
                                      Reservoi_releaseProvis_m3);  // eq-5
  NumericVector Reservoi_ratioCapacityInflow_1 = Reservoi_capacity_m3 / Reservoi_meanInflow_m3; // eq-7
  NumericVector temp_inflowRatio_ = 4 * Reservoi_ratioCapacityInflow_1 * Reservoi_ratioCapacityInflow_1; // eq-7

  return ifelse(Reservoi_ratioCapacityInflow_1 > 0.5,
                Reservoi_releaseCoefficient_1 * Reservoi_releaseProvis_m3,
                temp_inflowRatio_ * Reservoi_releaseCoefficient_1 * Reservoi_releaseProvis_m3 + (1 - temp_inflowRatio_) * Reservoi_inflow_m3); //  eq-7
}

//' @rdname reservoi
//' @param param_Reservoi_han_alpha <0,1> 0.85 parameter for [reservoireleasCoefficent_Hanasaki()],
//' @return new Reservoi_releaseCoefficient_1
//' @export
// [[Rcpp::export]]
NumericVector reservoiReleasCoefficent_Hanasaki(
   NumericVector Reservoi_water_m3,
   NumericVector Reservoi_capacity_m3,
   NumericVector Reservoi_releaseCoefficient_1,
   LogicalVector Reservoi_isOperateStart_01,
   NumericVector param_Reservoi_han_alpha
)
{

 NumericVector Reservoi_releaseCoefficient_1_NEW = Reservoi_water_m3 / (param_Reservoi_han_alpha * Reservoi_capacity_m3); // eq-3 // WG3 0.85
 Reservoi_releaseCoefficient_1_NEW = pmax(Reservoi_releaseCoefficient_1_NEW, 0.1);
 return ifelse(Reservoi_isOperateStart_01, Reservoi_releaseCoefficient_1_NEW, Reservoi_releaseCoefficient_1); // eq-2
 
}











#ifndef CONFLUEN
#define CONFLUEN

#include <RcppArmadillo.h>

arma::vec riverout_LinearResorvoir(
    const arma::vec& RIVER_water_m3,
    const arma::vec& RIVER_inflow_m3,
    const arma::vec& RIVER_velocity_km,
    const arma::vec& RIVER_length_km);

arma::vec riverlakout_LinearResorvoir(
    const arma::vec& Riverlak_water_m3,
    const arma::vec& Riverlak_inflow_m3,
    const arma::vec& Riverlak_capacity_m3,
    const arma::vec& param_Riverlak_lin_storeFactor);

arma::vec reservoireleas_Hanasaki(
    arma::vec Reservoi_water_m3,
    const arma::vec& Reservoi_inflow_m3,
    const arma::vec& Reservoi_demand_m3,
    const arma::vec& Reservoi_capacity_m3,
    const arma::vec& Reservoi_meanInflow_m3,
    const arma::vec& Reservoi_meanDemand_m3,
    const arma::vec& Reservoi_releaseCoefficient_1,
    const arma::uvec& Reservoi_isIrrigate_01);

#endif

#include <RcppArmadillo.h>
#include "confluen.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

arma::uvec find_in(const arma::uvec& x, const arma::uvec& table) {
  arma::uvec result(x.n_elem);
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    arma::uvec found = arma::find(table == x(i));
    result(i) = found.empty() ? 0 : found(0);  // 0 = not found
  }
  return result;
}

arma::vec inflow_add(const arma::vec& num_Outflow_LastStep, const arma::imat& int_InflowCell) {
  arma::vec num_Inflow_m3(int_InflowCell.n_rows, arma::fill::zeros);
  
  for (arma::uword i = 0; i < int_InflowCell.n_rows; ++i) {
    arma::ivec row = int_InflowCell.row(i).t();
    arma::uvec valid_idx = arma::find(row > 0);
    
    if (!valid_idx.empty()) {
      arma::uvec outflow_idx = arma::conv_to<arma::uvec>::from(row(valid_idx)) - 1;
      num_Inflow_m3(i) = arma::sum(num_Outflow_LastStep.elem(outflow_idx));
    }
  }
  return num_Inflow_m3;
}

//' confluen with sequential routing map
//' @name confluenSequential
//' @export
// [[Rcpp::export]]
arma::vec confluen_WaterGAP3_H(
    arma::vec& RIVER_water_m3,
    const arma::vec& RIVER_length_km,
    const arma::vec& RIVER_velocity_km,
    const arma::vec& RIVER_inflow_m3,
    const arma::field<arma::uvec>& CELL_cellNumberStep_int,
    const arma::field<arma::imat>& CELL_inflowCellNumberStep_int)
{
  arma::vec RIVER_water_m3_TEMP = RIVER_water_m3;
  arma::vec RIVER_outflow_m3(RIVER_inflow_m3.n_elem, arma::fill::zeros);
  
  arma::uvec idx_Cell_Step = CELL_cellNumberStep_int(0) - 1;
  
  arma::vec step_RiverOutflow_m3 = riverout_LinearResorvoir(
    RIVER_water_m3_TEMP.elem(idx_Cell_Step),
    RIVER_inflow_m3.elem(idx_Cell_Step),
    RIVER_velocity_km.elem(idx_Cell_Step),
    RIVER_length_km.elem(idx_Cell_Step)
  );
  
  arma::vec step_RIVER_Water_New = RIVER_water_m3_TEMP.elem(idx_Cell_Step)
    + RIVER_inflow_m3.elem(idx_Cell_Step)
    - step_RiverOutflow_m3;
    
    RIVER_outflow_m3.elem(idx_Cell_Step) = step_RiverOutflow_m3;
    RIVER_water_m3_TEMP.elem(idx_Cell_Step) = step_RIVER_Water_New;
    
    for (arma::uword i_Step = 1; i_Step < CELL_cellNumberStep_int.n_elem; ++i_Step) {
      idx_Cell_Step = CELL_cellNumberStep_int(i_Step) - 1;
      
      arma::vec step_UpstreamInflow_m3 = inflow_add(
        RIVER_outflow_m3,
        CELL_inflowCellNumberStep_int(i_Step)
      );
      step_UpstreamInflow_m3 += RIVER_inflow_m3.elem(idx_Cell_Step);
      
      step_RiverOutflow_m3 = riverout_LinearResorvoir(
        RIVER_water_m3_TEMP.elem(idx_Cell_Step),
        step_UpstreamInflow_m3,
        RIVER_velocity_km.elem(idx_Cell_Step),
        RIVER_length_km.elem(idx_Cell_Step)
      );
      
      step_RIVER_Water_New = RIVER_water_m3_TEMP.elem(idx_Cell_Step)
        + step_UpstreamInflow_m3
      - step_RiverOutflow_m3;
      
      RIVER_outflow_m3.elem(idx_Cell_Step) = step_RiverOutflow_m3;
      RIVER_water_m3_TEMP.elem(idx_Cell_Step) = step_RIVER_Water_New;
    }
    
    RIVER_water_m3 = RIVER_water_m3_TEMP;
    return RIVER_outflow_m3;
}



//' @rdname confluenSequential
//' @export
// [[Rcpp::export]]
arma::vec confluen_WaterGAP3_N(
    arma::vec& RIVER_water_m3,
    const arma::vec& RIVER_length_km,
    const arma::vec& RIVER_velocity_km,
    const arma::vec& RIVER_inflow_m3,
    const arma::field<arma::uvec>& CELL_cellNumberStep_int,
    const arma::field<arma::imat>& CELL_inflowCellNumberStep_int,
    const arma::uvec& Riverlak_cellNumber_int,
    const arma::vec& Riverlak_capacity_m3,
    const arma::vec& param_Riverlak_lin_storeFactor)
{
  arma::vec RIVER_water_m3_TEMP = RIVER_water_m3;
  arma::vec RIVER_outflow_m3 = arma::vec(RIVER_inflow_m3.n_elem, arma::fill::zeros);
  
  // Riverlak storage
  arma::vec Riverlak_water_m3 = RIVER_water_m3_TEMP.elem(Riverlak_cellNumber_int - 1);
  
  int n_Step = CELL_cellNumberStep_int.n_elem;
  
  // Process first step
  arma::uvec idx_Cell_Step = CELL_cellNumberStep_int(0) - 1; // Convert to 0-based
  
  // River segment calculations
  arma::vec step_RiverOutflow_m3 = riverout_LinearResorvoir(
    RIVER_water_m3_TEMP.elem(idx_Cell_Step),
    RIVER_inflow_m3.elem(idx_Cell_Step),
    RIVER_velocity_km.elem(idx_Cell_Step),
    RIVER_length_km.elem(idx_Cell_Step)
  );
  
  arma::vec step_RIVER_Water_New = RIVER_water_m3_TEMP.elem(idx_Cell_Step) + 
    RIVER_inflow_m3.elem(idx_Cell_Step) - 
    step_RiverOutflow_m3;
  
  // Riverlak processing
  arma::uvec idx_Riverlak_Step = arma::find(arma::find(Riverlak_cellNumber_int == arma::conv_to<arma::uvec>::from(idx_Cell_Step + 1)));
  if (idx_Riverlak_Step.n_elem > 0) {
    arma::uvec idx_Step_Riverlak = arma::find(arma::find(idx_Cell_Step + 1 == Riverlak_cellNumber_int));
    
    arma::vec step_RiverlakOutflow_m3 = riverlakout_LinearResorvoir(
      Riverlak_water_m3.elem(idx_Riverlak_Step),
      RIVER_inflow_m3.elem(idx_Cell_Step.elem(idx_Step_Riverlak)),
      Riverlak_capacity_m3.elem(idx_Riverlak_Step),
      param_Riverlak_lin_storeFactor.elem(idx_Riverlak_Step)
    );
    
    step_RiverOutflow_m3.elem(idx_Step_Riverlak) = step_RiverlakOutflow_m3;
    step_RIVER_Water_New.elem(idx_Step_Riverlak) = 
      Riverlak_water_m3.elem(idx_Riverlak_Step) + 
      RIVER_inflow_m3.elem(idx_Cell_Step.elem(idx_Step_Riverlak)) - 
      step_RiverlakOutflow_m3;
  }
  
  RIVER_outflow_m3.elem(idx_Cell_Step) = step_RiverOutflow_m3;
  RIVER_water_m3_TEMP.elem(idx_Cell_Step) = step_RIVER_Water_New;
  
  // Process remaining steps
  for (int i_Step = 1; i_Step < n_Step; ++i_Step) {
    idx_Cell_Step = CELL_cellNumberStep_int(i_Step) - 1;
    
    // Calculate upstream inflow
    arma::vec step_UpstreamInflow_m3 = inflow_add(
      RIVER_outflow_m3,
      CELL_inflowCellNumberStep_int(i_Step)
    );
    step_UpstreamInflow_m3 += RIVER_inflow_m3.elem(idx_Cell_Step);
    
    // River segment calculations
    step_RiverOutflow_m3 = riverout_LinearResorvoir(
      RIVER_water_m3_TEMP.elem(idx_Cell_Step),
      step_UpstreamInflow_m3,
      RIVER_velocity_km.elem(idx_Cell_Step),
      RIVER_length_km.elem(idx_Cell_Step)
    );
    
    step_RIVER_Water_New = RIVER_water_m3_TEMP.elem(idx_Cell_Step) + 
      step_UpstreamInflow_m3 - 
      step_RiverOutflow_m3;
    
    // Riverlak processing
    idx_Riverlak_Step = arma::find(arma::find(Riverlak_cellNumber_int == arma::conv_to<arma::uvec>::from(idx_Cell_Step + 1)));
    if (idx_Riverlak_Step.n_elem > 0) {
      arma::uvec idx_Step_Riverlak = arma::find(arma::find(idx_Cell_Step + 1 == Riverlak_cellNumber_int));
      
      arma::vec step_RiverlakOutflow_m3 = riverlakout_LinearResorvoir(
        Riverlak_water_m3.elem(idx_Riverlak_Step),
        step_UpstreamInflow_m3.elem(idx_Step_Riverlak),
        Riverlak_capacity_m3.elem(idx_Riverlak_Step),
        param_Riverlak_lin_storeFactor.elem(idx_Riverlak_Step)
      );
      
      step_RiverOutflow_m3.elem(idx_Step_Riverlak) = step_RiverlakOutflow_m3;
      step_RIVER_Water_New.elem(idx_Step_Riverlak) = 
        Riverlak_water_m3.elem(idx_Riverlak_Step) + 
        step_UpstreamInflow_m3.elem(idx_Step_Riverlak) - 
        step_RiverlakOutflow_m3;
    }
    
    RIVER_outflow_m3.elem(idx_Cell_Step) = step_RiverOutflow_m3;
    RIVER_water_m3_TEMP.elem(idx_Cell_Step) = step_RIVER_Water_New;
  }
  
  RIVER_water_m3 = RIVER_water_m3_TEMP;
  return RIVER_outflow_m3;
}


//' @rdname confluenSequential
//' @export
// [[Rcpp::export]]
arma::vec confluen_WaterGAP3_U(
    arma::vec& RIVER_water_m3,
    const arma::vec& RIVER_length_km,
    const arma::vec& RIVER_velocity_km,
    const arma::vec& RIVER_inflow_m3,
    const arma::field<arma::uvec>& CELL_cellNumberStep_int,
    const arma::field<arma::imat>& CELL_inflowCellNumberStep_int,
    const arma::uvec& Riverlak_cellNumber_int,
    const arma::vec& Riverlak_capacity_m3,
    const arma::uvec& Reservoi_cellNumber_int,
    const arma::vec& Reservoi_demand_m3,
    const arma::vec& Reservoi_capacity_m3,
    const arma::vec& Reservoi_meanInflow_m3,
    const arma::vec& Reservoi_meanDemand_m3,
    const arma::vec& Reservoi_releaseCoefficient_1,
    const arma::uvec& Reservoi_isIrrigate_01,
    const arma::vec& param_Riverlak_lin_storeFactor)
{
  arma::vec RIVER_water_m3_TEMP = RIVER_water_m3;
  arma::vec RIVER_outflow_m3(RIVER_inflow_m3.n_elem, arma::fill::zeros);
  
  arma::vec Riverlak_water_m3 = RIVER_water_m3_TEMP.elem(Riverlak_cellNumber_int - 1);
  arma::vec Reservoi_water_m3 = RIVER_water_m3_TEMP.elem(Reservoi_cellNumber_int - 1);
  
  int n_Step = CELL_cellNumberStep_int.n_elem;
  
  arma::uvec idx_Cell_Step = CELL_cellNumberStep_int(0) - 1;
  
  arma::vec step_RiverOutflow_m3 = riverout_LinearResorvoir(
    RIVER_water_m3_TEMP.elem(idx_Cell_Step),
    RIVER_inflow_m3.elem(idx_Cell_Step),
    RIVER_velocity_km.elem(idx_Cell_Step),
    RIVER_length_km.elem(idx_Cell_Step)
  );
  
  arma::vec step_RIVER_Water_New = RIVER_water_m3_TEMP.elem(idx_Cell_Step) + 
    RIVER_inflow_m3.elem(idx_Cell_Step) - 
    step_RiverOutflow_m3;
  
  arma::uvec idx_Riverlak_Step = find_in(Riverlak_cellNumber_int, idx_Cell_Step + 1);
  if (idx_Riverlak_Step.n_elem > 0) {
    arma::uvec idx_Step_Riverlak = find_in(idx_Cell_Step + 1, Riverlak_cellNumber_int);
    
    arma::vec step_RiverlakOutflow_m3 = riverlakout_LinearResorvoir(
      Riverlak_water_m3.elem(idx_Riverlak_Step),
      RIVER_inflow_m3.elem(idx_Cell_Step.elem(idx_Step_Riverlak)),
      Riverlak_capacity_m3.elem(idx_Riverlak_Step),
      param_Riverlak_lin_storeFactor.elem(idx_Riverlak_Step)
    );
    
    step_RiverOutflow_m3.elem(idx_Step_Riverlak) = step_RiverlakOutflow_m3;
    step_RIVER_Water_New.elem(idx_Step_Riverlak) = 
      Riverlak_water_m3.elem(idx_Riverlak_Step) + 
      RIVER_inflow_m3.elem(idx_Cell_Step.elem(idx_Step_Riverlak)) - 
      step_RiverlakOutflow_m3;
  }
  
  arma::uvec idx_Reservoi_Step = find_in(Reservoi_cellNumber_int, idx_Cell_Step + 1);
  if (idx_Reservoi_Step.n_elem > 0) {
    arma::uvec idx_Step_Reservoi = find_in(idx_Cell_Step + 1, Reservoi_cellNumber_int);
    
    arma::vec step_ReservoiOutflow_m3 = reservoireleas_Hanasaki(
      Reservoi_water_m3.elem(idx_Reservoi_Step),
      RIVER_inflow_m3.elem(idx_Cell_Step.elem(idx_Step_Reservoi)),
      Reservoi_demand_m3.elem(idx_Reservoi_Step),
      Reservoi_capacity_m3.elem(idx_Reservoi_Step),
      Reservoi_meanInflow_m3.elem(idx_Reservoi_Step),
      Reservoi_meanDemand_m3.elem(idx_Reservoi_Step),
      Reservoi_releaseCoefficient_1.elem(idx_Reservoi_Step),
      Reservoi_isIrrigate_01.elem(idx_Reservoi_Step)
    );
    
    step_RiverOutflow_m3.elem(idx_Step_Reservoi) = step_ReservoiOutflow_m3;
    step_RIVER_Water_New.elem(idx_Step_Reservoi) = 
      Reservoi_water_m3.elem(idx_Reservoi_Step) + 
      RIVER_inflow_m3.elem(idx_Cell_Step.elem(idx_Step_Reservoi)) - 
      step_ReservoiOutflow_m3;
  }
  
  RIVER_outflow_m3.elem(idx_Cell_Step) = step_RiverOutflow_m3;
  RIVER_water_m3_TEMP.elem(idx_Cell_Step) = step_RIVER_Water_New;
  
  for (int i_Step = 1; i_Step < n_Step; ++i_Step) {
    idx_Cell_Step = CELL_cellNumberStep_int(i_Step) - 1;
    
    arma::vec step_UpstreamInflow_m3 = inflow_add(
      RIVER_outflow_m3,
      CELL_inflowCellNumberStep_int(i_Step)
    );
    step_UpstreamInflow_m3 += RIVER_inflow_m3.elem(idx_Cell_Step);
    
    step_RiverOutflow_m3 = riverout_LinearResorvoir(
      RIVER_water_m3_TEMP.elem(idx_Cell_Step),
      step_UpstreamInflow_m3,
      RIVER_velocity_km.elem(idx_Cell_Step),
      RIVER_length_km.elem(idx_Cell_Step)
    );
    
    step_RIVER_Water_New = RIVER_water_m3_TEMP.elem(idx_Cell_Step) + 
      step_UpstreamInflow_m3 - 
      step_RiverOutflow_m3;
    
    idx_Riverlak_Step = find_in(Riverlak_cellNumber_int, idx_Cell_Step + 1);
    if (idx_Riverlak_Step.n_elem > 0) {
      arma::uvec idx_Step_Riverlak = find_in(idx_Cell_Step + 1, Riverlak_cellNumber_int);
      
      arma::vec step_RiverlakOutflow_m3 = riverlakout_LinearResorvoir(
        Riverlak_water_m3.elem(idx_Riverlak_Step),
        step_UpstreamInflow_m3.elem(idx_Step_Riverlak),
        Riverlak_capacity_m3.elem(idx_Riverlak_Step),
        param_Riverlak_lin_storeFactor.elem(idx_Riverlak_Step)
      );
      
      step_RiverOutflow_m3.elem(idx_Step_Riverlak) = step_RiverlakOutflow_m3;
      step_RIVER_Water_New.elem(idx_Step_Riverlak) = 
        Riverlak_water_m3.elem(idx_Riverlak_Step) + 
        step_UpstreamInflow_m3.elem(idx_Step_Riverlak) - 
        step_RiverlakOutflow_m3;
    }
    
    idx_Reservoi_Step = find_in(Reservoi_cellNumber_int, idx_Cell_Step + 1);
    if (idx_Reservoi_Step.n_elem > 0) {
      arma::uvec idx_Step_Reservoi = find_in(idx_Cell_Step + 1, Reservoi_cellNumber_int);
      
      arma::vec step_ReservoiOutflow_m3 = reservoireleas_Hanasaki(
        Reservoi_water_m3.elem(idx_Reservoi_Step),
        step_UpstreamInflow_m3.elem(idx_Step_Reservoi),
        Reservoi_demand_m3.elem(idx_Reservoi_Step),
        Reservoi_capacity_m3.elem(idx_Reservoi_Step),
        Reservoi_meanInflow_m3.elem(idx_Reservoi_Step),
        Reservoi_meanDemand_m3.elem(idx_Reservoi_Step),
        Reservoi_releaseCoefficient_1.elem(idx_Reservoi_Step),
        Reservoi_isIrrigate_01.elem(idx_Reservoi_Step)
      );
      
      step_RiverOutflow_m3.elem(idx_Step_Reservoi) = step_ReservoiOutflow_m3;
      step_RIVER_Water_New.elem(idx_Step_Reservoi) = 
        Reservoi_water_m3.elem(idx_Reservoi_Step) + 
        step_UpstreamInflow_m3.elem(idx_Step_Reservoi) - 
        step_ReservoiOutflow_m3;
    }
    
    RIVER_outflow_m3.elem(idx_Cell_Step) = step_RiverOutflow_m3;
    RIVER_water_m3_TEMP.elem(idx_Cell_Step) = step_RIVER_Water_New;
  }
  
  RIVER_water_m3 = RIVER_water_m3_TEMP;
  return RIVER_outflow_m3;
}

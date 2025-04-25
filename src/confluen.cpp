#include "confluen.h"
// [[Rcpp::interfaces(r, cpp)]]


//' **confluence**
//' @description 
//' \loadmathjax
//' 
//' In hydrological modeling, routing (named as [confluen] in HydroGallery) refers to the process of simulating the movement of water through a river network or other drainage system. 
//' It allows the model to predict the flow of water in rivers and streams. 
//' In hydrological models, routing is typically performed using mathematical algorithms that account for the physical properties of the river network, 
//' such as its geometry, roughness, and discharge capacity. 
//' The parameters that govern routing, such as flow velocity and channel roughness, 
//' can have a significant impact on the accuracy of the model.
//' 
//' `confluence` is a calculation function that causes water to flow into the gauge point.
//' - `IUH`: IUH (Instant Unit Hydrograph) with one watercourse, 
//' - `IUH2S`: IUH with two water sources, each with a different IUH vector, 
//' - `IUH3S`: IUH with three water sources, each with a different IUH vector.
//' 
//' Under the concept of the conceptual HM, the water flux to the water flow will be calculated using the confluence process. 
//' This process does not calculate the water balance, but rather the time-varying nature of the water flow. 
//' The "Instant Unit Hydrograph" method is the most effective way to deal with time-varying flows. 
//' In the first stage, only [confluenIUH] will be supported.
//' 
//' So we can give the function:
//' 
//' \mjsdeqn{Q = f_{confluen}(F, u)}
//' 
//' 
//' 
//' where
//' - \mjseqn{Q} is stream flow, but still in mm/TS not m3/TS or m3/S
//' - \mjseqn{F} is flux that will into river conflen, e.g.`LAND_runoff_mm`, `SOIL_interflow_mm` or `GROUND_baseflow_mm`
//' - \mjseqn{u} is Instant Unit Hydrograph series
//' 
//' @references
//' \insertAllCited{}
//' @name confluen
//' @inheritParams all_vari
//' @return confluenced water (mm/m2)
//' @export
// [[Rcpp::export]]
NumericVector confluen_IUH(
    NumericVector CONFLUEN_inputWater_mm, 
    NumericVector CONFLUEN_iuh_1
)
{
  
  int n_iuh = CONFLUEN_iuh_1.size(), n_time = CONFLUEN_inputWater_mm.size();
  NumericVector CONFLUEN_outputWater_mm (n_time);
  for (int i = 0; i < n_iuh; i++) {
    for (int j = 0; j <= i; j++) {
      CONFLUEN_outputWater_mm[i] += CONFLUEN_inputWater_mm[i-j] * CONFLUEN_iuh_1[j];
    }
  }
  for (int i = n_iuh; i < n_time; i++) {
    for (int j = 0; j < n_iuh; j++) {
      CONFLUEN_outputWater_mm[i] += CONFLUEN_inputWater_mm[i-j] * CONFLUEN_iuh_1[j];
    }
  }
  
  return CONFLUEN_outputWater_mm;
  
}

//' @rdname confluen
//' @export
// [[Rcpp::export]]
NumericVector confluen_IUH2S(
    NumericVector LAND_runoff_mm,
    NumericVector GROUND_baseflow_mm, 
    NumericVector CONFLUEN_iuhLand_1,
    NumericVector CONFLUEN_iuhGround_1
)
{
  NumericVector CONFLUEN_runoff_mm (LAND_runoff_mm.size()), CONFLUEN_baseflow_mm (GROUND_baseflow_mm.size());
  CONFLUEN_runoff_mm = confluen_IUH(
    LAND_runoff_mm, 
    CONFLUEN_iuhLand_1
  );
  CONFLUEN_baseflow_mm = confluen_IUH(
    GROUND_baseflow_mm, 
    CONFLUEN_iuhGround_1
  );
  

  return CONFLUEN_runoff_mm + CONFLUEN_baseflow_mm;
  
}

//' @rdname confluen
//' @export
// [[Rcpp::export]]
NumericVector confluen_IUH3S(
    NumericVector LAND_runoff_mm,
    NumericVector SOIL_interflow_mm, 
    NumericVector GROUND_baseflow_mm, 
    NumericVector CONFLUEN_iuhLand_1,
    NumericVector CONFLUEN_iuhSoil_1,
    NumericVector CONFLUEN_iuhGround_1
)
{
  NumericVector CONFLUEN_runoff_mm (LAND_runoff_mm.size()), CONFLUEN_interflow_mm (SOIL_interflow_mm.size()), CONFLUEN_baseflow_mm (GROUND_baseflow_mm.size());
  CONFLUEN_runoff_mm = confluen_IUH(
    LAND_runoff_mm, 
    CONFLUEN_iuhLand_1
  );
  CONFLUEN_interflow_mm = confluen_IUH(
    SOIL_interflow_mm, 
    CONFLUEN_iuhSoil_1
  );
  CONFLUEN_baseflow_mm = confluen_IUH(
    GROUND_baseflow_mm, 
    CONFLUEN_iuhGround_1
  );
  
  
  return CONFLUEN_runoff_mm + CONFLUEN_interflow_mm + CONFLUEN_baseflow_mm;
  
}


//' create **IUH** (Instant Unit Hydrograph)
//' @name confluenIUH
//' @inheritParams all_vari
//' @description
//' \loadmathjax
//' 
//' The process `confluenIUH` return a series of portions, that means how many flux water will
//' in those moment into the river.
//' The sum of this series will always in 1.
//' 
//' So we can give the function:
//' 
//' \mjsdeqn{u = f_{confluenIUH}(t_r, ...)}
//' 
//' 
//' 
//' where
//' - \mjseqn{u} is series of portions
//' - \mjseqn{t_r} is  `CONFLUEN_responseTime_TS`
//' 
//' @references
//' \insertAllCited{}
//' @return IUH (list of num vector) 
//' @details
//' # **_GR4J1** \insertCite{GR4J_Perrin_2003}{HydroGallery}: 
//'
//' 
//' \mjsdeqn{u(i) = S(i) - S(i-1)}
//' \mjsdeqn{S(i) = \left( \frac{i}{t_r} \right)^{2.5}, \quad 0 \leq i \leq t_r}
//' where
//'   - \mjseqn{u} is IUH series
//'   - \mjseqn{i} is index
//' @export
// [[Rcpp::export]]
NumericVector confluenIUH_GR4J1(
    double CONFLUEN_responseTime_TS
)
{
  double t_max = ceil(CONFLUEN_responseTime_TS);
  IntegerVector seq_t = seq(1, t_max);
  NumericVector seq_t2 = as<NumericVector>(seq_t);
  NumericVector SH_1 = pow(( seq_t2/ CONFLUEN_responseTime_TS), 2.5);
  SH_1(t_max - 1) = 1;
  SH_1[Range(1, t_max - 1)] = diff(SH_1);
  return SH_1;
}



//' @rdname confluenIUH
//' @details
//' # **_GR4J2** \insertCite{GR4J_Perrin_2003}{HydroGallery}: 
//'
//' 
//' \mjsdeqn{u(i) = S(i) - S(i-1)}
//' \mjsdeqn{S(i) = 0.5\left( \frac{i}{t_r} \right)^{2.5}, \quad 0 \leq i \leq t_r}
//' \mjsdeqn{S(i) = 1 - 0.5\left(2 - \frac{i}{t_r} \right)^{2.5}, \quad t_r < i < 2t_r}
//' \mjsdeqn{S(i) = 0, \quad i = 2t_r}
//' where
//'   - \mjseqn{u} is IUH series
//'   - \mjseqn{i} is index
//' @export
// [[Rcpp::export]]
NumericVector confluenIUH_GR4J2(
    double CONFLUEN_responseTime_TS
)
{
  double t_max_1 = ceil(CONFLUEN_responseTime_TS);
  double t_max_2 = ceil(2 * CONFLUEN_responseTime_TS);
  IntegerVector seq_t1 = seq(1, t_max_1 - 1);
  NumericVector seq_t12 = as<NumericVector>(seq_t1);
  IntegerVector seq_t2 = seq(t_max_1, (t_max_2 - 1));
  NumericVector seq_t22 = as<NumericVector>(seq_t2);
  
  NumericVector SH_2_1 = .5 * pow((seq_t12 / CONFLUEN_responseTime_TS),2.5);
  NumericVector SH_2_2 = 1 - .5 * pow((2 - seq_t22 / CONFLUEN_responseTime_TS),2.5);
  NumericVector SH_2(t_max_2, 1);
  SH_2[Range(0, t_max_1 - 2)] = SH_2_1;
  SH_2[Range(t_max_1 - 1, t_max_2 - 2)] = SH_2_2;
  SH_2[Range(1, t_max_2 - 1)] = diff(SH_2);
  
  return SH_2;
}


//' @rdname confluenIUH
//' @details
//' # **_Kelly** \insertCite{iuh_Kelly_1955}{HydroGallery}: 
//'
//' 
//' \mjsdeqn{u(i) = \frac{4}{t_r^2} \left( i + k \left( e^{-i/k} \right) \right), \quad i \leq t_r / 2 }
//' \mjsdeqn{u(i) = - \frac{4}{t_r^2}(i - k - t_r) + \frac{4ke^{-i/k}}{t_r^2} (1 - 2 e^{t_r/(2k)}), \quad t_r / 2 < i \leq t_r }
//' \mjsdeqn{u(i) =  \frac{4ke^{-i/k}}{t_r^2} (1 - 2 e^{t_r/(2k)} +  e^{t_r/k}), \quad i > t_r }
//' where
//'   - \mjseqn{k} is `param_CONFLUEN_kel_k`
//' @param param_CONFLUEN_kel_k <1, 4> parameter for[confluenIUH_Kelly()]
//' @export
// [[Rcpp::export]]
NumericVector confluenIUH_Kelly(
    double CONFLUEN_responseTime_TS,
    double param_CONFLUEN_kel_k
)
{
  double CONFLUEN_concentratTime_TS = CONFLUEN_responseTime_TS * param_CONFLUEN_kel_k;
  double num_temp_tc2 = (CONFLUEN_concentratTime_TS * CONFLUEN_concentratTime_TS);
  double num_temp_12_34 = 4 * CONFLUEN_responseTime_TS  / num_temp_tc2 * 
    (1 - 2 * exp(CONFLUEN_concentratTime_TS / CONFLUEN_responseTime_TS * 0.5));
  double num_temp_12_35 = 4 * CONFLUEN_responseTime_TS  / num_temp_tc2 * 
    (1 - 2 * exp(CONFLUEN_concentratTime_TS / CONFLUEN_responseTime_TS * 0.5) + exp(CONFLUEN_concentratTime_TS / CONFLUEN_responseTime_TS));
  double t_max = ceil(std::max(CONFLUEN_concentratTime_TS, - CONFLUEN_responseTime_TS * log(0.002 / num_temp_12_35)));
  NumericVector iuh_1, iuh_2, iuh_3, iuh_, vct_iuh, temp_etK;
  IntegerVector seq_t = seq(1, 20 * t_max);
  NumericVector seq_t2 = as<NumericVector>(seq_t) / 20.0;
  temp_etK = exp(- seq_t2 / CONFLUEN_responseTime_TS);
  iuh_1 = 4 / num_temp_tc2 * (seq_t2 + CONFLUEN_responseTime_TS * (temp_etK - 1));
  iuh_2 = num_temp_12_34 * temp_etK - 4 / num_temp_tc2 * (seq_t2 - CONFLUEN_responseTime_TS - CONFLUEN_concentratTime_TS);
  iuh_3 = num_temp_12_35 * temp_etK;
  iuh_ = ifelse(seq_t2 > CONFLUEN_concentratTime_TS * 0.5, iuh_2, iuh_1);
  iuh_ = ifelse(seq_t2 > CONFLUEN_concentratTime_TS, iuh_3, iuh_);
  NumericMatrix mat_iuh = NumericMatrix(20, t_max, iuh_.begin());
  vct_iuh = colMeans(mat_iuh);
  return vct_iuh / sum(vct_iuh);
}


//' @rdname confluenIUH
//' @details
//' # **_Nash** \insertCite{iuh_Nash_1957}{HydroGallery}: 
//'
//' 
//' \mjsdeqn{u(i) = \frac{1}{t_r\Gamma(n)} \left(\frac{4}{t_r^2}\right)^{n -1}e^{-i/t_r}}
//' where
//'   - \mjseqn{n} is `param_CONFLUEN_nas_n`
//' @param param_CONFLUEN_nas_n <1, 8> parameter for[confluenIUH_Nash()]
//' @export
// [[Rcpp::export]]
NumericVector confluenIUH_Nash(
    double CONFLUEN_responseTime_TS,
    double param_CONFLUEN_nas_n
)
{
  double t_max = ceil(std::max(4.0, param_CONFLUEN_nas_n) * 3 * CONFLUEN_responseTime_TS);
  NumericVector iuh_, vct_iuh;
  IntegerVector seq_t = seq(1, 20 * t_max);
  NumericVector seq_t2 = as<NumericVector>(seq_t) / 20.0;
  iuh_ = pow(seq_t2 / CONFLUEN_responseTime_TS, param_CONFLUEN_nas_n - 1) * exp(- seq_t2 / CONFLUEN_responseTime_TS) / 
    CONFLUEN_responseTime_TS / tgamma(param_CONFLUEN_nas_n);
  NumericMatrix mat_iuh = NumericMatrix(20, t_max, iuh_.begin());
  vct_iuh = colMeans(mat_iuh);
  return vct_iuh / sum(vct_iuh);
}

//' @rdname confluenIUH
//' @details
//' # **_Clark** \insertCite{iuh_Clark_1945}{HydroGallery}: 
//'
//' 
//' \mjsdeqn{u(i) = \frac{1}{t_r} e^{-i/t_r} }
//' where
//'   - \mjseqn{t_r} is `CONFLUEN_responseTime_TS`
//' @export
// [[Rcpp::export]]
NumericVector confluenIUH_Clark(
    double CONFLUEN_responseTime_TS
)
{
  double t_max = ceil(- CONFLUEN_responseTime_TS * log(CONFLUEN_responseTime_TS * 0.005));
  IntegerVector seq_t = seq(1, 20 * t_max);
  NumericVector seq_t2 = as<NumericVector>(seq_t) / 20.0;
  NumericVector iuh_ = 1 / CONFLUEN_responseTime_TS * exp(- seq_t2 / CONFLUEN_responseTime_TS);
  NumericMatrix mat_iuh = NumericMatrix(20, t_max, iuh_.begin());
  NumericVector vct_iuh = colMeans(mat_iuh);
  return vct_iuh / sum(vct_iuh);
}





NumericVector inflow_add(
    NumericVector num_Outflow_LastStep,
    IntegerMatrix int_InflowCell
)
{
  
  int n_Cell = int_InflowCell.nrow();  // Number of rows (cells)
  NumericVector num_Inflow_m3(n_Cell, 0.); // Initialize the result vector with NA
  
  for (int i_Cell = 0; i_Cell < n_Cell; i_Cell++) {
    double inflow_sum = 0.0;
    
    for (int i_InflowCell = 0; i_InflowCell < int_InflowCell.ncol(); i_InflowCell++) {
      int index = int_InflowCell(i_Cell, i_InflowCell);
      if (index == NA_INTEGER) break;
      inflow_sum += num_Outflow_LastStep[index - 1];
    }
    num_Inflow_m3[i_Cell] = inflow_sum;
  }
  
  return num_Inflow_m3;
}



IntegerVector get_idx_step(IntegerVector int_Cell, IntegerVector int_Step) {
  // Find the match between int_Step and int_Cell
  IntegerVector int_Match = match(int_Cell, int_Step);
  
  return na_omit(int_Match);
}

IntegerVector get_idx_cell(IntegerVector int_Cell, IntegerVector int_Step) {
  // Find the match between int_Step and int_Cell
  IntegerVector int_Match = match(int_Step, int_Cell);
  
  return na_omit(int_Match);
}


//' @rdname confluen
//' @export
// [[Rcpp::export]]
NumericVector confluen_WaterGAP3_H(
   NumericVector &RIVER_water_m3,
   NumericVector RIVER_length_km,
   NumericVector RIVER_velocity_km,
   NumericVector RIVER_outflow_m3,
   List CELL_cellNumberStep_int,
   List CELL_inflowCellNumberStep_int
)
{
 
 NumericVector step_RIVER_Water_New, step_RIVER_Outflow,
 step_RiverOutflow_m3;
 
 IntegerVector idx_Cell_Step,
 idx_Riverlak_Step, idx_Step_Riverlak;
 int n_Step = CELL_cellNumberStep_int.size();
 
 // Step i later with Inflow
 for (int i_Step = 1; i_Step < n_Step; i_Step++)
 {
   
   idx_Cell_Step = CELL_cellNumberStep_int[i_Step];
   NumericVector step_RiverWater = subset_get(RIVER_water_m3, idx_Cell_Step);
   NumericVector step_UpstreamInflow_m3;
   
   
   
   // Inflow upstream
   step_UpstreamInflow_m3 = inflow_add(
     RIVER_outflow_m3,
     CELL_inflowCellNumberStep_int[i_Step]
   );
   
   step_RiverOutflow_m3 = riverout_LinearResorvoir(
     step_RiverWater,
     step_UpstreamInflow_m3,
     subset_get(RIVER_velocity_km, idx_Cell_Step),
     subset_get(RIVER_length_km, idx_Cell_Step)
   );
   step_RIVER_Water_New = pmax(step_RiverWater + step_UpstreamInflow_m3 - step_RiverOutflow_m3, 0.0);
   step_RIVER_Outflow = step_RiverWater + step_UpstreamInflow_m3 - step_RIVER_Water_New;
   
   
   
   subset_put(RIVER_outflow_m3, idx_Cell_Step, step_RIVER_Outflow);
   subset_put(RIVER_water_m3, idx_Cell_Step,  step_RIVER_Water_New);
   
   
 }
 
 
 return RIVER_outflow_m3;
 
}

//' @rdname confluen
//' @export
// [[Rcpp::export]]
NumericVector confluen_WaterGAP3_L(
   NumericVector &RIVER_water_m3,
   NumericVector RIVER_length_km,
   NumericVector RIVER_velocity_km,
   NumericVector RIVER_outflow_m3,
   List CELL_cellNumberStep_int,
   List CELL_inflowCellNumberStep_int,
   IntegerVector Riverlak_cellNumber_int,
   NumericVector &Riverlak_water_m3,
   NumericVector Riverlak_capacity_m3,
   NumericVector param_Riverlak_lin_storeFactor
)
{
 
 NumericVector step_RiverOutflow_m3, step_RiverlakOutflow_m3, step_UpstreamInflow_m3;
 
 IntegerVector idx_Cell_Step,
 idx_Step_Upstream,
 idx_Riverlak_Step, idx_Step_Riverlak;
 int n_Step = CELL_cellNumberStep_int.size();
 
 // Overflow for Riverlak
 NumericVector Riverlak_overflow_m3 = pmax(Riverlak_water_m3 - Riverlak_capacity_m3, 0);
 Riverlak_water_m3 += -Riverlak_overflow_m3;
 
 
 // Step i later with Inflow
 for (int i_Step = 1; i_Step < n_Step; i_Step++)
 {
   
   idx_Cell_Step = CELL_cellNumberStep_int[i_Step];
   
   
   // Inflow upstream
   step_UpstreamInflow_m3 = inflow_add(
     RIVER_outflow_m3,
     CELL_inflowCellNumberStep_int[i_Step]
   );
   
   
   // river segment
   NumericVector step_RiverWater = subset_get(RIVER_water_m3, idx_Cell_Step);
   
   step_RiverOutflow_m3 = riverout_LinearResorvoir(
     step_RiverWater,
     step_UpstreamInflow_m3,
     subset_get(RIVER_velocity_km, idx_Cell_Step),
     subset_get(RIVER_length_km, idx_Cell_Step)
   );
   NumericVector step_RIVER_Water_New = pmax(step_RiverWater + step_UpstreamInflow_m3 - step_RiverOutflow_m3, 0.0);
   NumericVector step_RIVER_Outflow = step_RiverWater + step_UpstreamInflow_m3 - step_RIVER_Water_New;
   
   // Riverlak
   idx_Riverlak_Step = get_idx_cell(Riverlak_cellNumber_int, idx_Cell_Step);
   idx_Step_Riverlak = get_idx_step(Riverlak_cellNumber_int, idx_Cell_Step);
   if (idx_Riverlak_Step.size() > 0) {
     NumericVector step_RiverlakWater= subset_get(Riverlak_water_m3, idx_Riverlak_Step),
       step_RiverlakInflow = subset_get(step_UpstreamInflow_m3, idx_Step_Riverlak);
     
     
     step_RiverlakOutflow_m3 = riverlakout_LinearResorvoir(
       step_RiverlakWater,
       step_RiverlakInflow,
       subset_get(Riverlak_capacity_m3, idx_Riverlak_Step),
       subset_get(param_Riverlak_lin_storeFactor, idx_Riverlak_Step)
     );
     subset_put(step_RIVER_Outflow, idx_Step_Riverlak, step_RiverlakOutflow_m3);
     subset_put(Riverlak_water_m3, idx_Riverlak_Step,  step_RiverlakWater + step_RiverlakInflow - step_RiverlakOutflow_m3);
   }
   
   
   subset_put(RIVER_outflow_m3, idx_Cell_Step, step_RIVER_Outflow);
   subset_put(RIVER_water_m3, idx_Cell_Step,  step_RIVER_Water_New);
   
   
   
 }
 
 subset_add(RIVER_outflow_m3, Riverlak_cellNumber_int, Riverlak_overflow_m3);
 return RIVER_outflow_m3;
 
}



//' @rdname confluen
//' @export
// [[Rcpp::export]]
NumericVector confluen_WaterGAP3_N(
   NumericVector &RIVER_water_m3,
   NumericVector RIVER_length_km,
   NumericVector RIVER_velocity_km,
   NumericVector RIVER_outflow_m3,
   List CELL_cellNumberStep_int,
   List CELL_inflowCellNumberStep_int,
   IntegerVector Riverlak_cellNumber_int,
   NumericVector Riverlak_capacity_m3,
   NumericVector param_Riverlak_lin_storeFactor
)
{
 
 NumericVector step_RiverOutflow_m3, step_RiverlakOutflow_m3, step_UpstreamInflow_m3;
 
 IntegerVector idx_Cell_Step,
 idx_Step_Upstream,
 idx_Riverlak_Step, idx_Step_Riverlak;
 int n_Step = CELL_cellNumberStep_int.size();
 
 // Riverlak
 NumericVector Riverlak_water_m3 = subset_get(RIVER_water_m3, Riverlak_cellNumber_int);
 // Overflow for Riverlak
 NumericVector Riverlak_overflow_m3 = pmax(Riverlak_water_m3 - Riverlak_capacity_m3, 0);
 Riverlak_water_m3 += -Riverlak_overflow_m3;
 subset_put(RIVER_water_m3, Riverlak_cellNumber_int,  Riverlak_overflow_m3);
 
 // Step i later with Inflow
 for (int i_Step = 1; i_Step < n_Step; i_Step++)
 {
   
   idx_Cell_Step = CELL_cellNumberStep_int[i_Step];
   
   
   // Inflow upstream
   step_UpstreamInflow_m3 = inflow_add(
     RIVER_outflow_m3,
     CELL_inflowCellNumberStep_int[i_Step]
   );
   
   
   // river segment
   NumericVector step_RiverWater = subset_get(RIVER_water_m3, idx_Cell_Step);
   
   step_RiverOutflow_m3 = riverout_LinearResorvoir(
     step_RiverWater,
     step_UpstreamInflow_m3,
     subset_get(RIVER_velocity_km, idx_Cell_Step),
     subset_get(RIVER_length_km, idx_Cell_Step)
   );
   NumericVector step_RIVER_Water_New = pmax(step_RiverWater + step_UpstreamInflow_m3 - step_RiverOutflow_m3, 0.0);
   NumericVector step_RIVER_Outflow = step_RiverWater + step_UpstreamInflow_m3 - step_RIVER_Water_New;
   
   // Riverlak
   idx_Riverlak_Step = get_idx_cell(Riverlak_cellNumber_int, idx_Cell_Step);
   idx_Step_Riverlak = get_idx_step(Riverlak_cellNumber_int, idx_Cell_Step);
   if (idx_Riverlak_Step.size() > 0) {
     NumericVector step_RiverlakWater= subset_get(Riverlak_water_m3, idx_Riverlak_Step),
       step_RiverlakInflow = subset_get(step_UpstreamInflow_m3, idx_Step_Riverlak);
     
     
     step_RiverlakOutflow_m3 = riverlakout_LinearResorvoir(
       step_RiverlakWater,
       step_RiverlakInflow,
       subset_get(Riverlak_capacity_m3, idx_Riverlak_Step),
       subset_get(param_Riverlak_lin_storeFactor, idx_Riverlak_Step)
     );
     subset_put(step_RIVER_Outflow, idx_Step_Riverlak, step_RiverlakOutflow_m3);
     subset_put(step_RIVER_Water_New, idx_Step_Riverlak,  step_RiverlakWater + step_RiverlakInflow - step_RiverlakOutflow_m3);
   }
   
   
   subset_put(RIVER_outflow_m3, idx_Cell_Step, step_RIVER_Outflow);
   subset_put(RIVER_water_m3, idx_Cell_Step,  step_RIVER_Water_New);
   
   
   
 }
 
 subset_add(RIVER_outflow_m3, Riverlak_cellNumber_int, Riverlak_overflow_m3);
 return RIVER_outflow_m3;
 
}


//' @rdname confluen
//' @export
// [[Rcpp::export]]
NumericVector confluen_WaterGAP3(
   NumericVector &RIVER_water_m3,
   NumericVector RIVER_length_km,
   NumericVector RIVER_velocity_km,
   NumericVector RIVER_outflow_m3,
   List CELL_cellNumberStep_int,
   List CELL_inflowCellNumberStep_int,
   IntegerVector Riverlak_cellNumber_int,
   NumericVector Riverlak_capacity_m3,
   IntegerVector Reservoi_cellNumber_int,
   NumericVector Reservoi_inflow_m3,
   NumericVector Reservoi_demand_m3,
   NumericVector Reservoi_capacity_m3,
   NumericVector Reservoi_meanInflow_m3,
   NumericVector Reservoi_meanDemand_m3,
   NumericVector Reservoi_releaseCoefficient_1,
   LogicalVector Reservoi_isIrrigate_01,
   NumericVector param_Riverlak_lin_storeFactor
)
{
 
 NumericVector step_RiverOutflow_m3, 
 step_RiverlakOutflow_m3, step_ReservoiOutflow_m3,
 step_UpstreamInflow_m3;
 
 IntegerVector idx_Cell_Step,
 idx_Step_Upstream,
 idx_Riverlak_Step, idx_Step_Riverlak,
 idx_Reservoi_Step, idx_Step_Reservoi;
 int n_Step = CELL_cellNumberStep_int.size();
 
 // Riverlak
 NumericVector Riverlak_water_m3 = subset_get(RIVER_water_m3, Riverlak_cellNumber_int);
 // Overflow for Riverlak
 NumericVector Riverlak_overflow_m3 = pmax(Riverlak_water_m3 - Riverlak_capacity_m3, 0);
 Riverlak_water_m3 += -Riverlak_overflow_m3;
 subset_put(RIVER_water_m3, Riverlak_cellNumber_int,  Riverlak_overflow_m3);
 
 // Reservoi
 NumericVector Reservoi_water_m3 = subset_get(RIVER_water_m3, Reservoi_cellNumber_int);
 // Overflow for Reservoi
 NumericVector Reservoi_overflow_m3 = pmax(Reservoi_water_m3 - Reservoi_capacity_m3, 0);
 Reservoi_water_m3 += -Reservoi_overflow_m3;
 subset_put(RIVER_water_m3, Reservoi_cellNumber_int,  Reservoi_overflow_m3);
 
 // Step i later with Inflow
 for (int i_Step = 1; i_Step < n_Step; i_Step++)
 {
   
   idx_Cell_Step = CELL_cellNumberStep_int[i_Step];
   
   
   // Inflow upstream
   step_UpstreamInflow_m3 = inflow_add(
     RIVER_outflow_m3,
     CELL_inflowCellNumberStep_int[i_Step]
   );
   
   
   // river segment
   NumericVector step_RiverWater = subset_get(RIVER_water_m3, idx_Cell_Step);
   
   step_RiverOutflow_m3 = riverout_LinearResorvoir(
     step_RiverWater,
     step_UpstreamInflow_m3,
     subset_get(RIVER_velocity_km, idx_Cell_Step),
     subset_get(RIVER_length_km, idx_Cell_Step)
   );
   NumericVector step_RIVER_Water_New = pmax(step_RiverWater + step_UpstreamInflow_m3 - step_RiverOutflow_m3, 0.0);
   NumericVector step_RIVER_Outflow = step_RiverWater + step_UpstreamInflow_m3 - step_RIVER_Water_New;
   
   // Riverlak
   idx_Riverlak_Step = get_idx_cell(Riverlak_cellNumber_int, idx_Cell_Step);
   idx_Step_Riverlak = get_idx_step(Riverlak_cellNumber_int, idx_Cell_Step);
   if (idx_Riverlak_Step.size() > 0) {
     NumericVector step_RiverlakWater= subset_get(Riverlak_water_m3, idx_Riverlak_Step),
       step_RiverlakInflow = subset_get(step_UpstreamInflow_m3, idx_Step_Riverlak);
     
     
     step_RiverlakOutflow_m3 = riverlakout_LinearResorvoir(
       step_RiverlakWater,
       step_RiverlakInflow,
       subset_get(Riverlak_capacity_m3, idx_Riverlak_Step),
       subset_get(param_Riverlak_lin_storeFactor, idx_Riverlak_Step)
     );
     subset_put(step_RIVER_Outflow, idx_Step_Riverlak, step_RiverlakOutflow_m3);
     subset_put(step_RIVER_Water_New, idx_Step_Riverlak,  step_RiverlakWater + step_RiverlakInflow - step_RiverlakOutflow_m3);
   }
   
   
   subset_put(RIVER_outflow_m3, idx_Cell_Step, step_RIVER_Outflow);
   subset_put(RIVER_water_m3, idx_Cell_Step,  step_RIVER_Water_New);
   
   // Reservoi
   idx_Reservoi_Step = get_idx_cell(Reservoi_cellNumber_int, idx_Cell_Step);
   idx_Step_Reservoi = get_idx_step(Reservoi_cellNumber_int, idx_Cell_Step);
   if (idx_Reservoi_Step.size() > 0) {
     NumericVector step_ReservoiWater= subset_get(Reservoi_water_m3, idx_Reservoi_Step),
       step_ReservoiInflow = subset_get(step_UpstreamInflow_m3, idx_Step_Reservoi),
       step_ReservoiDemand = subset_get(Reservoi_demand_m3, idx_Reservoi_Step);
     
     
     step_ReservoiOutflow_m3 = reservoireleas_Hanasaki(
       step_ReservoiWater,
       step_ReservoiInflow,
       step_ReservoiDemand,
       subset_get(Reservoi_capacity_m3, idx_Reservoi_Step),
       subset_get(Reservoi_meanInflow_m3, idx_Reservoi_Step),
       subset_get(Reservoi_meanDemand_m3, idx_Reservoi_Step),
       subset_get(Reservoi_releaseCoefficient_1, idx_Reservoi_Step), 
       subset_get_logical(Reservoi_isIrrigate_01, idx_Reservoi_Step)
     );
     subset_put(step_RIVER_Outflow, idx_Step_Reservoi, step_ReservoiOutflow_m3);
     subset_put(step_RIVER_Water_New, idx_Step_Reservoi,  step_ReservoiWater + step_ReservoiInflow - step_ReservoiOutflow_m3);
   }
   
   
   subset_put(RIVER_outflow_m3, idx_Cell_Step, step_RIVER_Outflow);
   subset_put(RIVER_water_m3, idx_Cell_Step,  step_RIVER_Water_New);
   
   
 }
 
 subset_add(RIVER_outflow_m3, Riverlak_cellNumber_int, Riverlak_overflow_m3);
 return RIVER_outflow_m3;
 
}


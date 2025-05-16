#include "utils.h"
// [[Rcpp::depends(RcppArmadillo)]]
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
arma::vec confluen_IUH(
    const arma::vec& CONFLUEN_inputWater_mm, 
    const arma::vec& CONFLUEN_iuh_1) {
  
  int n_iuh = CONFLUEN_iuh_1.n_elem;
  int n_time = CONFLUEN_inputWater_mm.n_elem;
  arma::vec CONFLUEN_outputWater_mm(n_time, arma::fill::zeros);
  
  for (int i = 0; i < n_iuh; i++) {
    for (int j = 0; j <= i; j++) {
      CONFLUEN_outputWater_mm[i] += CONFLUEN_inputWater_mm[i - j] * CONFLUEN_iuh_1[j];
    }
  }
  
  for (int i = n_iuh; i < n_time; i++) {
    for (int j = 0; j < n_iuh; j++) {
      CONFLUEN_outputWater_mm[i] += CONFLUEN_inputWater_mm[i - j] * CONFLUEN_iuh_1[j];
    }
  }
  
  return CONFLUEN_outputWater_mm;
}

//' @rdname confluen
//' @export
// [[Rcpp::export]]
arma::vec confluen_IUH2S(
    const arma::vec& LAND_runoff_mm,
    const arma::vec& GROUND_baseflow_mm, 
    const arma::vec& CONFLUEN_iuhLand_1,
    const arma::vec& CONFLUEN_iuhGround_1) {
  
  arma::vec CONFLUEN_runoff_mm = confluen_IUH(LAND_runoff_mm, CONFLUEN_iuhLand_1);
  arma::vec CONFLUEN_baseflow_mm = confluen_IUH(GROUND_baseflow_mm, CONFLUEN_iuhGround_1);
  
  return CONFLUEN_runoff_mm + CONFLUEN_baseflow_mm;
}

//' @rdname confluen
//' @export
// [[Rcpp::export]]
arma::vec confluen_IUH3S(
    const arma::vec& LAND_runoff_mm,
    const arma::vec& SOIL_interflow_mm, 
    const arma::vec& GROUND_baseflow_mm, 
    const arma::vec& CONFLUEN_iuhLand_1,
    const arma::vec& CONFLUEN_iuhSoil_1,
    const arma::vec& CONFLUEN_iuhGround_1) {
  
  arma::vec CONFLUEN_runoff_mm = confluen_IUH(LAND_runoff_mm, CONFLUEN_iuhLand_1);
  arma::vec CONFLUEN_interflow_mm = confluen_IUH(SOIL_interflow_mm, CONFLUEN_iuhSoil_1);
  arma::vec CONFLUEN_baseflow_mm = confluen_IUH(GROUND_baseflow_mm, CONFLUEN_iuhGround_1);
  
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
//' \mjsdeqn{u(i) = S(i) - S(i-1)}
//' \mjsdeqn{S(i) = \left( \frac{i}{t_r} \right)^{2.5}, \quad 0 \leq i \leq t_r}
//' where
//'   - \mjseqn{u} is IUH series
//'   - \mjseqn{i} is index
//' @export
// [[Rcpp::export]]
arma::vec confluenIUH_GR4J1(double CONFLUEN_responseTime_TS) {
  int t_max = std::ceil(CONFLUEN_responseTime_TS);
  arma::vec seq_t = arma::regspace(1, t_max);
  arma::vec SH_1 = arma::pow(seq_t / CONFLUEN_responseTime_TS, 2.5);
  SH_1(t_max - 1) = 1.0;
  SH_1.subvec(1, t_max - 1) = arma::diff(SH_1);
  return SH_1;
}

//' @rdname confluenIUH
//' @details
//' # **_GR4J2** \insertCite{GR4J_Perrin_2003}{HydroGallery}: 
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
arma::vec confluenIUH_GR4J2(double CONFLUEN_responseTime_TS) {
  int t_max_1 = std::ceil(CONFLUEN_responseTime_TS);
  int t_max_2 = std::ceil(2 * CONFLUEN_responseTime_TS);
  
  arma::vec seq_t1 = arma::regspace(1, t_max_1 - 1);
  arma::vec seq_t2 = arma::regspace(t_max_1, t_max_2 - 1);
  
  arma::vec SH_2_1 = 0.5 * arma::pow(seq_t1 / CONFLUEN_responseTime_TS, 2.5);
  arma::vec SH_2_2 = 1.0 - 0.5 * arma::pow(2.0 - seq_t2 / CONFLUEN_responseTime_TS, 2.5);
  
  arma::vec SH_2(t_max_2, arma::fill::ones);
  SH_2.subvec(0, t_max_1 - 2) = SH_2_1;
  SH_2.subvec(t_max_1 - 1, t_max_2 - 2) = SH_2_2;
  SH_2.subvec(1, t_max_2 - 1) = arma::diff(SH_2);
  
  return SH_2;
}

//' @rdname confluenIUH
//' @details
//' # **_Kelly** \insertCite{iuh_Kelly_1955}{HydroGallery}: 
//'
//' \mjsdeqn{u(i) = \frac{4}{t_r^2} \left( i + k \left( e^{-i/k} \right) \right), \quad i \leq t_r / 2 }
//' \mjsdeqn{u(i) = - \frac{4}{t_r^2}(i - k - t_r) + \frac{4ke^{-i/k}}{t_r^2} (1 - 2 e^{t_r/(2k)}), \quad t_r / 2 < i \leq t_r }
//' \mjsdeqn{u(i) =  \frac{4ke^{-i/k}}{t_r^2} (1 - 2 e^{t_r/(2k)} +  e^{t_r/k}), \quad i > t_r }
//' where
//'   - \mjseqn{k} is `param_CONFLUEN_kel_k`
//' @param param_CONFLUEN_kel_k <1, 4> parameter for[confluenIUH_Kelly()]
//' @export
// [[Rcpp::export]]
arma::vec confluenIUH_Kelly(double CONFLUEN_responseTime_TS, double param_CONFLUEN_kel_k) {
  double tc = CONFLUEN_responseTime_TS * param_CONFLUEN_kel_k;
  double tc2 = tc * tc;

  double part_34 = 4 * CONFLUEN_responseTime_TS / tc2 *
    (1 - 2 * std::exp(tc / (2 * CONFLUEN_responseTime_TS)));
  double part_35 = 4 * CONFLUEN_responseTime_TS / tc2 *
    (1 - 2 * std::exp(tc / (2 * CONFLUEN_responseTime_TS)) + std::exp(tc / CONFLUEN_responseTime_TS));
  
  int t_max = std::ceil(std::max(tc, -CONFLUEN_responseTime_TS * std::log(0.002 / part_35)));

  arma::vec seq_t = arma::regspace(1, 20 * t_max) / 20.0;
  arma::vec etK = arma::exp(-seq_t / CONFLUEN_responseTime_TS);

  arma::vec iuh_1 = 4.0 / tc2 * (seq_t + CONFLUEN_responseTime_TS * (etK - 1.0));
  arma::vec iuh_2 = part_34 * etK - 4.0 / tc2 * (seq_t - CONFLUEN_responseTime_TS - tc);
  arma::vec iuh_3 = part_35 * etK;
  
  arma::vec iuh(seq_t.n_elem, arma::fill::zeros);
  arma::uvec idx1 = arma::find(seq_t <= 0.5 * tc);
  arma::uvec idx2 = arma::find(seq_t > 0.5 * tc && seq_t <= tc);
  arma::uvec idx3 = arma::find(seq_t > tc);
  
  iuh.elem(idx1) = iuh_1.elem(idx1);
  iuh.elem(idx2) = iuh_2.elem(idx2);
  iuh.elem(idx3) = iuh_3.elem(idx3);
  
  arma::mat mat_iuh = arma::reshape(iuh, 20, t_max);
  arma::vec vct_iuh = arma::mean(mat_iuh, 0).t();
  
  return vct_iuh / arma::sum(vct_iuh);
}

//' @rdname confluenIUH
//' @details
//' # **_Nash** \insertCite{iuh_Nash_1957}{HydroGallery}: 
//'
//' \mjsdeqn{u(i) = \frac{1}{t_r\Gamma(n)} \left(\frac{4}{t_r^2}\right)^{n -1}e^{-i/t_r}}
//' where
//'   - \mjseqn{n} is `param_CONFLUEN_nas_n`
//' @param param_CONFLUEN_nas_n <1, 8> parameter for[confluenIUH_Nash()]
//' @export
// [[Rcpp::export]]
arma::vec confluenIUH_Nash(double CONFLUEN_responseTime_TS, double param_CONFLUEN_nas_n) {
  int t_max = std::ceil(std::max(4.0, param_CONFLUEN_nas_n) * 3 * CONFLUEN_responseTime_TS);
  arma::vec seq_t = arma::regspace(1, 20 * t_max) / 20.0;
  
  arma::vec iuh = arma::pow(seq_t / CONFLUEN_responseTime_TS, param_CONFLUEN_nas_n - 1) %
                  arma::exp(-seq_t / CONFLUEN_responseTime_TS) /
                  (CONFLUEN_responseTime_TS * tgamma(param_CONFLUEN_nas_n));
  
  arma::mat mat_iuh = arma::reshape(iuh, 20, t_max);
  arma::vec vct_iuh = arma::mean(mat_iuh, 0).t();
  
  return vct_iuh / arma::sum(vct_iuh);
}

//' @rdname confluenIUH
//' @details
//' # **_Clark** \insertCite{iuh_Clark_1945}{HydroGallery}: 
//'
//' \mjsdeqn{u(i) = \frac{1}{t_r} e^{-i/t_r} }
//' where
//'   - \mjseqn{t_r} is `CONFLUEN_responseTime_TS`
//' @export
// [[Rcpp::export]]
arma::vec confluenIUH_Clark(double CONFLUEN_responseTime_TS) {
  int t_max = std::ceil(-CONFLUEN_responseTime_TS * std::log(CONFLUEN_responseTime_TS * 0.005));
  arma::vec seq_t = arma::regspace(1, 20 * t_max) / 20.0;
  
  arma::vec iuh = arma::exp(-seq_t / CONFLUEN_responseTime_TS) / CONFLUEN_responseTime_TS;
  
  arma::mat mat_iuh = arma::reshape(iuh, 20, t_max);
  arma::vec vct_iuh = arma::mean(mat_iuh, 0).t();
  
  return vct_iuh / arma::sum(vct_iuh);
}
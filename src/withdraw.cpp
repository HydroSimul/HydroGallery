#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

//' Withdraw water
//' @name withdraw
//' @description
//' [withdraw_SingleCell] This function handles water withdrawal from a single cell's water storage.
//' It updates the withdrawal deficit and the remaining water volume.
//' @param CELL_withdrawal_m3 Vector of withdrawal requirements in cubic meters.
//'   This represents water demand and will be updated to reflect remaining deficit.
//' @param CELL_water_m3 Vector of available water volumes in cubic meters.
//'   This will be updated to reflect remaining water after withdrawal.
//' @return None. Parameters are updated by reference.
//' @export
// [[Rcpp::export]]
void withdraw_SingleCell(arma::vec& CELL_withdrawal_m3,
                         arma::vec& CELL_water_m3) {

  arma::vec withdrawActual = arma::min(CELL_withdrawal_m3, CELL_water_m3);

  CELL_withdrawal_m3 -= withdrawActual;
  CELL_water_m3 -= withdrawActual;
}

//' [withdrawSurface_AroundMax] This function identifies the surrounding cell with the maximum water availability
//' and withdraws water from it to satisfy the demand.
//' @rdname withdraw
//' @param CELL_cellNumberAround_int Matrix of cell indices that define surrounding cells.
//' @export
// [[Rcpp::export]]
void withdrawSurface_AroundMax(arma::vec& CELL_withdrawal_m3,
                               arma::vec& RIVER_water_m3,
                               arma::vec& Lake_water_m3,
                               const arma::imat& CELL_cellNumberAround_int) {

  arma::uword n_Spat = CELL_withdrawal_m3.n_elem;

  for (arma::uword i_Spat = 0; i_Spat < n_Spat; i_Spat++) {
    if (CELL_withdrawal_m3[i_Spat] <= 0) continue;

    arma::ivec idx_Around = CELL_cellNumberAround_int.col(i_Spat) - 1;
    idx_Around = idx_Around(find(idx_Around >= 0));
    arma::uword n_Around = idx_Around.n_elem;
    if (n_Around == 0) continue;

    arma::vec allWater_Around(n_Around);
    for (arma::uword j = 0; j < n_Around; j++) {
      int cell_idx = idx_Around[j];
      allWater_Around[j] = RIVER_water_m3[cell_idx] + Lake_water_m3[cell_idx];
    }

    arma::uword idx_MaxAround = allWater_Around.index_max();
    int target_cell = idx_Around[idx_MaxAround];
    double withdrawal = CELL_withdrawal_m3[i_Spat];
    double river = RIVER_water_m3[target_cell];
    double lake = Lake_water_m3[target_cell];
    double total_water = river + lake;

    if (withdrawal > total_water) {
      RIVER_water_m3[target_cell] = 0;
      Lake_water_m3[target_cell] = 0;
      CELL_withdrawal_m3[i_Spat] -= total_water;
    } else if (withdrawal > river) {
      RIVER_water_m3[target_cell] = 0;
      Lake_water_m3[target_cell] = total_water - withdrawal;
      CELL_withdrawal_m3[i_Spat] = 0;
    } else {
      RIVER_water_m3[target_cell] -= withdrawal;
      CELL_withdrawal_m3[i_Spat] = 0;
    }
  }
}

//' [withdrawSurface_Around] This function withdraws water from all surrounding cells proportionally
//' based on their water availability.
//' @rdname withdraw
//' @export
// [[Rcpp::export]]
void withdrawSurface_Around(arma::vec& CELL_withdrawal_m3,
                            arma::vec& RIVER_water_m3,
                            arma::vec& Lake_water_m3,
                            const arma::imat& CELL_cellNumberAround_int) {

  arma::uword n_Spat = CELL_withdrawal_m3.n_elem;

  for (arma::uword i_Spat = 0; i_Spat < n_Spat; i_Spat++) {
    if (CELL_withdrawal_m3[i_Spat] <= 0) continue;

    arma::ivec idx_Around = CELL_cellNumberAround_int.col(i_Spat) - 1;
    idx_Around = idx_Around(find(idx_Around >= 0));
    arma::uword n_Around = idx_Around.n_elem;
    if (n_Around == 0) continue;

    arma::ivec valid_indices = idx_Around(find(idx_Around >= 0));
    if (valid_indices.n_elem == 0) continue;

    double total_available = 0;
    arma::vec totalWater_Around(valid_indices.n_elem);

    for (arma::uword j = 0; j < valid_indices.n_elem; j++) {
      int idx = valid_indices[j];
      double water = RIVER_water_m3[idx] + Lake_water_m3[idx];
      totalWater_Around[j] = water;
      total_available += water;
    }

    double withdrawal = CELL_withdrawal_m3[i_Spat];

    if (withdrawal > total_available) {
      for (arma::uword j = 0; j < valid_indices.n_elem; j++) {
        int idx = valid_indices[j];
        RIVER_water_m3[idx] = 0;
        Lake_water_m3[idx] = 0;
      }
      CELL_withdrawal_m3[i_Spat] -= total_available;
    } else {
      arma::vec weights = totalWater_Around / total_available;

      for (arma::uword j = 0; j < valid_indices.n_elem; j++) {
        int idx = valid_indices[j];
        double cell_withdrawal = withdrawal * weights[j];

        double river = RIVER_water_m3[idx];
        double lake = Lake_water_m3[idx];
        double total_cell_water = river + lake;

        if (cell_withdrawal >= total_cell_water) {
          RIVER_water_m3[idx] = 0;
          Lake_water_m3[idx] = 0;
        } else if (cell_withdrawal > river) {
          RIVER_water_m3[idx] = 0;
          Lake_water_m3[idx] = total_cell_water - cell_withdrawal;
        } else {
          RIVER_water_m3[idx] -= cell_withdrawal;
        }
      }

      CELL_withdrawal_m3[i_Spat] = 0;
    }
  }
}

//' [withdrawSurface_WithdrawNet] This function withdraws water by following a predefined network of cells
//' from which water can be withdrawn.
//' @name withdraw
//' @param CELL_cellNumberWithdrawNet_int Matrix defining the withdrawal network.
//' @export
// [[Rcpp::export]]
void withdrawSurface_WithdrawNet(arma::vec& CELL_withdrawal_m3,
                                 arma::vec& RIVER_water_m3,
                                 arma::vec& Lake_water_m3,
                                 const arma::imat& CELL_cellNumberWithdrawNet_int) {

  arma::uword n_Spat = CELL_withdrawal_m3.n_elem;

  for (arma::uword i_Spat = 0; i_Spat < n_Spat; i_Spat++) {
    if (CELL_withdrawal_m3[i_Spat] <= 0) continue;

    arma::ivec idx_WithdrawNet = CELL_cellNumberWithdrawNet_int.col(i_Spat) - 1;

    arma::uword idx_WD = 0;
    while (CELL_withdrawal_m3[i_Spat] > 0 && idx_WD < static_cast<arma::uword>(idx_WithdrawNet.n_elem)) {
      int idx_Spat_WD = idx_WithdrawNet[idx_WD];

      if (idx_Spat_WD < 0) {
        idx_WD++;
        continue;
      }

      double river = RIVER_water_m3[idx_Spat_WD];
      double lake = Lake_water_m3[idx_Spat_WD];
      double total_Water = river + lake;
      double withdrawal = CELL_withdrawal_m3[i_Spat];

      if (withdrawal > total_Water) {
        RIVER_water_m3[idx_Spat_WD] = 0;
        Lake_water_m3[idx_Spat_WD] = 0;
        CELL_withdrawal_m3[i_Spat] -= total_Water;
      } else if (withdrawal > river) {
        RIVER_water_m3[idx_Spat_WD] = 0;
        Lake_water_m3[idx_Spat_WD] = total_Water - withdrawal;
        CELL_withdrawal_m3[i_Spat] = 0;
        break;
      } else {
        RIVER_water_m3[idx_Spat_WD] -= withdrawal;
        CELL_withdrawal_m3[i_Spat] = 0;
        break;
      }

      idx_WD++;
    }
  }
}

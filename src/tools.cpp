#include <iostream>
#include <cmath>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  int estimations_size = estimations.size();

  if (estimations_size != ground_truth.size()) {
    std::cout << "Estimations and ground truth must be equal size\n";
    return rmse;
  }

  if (estimations_size == 0) {
    std::cout << "Cannot calculate RMSE of 0 estimations\n";
    return rmse;
  }

  for (int i = 0; i < estimations_size; ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  rmse = rmse / estimations_size;
  rmse = rmse.array().sqrt();
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
}

VectorXd Tools::CalculateZpred(const VectorXd &x_state) {
  VectorXd z_pred = VectorXd(3);
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(2);
  double px2 = px * px;
  double py2 = py * py;
  double sqrt_px2_py2 = std::sqrt(px2 + py2);
  z_pred(0) = sqrt_px2_py2;
  z_pred(1) = std::atan2(px, py);
  z_pred(2) = (px * vx + py * vy) / sqrt_px2_py2;
  return z_pred;
}

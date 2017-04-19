#include <iostream>
#include <cmath>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

const double delta = 0.001;

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
  const double px = x_state(0);
  const double py = x_state(1);
  const double vx = x_state(2);
  const double vy = x_state(2);

  const double c1 = px * px + py * py;
  const double c2 = std::sqrt(c1);
  const double c3 = c1 * c2;

  MatrixXd Hj = MatrixXd(3, 4);
  if (std::fabs(c1) < delta) {
    std::cout << "Division by zero in CalculateJacobian() - Returning zero Hj\n";
    Hj << 0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0;
    return Hj;
  }

  Hj << (px/c2), (py/c2), 0, 0,
        -(py/c1), (px/c1), 0, 0,
        py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

  return Hj;
}

VectorXd Tools::CalculateZpred(const VectorXd &x_state) {
  const double px = x_state(0);
  const double py = x_state(1);
  const double vx = x_state(2);
  const double vy = x_state(3);

  const double c1 = px * px + py * py;
  const double c2 = std::sqrt(c1);

  VectorXd z_pred = VectorXd(3);
  if (std::fabs(c1) < delta || std::fabs(px) < delta) {
    std::cout << "Division by zero in CalculateZpred() - Returning zero z_pred\n";
    z_pred << 0, 0, 0;
    return z_pred;
  }

  z_pred(0) = c2;
  z_pred(1) = std::atan2(py, px);
  z_pred(2) = (px*vx + py*vy)/c2;

  return z_pred;
}

#include <math.h>
#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z, const MatrixXd &R) {
  H_ = MatrixXd(2, 4);
  H_ << 1, 0, 0, 0,
        0, 1, 0, 0;
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;

  FinishUpdate(y, R);
}

void KalmanFilter::UpdateEKF(const VectorXd &z, const MatrixXd &R) {
  VectorXd z_pred = tools.CalculateZpred(x_);
  VectorXd y = z - z_pred;
  // normalize phi between -pi and +pi
  while (y(1) < -M_PI) {
      y(1) += 2 * M_PI;
  }
  while (y(1) > M_PI) {
      y(1) -= 2 * M_PI;
  }
  H_ = tools.CalculateJacobian(x_);

  FinishUpdate(y, R);
}

void KalmanFilter::FinishUpdate(const VectorXd &y, const MatrixXd &R) {
  MatrixXd Ht = H_.transpose();
  MatrixXd P_Ht = P_ * Ht;
  MatrixXd S = H_ * P_Ht + R;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_Ht * Si;

  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1000, 0,
               0, 0, 0, 1000;

    VectorXd raw_meas = measurement_pack.raw_measurements_;

    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1, 0,
               0, 0, 0, 1;

    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << 0, 0, 0, 0,
               0, 0, 0, 0,
               0, 0, 0, 0,
               0, 0, 0, 0;

    previous_timestamp_ = measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double x_proj = cos(raw_meas[1]);
      double y_proj = sin(raw_meas[1]);
      ekf_.x_[0] = raw_meas[0] * x_proj;
      ekf_.x_[1] = raw_meas[0] * y_proj;
      ekf_.x_[2] = raw_meas[2] * x_proj;
      ekf_.x_[3] = raw_meas[2] * y_proj;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_[0] = raw_meas[0];
      ekf_.x_[1] = raw_meas[1];
      ekf_.x_[2] = 0;
      ekf_.x_[3] = 0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  long long current_timestamp_ = measurement_pack.timestamp_;

  int noise_ax = 9;
  int noise_ay = 9;
  // dt expressed in seconds
  float dt = (current_timestamp_ - previous_timestamp_) / 1000000.0;
  float dt_2 = dt * dt;
  float dt_3 = dt * dt_2;
  float dt_4 = dt * dt_3;

  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  ekf_.Q_(0, 0) = dt_4 / 4 * noise_ax;
  ekf_.Q_(1, 1) = dt_4 / 4 * noise_ay;
  ekf_.Q_(2, 2) = dt_2 * noise_ax;
  ekf_.Q_(3, 3) = dt_2 * noise_ay;
  ekf_.Q_(0, 2) = dt_3 / 2 * noise_ax;
  ekf_.Q_(2, 0) = dt_3 / 2 * noise_ax;
  ekf_.Q_(1, 3) = dt_3 / 2 * noise_ay;
  ekf_.Q_(3, 1) = dt_3 / 2 * noise_ay;

  previous_timestamp_ = current_timestamp_;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    VectorXd z = VectorXd(3);
    z << measurement_pack.raw_measurements_[0],
         measurement_pack.raw_measurements_[1],
         measurement_pack.raw_measurements_[2];
    ekf_.UpdateEKF(z, R_radar_);
  } else {
    // Laser updates
    VectorXd z = VectorXd(2);
    z << measurement_pack.raw_measurements_[0],
         measurement_pack.raw_measurements_[1];
    ekf_.Update(z, R_laser_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

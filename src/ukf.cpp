

#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
//  std_a_ = 30;
  std_a_ = 4.0;     

  // Process noise standard deviation yaw acceleration in rad/s^2
//  std_yawdd_ = 30;
  std_yawdd_ = 4.0;

  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

   // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
   * End DO NOT MODIFY section for measurement noise values
   */

   /**
    * TODO: Complete the initialization. See ukf.h for other member properties.
    * Hint: one or more values initialized above might be wildly off...
    */

  // The first measurement initialises the filter.
  is_initialized_ = false;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // set radar measurement dimension (r, phi and r_dot)
  n_z_R_ = 3;

  // set Laser measurement dimension (rx, ry)
  n_z_L_ = 2;

  lambda_ = 3. - n_aug_;

  weights_ = VectorXd(2 * n_aug_ + 1);

  // Compute the wights for calculating the sigma points
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  double weight = 0.5 / (lambda_ + n_aug_);
  weights_(0) = weight_0;

  for (int i = 1; i < 2 * n_aug_ + 1; ++i) {
    weights_(i) = weight;
  }

  // predicted sigma points.
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // radar mean predicted measurement
  z_pred_r_ = VectorXd(n_z_R_);

  // laser mean predicted measurement
  z_pred_l_ = VectorXd(n_z_L_);

  // radar measurement space sigma points
  Zsig_r_ = MatrixXd(n_z_R_, 2 * n_aug_ + 1);

  // laser measurement space sigma points
  Zsig_l_ = MatrixXd(n_z_L_, 2 * n_aug_ + 1);;

  // radar measurement covariance matrix
  S_r_ = MatrixXd(n_z_R_, n_z_R_);

  // laser measurement covariance matrix
  S_l_ = MatrixXd(n_z_L_, n_z_L_);

  // the laser measurement noise matrix
  R_l = MatrixXd(n_z_L_, n_z_L_);
  R_l << std::pow(std_laspx_, 2), 0,
    0, std::pow(std_laspy_, 2);

  // the radar measurement noise matrix
  R_r = MatrixXd(n_z_R_, n_z_R_);
  R_r << std::pow(std_radr_, 2), 0, 0,
    0, std::pow(std_radphi_, 2), 0,
    0, 0, std::pow(std_radrd_, 2);

  proc_noise_ = MatrixXd(2, 2);
  proc_noise_ << std::pow(std_a_, 2), 0,
    0, std::pow(std_yawdd_, 2);

  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
  if (!is_initialized_)
  {
    /// Check if the measurement comes from the laser
    if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      // Initialise the state from the laser
      // set the state with the initial location and zero velocity
      const double px = meas_package.raw_measurements_[0];
      const double py = meas_package.raw_measurements_[1];
      x_ << px, py, 0, 0, 0;

      // Initiase the state covariance matrix
      P_ << std::pow(std_laspx_, 2), 0, 0, 0, 0,
        0, std::pow(std_laspy_, 2), 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;      

    }
    else // Process a radar measurement
    {
      // Initialise the state from the radar
      const double r = meas_package.raw_measurements_[0];
      const double phi = meas_package.raw_measurements_[1];
      const double r_dot = meas_package.raw_measurements_[2];

      const double px = r * cos(phi);
      const double py = r * sin(phi);

      x_ << px, py, 0, 0, 0;

      // Initiase the state covariance matrix
      P_ << std::pow(std_radr_, 2), 0, 0, 0, 0,
        0, std::pow(std_radr_, 2), 0, 0, 0,
        0, 0, std::pow(std_radrd_, 2), 0, 0,
        0, 0, 0, std::pow(std_radphi_, 2), 0,
        0, 0, 0, 0, std::pow(std_radphi_, 2);   

    }

    is_initialized_ = true;

    // previous timestamp
    time_us_ = meas_package.timestamp_;

    return;
  }

  // Compute the time delta
  double delta_t = (meas_package.timestamp_ - time_us_) / 1e6;    // convert to seconds
  time_us_ = meas_package.timestamp_; // Store for the next update

  // Do the predition step
  Prediction(delta_t);

  // And now do the update step
  /// Check if the measurement comes from the laser
  if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    UpdateLidar(meas_package);
  }
  else
  {
    UpdateRadar(meas_package);
  }
}

void UKF::Prediction(double delta_t) 
{
  /**
   * TODO: Complete this function! Estimate the object's location.
   * Modify the state vector, x_. Predict sigma points, the state,
   * and the state covariance matrix.
   */

   // Create the Augmented sigma points
  AugmentedSigmaPoints();

  // now predict them
  SigmaPointPrediction(delta_t);

  // Predict the new mean and covariance
  PredictMeanAndCovariance();
}

void UKF::UpdateLidar(MeasurementPackage meas_package) 
{
  /**
   * TODO: Complete this function! Use lidar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

  PredictLidarMeasurement();

  UpdateStateLidar(meas_package);
}

void UKF::UpdateRadar(MeasurementPackage meas_package) 
{
  /**
   * TODO: Complete this function! Use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  PredictRadarMeasurement();

  UpdateStateRadar(meas_package);
}

// Assignment functions

void UKF::AugmentedSigmaPoints()
{
  // create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // create augmented mean state
  x_aug.head(n_x_) = x_;
  for (int i = n_x_; i < n_aug_; ++i)
  {
    x_aug[i] = 0;
  }

  // create augmented covariance matrix
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  for (int i = n_x_; i < n_aug_; ++i)
  {
    P_aug.row(i).setZero();
    P_aug.col(i).setZero();
  }

  P_aug.block(n_x_, n_x_, 2, 2) = proc_noise_;

  // The first block 
  MatrixXd lambdaTerm(n_aug_, n_aug_);
  lambdaTerm = (lambda_ + n_aug_) * P_aug;
  lambdaTerm = lambdaTerm.llt().matrixL();
  MatrixXd first(n_aug_, n_aug_);
  for (int r = 0; r < n_aug_; ++r)
  {
    for (int c = 0; c < n_aug_; ++c)
    {
      first(r, c) = x_aug(r) + lambdaTerm(r, c);
    }
  }

  // The second block
  MatrixXd second(n_aug_, n_aug_);
  for (int r = 0; r < n_aug_; ++r)
  {
    for (int c = 0; c < n_aug_; ++c)
    {
      second(r, c) = x_aug(r) - lambdaTerm(r, c);
    }
  }

  // set sigma points as columns of matrix Xsig
  Xsig_aug.block(0, 0, n_aug_, 1) = x_aug;
  Xsig_aug.block(0, 1, n_aug_, n_aug_) = first;
  Xsig_aug.block(0, 1 + n_aug_, n_aug_, n_aug_) = second;

  // write result
  Xsig_aug_ = Xsig_aug;
}

void UKF::SigmaPointPrediction(double delta_t)
{
    // predict sigma points
  for (int c = 0; c < 2 * n_aug_ + 1; ++c)
  {
    MatrixXd col(n_aug_, 1);
    MatrixXd colOut(n_x_, 1);

    col = Xsig_aug_.col(c);

    const double px = col(0, 0);
    const double py = col(1, 0);
    const double v = col(2, 0);
    const double psi = col(3, 0);
    const double psi_dot = col(4, 0);
    const double nu_a = col(5, 0);
    const double nu_psi_dot_dot = col(6, 0);

    double px_pred;
    double py_pred;
    double v_pred;
    double psi_pred;
    double psi_dot_pred;

    if (fabs(psi_dot) > 0.0001)
    {
      px_pred = px + ((v / psi_dot) * (sin(psi + psi_dot * delta_t) - sin(psi))) + (0.5 * std::pow(delta_t, 2) * cos(psi) * nu_a);
      py_pred = py + ((v / psi_dot) * (-cos(psi + psi_dot * delta_t) + cos(psi))) + (0.5 * std::pow(delta_t, 2) * sin(psi) * nu_a);
      v_pred = v + (delta_t * nu_a);
      psi_pred = psi + (psi_dot * delta_t) + (0.5 * std::pow(delta_t, 2) * nu_psi_dot_dot);
      psi_dot_pred = psi_dot + (delta_t * nu_psi_dot_dot);
    }
    else
    {
      px_pred = px + (v * cos(psi) * delta_t) + (0.5 * std::pow(delta_t, 2) * cos(psi) * nu_a);
      py_pred = py + (v * sin(psi) * delta_t) + (0.5 * std::pow(delta_t, 2) * sin(psi) * nu_a);
      v_pred = v + delta_t * nu_a;
      psi_pred = psi + (psi_dot * delta_t) + (0.5 * std::pow(delta_t, 2) * nu_psi_dot_dot);
      psi_dot_pred = psi_dot + delta_t * nu_psi_dot_dot;
    }

    colOut << px_pred, py_pred, v_pred, psi_pred, psi_dot_pred;
    Xsig_pred_.col(c) = colOut;
  }
}

void UKF::PredictMeanAndCovariance()
{
  // predict state mean
  x_.setZero();
  for (int i = 0; i < n_x_; ++i)
  {
    for (int j = 0; j < 2 * n_aug_ + 1; ++j)
    {
      x_(i) = x_(i) + (weights_(j) * Xsig_pred_(i, j));
    }
  }

  // predict state covariance matrix
  P_.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    VectorXd delta = Xsig_pred_.col(i) - x_;

    // angle wrap around.
    while (delta(3) > M_PI)
      delta(3) -= 2.0 * M_PI;
    while (delta(3) < -M_PI)
      delta(3) += 2.0 * M_PI;

    P_ = P_ + weights_(i) * delta * delta.transpose();
  }
}

void UKF::PredictRadarMeasurement()
{
  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    const double px = Xsig_pred_(0, i);
    const double py = Xsig_pred_(1, i);
    const double v = Xsig_pred_(2, i);
    const double phi = Xsig_pred_(3, i);
    const double phi_dot = Xsig_pred_(4, i);

    const double range = sqrt(std::pow(px, 2) + std::pow(py, 2));
    const double psi = atan2(py, px);
    const double r_dot = (px * cos(phi) * v + py * sin(phi) * v) / range;

    Zsig_r_.col(i) << range, psi, r_dot;
  }

  // predict state mean (radar)
  z_pred_r_.setZero();
  for (int i = 0; i < n_z_R_; ++i)
  {
    for (int j = 0; j < 2 * n_aug_ + 1; ++j)
    {
      z_pred_r_(i) = z_pred_r_(i) + (weights_(j) * Zsig_r_(i, j));
    }
  }

  // calculate innovation covariance matrix S (radar)
  S_r_.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    VectorXd delta = Zsig_r_.col(i) - z_pred_r_;
    while (delta(1) > M_PI)
    {
      delta(1) -= 2.0 * M_PI;
    }
    while (delta(1) < -M_PI)
    {
      delta(1) += 2.0 * M_PI;
    }

    S_r_ = S_r_ + weights_(i) * delta * delta.transpose();
  }

  S_r_ = S_r_ + R_r;
}

void UKF::PredictLidarMeasurement()
{
  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  { 
    Zsig_l_(0, i) = Xsig_pred_(0, i);
    Zsig_l_(1, i) = Xsig_pred_(1, i);
  }

  // predict state mean (lidar)
  z_pred_l_.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) 
  {
    z_pred_l_ = z_pred_l_ + weights_(i) * Zsig_l_.col(i);
  }

  // calculate innovation covariance matrix S (lidar)
  S_l_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) 
  {                                               
    VectorXd delta = Zsig_l_.col(i) - z_pred_l_;
    S_l_ = S_l_ + weights_(i) * delta * delta.transpose();
  }

  // add measurement noise covariance matrix
  S_l_ = S_l_ + R_l;
}

void UKF::UpdateStateRadar(const MeasurementPackage &meas_package)
{
  // Extract the measurements
  VectorXd z = VectorXd(n_z_R_);
  z(0) = meas_package.raw_measurements_[0];
  z(1) = meas_package.raw_measurements_[1];
  z(2) = meas_package.raw_measurements_[2];

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_R_);

  // calculate cross correlation matrix
  VectorXd xTerm = VectorXd(n_x_);
  xTerm.setZero();
  VectorXd zTerm = VectorXd(n_z_R_);
  zTerm.setZero();
  Tc.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    zTerm = Zsig_r_.col(i) - z_pred_r_;
    xTerm = Xsig_pred_.col(i) - x_;

    while (zTerm(1) > M_PI) zTerm(1) -= 2.*M_PI;
    while (zTerm(1) < -M_PI) zTerm(1) += 2.*M_PI;

    while (xTerm(3) > M_PI) xTerm(3) -= 2.*M_PI;
    while (xTerm(3) < -M_PI) xTerm(3) += 2.*M_PI;

    Tc = Tc + weights_(i) * xTerm * zTerm.transpose();
  }

  // calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_x_, n_z_R_);
  K.setZero();
  K = Tc * S_r_.inverse();

  // update state mean and covariance matrix
  zTerm = z - z_pred_r_;

  // angle normalization
  while (zTerm(1) > M_PI) zTerm(1) -= 2.*M_PI;
  while (zTerm(1) < -M_PI) zTerm(1) += 2.*M_PI;

  x_ = x_ + K * zTerm;
  P_ = P_ - K * S_r_ * K.transpose();

  nis_radar = zTerm.transpose() * S_r_.inverse() * zTerm;
  std::ofstream nis_radar_file; 
  static bool first_radar = true;
  if (first_radar == true)
  {
    nis_radar_file.open("nis_radar.csv", std::ios::out | std::ios::trunc);
    first_radar = false;
  }
  else
    nis_radar_file.open("nis_radar.csv", std::ios::out | std::ios::ate | std::ios::app);

  nis_radar_file << nis_radar << std::endl;
  nis_radar_file.close();
}

void UKF::UpdateStateLidar(const MeasurementPackage &meas_package) 
{
  // Extract the measurements
  VectorXd z = VectorXd(n_z_L_);
  z(0) = meas_package.raw_measurements_[0];
  z(1) = meas_package.raw_measurements_[1];

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_L_);

  // calculate cross correlation matrix
  VectorXd xTerm = VectorXd(n_x_);
  xTerm.setZero();
  VectorXd zTerm = VectorXd(n_z_R_);
  zTerm.setZero();
  Tc.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {                                                
    zTerm = Zsig_l_.col(i) - z_pred_l_;

    // state difference
    xTerm = Xsig_pred_.col(i) - x_;

    // angle normalization
    while (xTerm(3)> M_PI) xTerm(3) -= 2.*M_PI;
    while (xTerm(3)<-M_PI) xTerm(3) += 2.*M_PI;

    Tc = Tc + weights_(i) * xTerm * zTerm.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S_l_.inverse();

  zTerm = z - z_pred_l_;

  // update state mean and covariance matrix
  x_ = x_ + K * zTerm;
  P_ = P_ - K*S_l_*K.transpose();

  nis_lidar = zTerm.transpose() * S_l_.inverse() * zTerm;
  static bool first_lidar = true;
  std::ofstream nis_lidar_file; 
  if (first_lidar == true)
  {
    nis_lidar_file.open("nis_lidar.csv", std::ios::out | std::ios::trunc);
    first_lidar = false;
  }
  else
    nis_lidar_file.open("nis_lidar.csv", std::ios::out | std::ios::ate | std::ios::app);
  
  nis_lidar_file << nis_lidar << std::endl;
  nis_lidar_file.close();
}

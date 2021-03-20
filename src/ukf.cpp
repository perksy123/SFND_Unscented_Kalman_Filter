
#include "stdafx.h"

#include "ukf.h"
#include "Eigen/Dense"

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
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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

    // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // set radar measurement dimension (r, phi and r_dot)
  n_z_R_ = 3;

  // set Laser measurement dimension (rx, ry)
  n_z_L_ = 2;

  lambda_ = 3. - n_aug_;

  // Compute the wights for calculating the sigma points
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  double weight = 0.5 / (lambda_ + n_aug_);
  weights_(0) = weight_0;

  for (int i = 1; i < 2 * n_aug_ + 1; ++i) {
    weights_(i) = weight;
  }

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
    }
    else // Process a rader measurement
    {
      // Initialise the state from the radar
      const double r = meas_package.raw_measurements_[0];
      const double phi = meas_package.raw_measurements_[1];
      const double r_dot = meas_package.raw_measurements_[2];


      const double px = r * cos(phi);
      const double py = r * sin(phi);
      const double vx = r_dot * cos(phi);
      const double vy = r_dot * sin(phi);
      const double v = sqrt(std::pow(vx, 2) + std::pow(vy, 2));
      const double yaw = atan2(vy, vx);

      x_ << px, py, v, yaw, 0;

    }

    is_initialized_ = true;

    // previous timestamp
    time_us_ = meas_package.timestamp_;

    // Initiase the state covariance matrix
    P_.setZero();
    P_ << std::pow(std_laspx_, 2), 0, 0, 0, 0,
      0, std::pow(std_laspy_, 2), 0, 0, 0,
      0, 0, std::pow(std_radrd_, 2), 0, 0,
      0, 0, 0, std::pow(std_radphi_, 2), 0,
      0, 0, 0, 0, std::pow(std_radrd_, 2);      // Not sure about the phi_dot noise....

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
    UpdateRadar(meas_package);
  }
  else
  {
    UpdateRadar(meas_package);
  }
}

void UKF::Prediction(double delta_t) {
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

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

  PredictLidarMeasurement();

  UpdateStateLidar(meas_package);

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
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
  /*    // set state dimension
      int n_x = 5;

      // set augmented dimension
      int n_aug = 7;

      // Process noise standard deviation longitudinal acceleration in m/s^2
      double std_a = 0.2;

      // Process noise standard deviation yaw acceleration in rad/s^2
      double std_yawdd = 0.2;

      // define spreading parameter
      double lambda = 3 - n_aug;

      // set example state
      VectorXd x = VectorXd(n_x);
      x << 5.7441,
          1.3800,
          2.2049,
          0.5015,
          0.3528;

      // set example covariance matrix
      MatrixXd P = MatrixXd(n_x, n_x);
      P << 0.0043, -0.0013, 0.0030, -0.0022, -0.0020,
          -0.0013, 0.0077, 0.0011, 0.0071, 0.0060,
          0.0030, 0.0011, 0.0054, 0.0007, 0.0008,
          -0.0022, 0.0071, 0.0007, 0.0098, 0.0100,
          -0.0020, 0.0060, 0.0008, 0.0100, 0.0123;

      // create augmented mean vector
      VectorXd x_aug = VectorXd(n_aug);

      // create augmented state covariance
      MatrixXd P_aug = MatrixXd(n_aug, n_aug);

      // create sigma point matrix
      MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);*/

      /**
      * Student part begin
      */

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

  // The first point 
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

  // print result
//    std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

  // write result
  Xsig_aug_ = Xsig_aug;
}

void UKF::SigmaPointPrediction(double delta_t)
{
  /*
      // set state dimension
      int n_x = 5;

      // set augmented dimension
      int n_aug = 7;

      // create example sigma point matrix
      MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
      Xsig_aug <<
          5.7441, 5.85768, 5.7441, 5.7441, 5.7441, 5.7441, 5.7441, 5.7441, 5.63052, 5.7441, 5.7441, 5.7441, 5.7441, 5.7441, 5.7441,
          1.38, 1.34566, 1.52806, 1.38, 1.38, 1.38, 1.38, 1.38, 1.41434, 1.23194, 1.38, 1.38, 1.38, 1.38, 1.38,
          2.2049, 2.28414, 2.24557, 2.29582, 2.2049, 2.2049, 2.2049, 2.2049, 2.12566, 2.16423, 2.11398, 2.2049, 2.2049, 2.2049, 2.2049,
          0.5015, 0.44339, 0.631886, 0.516923, 0.595227, 0.5015, 0.5015, 0.5015, 0.55961, 0.371114, 0.486077, 0.407773, 0.5015, 0.5015, 0.5015,
          0.3528, 0.299973, 0.462123, 0.376339, 0.48417, 0.418721, 0.3528, 0.3528, 0.405627, 0.243477, 0.329261, 0.22143, 0.286879, 0.3528, 0.3528,
          0, 0, 0, 0, 0, 0, 0.34641, 0, 0, 0, 0, 0, 0, -0.34641, 0,
          0, 0, 0, 0, 0, 0, 0, 0.34641, 0, 0, 0, 0, 0, 0, -0.34641;

      // create matrix with predicted sigma points as columns
      MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

      double delta_t = 0.1; // time diff in sec

    */

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
  /*
      // set state dimension
      int n_x = 5;

      // set augmented dimension
      int n_aug = 7;

      // define spreading parameter
      double lambda = 3 - n_aug;

      // create example matrix with predicted sigma points
      MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
      Xsig_pred <<
          5.9374, 6.0640, 5.925, 5.9436, 5.9266, 5.9374, 5.9389, 5.9374, 5.8106, 5.9457, 5.9310, 5.9465, 5.9374, 5.9359, 5.93744,
          1.48, 1.4436, 1.660, 1.4934, 1.5036, 1.48, 1.4868, 1.48, 1.5271, 1.3104, 1.4787, 1.4674, 1.48, 1.4851, 1.486,
          2.204, 2.2841, 2.2455, 2.2958, 2.204, 2.204, 2.2395, 2.204, 2.1256, 2.1642, 2.1139, 2.204, 2.204, 2.1702, 2.2049,
          0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337, 0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188, 0.5367, 0.535048,
          0.352, 0.29997, 0.46212, 0.37633, 0.4841, 0.41872, 0.352, 0.38744, 0.40562, 0.24347, 0.32926, 0.2214, 0.28687, 0.352, 0.318159;

      // create vector for weights
      VectorXd weights = VectorXd(2 * n_aug + 1);
  */
  // create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  // create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  /**
  * Student part begin
  */

  // predict state mean
  x.setZero();
  for (int i = 0; i < n_x_; ++i)
  {
    for (int j = 0; j < 2 * n_aug_ + 1; ++j)
    {
      x(i) = x(i) + (weights_(j) * Xsig_pred_(i, j));
    }
  }

  // predict state covariance matrix
  P.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    VectorXd delta = Xsig_pred_.col(i) - x;
    while (delta(3) > M_PI)
    {
      delta(3) -= 2.0 * M_PI;
    }
    while (delta(3) < -M_PI)
    {
      delta(3) += 2.0 * M_PI;
    }

    P = P + weights_(i) * delta * delta.transpose();
  }

  // write result
  x_ = x;
  P_ = P;
}

void UKF::PredictRadarMeasurement()
{
  /*
    // set state dimension
    int n_x = 5;

    // set augmented dimension
    int n_aug = 7;

    // set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;

    // define spreading parameter
    double lambda = 3 - n_aug;

    // set vector for weights
    VectorXd weights = VectorXd(2 * n_aug + 1);
    double weight_0 = lambda / (lambda + n_aug);
    double weight = 0.5 / (lambda + n_aug);
    weights(0) = weight_0;

    for (int i = 1; i<2 * n_aug + 1; ++i) {
        weights(i) = weight;
    }

    // radar measurement noise standard deviation radius in m
    double std_radr = 0.3;

    // radar measurement noise standard deviation angle in rad
    double std_radphi = 0.0175;

    // radar measurement noise standard deviation radius change in m/s
    double std_radrd = 0.1;

    // create example matrix with predicted sigma points
    MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
    Xsig_pred <<
        5.9374, 6.0640, 5.925, 5.9436, 5.9266, 5.9374, 5.9389, 5.9374, 5.8106, 5.9457, 5.9310, 5.9465, 5.9374, 5.9359, 5.93744,
        1.48, 1.4436, 1.660, 1.4934, 1.5036, 1.48, 1.4868, 1.48, 1.5271, 1.3104, 1.4787, 1.4674, 1.48, 1.4851, 1.486,
        2.204, 2.2841, 2.2455, 2.2958, 2.204, 2.204, 2.2395, 2.204, 2.1256, 2.1642, 2.1139, 2.204, 2.204, 2.1702, 2.2049,
        0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337, 0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188, 0.5367, 0.535048,
        0.352, 0.29997, 0.46212, 0.37633, 0.4841, 0.41872, 0.352, 0.38744, 0.40562, 0.24347, 0.32926, 0.2214, 0.28687, 0.352, 0.318159;*/

        /**
        * Student part begin
        */

        // transform sigma points into measurement space
  Zsig_r_.setZero();
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

  // predict state mean
  z_pred_r_.setZero();
  for (int i = 0; i < n_z_R_; ++i)
  {
    for (int j = 0; j < 2 * n_aug_ + 1; ++j)
    {
      z_pred_r_(i) = z_pred_r_(i) + (weights_(j) * Zsig_r_(i, j));
    }
  }


  // calculate innovation covariance matrix S
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
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    Zsig_l_(0, i) = p_x;
    Zsig_l_(1, i) = p_y;
  }

  // calculate the mean predicted measurement
  z_pred_l_.fill(0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) 
  {
    z_pred_l_ = z_pred_l_ + weights_(i) * Zsig_l_.col(i);
  }

  S_l_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) 
  {  
                                              
    VectorXd z_diff = Zsig_l_.col(i) - z_pred_l_;

    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

    S_l_ = S_l_ + weights_(i) * z_diff * z_diff.transpose();

  }
  // add measurement noise covariance matrix
  S_l_ = S_l_ + R_l;
}


void UKF::UpdateStateRadar(const MeasurementPackage &meas_package)
{
  /*
    // set state dimension
    int n_x = 5;

    // set augmented dimension
    int n_aug = 7;

    // set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;

    // define spreading parameter
    double lambda = 3 - n_aug;

    // set vector for weights
    VectorXd weights = VectorXd(2 * n_aug + 1);
    double weight_0 = lambda / (lambda + n_aug);
    double weight = 0.5 / (lambda + n_aug);
    weights(0) = weight_0;

    for (int i = 1; i<2 * n_aug + 1; ++i) {
        weights(i) = weight;
    }

    // create example matrix with predicted sigma points in state space
    MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
    Xsig_pred <<
        5.9374, 6.0640, 5.925, 5.9436, 5.9266, 5.9374, 5.9389, 5.9374, 5.8106, 5.9457, 5.9310, 5.9465, 5.9374, 5.9359, 5.93744,
        1.48, 1.4436, 1.660, 1.4934, 1.5036, 1.48, 1.4868, 1.48, 1.5271, 1.3104, 1.4787, 1.4674, 1.48, 1.4851, 1.486,
        2.204, 2.2841, 2.2455, 2.2958, 2.204, 2.204, 2.2395, 2.204, 2.1256, 2.1642, 2.1139, 2.204, 2.204, 2.1702, 2.2049,
        0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337, 0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188, 0.5367, 0.535048,
        0.352, 0.29997, 0.46212, 0.37633, 0.4841, 0.41872, 0.352, 0.38744, 0.40562, 0.24347, 0.32926, 0.2214, 0.28687, 0.352, 0.318159;

    // create example vector for predicted state mean
    VectorXd x = VectorXd(n_x);
    x <<
        5.93637,
        1.49035,
        2.20528,
        0.536853,
        0.353577;

    // create example matrix for predicted state covariance
    MatrixXd P = MatrixXd(n_x, n_x);
    P <<
        0.0054342, -0.002405, 0.0034157, -0.0034819, -0.00299378,
        -0.002405, 0.01084, 0.001492, 0.0098018, 0.00791091,
        0.0034157, 0.001492, 0.0058012, 0.00077863, 0.000792973,
        -0.0034819, 0.0098018, 0.00077863, 0.011923, 0.0112491,
        -0.0029937, 0.0079109, 0.00079297, 0.011249, 0.0126972;

    // create example matrix with sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);
    Zsig <<
        6.1190, 6.2334, 6.1531, 6.1283, 6.1143, 6.1190, 6.1221, 6.1190, 6.0079, 6.0883, 6.1125, 6.1248, 6.1190, 6.1188, 6.12057,
        0.24428, 0.2337, 0.27316, 0.24616, 0.24846, 0.24428, 0.24530, 0.24428, 0.25700, 0.21692, 0.24433, 0.24193, 0.24428, 0.24515, 0.245239,
        2.1104, 2.2188, 2.0639, 2.187, 2.0341, 2.1061, 2.1450, 2.1092, 2.0016, 2.129, 2.0346, 2.1651, 2.1145, 2.0786, 2.11295;

    // create example vector for mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred <<
        6.12155,
        0.245993,
        2.10313;

    // create example matrix for predicted measurement covariance
    MatrixXd S = MatrixXd(n_z, n_z);
    S <<
        0.0946171, -0.000139448, 0.00407016,
        -0.000139448, 0.000617548, -0.000770652,
        0.00407016, -0.000770652, 0.0180917;

    // create example vector for incoming radar measurement
    VectorXd z = VectorXd(n_z);
    z <<
        5.9214,   // rho in m
        0.2187,   // phi in rad
        2.0062;   // rho_dot in m/s
  */

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

    Tc = Tc + weights_(i) * xTerm * xTerm.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S_l_.inverse();

  zTerm = z - z_pred_l_;

//  double NIS = z_diff.transpose()*S_l_.inverse()*z_diff;
  //std::cout << "\nLASER NIS: " << NIS << std::endl;

  // update state mean and covariance matrix
  x_ = x_ + K * zTerm;
  P_ = P_ - K*S_l_*K.transpose();

}

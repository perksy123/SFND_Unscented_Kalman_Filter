#ifndef UKF_H
#define UKF_H

#include "Eigen/Dense"
#include "measurement_package.h"

class UKF {
 public:
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  /**
  * Student assignment functions
  */

  void AugmentedSigmaPoints();
  void SigmaPointPrediction(double delta_t);
  void PredictMeanAndCovariance();

  void PredictRadarMeasurement();
  void PredictLidarMeasurement();

  void UpdateStateRadar(const MeasurementPackage &meas_package);
  void UpdateStateLidar(const MeasurementPackage &meas_package);

  // initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // predicted sigma points matrix
  Eigen::MatrixXd Xsig_pred_;

  // time when the state is true, in us
  long long time_us_;

  // radar mean predicted measurement
  Eigen::VectorXd z_pred_r_;

  // laser mean predicted measurement
  Eigen::MatrixXd z_pred_l_;

  // radar measurement space sigma points
  Eigen::MatrixXd Zsig_r_;

  // laser measurement space sigma points
  Eigen::MatrixXd Zsig_l_;

  // radar measurement covariance matrix
  Eigen::MatrixXd S_r_;

  // laser measurement covariance matrix
  Eigen::MatrixXd S_l_;

  // Laser measurement noise matrix
  Eigen::MatrixXd R_l;

  // Radar measurement noise matrix
  Eigen::MatrixXd R_r;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  // Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  // Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  // Radar measurement noise standard deviation radius in m
  double std_radr_;

  // Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  // Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  // Weights of sigma points
  Eigen::VectorXd weights_;

  // State dimension
  int n_x_;

  // Augmented state dimension
  int n_aug_;

  // Radar measurement dimention
  int n_z_R_;

  // Laser measurement dimention
  int n_z_L_;

  // Sigma point spreading parameter
  double lambda_;

  // The process noise matrix
  Eigen::MatrixXd proc_noise_;

  // Augmented sigma point matrix
  Eigen::MatrixXd Xsig_aug_;

};

#endif  // UKF_H
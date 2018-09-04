#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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
  std_a_ = 1.02;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.52;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:
  Complete the initialization. See ukf.h for other member properties.
  Hint: one or more values initialized above might be wildly off...
  */

  //set state dimension and augmented dimension
  n_x_ = 5;
  n_aug_ = 7;

  //define spreading parameter (rule of thumb)
  lambda_ = 3 - n_aug_;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Initialization
  is_initialized_ = false;

  // Weights
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_/(lambda_+n_aug_);
  for (int ii=1; ii<2*n_aug_+1; ii++) {
    weights_(ii) = 0.5/(n_aug_+lambda_);
  }
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */


  if(not is_initialized_) {
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;
    x_ << 0.5, 0.5, 0.01, 0.01, 0.01;

    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      cout << "INITIALIZED BY LASER" << endl;
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);

    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // Extract radar measurements
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float slope = tan(phi);
      float speed = meas_package.raw_measurements_[2];

      // Convert to cartesian basis
      float px = abs(sqrt((rho*rho) / (1.0 + slope*slope)));
      float py = abs(px * slope);
      float vx = abs(sqrt((speed*speed) / (1.0 + slope*slope)));
      float vy = abs(vx * slope);

      // Sign correction (because of the quadratic equations above)
      if (phi < 0) {
        py = -py;
        vy = -vy;
      }
      if (phi > -M_PI/2.0 and phi < M_PI/2.0) {
        px = -px;
        vx = -vx;
      }
      float v = sqrt(vx*vx + vy*vy);
      x_(0) = px;
      x_(1) = py;
      x_(2) = v;
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    // Switch off radar/lidar
    //use_radar_ = false;
    return;
  }

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    double delta_t = (meas_package.timestamp_ - time_us_)/double(1e6);
    time_us_ = meas_package.timestamp_;
    Prediction(delta_t);
    UpdateRadar(meas_package);
  } else if((meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)) {
    double delta_t = (meas_package.timestamp_ - time_us_)/double(1e6);
    time_us_ = meas_package.timestamp_;
    Prediction(delta_t);
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  //************************* Generate sigma points ****************************
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);
  // Calculate square root of P
  MatrixXd A = P_.llt().matrixL();
  // Set first column (mean)
  Xsig.col(0) = x_;
  // Set the remaining sigma points
  double lambda = 3 - n_x_;
  for (int ii = 0; ii < n_x_; ii++)
  {
    Xsig.col(ii+1)      = x_ + sqrt(lambda+n_x_) * A.col(ii);
    Xsig.col(ii+1+n_x_) = x_ - sqrt(lambda+n_x_) * A.col(ii);
  }

  //************************* UKF Augmentation *********************************
  // Create augmented mean state, state covarince and sigma point matrices
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_+1) = 0;

  // Create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_+1, n_x_+1) = std_yawdd_ * std_yawdd_;

  // Create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // Create augmented sigma point matrix
  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int ii = 0; ii< n_aug_; ii++)
  {
    Xsig_aug.col(ii+1)        = x_aug + sqrt(lambda_+n_aug_) * L.col(ii);
    Xsig_aug.col(ii+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(ii);
  }

  //*************************** Predict sigma points ***************************
  for (int ii = 0; ii< 2*n_aug_+1; ii++)
  {
    // Extract values for better readability
    double p_x = Xsig_aug(0,ii);
    double p_y = Xsig_aug(1,ii);
    double v = Xsig_aug(2,ii);
    double yaw = Xsig_aug(3,ii);
    double yawd = Xsig_aug(4,ii);
    double nu_a = Xsig_aug(5,ii);
    double nu_yawdd = Xsig_aug(6,ii);

    // Predicted state values
    double px_p, py_p;

    // Avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // Add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;
    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // Write predicted sigma point into right column
    Xsig_pred_(0,ii) = px_p;
    Xsig_pred_(1,ii) = py_p;
    Xsig_pred_(2,ii) = v_p;
    Xsig_pred_(3,ii) = yaw_p;
    Xsig_pred_(4,ii) = yawd_p;
  }

  // ******************** Predict mean (state) and covariance ******************
  // Predict mean
  x_.fill(0.0);
  for (int ii = 0; ii < 2 * n_aug_ + 1; ii++) {  //iterate over sigma points
    x_ = x_ + weights_(ii) * Xsig_pred_.col(ii);
  }

  // Predicted state covariance matrix
  P_.fill(0.0);
  for (int ii = 0; ii < 2 * n_aug_ + 1; ii++) {  //iterate over sigma points

    // State difference
    VectorXd x_diff = Xsig_pred_.col(ii) - x_;
    // Angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(ii) * x_diff * x_diff.transpose() ;
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:
  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.
  You'll also need to calculate the lidar NIS.
  */

  int n_z = 2;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // Transform sigma points into measurement space
  for (int ii = 0; ii < 2 * n_aug_ + 1; ii++) {

    // Extract values for better readibility
    double p_x = Xsig_pred_(0,ii);
    double p_y = Xsig_pred_(1,ii);

    // Measurement model
    Zsig(0,ii) = p_x;                        //x
    Zsig(1,ii) = p_y;                        //y
  }

  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int ii=0; ii < 2*n_aug_+1; ii++) {
      z_pred = z_pred + weights_(ii) * Zsig.col(ii);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int ii = 0; ii < 2 * n_aug_ + 1; ii++) {
    VectorXd z_diff = Zsig.col(ii) - z_pred;

    // Angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(ii) * z_diff * z_diff.transpose();
  }

  // Add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;
  S = S + R;

  // Create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // Calculate cross correlation matrix
  Tc.fill(0.0);
  for (int ii = 0; ii < 2 * n_aug_ + 1; ii++) {
    VectorXd z_diff = Zsig.col(ii) - z_pred;

    // State difference
    VectorXd x_diff = Xsig_pred_.col(ii) - x_;
    // Angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(ii) * x_diff * z_diff.transpose();
  }

  // Kalman gain K
  MatrixXd K = Tc * S.inverse();

  // Residual
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;

  // Update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  // Normalized Innovation Squared
  // cout << "NIS RADAR: " << z_diff.transpose() * S.inverse() * z_diff << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.
  You'll also need to calculate the radar NIS.
  */

  // Create matrices for sigma points in measurement space
  int n_z = 3;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // Transform sigma points into measurement space
  for (int ii = 0; ii < 2 * n_aug_ + 1; ii++) {

    // Extract values for better readibility
    double p_x = Xsig_pred_(0,ii);
    double p_y = Xsig_pred_(1,ii);
    double v   = Xsig_pred_(2,ii);
    double yaw = Xsig_pred_(3,ii);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // Measurement model
    Zsig(0,ii) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,ii) = atan2(p_y,p_x);                                 //phi
    Zsig(2,ii) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int ii=0; ii < 2*n_aug_+1; ii++) {
      z_pred = z_pred + weights_(ii) * Zsig.col(ii);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int ii = 0; ii < 2 * n_aug_ + 1; ii++) {
    VectorXd z_diff = Zsig.col(ii) - z_pred;

    // Angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(ii) * z_diff * z_diff.transpose();
  }

  // Add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R << std_radr_*std_radr_, 0, 0,
       0, std_radphi_*std_radphi_, 0,
       0, 0,std_radrd_*std_radrd_;
  S = S + R;


  // Create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // Calculate cross correlation matrix
  Tc.fill(0.0);
  for (int ii = 0; ii < 2 * n_aug_ + 1; ii++) {
    VectorXd z_diff = Zsig.col(ii) - z_pred;
    // Angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // State difference
    VectorXd x_diff = Xsig_pred_.col(ii) - x_;
    // Angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(ii) * x_diff * z_diff.transpose();
  }

  // Kalman gain K
  MatrixXd K = Tc * S.inverse();

  // Residual
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;

  // Angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // Update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  // Normalized Innovation Squared
  // cout << "NIS RADAR: " << z_diff.transpose() * S.inverse() * z_diff << endl;
}

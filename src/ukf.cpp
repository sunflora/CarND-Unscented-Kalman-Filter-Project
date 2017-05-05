#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  is_debug_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  // std_a_ = 30;
   std_a_ = 1;  // in quiz: 0.2


  // Process noise standard deviation yaw acceleration in rad/s^2
  //std_yawdd_ = 30;
   std_yawdd_ = 1; // in quiz: 0.2

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
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  is_initialized_ = false;
  time_us_ = 0;

  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  ///* Sigma point spreading parameter
  lambda_ = 3.0 - n_aug_;

  ///* Weights of sigma points
  //set vector for weights
  weights_ = VectorXd(2*n_aug_+1);

  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  ///* the current NIS for radar
  NIS_radar_ = 0.0;

  ///* the current NIS for laser
  NIS_laser_ = 0.0;

  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

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

    if (is_debug_) std::cout << "In ProcessMeasurement" << std::endl;

    /*****************************************************************************
    *  Initialization
    ****************************************************************************/
    if (!is_initialized_) {

        if (is_debug_) std::cout << "Begin Initialization." << std::endl;

        P_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0,
              0, 0, 1, 0, 0,
              0, 0, 0, 1, 0,
              0, 0, 0, 0, 1;

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */

            //measurement matrix

            double rho = meas_package.raw_measurements_[0];
            double phi = meas_package.raw_measurements_[1];
            double rho_dot = meas_package.raw_measurements_[2];

            double px = rho * cos(phi);
            double py = rho * sin(phi);

            double vx = rho_dot * cos(phi); // 0.0;
            double vy = rho_dot * sin(phi);

            x_ << px, py, vx, vy, 0.0;

        } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            /**
            Initialize state.
            */

            double px = meas_package.raw_measurements_[0];
            double py = meas_package.raw_measurements_[1];
            double vx = 0.0;
            double vy = 0.0;

            x_ << px, py, vx, vy, 0.0;

        }

        time_us_ = meas_package.timestamp_;

        // done initializing, no need to predict or update
        is_initialized_ = true;

        if (is_debug_) cout << "Initialized!" << endl;
        if (is_debug_) cout << "x_: " << endl;
        if (is_debug_) cout << x_ << endl;

        return;
    }

    if ((use_laser_ == false) && (meas_package.sensor_type_ == MeasurementPackage::LASER)) {
        return;
    }

    if ((use_radar_ == false) && (meas_package.sensor_type_ == MeasurementPackage::RADAR)) {
        return;
    }

    double dt = (double) ((double) meas_package.timestamp_ - (double) time_us_) / 1000000.0;

    time_us_ = meas_package.timestamp_;


    if (dt > 0.1) {
        dt = dt - round(dt / 0.05);
    }

    Prediction(dt);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        UpdateRadar(meas_package);
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
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

    //create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

    AugmentedSigmaPoints(Xsig_aug);
    SigmaPointPrediction(Xsig_aug, delta_t);
    PredictMeanAndCovariance();

}

void UKF::AugmentedSigmaPoints(MatrixXd &Xsig_aug) {

    //create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_);

    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

    //create augmented mean state
    x_aug.head(n_x_) = x_;
    x_aug(n_x_) = 0;
    x_aug(n_x_+1) = 0;

    //create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug(n_x_, n_x_) = std_a_ * std_a_;
    P_aug(n_x_+1, n_x_+1) = std_yawdd_ * std_yawdd_;

    //create square root matrix
    MatrixXd L = P_aug.llt().matrixL();

    //create augmented sigma points
    Xsig_aug.col(0) = x_aug;
    for (int i = 0; i < n_aug_; i++) {
        Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
        Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
    }

    //print result
    if (is_debug_) std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

}

void UKF::SigmaPointPrediction(MatrixXd &Xsig_aug, double delta_t) {

    //predict sigma points
    for (int i = 0; i< 2*n_aug_+1; i++)
    {
        //extract values for better readability
        double p_x = Xsig_aug(0,i);
        double p_y = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double yawd = Xsig_aug(4,i);
        double nu_a = Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);

        //predicted state values
        double px_p, py_p;

        //avoid division by zero
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

        //add noise
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;

        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;

        //write predicted sigma point into right column
        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;
        Xsig_pred_(2,i) = v_p;
        Xsig_pred_(3,i) = yaw_p;
        Xsig_pred_(4,i) = yawd_p;
    }

    //print result
    if (is_debug_) std::cout << "Xsig_pred_ = " << std::endl << Xsig_pred_ << std::endl;

}

double UKF::NormalizeAngle(double pValue) {
    if (fabs(pValue) > M_PI) {
        pValue -= round(pValue / (2.0 * M_PI)) * (2.0 * M_PI);
    }
    return pValue;
}

void UKF::PredictMeanAndCovariance() {


    //predicted state mean
    x_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        x_ = x_+ weights_(i) * Xsig_pred_.col(i);
    }
    if (is_debug_) std::cout << "predicted state mean x_ is: " << std::endl;
    if (is_debug_) std::cout << x_ << std::endl;

    //predicted state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        if (is_debug_) std::cout << "x_diff(3): " << x_diff(3) <<  std::endl;

        //while (x_diff(3)>  M_PI) x_diff(3)-=2.*M_PI;
        //while (x_diff(3)< -M_PI) x_diff(3)+=2.*M_PI;

        NormalizeAngle(x_diff(3));

        P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
    }

    if (is_debug_) std::cout << "state covariance matrix P_ is " << std::endl;
    if (is_debug_) std::cout << P_ << std::endl;

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

    //set measurement dimension, radar can measure r, phi, and r_dot
    // int n_z = 3;

    int n_z = meas_package.raw_measurements_.rows();


    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

        // extract values for better readibility
        double p_x = Xsig_pred_(0, i);
        double p_y = Xsig_pred_(1, i);
        double v = Xsig_pred_(2, i);
        double yaw = Xsig_pred_(3, i);

        double v1 = cos(yaw) * v;
        double v2 = sin(yaw) * v;

        // measurement model
        Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);                            //r
        Zsig(1, i) = atan2(p_y, p_x);                                        //phi

        if (Zsig(0, i) == 0) {
            Zsig(2, i) = 0;
        } else {
            Zsig(2, i) = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y);    //r_dot
        }
    }

    if (is_debug_) std::cout << "Zsig: " << std::endl << Zsig << std::endl;

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        //angle normalization
        //while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
        //while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

        NormalizeAngle(z_diff(1));

        S = S + weights_(i) * z_diff * z_diff.transpose();
    }

    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z, n_z);
    R << std_radr_ * std_radr_, 0, 0,
            0, std_radphi_ * std_radphi_, 0,
            0, 0, std_radrd_ * std_radrd_;
    S = S + R;

    UpdateState(z_pred,Zsig,S,meas_package.raw_measurements_);

    VectorXd z_diff = meas_package.raw_measurements_-z_pred;

    //while (z_diff(2)> M_PI) z_diff(2)-=2.*M_PI;
    //while (z_diff(2)<-M_PI) z_diff(2)+=2.*M_PI;

    NormalizeAngle(z_diff(2));

    NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;



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

    //set measurement dimension, lidar can measure px and py
    // int n_z = 2;
    int n_z = meas_package.raw_measurements_.rows();

    VectorXd z_pred = VectorXd(n_z);
    MatrixXd S = MatrixXd(n_z, n_z);

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

        // extract values for better readibility
        double p_x = Xsig_pred_(0, i);
        double p_y = Xsig_pred_(1, i);

        // measurement model
        Zsig(0, i) = p_x;
        Zsig(1, i) = p_y;

    }

    //mean predicted measurement

    z_pred.fill(0.0);

    for (int i=0; i < 2*n_aug_ + 1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    //measurement covariance matrix S

    S.fill(0.0);

    for (int i = 0; i < 2*n_aug_ + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        //angle normalization
        //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        NormalizeAngle(z_diff(1));

        S = S + weights_(i) * z_diff * z_diff.transpose();
    }

    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z, n_z);
    R << std_laspx_ * std_laspx_, 0,
            0, std_laspy_ * std_laspy_;

    S = S + R;

    UpdateState(z_pred,Zsig,S,meas_package.raw_measurements_);

    // calculate NIS
    VectorXd z_diff = meas_package.raw_measurements_-z_pred;
    NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

}

void UKF::UpdateState(VectorXd &z_pred,  MatrixXd &Zsig, MatrixXd &S, VectorXd &z) {

    //set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = z_pred.rows();

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        NormalizeAngle(z_diff(1));

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        //while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        //while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        NormalizeAngle(x_diff(3));

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    //residual
    VectorXd z_diff = z - z_pred;

    //angle normalization
    //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    NormalizeAngle(z_diff(1));

    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();

};



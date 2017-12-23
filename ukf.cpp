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
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

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
  n_x_ = 5;
  n_aug_ = 7;
  P_ = MatrixXd(n_x_, n_x_);
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);
  time_us_ = 1;
  NIS_radar_ = 0.5;
  NIS_laser_ = 0.5;
  is_initialized_ = false;


  lambda_ = 3 - n_aug_;
  weights_ = VectorXd(2*n_aug_ + 1);
  // set weights
  double weight_0_ = lambda_/(lambda_ + n_aug_);
  weights_(0) = weight_0_;
  for (int i=1;  i < 2*n_aug_ + 1; i++) {  //2n+1 weights
    double weight_ = 0.5/(n_aug_ + lambda_);
    weights_(i) = weight_;
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
  omplete this function! Make sure you switch between lidar and radar
  measurements.
  */

  float p_x_, p_y_, rho_, theta_;
  // double t0 = 1.47701e15; // beginning of time to make time calculations easier
  float delta_t_;

  if (!is_initialized_) {

   // Initialization of P
    P_ << 1, 0, 0, 0, 0,
    	  0, 1, 0, 0, 0,
		  0, 0, 1, 0, 0,
		  0, 0, 0, 1, 0,
		  0, 0, 0, 0, 0.2;

	// Initialization of the state
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
    	rho_ = meas_package.raw_measurements_(0);
    	theta_ = meas_package.raw_measurements_(1);

    	p_x_ = rho_*cos(theta_);
    	p_y_ = rho_*sin(theta_);

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
    	p_x_ = meas_package.raw_measurements_(0);
    	p_y_ = meas_package.raw_measurements_(1);
    }
    // done initializing, no need to predict or update
    x_ << p_x_, p_y_, 0, 0, 0;
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
  }

  delta_t_ = (meas_package.timestamp_ - time_us_)/1.0e6;	//dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  // Prediction Step
  Prediction(delta_t_);

  // Update Step
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR and use_radar_){
	  UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER and use_laser_){
	  UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t_) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // *******************************************//
  // Generating Sigma points
  MatrixXd Xsig_ = MatrixXd(n_x_, 2 *n_x_ + 1);
  MatrixXd A_ = P_.llt().matrixL();
  lambda_ = 3 - n_x_;
  Xsig_.col(0)  = x_;
  //set remaining sigma points
  for (int i = 0; i < n_x_; i++)
  {
  Xsig_.col(i+1)      = x_ + sqrt(lambda_ + n_x_) * A_.col(i);
  Xsig_.col(i+1+n_x_) = x_ - sqrt(lambda_ + n_x_) * A_.col(i);
  }

  // *******************************************//
  // Augmentation matrix
  VectorXd x_aug_ = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  lambda_ = 3 - n_aug_;
  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;


  //Create augmented covariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = std_a_ * std_a_;
  P_aug_(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug_.llt().matrixL();

  //create augmented sigma points
  Xsig_aug_.col(0)  = x_aug_;

  float sqrt_lambda_naug = sqrt(lambda_ + n_aug_);

  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug_.col(i+1)        = x_aug_ + sqrt_lambda_naug*L.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt_lambda_naug*L.col(i);
  }

  // *******************************************//
  // Predict Sigma Points
//  MatrixXd Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);

  for (int i = 0; i< 2*n_aug_ + 1; i++)
   {
     //extract values for better readability
     const double p_x_ 	= Xsig_aug_(0,i);
     const double p_y_ 	= Xsig_aug_(1,i);
     const double v_ 	= Xsig_aug_(2,i);
     const double yaw_ 	= Xsig_aug_(3,i);
     const double yawd_ = Xsig_aug_(4,i);
     const double nu_a_ = Xsig_aug_(5,i);
     const double nu_yawdd_ = Xsig_aug_(6,i);

     //predicted state values
     double px_p_, py_p_;

     //avoid division by zero
     if (fabs(yawd_) > 0.001) {
         px_p_ = p_x_ + v_/yawd_*(sin(yaw_ + yawd_*delta_t_) - sin(yaw_));
         py_p_ = p_y_ + v_/yawd_*(cos(yaw_) - cos(yaw_ + yawd_*delta_t_));
     }
     else {
         px_p_ = p_x_ + v_*delta_t_*cos(yaw_);
         py_p_ = p_y_ + v_*delta_t_*sin(yaw_);
     }

     double v_p_ = v_;
     double yaw_p_ = yaw_ + yawd_*delta_t_;
     double yawd_p_ = yawd_;

     //add noise
     px_p_ = px_p_ + 0.5*nu_a_*delta_t_*delta_t_ * cos(yaw_);
     py_p_ = py_p_ + 0.5*nu_a_*delta_t_*delta_t_ * sin(yaw_);
     v_p_ = v_p_ + nu_a_*delta_t_;

     yaw_p_ = yaw_p_ + 0.5*nu_yawdd_*delta_t_*delta_t_;
     yawd_p_ = yawd_p_ + nu_yawdd_*delta_t_;

     //write predicted sigma point into right column
     Xsig_pred_(0,i) = px_p_;
     Xsig_pred_(1,i) = py_p_;
     Xsig_pred_(2,i) = v_p_;
     Xsig_pred_(3,i) = yaw_p_;
     Xsig_pred_(4,i) = yawd_p_;

   }

  // *******************************************//
  // Predicting mean and covariance

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_+ weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff_ = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff_(3)> M_PI) x_diff_(3)-=2.*M_PI;
    while (x_diff_(3)<-M_PI) x_diff_(3)+=2.*M_PI;
    P_ = P_ + weights_(i) * x_diff_ * x_diff_.transpose() ;
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
	// update the state by using Kalman Filter equations

	int n_z_ = 2;

	// read the laser measurement data
	VectorXd z_(n_z_);
	z_ << meas_package.raw_measurements_(0), //p_x
			meas_package.raw_measurements_(1); // p_y

	MatrixXd H_ = MatrixXd(2,5);
	H_ << 1,0,0,0,0,
		  0,1,0,0,0;

	MatrixXd R_ = MatrixXd(n_z_, n_z_);
	R_ << std_laspx_*std_laspx_, 0,
		 0, std_laspy_*std_laspy_;

	MatrixXd I_ = MatrixXd::Identity(n_x_, n_x_);
	MatrixXd Ht_ = H_.transpose();

	VectorXd z_pred_ = H_ * x_;
	VectorXd y = z_ - z_pred_;
	MatrixXd S_ = H_ * P_ * Ht_ + R_;
	MatrixXd Si_ = S_.inverse();
	MatrixXd PHt_ = P_ * Ht_;
	MatrixXd K_ = PHt_ * Si_;

	//new estimate
	x_ = x_ + (K_ * y);
	P_ = (I_ - K_ * H_) * P_;
    // Calculate the NIS radar
    NIS_laser_ = (z_-z_pred_).transpose()*S_.inverse()*(z_-z_pred_);
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

  int n_z_ = 3;

  // read the radar measurement data
  VectorXd z_(n_z_);

  // read the measurement data
  z_ << meas_package.raw_measurements_(0), //r
		meas_package.raw_measurements_(1), //phi
		meas_package.raw_measurements_(2); // r_dot

  //transform sigma points into measurement space

  //create matrix for sigma points in measurement space
  MatrixXd Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred_ = VectorXd(n_z_);

  //measurement covariance matrix S
  MatrixXd S_ = MatrixXd(n_z_, n_z_);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

  // extract values for better readibility
    double p_x_ = Xsig_pred_(0,i);
    double p_y_ = Xsig_pred_(1,i);
    double v_ = Xsig_pred_(2,i);
    double yaw_ = Xsig_pred_(3,i);

    double v1_ = cos(yaw_)*v_;
    double v2_ = sin(yaw_)*v_;

    // measurement model
    Zsig_(0,i) = sqrt(p_x_*p_x_ + p_y_*p_y_);                        //r
	if (abs(p_x_) < 0.001){
		Zsig_(1,i) = 0;}
	else {
		Zsig_(1,i) = atan2(p_y_,p_x_);
	}	                                  //phi
    Zsig_(2,i) = (p_x_*v1_ + p_y_*v2_) / max(0.001, sqrt(p_x_*p_x_ + p_y_*p_y_));//r_dot
    }
    // mean predicted measurement
    z_pred_.fill(0.0);
    for (int i=0; i < 2*n_aug_ + 1; i++) {
      z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
    }

    //measurement covariance matrix S
    S_.fill(0.0);
    for (int i = 0; i < 2*n_aug_ + 1; i++) {  //2n+1 simga points
      //residual
      VectorXd z_diff_ = Zsig_.col(i) - z_pred_;

      //angle normalization
      z_diff_(1) = fmod(z_diff_(1),2.*M_PI );
      if(z_diff_(1) > M_PI){
      	z_diff_(1) -= 2.*M_PI;
      }else if(z_diff_(1) < -M_PI){
      	z_diff_(1) += 2.*M_PI;
      }

      S_ = S_ + weights_(i)*z_diff_*z_diff_.transpose();
    }
    //add measurement noise covariance matrix
    MatrixXd R_ = MatrixXd(n_z_, n_z_);
    R_ << std_radr_*std_radr_, 0, 0,
         0, std_radphi_*std_radphi_, 0,
         0, 0,std_radrd_*std_radrd_;
    S_ = S_ + R_;

  // UKF Update
  MatrixXd Tc_ = MatrixXd(n_x_, n_z_);
  Tc_.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff_ = Zsig_.col(i) - z_pred_;
    //angle normalization
    z_diff_(1) = fmod(z_diff_(1),2.*M_PI );
    if(z_diff_(1) > M_PI){
    	z_diff_(1) -= 2.*M_PI;
    }else if(z_diff_(1) < -M_PI){
    	z_diff_(1) += 2.*M_PI;
    }

    // state difference
    VectorXd x_diff_ = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff_(3) = fmod(x_diff_(3),2.*M_PI );
    if(x_diff_(3) > M_PI){
    	x_diff_(3) -= 2.*M_PI;
    }else if(z_diff_(1) < -M_PI){
    	x_diff_(3) += 2.*M_PI;
    }

    Tc_ = Tc_ + weights_(i)*x_diff_*z_diff_.transpose();
    }
    //Kalman gain K;
    MatrixXd K_ = Tc_*S_.inverse();

    //residual
    VectorXd z_diff_ = z_ - z_pred_;
    //angle normalization

    z_diff_(1) = fmod(z_diff_(1),2.*M_PI );
    if(z_diff_(1) > M_PI){
    	z_diff_(1) -= 2.*M_PI;
    }else if(z_diff_(1) < -M_PI){
    	z_diff_(1) += 2.*M_PI;
    }

    //update state mean and covariance matrix
    x_ = x_ + K_*z_diff_;
    P_ = P_ - K_*S_*K_.transpose();

    // Calculate the NIS radar
    NIS_radar_ = (z_-z_pred_).transpose()*S_.inverse()*(z_-z_pred_);
}




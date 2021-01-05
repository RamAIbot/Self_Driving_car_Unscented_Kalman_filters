#include "ukf.h"
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
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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
 n_x_= 5;
 n_aug_ = 7;
 lambda_ = 3 - n_x_;
 is_initialized_ = false;
 time_us_ = 0.0;
 
 weights_ = VectorXd(2*n_aug_+1);
 
 
 weights_(0) = lambda_ / (lambda_ + n_aug_);
 for(int i=1;i<=2*n_aug_;i++)
 {
 	weights_(i) = 1 / (2 * (lambda_ + n_aug_));
 }
   
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
  //Augmented sigma points matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  
  //Augmented mean vector
  VectorXd x_aug = VectorXd(7);
  //Augmented state Covariance
  MatrixXd P_aug = MatrixXd(7,7);
  
  
}

UKF::~UKF() {}



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
  //Creating augmented mean state
  VectorXd x_aug = VectorXd(7);
  MatrixXd P_aug = MatrixXd(7,7);
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
  //Augmented sigma points matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  weights_ = VectorXd(2*n_aug_+1);
 
 
 weights_(0) = lambda_ / (lambda_ + n_aug_);
 for(int i=1;i<=2*n_aug_;i++)
 {
 	weights_(i) = 1 / (2 * (lambda_ + n_aug_));
 }
 
 
 
  
  cout<<"inside prediction"<<endl;
  x_aug.fill(0.0);
  x_aug.head(5) = x_;
  
  		   
  //cout<<"2"<<endl;
  x_aug(5) = 0;
  
  x_aug(6) = 0;
  
  //Creating augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;
  //cout<<"first part -2"<<endl;
  //Square root of P_aug matrix
  MatrixXd L = P_aug.llt().matrixL();
  
  //Creating Augmented Sigma points
  Xsig_aug.col(0) = x_aug; 
  
  for(int i=0;i<n_aug_;i++)
  {
  	Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
  	Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  
  //Predicting Sigma Points
  
  for(int i=0;i<2*n_aug_+1;i++)
  {
  	double p_x = Xsig_aug(0,i);
  	double p_y = Xsig_aug(1,i);
  	double v = Xsig_aug(2,i);
  	double yaw = Xsig_aug(3,i);
  	double yawd = Xsig_aug(4,i);
  	double nu_a = Xsig_aug(5,i);
  	double nu_yawdd = Xsig_aug(6,i);
  	
  	double px_p,py_p;
  	if(fabs(yawd) > 0.001)
  	{
  		px_p = p_x + v/yawd * (sin (yaw + yawd*delta_t) - sin(yaw));
  		py_p = p_y + v/yawd * (cos(yaw) - cos(yaw + yawd*delta_t));
	}
	else
	{
		px_p = p_x + v * cos(yaw) * delta_t;
		py_p = p_y + v * sin(yaw) * delta_t;	
	}
	double v_p = v;
	double yaw_p = yaw + yawd * delta_t;
	double yawd_p = yawd;
	
	px_p = px_p + 0.5 * delta_t * delta_t * nu_a * cos(yaw);
	py_p = py_p + 0.5 *delta_t * delta_t * nu_a * sin(yaw);
	
	v_p = v_p + delta_t * nu_a;
	
	yaw_p = yaw_p + 0.5 * delta_t * delta_t * nu_yawdd;
	
	yawd_p = yawd_p + delta_t * nu_yawdd;
	
	Xsig_pred_(0,i) = px_p;
	Xsig_pred_(1,i) = py_p;
	Xsig_pred_(2,i) = v_p;
	Xsig_pred_(3,i) = yaw_p;
	Xsig_pred_(4,i) = yawd_p;
  }
  
  //Predict State mean
  x_.fill(0.0);
  
  for(int i=0;i< 2*n_aug_+1;i++)
  {
  	x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
  
  //Predicting Covariance matrix
  P_.fill(0.0);
  for(int i=0;i< 2*n_aug_+1;i++)
  {
  	VectorXd x_diff = Xsig_pred_.col(i) - x_;
  	
  	while(x_diff(3) > M_PI) x_diff(3) -= 2. *M_PI;
  	while(x_diff(3) < -M_PI) x_diff(3) += 2. *M_PI;
  	
  	P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
  cout<<"inside prediction ends"<<endl;
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
  cout<<"Inside Update Lidar"<<endl;
   int n_z_ = 2;
  //Matrix of Sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_ , 2*n_aug_+1);
  //Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  //Measurement Covariance Matrix
  MatrixXd S = MatrixXd(n_z_,n_z_);
  
  VectorXd z = VectorXd(n_z_);
  
  z << meas_package.raw_measurements_(0),
  	   meas_package.raw_measurements_(1);
 //cout<<"Lidar update 1"<<endl;
  for(int i=0;i< 2*n_aug_+1 ; i++)
  {
  	double p_x = Xsig_pred_(0,i);
	double p_y = Xsig_pred_(1,i);
	
	
	
	
	Zsig(0,i) = p_x;
	Zsig(1,i) = p_y;
	
	
  }
  //cout<<"Lidar update 2"<<endl;
  z_pred.fill(0.0);
  
  for(int i=0;i <2*n_aug_+1 ; i++)
  {
  	z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  S.fill(0.0);
  for(int i=0;i<2*n_aug_+1 ;i++)
  {
  	VectorXd z_diff = Zsig.col(i) - z_pred;
//  	while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
//    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
  	S = S + weights_(i) * z_diff * z_diff.transpose();
  }  
  //cout<<"Lidar update 3"<<endl;
  MatrixXd R = MatrixXd(n_z_, n_z_);
  R << std_laspx_*std_laspx_ , 0 ,
  		0 , std_laspy_*std_laspy_ ;
  		
  		
  S = S + R;
 // cout<<"Lidar update 3"<<endl;
  MatrixXd Tc = MatrixXd(n_x_,n_z_);
  
  Tc.fill(0.0);
  
  for(int i=0;i< 2*n_aug_+1 ; i++)
  {
  	VectorXd z_diff = Zsig.col(i) - z_pred;
  	
//  	while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
//    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  	
  	VectorXd x_diff = Xsig_pred_.col(i) - x_;
  	
  	while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    
    
  	Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  
  MatrixXd K = Tc * S.inverse();
  
  VectorXd z_diff = z - z_pred;
  
//  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
//  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
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
  cout<<"Inside Update RADAR"<<endl;
  int n_z_ = 3;
  //Matrix of Sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_ , 2*n_aug_+1);
  //Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  //Measurement Covariance Matrix
  MatrixXd S = MatrixXd(n_z_,n_z_);
  
  VectorXd z = VectorXd(n_z_);
  
  z << meas_package.raw_measurements_(0),
  	 meas_package.raw_measurements_(1),
  	 meas_package.raw_measurements_(2);
 
  for(int i=0;i< 2*n_aug_+1 ; i++)
  {
  	double p_x = Xsig_pred_(0,i);
	double p_y = Xsig_pred_(1,i);
	double v = Xsig_pred_(2,i);
	double yaw = Xsig_pred_(3,i);
	
	double v1 = v * cos(yaw);
	double v2 = v * sin(yaw);
	
	double rho = sqrt(p_x * p_x + p_y* p_y);
	double phi = atan2(p_y,p_x);
	double rho_dot = (p_x*v1 + p_y*v2) / (sqrt(p_x * p_x + p_y* p_y));
	
	Zsig(0,i) = rho;
	Zsig(1,i) = phi;
	Zsig(2,i) = rho_dot;
	
  }
  
  z_pred.fill(0.0);
  
  for(int i=0;i <2*n_aug_+1 ; i++)
  {
  	z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  S.fill(0.0);
  for(int i=0;i<2*n_aug_+1 ;i++)
  {
  	VectorXd z_diff = Zsig.col(i) - z_pred;
  	while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
  	S = S + weights_(i) * z_diff * z_diff.transpose();
  }  
  
  MatrixXd R = MatrixXd(n_z_, n_z_);
  R << std_radr_*std_radr_ , 0 , 0,
  		0 , std_radphi_*std_radphi_ , 0,
  		0 , 0 , std_radrd_*std_radrd_;
  		
  S = S + R;
  
  MatrixXd Tc = MatrixXd(n_x_,n_z_);
  
  Tc.fill(0.0);
  
  for(int i=0;i< 2*n_aug_+1 ; i++)
  {
  	VectorXd z_diff = Zsig.col(i) - z_pred;
  	
  	while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  	
  	VectorXd x_diff = Xsig_pred_.col(i) - x_;
  	
  	while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    
    
  	Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  
  MatrixXd K = Tc * S.inverse();
  
  VectorXd z_diff = z - z_pred;
  
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
	
  
}



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
  
  if(((meas_package.sensor_type_ == MeasurementPackage::RADAR) && (use_radar_)) || ((meas_package.sensor_type_ == MeasurementPackage::LASER) && (use_laser_)))
  {
  	if(!is_initialized_)
  	{
  			cout<<"Inside is_initialized"<<endl;
  			x_ << 1, 1, 1, 1, 0.1;
			  	  
  			
  			P_ << 0.15,0,0,0,0,
  				  0,0.15,0,0,0,
  				  0,0,0.3,0,0,
  				  0,0,0,0.03,0,
  				  0,0,0,0,0.3;
  				
  			
  			time_us_ = meas_package.timestamp_;
  			
  			if((meas_package.sensor_type_ == MeasurementPackage::LASER) && (use_laser_))
  			{
  				cout<<"Lidar initialization"<<endl;
  				
  				x_(0) = meas_package.raw_measurements_(0);
  				x_(1) = meas_package.raw_measurements_(1);
  				x_(2) = 4;
  				x_(3) = 0.4;
  				x_(4) = 0.0;
  				P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
            		  0, std_laspy_*std_laspy_, 0, 0, 0,
            		  0, 0, 1, 0, 0,
            		  0, 0, 0, 1, 0,
            		  0, 0, 0, 0, 1;
  				 
  				
  				cout<<"Lidar initialization done"<<endl;
			}
			else if((meas_package.sensor_type_ == MeasurementPackage::RADAR) && (use_radar_))
			{
				cout<<"Radar initialization"<<endl;
				float ro = meas_package.raw_measurements_(0);
				float phi = meas_package.raw_measurements_(1);
				float ro_dot =  meas_package.raw_measurements_(2);
				
				x_(0) = ro * cos(phi);
				x_(1) = ro * sin(phi);
				x_(2) = 4;
  				x_(3) = ro_dot*cos(phi);
  				x_(4) = ro_dot*sin(phi);
  				
  				P_ << std_radr_*std_radr_, 0, 0, 0, 0,
            		   0, std_radr_*std_radr_, 0, 0, 0,
            		   0, 0, 1, 0, 0,
            		   0, 0, 0, std_radphi_*std_radphi_, 0,
            	 	   0, 0, 0, 0, std_radphi_*std_radphi_;
			}
			cout<<"Initialization done"<<endl;
			is_initialized_ = true;
			return;
			
	}
	//Prediction
	double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
	time_us_ = meas_package.timestamp_;
	
	Prediction(dt);
	
	//Update
	
	if(meas_package.sensor_type_ == MeasurementPackage::LASER)
	{
		UpdateLidar(meas_package);
	}
	else if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
	{
		UpdateRadar(meas_package);
	}
  }
}

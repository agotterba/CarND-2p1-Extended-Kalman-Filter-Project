#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  //cout <<"setting previous_timestamp_\n";
  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  F_ = MatrixXd(4, 4);
  Q_ = MatrixXd(4, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225,      0,
                   0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09,      0,    0,
                 0, 0.0009,    0,
                 0,      0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
    x, F, H_laser, H_jacobian, P, etc.  Q?
  */
  //cout <<"setting H_ inits\n";
  H_laser_ << 1,0,0,0, //constant
              0,1,0,0;

  Hj_ <<  1,1,0,0, // depends on position;
         -1,1,0,0, // just filling active fields with 1's for now.
          1,1,1,1;
  //cout <<"setting F_ \n";
  F_ << 1,0,1,0,
        0,1,0,1,
        0,0,1,0,
        0,0,0,1;
  //cout <<"setting Q_ \n";

  Q_ << 1,0,1,0,
        0,1,0,1,
        1,0,1,0,
        0,1,0,1;
  //cout <<"FusionEKF finished initializing\n";

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  cout <<"\nNOTE: ProcessMeasurement: Processing Measurement "<<meas_count<<"\n";
  meas_count += 1;

  //cout <<"setting timestamp\n";
  long long current_timestamp = measurement_pack.timestamp_;
  //cout <<"  set timestamp to "<<current_timestamp<<"\n";
  double ro;
  double theta;
  double ro_dot;
  VectorXd rothetanew = VectorXd(3);

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

    //convert to cartesean coordinates
    ro = measurement_pack.raw_measurements_(0);
    theta = measurement_pack.raw_measurements_(1);
    ro_dot = measurement_pack.raw_measurements_(2);
    //cout <<"  measurement is radar with ro,theta,ro_dot "<<ro<<","<<theta<<","<<ro_dot<<"\n";

    while (theta < -1 * M_PI){
      //cout << "  adding 2pi to theta of "<<theta<<"\n";
      theta += 2* M_PI;
      //cout <<"  theta is now "<<theta<<"\n";
    }
    while (theta > M_PI){
      //cout << "  subtracting 2pi from theta of "<<theta<<"\n";
      theta -= 2 * M_PI;
      //cout <<"  theta is now "<<theta<<"\n";
    }
  }
  rothetanew << ro,theta,ro_dot;
    
 
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    //cout << "EKF: " << endl;
    //ekf_.x_ = VectorXd(4);
    //ekf_.x_ << 1, 1, 1, 1;
    //cout << "ekf_ initialized with ekf_.x_ as: "<< endl;
    //cout << ekf_.x_;
    //cout <<"\n finished writing ekf_.x_\n";
    //cout<<"intializing x\n";
    VectorXd xynew = VectorXd(2);
    xynew << 0,0;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      xynew(0) = ro * cos(theta);
      xynew(1) = ro * sin(theta);
      //cout <<"  first measurement is radar with x,y "<<xynew(0)<<","<<xynew(1)<<"\n";
    } else {
      xynew(0) = measurement_pack.raw_measurements_(0);
      xynew(1) = measurement_pack.raw_measurements_(1);
      //cout <<"  first measurement is laser with x,y "<<xynew(0)<<","<<xynew(1)<<"\n";
    }

    ekf_.x_ = VectorXd(4);
    ekf_.x_ << xynew[0],xynew[1],0,0;
    //cout<<"  intializing P\n";
    ekf_.P_ = MatrixXd(4,4);
    ekf_.P_ << 1,  0,  0,  0,
               0,  1,  0,  0,
               0,  0,1e3,  0,
               0,  0,  0,1e3;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    previous_timestamp_ = current_timestamp;
    //cout <<"  ProcessMeasurement: Finished initialization\n";
    //cout <<"  initial x_ = \n" << ekf_.x_ <<"\n";
    //cout <<"  initial P_ = \n" << ekf_.P_ <<"\n";
    //cout <<"  FusionEKF finished initalization\n";

    cout << "x_ = \n" << ekf_.x_ << endl;
    cout << "P_ = \n" << ekf_.P_ << endl;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  double dt = (current_timestamp - previous_timestamp_) / 1e6;
  //cout <<"  set dt to "<<dt<<", as timestamp changed from "<<previous_timestamp_<<" to "<<current_timestamp<<"\n";
  //cout <<"setting F\n";
  F_ << 1,0,dt,0,
        0,1,0,dt,
        0,0,1,0,
        0,0,0,1;
  //cout <<"setting Q\n";
  double noise_ax = 9.0;
  double noise_ay = 9.0;
  double ax_t44 = noise_ax * pow(dt,4.0) / 4.0;
  double ax_t32 = noise_ax * pow(dt,3.0) / 2.0;
  double ax_t2  = noise_ax * pow(dt,2.0);
  double ay_t44 = noise_ay * pow(dt,4.0) / 4.0;
  double ay_t32 = noise_ay * pow(dt,3.0) / 2.0;
  double ay_t2  = noise_ay * pow(dt,2.0);
  Q_ <<   ax_t44,      0, ax_t32,      0,
               0, ay_t44,      0, ay_t32,
          ax_t32,      0,  ax_t2,      0,
               0, ay_t32,      0,  ay_t2;
  //cout <<"set Q; initalizing ekf\n";
  ekf_.Init(ekf_.x_,ekf_.P_,F_,H_laser_,R_laser_,Q_); //default to laser; it doesn't matter yet
  //ekf_.F_ = F_;
  //ekf_.Q_ = Q_;
  //cout <<"initalized ekf; calling predict\n";
  ekf_.Predict();
  
  //cout <<"FusionEKF finished predict\n";

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

    if (0){
      cout <<"ignoring radar update for laser debug\n";
      previous_timestamp_ = current_timestamp;
      cout << "x_ = \n" << ekf_.x_ << endl;
      cout << "P_ = \n" << ekf_.P_ << endl;
      return; //ignore radar to debug laser
    }

    //cout <<"  calculating jacobian based on x = \n"<<ekf_.x_<<"\n";
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    //cout <<"  calculated jacobian = \n"<<Hj_<<"\n";;

    ekf_.Init(ekf_.x_,ekf_.P_,F_,Hj_,R_radar_,Q_);
    //ekf_.H_ = Hj_;
    //ekf_.R_ = R_radar_;
    //cout << "radar init finished\n";
    //cout << "FusionEKF calling UpdateEKF for radar\n";
    ekf_.UpdateEKF(rothetanew);
 
  } else {
    // Laser updates

    if (0) {
      cout <<"ignoring laser update for radar debug\n";
      previous_timestamp_ = current_timestamp;
      cout << "x_ = \n" << ekf_.x_ << endl;
      cout << "P_ = \n" << ekf_.P_ << endl;
      return; //ignore laser to debug radar
    }

    //cout <<"initializing ekf for laser update\n";
    ekf_.Init(ekf_.x_,ekf_.P_,F_,H_laser_,R_laser_,Q_);
    //ekf_.H_ = H_laser_;
    //ekf_.R_ = R_laser_;
    //cout<<"calling ekf_.Update\n";
    ekf_.Update(measurement_pack.raw_measurements_);
    //cout <<"FusionEKF returned from Update\n";
  }

  previous_timestamp_ = current_timestamp;
  // print the output
  cout << "x_ = \n" << ekf_.x_ << endl;
  cout << "P_ = \n" << ekf_.P_ << endl;

  return;
}

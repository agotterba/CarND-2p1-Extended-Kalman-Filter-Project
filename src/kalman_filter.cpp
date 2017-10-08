#include "kalman_filter.h"
#include <iostream>
#include <math.h>
using namespace std;

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
  /**
  TODO:
    * predict the state
  */

  MatrixXd Ft = F_.transpose();
  x_ = F_ * x_;
  P_ = F_ * P_ * Ft + Q_;

  //cout <<"F_ = \n" << F_ <<"\n";
  //cout <<"P_ = \n" << P_ <<"\n";
  //cout <<"Ft = \n" << Ft <<"\n";
  //cout <<"Q_ = \n" << Q_ <<"\n";
  //cout <<"kalman_filter predicts x = \n"<<x_<<"\n";
  //cout <<"kalman_filter predicts P = \n"<<P_<<"\n";
  
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

  VectorXd z_pred = H_ * x_;

  //cout<<"z = \n"<<z<<"\n";
  //cout <<"H_ = \n" << H_ <<"\n";
  //cout <<"z_pred = \n" << z_pred <<"\n";

  calcUpdate(z,z_pred);

  return;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  //cout<<"z = \n"<<z<<"\n";
  //cout <<"H_ = \n" << H_ <<"\n";
  double px = x_[0];
  double py = x_[1];
  double vx = x_[2];
  double vy = x_[3];
  //double x_ro = pow((x_(0) * x_(0)) + (x_(1) * x_(1)),0.5);
  double x_ro = sqrt(pow(px,2) + pow(py,2));
  double x_theta = atan2(py,px);
  double diff = z(1) - x_theta;
  while (diff > M_PI){
    cout << "NOTE: UpdateEKF: adding 2pi to x_theta of "<<x_theta<<" as diff is "<<diff<<"\n";
    x_theta += 2* M_PI;
    diff = z(1) - x_theta;
    cout << "  now x_theta is "<<x_theta<<" and diff is "<<diff<<"\n";
    }
  while (diff < -1 * M_PI){
    cout << "NOTE: UpdateEKF: subtracting 2pi from x_theta of "<<x_theta<<" as diff is "<<diff<<"\n";
    x_theta -= 2 * M_PI;
    diff = z(1) - x_theta;
    cout << "  now x_theta is "<<x_theta<<" and diff is "<<diff<<"\n";
  }

  double x_ro_dot;
  //  if (x_ro > 1e-6){ //has to be positive, as an rss
  //    x_ro_dot = ( (px*vx)+(py*vy) )/x_ro;
  //  }else{
    //oct7 x_ro_dot = ( (x_(0)*x_(2))+(x_(1)*x_(3)) )/1e-4;
  //    x_ro_dot = 0;
  //  }
  if (x_ro < 1e-6){
    cout << "WARNING: UpdateEKF: x_ro less than 1e-6.  Setting to 1e-6\n";
    x_ro = 1e-6;
  }
  x_ro_dot = ( (px*vx)+(py*vy) )/x_ro;
  
  VectorXd z_pred = VectorXd(3);
  z_pred << x_ro,x_theta,x_ro_dot;
  //cout<<"z_pred = \n"<<z_pred<<"\n";
  calcUpdate(z,z_pred);
  return;

  /* Original, before refactoring
  VectorXd y = z - z_pred;
  cout<<"y = \n"<<y<<"\n";
  MatrixXd Ht = H_.transpose();
  cout <<"Ht = \n"<<Ht <<"\n";
  cout <<"P_ = \n" << P_ <<"\n";
  cout <<"R_ = \n" << R_ <<"\n";
  MatrixXd S = H_ * P_ * Ht + R_;
  cout <<"S = \n" << S <<"\n";
  MatrixXd Si = S.inverse();
  cout <<"Si = \n" << Si <<"\n";
  MatrixXd PHt = P_ * Ht;
  cout <<"PHt = \n" << PHt <<"\n";
  MatrixXd K = PHt * Si;
  cout <<"K = \n" << K <<"\n";
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size,x_size);
  cout <<"I = \n" << I <<"\n";
  P_ = (I - K * H_) * P_;
  */
}

void KalmanFilter::calcUpdate(const VectorXd &z,const VectorXd &z_pred) {
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

  //cout <<"y = \n"<<y<<"\n";
  //cout <<"Ht = \n"<<Ht <<"\n";
  //cout <<"R_ = \n" << R_ <<"\n";
  //cout <<"S = \n" << S <<"\n";
  //cout <<"Si = \n" << Si <<"\n";
  //cout <<"PHt = \n" << PHt <<"\n";
  //cout <<"K = \n" << K <<"\n";
  //cout <<"I = \n" << I <<"\n";
  //cout <<"x = \n" << x_ <<"\n";
  //cout <<"P = \n" << P_ <<"\n";
  
  return;
}

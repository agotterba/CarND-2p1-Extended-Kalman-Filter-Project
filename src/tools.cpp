#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  Eigen::MatrixXd rmse = VectorXd(4);
  rmse << 0,0,0,0;
  if (estimations.size() == 0){
    cout<<"NOTE: CalculateRMSE: called with estimations list of size 0\n";
    return rmse;
  }
  if (ground_truth.size() != estimations.size()){
    cout<<"ERROR: CalculateRMSE: estimations vector has size "<<estimations.size()
        <<" but ground_truth vector has size " <<ground_truth.size()<<" \n";
    return rmse;
  }

  //VectorXd diff = estimations.array() - ground_truth.array();
  //VectorXd diff2 = diff.array() * diff.array();
  //accumulate squared residuals

  for(int i=0; i < estimations.size(); ++i){
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }
  
  //calculate the mean
  // ... your code here
  VectorXd mean;
  rmse = rmse / estimations.size();
  
  //calculate the squared root
  // ... your code here
  rmse = rmse.array().sqrt();
  
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
  Hj << 0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0;
  //recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  double x2y2 = (px * px) + (py * py);
  if (x2y2 < 1e-6){
    cout<<"WARNING: CalculateJacobian: Division by zero as px is 0 and py is 0\n  setting x2y2 to 1e-6\n";
    x2y2 = 1e-6;
  }
  double sx2y2 = pow(x2y2,0.5);
  double x2y21p5  = pow(x2y2,1.5);
  double vx1  = vx*py - vy*px;
  double vy1  = vy*px - vx*py;

  //check division by zero
  //  if (x2y2 < 1e-6){
    
    //Hj <<        px*1e6,       py*1e6,        0,        0,
    //          -1*py*1e6,       px*1e6,        0,        0,
    //         py*vx1*1e6,   px*vy1*1e6,   px*1e6,   py*1e6;
  //return Hj;
  //  }
	
  //compute the Jacobian matrix
    
  Hj <<       px/sx2y2,       py/sx2y2,        0,        0,
            -1*py/x2y2,        px/x2y2,        0,        0,
        py*vx1/x2y21p5, px*vy1/x2y21p5, px/sx2y2, py/sx2y2;
          
  return Hj;
}

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

  VectorXd rmse = VectorXd(4);
  rmse << 0,0,0,0;

  int num_samples = estimations.size();
  VectorXd diff = VectorXd(4);

  for (int i = 0; i < num_samples; ++i) {
    diff = estimations[i] - ground_truth[i];
    diff = diff.array() * diff.array();
    rmse = rmse + diff;
  }

  rmse = rmse / num_samples;
  rmse = rmse.array().sqrt();

  /* Automatic file saving doesn't seem to work for some reason so I'll just
  collect the outputs from the terminal */
  VectorXd gt = ground_truth[num_samples-1];
  VectorXd es = estimations[num_samples-1];
  cout << "GT  : " << gt[0] << " " << gt[1] << " " << gt[2] << " " << gt[3] << endl;
  cout << "ES  : " << es[0] << " " << es[1] << " " << es[2] << " " << es[3] << endl;
  cout << "RMSE: " << rmse[0] << " " << rmse[1] << " " << rmse[2] << " " << rmse[3] << endl;

  return rmse;

}

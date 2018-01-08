#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {
}

Tools::~Tools() {
}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
		const vector<VectorXd> &ground_truth) {
	/**
	 TODO:
	 * Calculate the RMSE here.
	 */

	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for (unsigned int i = 0; i < estimations.size(); ++i) {

		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array() * residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse / estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	cout << "RMSE = " << rmse.transpose() << endl;
	return rmse;
}

void Tools::ConvertPolar2Cartesian(double rho, double phi, double& px,
		double& py) {
	double rho2 = rho * rho;
	double tan_phi = tan(phi);
	double px2 = rho2 / (1.0 + tan_phi * tan_phi);
	px = sqrt(px2);
	py = sqrt(rho2 - px2);
}

double Tools::NormalizeAngle(double rad) {
	const double TWO_PI = 2.0 * M_PI;

	if (rad < -M_PI) {
		int add_2pi = (int) ((-rad + M_PI) / TWO_PI);
		rad += add_2pi * TWO_PI;
	} else if (rad > M_PI) {
		int sub_2pi = (int) ((rad + M_PI) / TWO_PI);
		rad -= sub_2pi * TWO_PI;
	}

	return rad;
}

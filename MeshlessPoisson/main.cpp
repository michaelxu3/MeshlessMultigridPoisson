#include "grid.h"
#include <iostream>
#include <ctime>
#define pi 3.141592653589793238462643383279
using std::cout;
using std::endl;
int main() {
	std::vector<std::tuple<double, double, double>> points;
	vector<std::tuple<vector<int>, int, vector<double>>> boundaries;
	boundaries.resize(1);
	std::tuple<vector<int>, int, vector<double>> & bound = boundaries[0];
	vector<double> actualSoln_homog;
	vector<double> actualSoln_inhomog;

	int nx = 41;
	int ny = 41;
	Eigen::VectorXd source = Eigen::VectorXd(nx*ny);
	int pointIndex = 0;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			points.push_back(std::tuple<double, double, double>(i*1.0/(nx-1), j*1.0/(ny-1), 0));
			if (i == nx-1 || i == 0 || j == 0 || j == ny-1) {
				std::get<0>(bound).push_back(pointIndex);
			}
			source(pointIndex) = -2 * pi*pi*std::sin(pi*i / (double)(nx-1))*std::sin(pi*j / (double)(ny-1));
			actualSoln_homog.push_back(std::sin(pi*i / (double)(nx-1))*std::sin(pi*j / (double)(ny-1)));
			pointIndex++;
		}
	}
	Grid test = Grid(points, boundaries, 3, 7, 40, 1.5, source, 1000);
	vector<double> bcValues = vector<double>(4*(nx-1), 0.0);
	test.setBCFlag(0, std::string("dirichlet"), bcValues);
	test.build_laplacian();

	test.sor();
	Eigen::VectorXd* soln = test.values_;
	double l1_error = 0;
	for (int i = 0; i < nx*nx; i++) {
		//cout << "soln: " << soln->coeff(i) << " actual soln:" << actualSoln_homog.at(i) << endl;
		l1_error += std::abs(soln->coeff(i) - actualSoln_homog.at(i))/nx/nx;
	}
	cout << " L1 error: " << l1_error << endl;
}
#include "grid.h"
#include <iostream>
#define pi 3.141592653589793238462643383279
using std::cout;
using std::endl;
int main() {
	std::vector<std::tuple<double, double, double>> points;
	vector<std::tuple<vector<int>, int, vector<double>>> boundaries;
	boundaries.resize(1);
	std::tuple<vector<int>, int, vector<double>> & bound = boundaries[0];
	Eigen::VectorXd source = Eigen::VectorXd(441);
	int pointIndex = 0;
	for (int i = 0; i <= 20; i++) {
		for (int j = 0; j <= 20; j++) {
			points.push_back(std::tuple<double, double, double>(i/20.0, j/20.0, 0));
			if (i == 20 || i == 0 || j == 0 || j ==20) {
				std::get<0>(bound).push_back(pointIndex);
			}
			source(pointIndex) = -2 * pi*pi*std::sin(pi*i / 20.0)*std::sin(pi*j / 20.0);
			pointIndex++;
		}
	}
	Grid test = Grid(points, boundaries, 3, 3, 40, 1.0);
	vector<double> bcValues = vector<double>(80, 0.0);
	test.setBCFlag(0, std::string("dirichlet"), bcValues);
	test.build_laplacian();
	test.sor(1, source);

}
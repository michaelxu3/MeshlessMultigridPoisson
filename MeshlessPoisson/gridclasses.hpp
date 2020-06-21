#pragma once
#include <vector>
class GridProperties {
public:
	int rbfExp;
	int polyDeg;
	int laplaceMatSize;
	int stencilSize;
	double omega;
	int iters;
};
class Boundary {
public:
	int type;
	std::vector<int> bcPoints;
	std::vector<double> values;
};
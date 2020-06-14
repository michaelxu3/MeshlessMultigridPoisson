#ifndef GRID_H
#define GRID_H

#include <vector>
#include <tuple>
#include <string>
#include <Eigen/Dense>
using std::vector;
class Grid {
public:
	Grid(vector<std::tuple<double, double, double>> points, vector<std::tuple<vector<int>, int, double>> boundaries, int rbfExp, int polyDeg);
	~Grid();
	void setBC(int boundary, std::string type, double boundValue);
	void build_laplacian();
	vector <int> kNearestNeighbors(int pointNumber, int k);
private:
	vector<std::tuple<double, double, double>> points_;
	vector<std::tuple<vector<int>, int, double>> boundaries_;
	int rbfExp_;
	int polyDeg_;
	int laplaceMatSize_;
	Eigen::VectorXd* values_;
	Eigen::MatrixXd* laplaceMat_;

	std::pair<Eigen::MatrixXd, vector<int>> buildCoeffMatrix(int pointNum, int k);
	std::pair<Eigen::VectorXd, vector<int>> laplaceWeights(int pointID, int k);
};
#endif

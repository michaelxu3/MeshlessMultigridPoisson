#ifndef GRID_H
#define GRID_H

#include <vector>
#include <tuple>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
using std::vector;
class Grid {
public:
	Grid(vector<std::tuple<double, double, double>> points, vector<std::tuple<vector<int>, int, vector<double>>> boundaries, int rbfExp, int polyDeg, int stencilSize, double omega);
	~Grid();
	void setBCFlag(int boundary, std::string type, vector<double> boundValue);
	void boundaryOp();
	void build_laplacian(int k);
	void sor(int numIts, Eigen::VectorXd source);
	vector <int> kNearestNeighbors(int pointNumber);
private:
	vector<std::tuple<double, double, double>> points_;
	vector<std::tuple<vector<int>, int, vector<double>>> boundaries_;
	int rbfExp_;
	int polyDeg_;
	int laplaceMatSize_;
	int stencilSize_;
	vector<int> bcFlags_;
	double omega_;
	Eigen::VectorXd* values_;
	Eigen::SparseMatrix<double>* laplaceMat_;

	std::pair<Eigen::MatrixXd, vector<int>> buildCoeffMatrix(int pointNum);
	std::pair<Eigen::VectorXd, vector<int>> laplaceWeights(int pointID);
};
#endif

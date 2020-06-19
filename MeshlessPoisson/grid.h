#ifndef GRID_H
#define GRID_H

#include <vector>
#include <tuple>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "gridclasses.hpp"
using std::vector;
class Grid {
public:
	Eigen::VectorXd* values_;

	Grid(vector<std::tuple<double, double, double>> points, vector<std::tuple<vector<int>, int, vector<double>>> boundaries, 
										int rbfExp, int polyDeg, int stencilSize, double omega, Eigen::VectorXd source, int sorIters);
	~Grid();
	void setBCFlag(int boundary, std::string type, vector<double> boundValue);
	void build_laplacian();
	void sor();
	vector<int> kNearestNeighbors(std::tuple<double, double, double> point);
	std::pair<Eigen::MatrixXd, vector<int>> buildCoeffMatrix(std::tuple<double, double, double> point);
	std::pair<Eigen::VectorXd, vector<int>> pointInterpWeights(std::tuple<double, double, double> point);
	
	int getSize();
	int getStencilSize();
	int getPolyDeg();
	std::tuple<double, double, double> getPoint(int index);
private:
	vector<std::tuple<double, double, double>> points_;
	vector<std::tuple<vector<int>, int, vector<double>>> boundaries_;
	GridProperties properties_;
	int laplaceMatSize_;
	Eigen::SparseMatrix<double>* laplaceMat_;
	vector<int> bcFlags_;
	Eigen::VectorXd source_;
	int sorIters_;

	void boundaryOp();
	vector <int> kNearestNeighbors(int pointNumber);
	std::pair<Eigen::MatrixXd, vector<int>> buildCoeffMatrix(int pointNum);
	std::pair<Eigen::VectorXd, vector<int>> laplaceWeights(int pointID);
};
#endif

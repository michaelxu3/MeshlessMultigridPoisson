#ifndef GRID_H
#define GRID_H

#include <vector>
#include <tuple>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "gridclasses.hpp"
typedef std::tuple<double, double, double> Point;
using std::vector;
class Grid {
public:
	Eigen::VectorXd* values_;
	Eigen::VectorXd* residuals_;
	Eigen::VectorXd source_;

	Grid(vector<std::tuple<double, double, double>> points, vector<Boundary> boundaries, 
				GridProperties properties, Eigen::VectorXd source);
	~Grid();
	void setBCFlag(int boundary, std::string type, vector<double> boundValue);
	void build_laplacian();
	void sor(Eigen::SparseMatrix<double, 1>* matrix, Eigen::VectorXd* values, Eigen::VectorXd* rhs);
	void modifyCoeffDirichlet();
	void directSolve();
	std::pair<double, double> minMaxCoord(vector<int> pointIDs, char coord);
	std::vector<std::tuple<double, double, double>> shifting_scaling(vector<int> pointIDs, Point evalPoint);
	Eigen::VectorXd residual();

	 
	vector<int> kNearestNeighbors(Point point);
	std::tuple<Eigen::MatrixXd, vector<int>, vector<Point>> buildCoeffMatrix(Point point);
	std::tuple<Eigen::MatrixXd, vector<int>, vector<Point>> buildCoeffMatrix(int pointNum);
	std::pair<Eigen::VectorXd, vector<int>> pointInterpWeights(std::tuple<double, double, double> point);
	
	int getSize();
	int getStencilSize();
	int getPolyDeg();
	std::tuple<double, double, double> getPoint(int index);
	vector<std::tuple<double, double, double>> points_;
	vector<Boundary> boundaries_;
	GridProperties properties_;
	int laplaceMatSize_;
	Eigen::SparseMatrix<double, Eigen::RowMajor>* laplaceMat_;
	vector<int> bcFlags_;
	void boundaryOp(std::string coarse);
	vector <int> kNearestNeighbors(int pointNumber);
	std::pair<Eigen::VectorXd, vector<int>> laplaceWeights(int pointID);
};
#endif

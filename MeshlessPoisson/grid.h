#ifndef GRID_H
#define GRID_H

#include <vector>
#include <tuple>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "gridclasses.hpp"
#include "fileReadingFunctions.h"
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
	void fix_vector_bound_coarse (Eigen::VectorXd* vec);
	void print_bc_values(Eigen::VectorXd vec);
	std::pair<double, double> minMaxCoord(vector<int> pointIDs, char coord);
	std::vector<std::tuple<double, double, double>> shifting_scaling(vector<int> pointIDs, Point evalPoint);
	Eigen::VectorXd residual();
	void diagonal_scaling(Eigen::SparseMatrix<double, 1>* matrix, Eigen::VectorXd* rhs);
	void cuthill_mckee_ordering(vector<vector<int>> &adjacency, vector<int> &order);
	void reverse_cuthill_mckee_ordering(vector<vector<int>> &adjacency, vector<int> &order);
	void rcm_order_points();
	vector<int> kNearestNeighbors(Point point);
	vector <int> kNearestNeighbors(int pointNumber);
	std::tuple<Eigen::MatrixXd, vector<int>, vector<Point>> buildCoeffMatrix(Point point);
	std::tuple<Eigen::MatrixXd, vector<int>, vector<Point>> buildCoeffMatrix(int pointNum);
	std::pair<Eigen::VectorXd, vector<int>> laplaceWeights(int pointID);
	std::pair<Eigen::VectorXd, vector<int>> pointInterpWeights(std::tuple<double, double, double> point);
	
	int getSize();
	int getStencilSize();
	int getPolyDeg();



	vector<std::tuple<double, double, double>> points_;
	vector<Boundary> boundaries_;
	GridProperties properties_;
	int laplaceMatSize_;
	Eigen::SparseMatrix<double, Eigen::RowMajor>* laplaceMat_;
	vector<int> bcFlags_;
	void boundaryOp(std::string coarse);
	vector<double> cond_rbf;
	//Eigen::SparseMatrix<double, Eigen::RowMajor>* adjacencyMat_;
};
#endif

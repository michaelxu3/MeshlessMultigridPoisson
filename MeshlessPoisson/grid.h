#ifndef GRID_H
#define GRID_H
#include "gridclasses.hpp"
#include <tuple>
#include <string>
#include "fileReadingFunctions.h"
#include "general_computation_functions.h"
#include <stdexcept>
#include "math.h"
#include <iostream>
#include <algorithm>
#include <queue>
#include <time.h>
#include <unordered_map>  

typedef std::tuple<double, double, double> Point;
using std::vector;
class Grid {
public:
	//Instance Variables
	Eigen::VectorXd* values_;
	Eigen::VectorXd* residuals_;
	Eigen::VectorXd source_;
	vector<Point> points_;
	vector<Boundary> boundaries_;
	//first = bc point id, second = weights, third = corresponding neighbors.
	vector<deriv_normal_bc> deriv_normal_coeffs_;
	GridProperties properties_;
	int laplaceMatSize_;
	Eigen::SparseMatrix<double, Eigen::RowMajor>* laplaceMat_;
	vector<int> bcFlags_;
	bool neumannFlag_;
	bool implicitFlag_;

	Grid(vector<std::tuple<double, double, double>> points, vector<Boundary> boundaries, 
				GridProperties properties, Eigen::VectorXd source);
	~Grid();

	void boundaryOp(std::string coarse);
	void setBCFlag(int boundary, std::string type, vector<double> boundValue);
	void build_laplacian();
	void build_deriv_normal_bound();
	void sor(Eigen::SparseMatrix<double, 1>* matrix, Eigen::VectorXd* values, Eigen::VectorXd* rhs);

	void setNeumannFlag();
	void modify_coeff_neumann();
	void bound_eval_neumann();
	void fix_vector_bound_coarse (Eigen::VectorXd* vec);
	void print_bc_values(Eigen::VectorXd vec);
	void print_check_bc_normal_derivs();

	double cond_L();
	Eigen::VectorXd residual();

	void rcm_order_points();
	vector<Point> pointIDs_to_vector(const vector<int> &pointIDs);
	vector<int> kNearestNeighbors(Point point, bool neumannFlag, bool pointBCFlag);
	vector <int> kNearestNeighbors(int pointNumber, bool neumannFlag);

	std::tuple<Eigen::MatrixXd, vector<int>, vector<Point>> buildCoeffMatrix(Point point, bool neumann, bool pointBCFlag);
	std::tuple<Eigen::MatrixXd, vector<int>, vector<Point>> buildCoeffMatrix(int pointNum, bool neumann);
	std::pair<Eigen::VectorXd, vector<int>> laplaceWeights(int pointID);
	std::pair<Eigen::VectorXd, vector<int>> derivx_weights(int pointID);
	std::pair<Eigen::VectorXd, vector<int>> derivy_weights(int pointID);
	std::pair<Eigen::VectorXd, vector<int>> pointInterpWeights(Point point);
	
	int getSize();
	int getStencilSize();
	int getPolyDeg();
	vector<double> cond_rbf;
};
#endif

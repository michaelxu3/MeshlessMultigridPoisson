#pragma once
#include "multigrid.h"
#include "fileReadingFunctions.h"
#include <iostream>
#include <ctime>
#include <fstream>
#include <iomanip>
#define pi 3.141592653589793238462643383279
using std::cout;
using std::endl;

Grid* generateHomogDirichletGrid(int nx, int ny, int iters) {
	std::vector<std::tuple<double, double, double>> points;
	vector<int> bPts;
	vector<double> bValues;
	Eigen::VectorXd source = Eigen::VectorXd(nx*ny);
	double x, y;
	int pointIndex = 0;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			x = i * 1.0 / (nx - 1);
			y = j * 1.0 / (ny - 1);
			points.push_back(std::tuple<double, double, double>(x, y, 0));
			if (x == 0 || x == 1 || y == 0 || y == 1) {
				bPts.push_back(pointIndex);
				bValues.push_back(0.0);
			}
			source(pointIndex) = -2 * pi*pi*std::sin(pi*x)*std::sin(pi*y);
			pointIndex++;
		}
	}
	Boundary boundary;
	boundary.bcPoints = bPts;
	boundary.type = 1;
	boundary.values = bValues;
	vector<Boundary> bcs;
	bcs.push_back(boundary);
	GridProperties props;
	props.rbfExp = 3;
	props.iters = iters;
	props.polyDeg = 5;
	// cloud size
	props.stencilSize = (int)(1*(props.polyDeg + 1) * (props.polyDeg + 2));
	props.omega = 1.4;
	Grid* grid = new Grid(points, bcs, props, source);;
	grid->setBCFlag(0, std::string("dirichlet"), bValues);
	grid->build_laplacian();
	return grid;
}

void writeVectorToTxt(vector<double> vec, const char* filename)
{
	std::ofstream file;
	file.open(filename);
	//f << fixed << setprecision(2) << endl;
	for (size_t i = 0; i < vec.size(); i++) {
		file << vec[i] << "\n";
	}
	file.close();
}

double calc_l1_error(Grid* grid) {
	vector<Point> points = grid->points_;
	Eigen::VectorXd actual(points.size());
	double x, y;
	for (size_t i = 0; i < points.size(); i++) {
		x = std::get<0>(points[i]);
		y = std::get<1>(points[i]);
		actual(i) = std::sin(pi*x)*std::sin(pi*y);
	}
	return (*grid->values_ - actual).lpNorm<1>() / grid->source_.lpNorm<1>();
}
void testCartesianSingleGrid() {

	//10
	Grid* testGrid = generateHomogDirichletGrid(100, 100, 5);
	vector<double> res;
	/*
	vector<Point> points = testGrid->points_;
	Eigen::VectorXd actual(testGrid->laplaceMatSize_);
	testGrid->boundaryOp("fine");
	double x, y;
	for (size_t i = 0; i < points.size(); i++) {
		x = std::get<0>(points[i]);
		y = std::get<1>(points[i]);
		actual(i) = std::sin(pi*x)*std::sin(pi*y);
	}
	*testGrid->values_ = actual;
	//testGrid->modifyCoeffDirichlet();
	cout << "l1 error: " << calc_l1_error(testGrid) << endl;
	cout << "residual: " << testGrid->residual().norm() / testGrid->source_.norm() << endl;
	*/
	for (int i = 0; i < 50; i++) {
		res.push_back(testGrid->residual().norm() / testGrid->source_.norm());
		cout << "residual: " << testGrid->residual().norm() / testGrid->source_.norm() << endl;
		testGrid->sor(testGrid->laplaceMat_, testGrid->values_, &testGrid->source_);
	}
}
void testCartesianMultigrid(){
	
	Multigrid mg = Multigrid();
	clock_t start = std::clock();
	//mg.addGrid(generateHomogDirichletGrid(1000, 1000));
	//mg.addGrid(generateHomogDirichletGrid(500, 500));
	//mg.addGrid(generateHomogDirichletGrid(250, 250));
	//mg.addGrid(generateHomogDirichletGrid(125, 125));
	mg.addGrid(generateHomogDirichletGrid(100,100, 5));
	mg.addGrid(generateHomogDirichletGrid(50, 50, 5));
	//mg.addGrid(generateHomogDirichletGrid(25, 25));
	//mg.addGrid(generateHomogDirichletGrid(13, 13));
	mg.buildMatrices();
	clock_t build = std::clock();
	cout << (build - start) / ((double)CLOCKS_PER_SEC) << endl;
	vector<double> res;
	for (int i = 0; i < 30; i++) {
		res.push_back(mg.residual());
		mg.vCycle();
	}
	cout << (std::clock() - build) / ((double)CLOCKS_PER_SEC) << endl;
	
}
Grid* genGmshGrid(const char* filename, int polydeg, int iters, std::string filetype) {
	vector<std::tuple<double, double, double>> points;
	if (filetype.compare("txt") == 0) {
		points = pointsFromTxts(filename);
	}
	else {
		points = pointsFromMshFile(filename);
	}
	double x, y;
	vector<int> bPts;
	vector<double> bValues;

	Eigen::VectorXd source(points.size());
	for (size_t i = 0; i < points.size(); i++) {
		x = std::get<0>(points[i]);
		y = std::get<1>(points[i]);
		source(i) = -2 * pi*pi*std::sin(pi*x)*std::sin(pi*y);
		if (x == 0 || x == 1 || y == 0 || y == 1) {
			bPts.push_back(i);
			bValues.push_back(0.0);
		}
	}
	cout << points.size() << endl;
	Boundary boundary;
	boundary.bcPoints = bPts;
	boundary.type = 1;
	boundary.values = bValues;
	vector<Boundary> bcs;
	bcs.push_back(boundary);
	GridProperties props;
	props.rbfExp = 3;
	props.iters = iters;
	props.polyDeg = polydeg;
	// cloud size
	props.stencilSize = (int)(0.75*(props.polyDeg + 1) * (props.polyDeg + 2));
	props.omega = 1.4;
	Grid* grid = new Grid(points, bcs, props, source);
	grid->setBCFlag(0, std::string("dirichlet"), bValues);
	grid->rcm_order_points();
	grid->build_laplacian();
	// cout << bPts.size() << endl;
	return grid;
}
void testGmshSingleGrid() {
	Grid* testGrid = genGmshGrid("gmshtest0/square_11825.msh", 5, 5, "msh");

	vector<double> res;
	testGrid->boundaryOp("fine");

	
	vector<Point> points = testGrid->points_;
	Eigen::VectorXd actual(points.size());

	double x, y;
	vector<double> xv, yv;
	for (size_t i = 0; i < points.size(); i++) {
		x = std::get<0>(points[i]);
		y = std::get<1>(points[i]);
		xv.push_back(x);
		yv.push_back(y);
	}
	writeVectorToTxt(xv, "x.txt");
	writeVectorToTxt(yv, "y.txt");

	//testGrid->diagonal_scaling(testGrid->laplaceMat_, &testGrid->source_);
	//cout << "num bc points: " << testGrid->boundaries_[0].bcPoints.size() << endl;
	for (int i = 0; i < 8000; i++) {
		cout << "residual: " << testGrid->residual().lpNorm<1>() / testGrid->source_.lpNorm<1>()<< endl;
		res.push_back(testGrid->residual().lpNorm<1>() / testGrid->source_.lpNorm<1>());

		//testGrid->directSolve();
		testGrid->sor(testGrid->laplaceMat_, testGrid->values_, &testGrid->source_);
	}
	vector<double> rawRes;
	Eigen::VectorXd resid = testGrid->residual();
	for (int i = 0; i < resid.rows(); i++) {
		rawRes.push_back(resid.coeff(i));
	}
	writeVectorToTxt(rawRes, "raw_res.txt");
	cout << " values norm: " << testGrid->values_->lpNorm<1>() << endl;
	cout << testGrid->values_->maxCoeff() << endl;
	cout << testGrid->values_->minCoeff() << endl;
	cout << calc_l1_error(testGrid) << endl;
	
}
void testGmshMultigrid() {
	Multigrid mg;

	mg.addGrid(genGmshGrid("finer_mesh.txt", 5, 7, "txt"));
	mg.addGrid(genGmshGrid("fine_mesh.txt", 3, 7, "txt"));
	mg.addGrid(genGmshGrid("inter_mesh.txt", 3, 7, "txt"));
	mg.addGrid(genGmshGrid("coar_mesh.txt", 3, 7, "txt"));

	//mg.addGrid(genGmshGrid("gmshtest0/square_11825.msh", 5, 5, "msh"));
	//mg.addGrid(genGmshGrid("gmshtest0/square_2555.msh", 5, 5, "msh"));
	//mg.addGrid(genGmshGrid("gmshtest0/square_675.msh", 5, 5, "msh"));
	//mg.addGrid(genGmshGrid("gmshtest0/square_171.msh", 5, 5, "msh"));

	//mg.addGrid(genGmshGrid("square_test1_9482.msh", 5, 7, "msh"));
	//mg.addGrid(genGmshGrid("square_test1_2221.msh", 3, 7, "msh"));
	//mg.addGrid(genGmshGrid("square_test1_552.msh", 3, 7, "msh"));
	//mg.addGrid(genGmshGrid("square_test1_149.msh", 3, 7, "msh"));
	mg.buildMatrices();

	vector<double> res;
	for (int i = 0; i < 150; i++) {
		res.push_back(mg.residual());
		mg.vCycle();
	}
	cout << calc_l1_error(mg.grids_[mg.grids_.size()-1].second) << endl;
	writeVectorToTxt(res, "residual.txt");


}
int main() {
	//testGmshSingleGrid();
	testGmshMultigrid();
	//testCartesianMultigrid();
	return 0;
}
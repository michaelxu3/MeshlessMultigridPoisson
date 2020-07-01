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

Grid* generateHomogDirichletGrid(int nx, int ny) {
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
	props.iters = 7;
	props.polyDeg = 5;
	// cloud size
	props.stencilSize = (int)(1.25*(props.polyDeg + 1) * (props.polyDeg + 2));
	props.omega = 1.5;
	Grid* grid = new Grid(points, bcs, props, source);;
	grid->setBCFlag(0, std::string("dirichlet"), bValues);
	grid->build_laplacian();
	return grid;
}

int writeRes(vector<double> res)
{
	std::ofstream file;
	file.open("residual.txt");
	//f << fixed << setprecision(2) << endl;
	for (size_t i = 0; i < res.size(); i++) {
		file << res[i] << "\n";
	}
	file.close();
	return 0;
}
 void testCartesianGrid() {
	 
	//10
	Grid* testGrid = generateHomogDirichletGrid(41, 41);
	vector<double> res;
	testGrid->boundaryOp("fine");
	for (int i = 0; i < 1000; i++) {
		testGrid->sor(testGrid->laplaceMat_, testGrid->values_, &testGrid->source_);
		res.push_back(testGrid->residual().norm() / testGrid->source_.norm());
		cout << "residual: " << testGrid->residual().norm() / testGrid->source_.norm() << endl;
	}
	//writeRes(res);
	
	/*
	Multigrid mg = Multigrid();
	clock_t start = std::clock();
	//mg.addGrid(generateHomogDirichletGrid(1000, 1000));
	//mg.addGrid(generateHomogDirichletGrid(500, 500));
	//mg.addGrid(generateHomogDirichletGrid(250, 250));
	//mg.addGrid(generateHomogDirichletGrid(125, 125));
	//mg.addGrid(generateHomogDirichletGrid(100,100));
	mg.addGrid(generateHomogDirichletGrid(50, 50));
	mg.addGrid(generateHomogDirichletGrid(25, 25));
	mg.addGrid(generateHomogDirichletGrid(13, 13));
	mg.buildMatrices();
	clock_t build = std::clock();
	cout << (build - start) / ((double)CLOCKS_PER_SEC) << endl;
	vector<double> res;
	for (int i = 0; i < 30; i++) {
		res.push_back(mg.residual());
		mg.vCycle();
	}
	writeRes(res);
	cout << (std::clock() - build) / ((double)CLOCKS_PER_SEC) << endl;
	*/
}
 Grid* genGmshGrid(const char* filename) {
	 vector<std::tuple<double, double, double>> points = pointsFromMshFile(filename);
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
	 Boundary boundary;
	 boundary.bcPoints = bPts;
	 boundary.type = 1;
	 boundary.values = bValues;
	 vector<Boundary> bcs;
	 bcs.push_back(boundary);
	 GridProperties props;
	 props.rbfExp = 3;
	 props.iters = 7;
	 props.polyDeg = 5;
	 // cloud size
	 props.stencilSize = (int)(1.25*(props.polyDeg + 1) * (props.polyDeg + 2));
	 props.omega = 1;
	 Grid* grid = new Grid(points, bcs, props, source);;
	 grid->setBCFlag(0, std::string("dirichlet"), bValues);
	 grid->build_laplacian();
	// cout << bPts.size() << endl;
	 return grid;
 }
 void testGmshGrid() {
	 Grid* testGrid = genGmshGrid("square.msh");
	 vector<double> res;
	 testGrid->boundaryOp("fine");
	 //testGrid->modifyCoeffDirichlet();
	 for (int i = 0; i < 1000; i++) {
		 cout << "residual: " << testGrid->residual().norm() / testGrid->source_.norm() << endl;
		 //testGrid->directSolve();
		 testGrid->sor(testGrid->laplaceMat_, testGrid->values_, &testGrid->source_);

		 res.push_back(testGrid->residual().norm() / testGrid->source_.norm());
	 }
	 cout << testGrid->values_->maxCoeff() << endl;
	 cout << testGrid->values_->minCoeff() << endl;

 }
 int main() {
	 testGmshGrid();
	 return 0;
 }

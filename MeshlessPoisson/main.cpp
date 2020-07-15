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
	props.omega = 1.3;
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

double calc_l1_error(Grid* grid, bool neumannFlag) {
	vector<Point> points = grid->points_;
	Eigen::VectorXd actual(points.size());
	double x, y;
	for (size_t i = 0; i < points.size(); i++) {
		x = std::get<0>(points[i]);
		y = std::get<1>(points[i]);
		actual(i) = neumannFlag? std::cos(pi*x)*std::cos(pi*y) : std::sin(pi*x)*std::sin(pi*y);
	}
	double solMean = 0;
	double manufacturedMean = 0;
	if (!neumannFlag) {
		return (*grid->values_ - actual).lpNorm<1>() / grid->laplaceMatSize_;
	}

	for (int i = 0; i < grid->laplaceMatSize_; i++) {
		solMean += grid->values_->coeff(i) / grid->laplaceMatSize_;
		manufacturedMean += actual.coeff(i) / grid->laplaceMatSize_;
	}

	for (int i = 0; i < grid->laplaceMatSize_; i++) {
		grid->values_->coeffRef(i) += (manufacturedMean - solMean);
	}

	double error = 0;
	for (int i = 0; i < points.size(); i++) {
		error += std::abs(grid->values_->coeff(i) - actual.coeff(i));
	}
	error /= points.size();
	return error;
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
	testGrid->boundaryOp("fine");
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
Grid* genGmshGridNeumann(const char* filename, int polydeg, int iters, std::string filetype) {
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
	Eigen::VectorXd actual(points.size() + 1);
	Eigen::VectorXd source(points.size() + 1);
	for (size_t i = 0; i < points.size(); i++) {
		x = std::get<0>(points[i]);
		y = std::get<1>(points[i]);
		actual(i) = std::cos(pi*x)*std::cos(pi*y);
		source(i) = -2*pi*pi*std::cos(pi*x)*std::cos(pi*y);
		if (x == 0 || x == 1 || y == 0 || y == 1) {
			//cout << "x: " << x << " y: " << y << endl;
			bPts.push_back(i);
			bValues.push_back(0.0);
		}
	}
	cout << bPts.size() << endl;
	actual(actual.rows() - 1) = 0;

	source(source.rows() - 1) = 0;
	cout << "Grid Size: " << points.size() << endl;
	Boundary boundary;
	boundary.bcPoints = bPts;
	boundary.type = 2;
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
	grid->setBCFlag(0, std::string("neumann"), bValues);
	//cout << grid->neumannFlag_ << endl;
	//grid->rcm_order_points();
	grid->build_deriv_normal_bound();
	grid->build_laplacian();
	grid->modify_coeff_neumann();
	//*grid->values_ = actual;
	// cout << bPts.size() << endl;
	return grid;
}
void testGmshSingleGrid() {
	Grid* testGrid = genGmshGridNeumann("inter_mesh.txt", 5, 5, "txt");
	//cout << testGrid->laplaceMat_->toDense() << endl;
	vector<double> res;
	testGrid->boundaryOp("fine");
	for (int i = 0; i < 10000; i++) {
		cout << "residual: " << testGrid->residual().lpNorm<1>() / testGrid->source_.lpNorm<1>() << endl;
		res.push_back(testGrid->residual().lpNorm<1>() / testGrid->source_.lpNorm<1>());
		testGrid->sor(testGrid->laplaceMat_, testGrid->values_, &testGrid->source_);
	}
	testGrid->print_check_bc_normal_derivs();
	//band matrix plotting.
	
	Eigen::SparseMatrix<double, 1> * matrix = testGrid->laplaceMat_;
	const double* laplaceValues = matrix->valuePtr();
	//column indices of values
	const int* innerValues = matrix->innerIndexPtr();
	//index in values of first nonzero entry of each row, has size of (laplaceMatSize + 1) with last value being the # of nonzeros/end flag.
	const int* outerValues = matrix->outerIndexPtr();
	int valueIdx, innerIdx, outerIdx, rowStartIdx, rowEndIdx;
	int numNonZeros = matrix->nonZeros();
	int inner;
	valueIdx = 0;
	innerIdx = 0;
	outerIdx = 0;

	vector<double> iv, jv;
	for (int i = 0; i < matrix->rows(); i++) {
		//diagCoeff = laplaceMat_->coeff(i, i);
		rowStartIdx = outerValues[outerIdx];
		rowEndIdx = outerValues[outerIdx + 1];
		//sum coeffs*x_i_old on the row i
		for (int j = rowStartIdx; j < rowEndIdx; j++) {
			iv.push_back(i);
			jv.push_back(innerValues[j]);
		}
		outerIdx++;
	}
	writeVectorToTxt(iv, "i.txt");
	writeVectorToTxt(jv, "j.txt");

	//Residual Plotting
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
	
	vector<double> rawRes;
	Eigen::VectorXd resid = testGrid->residual();
	for (int i = 0; i < resid.rows()- 1; i++) {
		//rawRes.push_back(resid.coeff(i));
		rawRes.push_back(testGrid->values_->coeff(i));
	}
	writeVectorToTxt(rawRes, "raw_res.txt");
	cout << " values norm: " << testGrid->values_->lpNorm<1>() << endl;
	

	cout << testGrid->values_->maxCoeff() << endl;
	cout << testGrid->values_->minCoeff() << endl;
	cout << calc_l1_error(testGrid, true) << endl;
	
}
void testGmshDirichletMultigrid() {
	Multigrid mg;

	mg.addGrid(genGmshGrid("finer_mesh.txt", 5, 7, "txt"));
	mg.addGrid(genGmshGrid("fine_mesh.txt", 3, 7, "txt"));
	mg.addGrid(genGmshGrid("inter_mesh.txt", 3, 7, "txt"));
	mg.addGrid(genGmshGrid("coar_mesh.txt", 3, 7, "txt"));
	mg.buildMatrices();
	cout << 1 << endl;
	vector<double> res;
	for (int i = 0; i < 100; i++) {
		res.push_back(mg.residual());
		mg.vCycle();
	}
	cout << calc_l1_error(mg.grids_[mg.grids_.size()-1].second, false) << endl;
	writeVectorToTxt(res, "residual.txt");


}
void testGmshNeumannMultigrid() {
	Multigrid mg;

	mg.addGrid(genGmshGridNeumann("finer_mesh.txt", 5, 7, "txt"));
	mg.addGrid(genGmshGridNeumann("fine_mesh.txt", 3, 7, "txt"));
	mg.addGrid(genGmshGridNeumann("inter_mesh.txt", 3, 7, "txt"));
	mg.addGrid(genGmshGridNeumann("coar_mesh.txt", 3, 7, "txt"));
	mg.buildMatrices();
	vector<double> res;
	for (int i = 0; i < 120; i++) {
		res.push_back(mg.residual());
		mg.vCycle();
	}
	cout << calc_l1_error(mg.grids_[mg.grids_.size() - 1].second, true) << endl;
	writeVectorToTxt(res, "residual.txt");
}
int main() {

	testGmshSingleGrid();
	//testGmshDirichletMultigrid();
	//testGmshNeumannMultigrid();
	return 0;
}
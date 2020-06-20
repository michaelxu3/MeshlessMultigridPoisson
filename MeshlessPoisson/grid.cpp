#include "grid.h"
#include <stdexcept>
#include "math.h"
#include <iostream>
#include <algorithm>
Grid::Grid(vector<std::tuple<double, double, double>> points, vector<std::tuple<vector<int>, int, vector<double>>> boundaries,
										int rbfExp, int polyDeg, int stencilSize, double omega, Eigen::VectorXd source, int sorIters) {
	points_ = points;
	boundaries_ = boundaries;
	properties_.rbfExp = rbfExp;
	properties_.polyDeg = polyDeg;
	laplaceMatSize_ = (int)(points.size());
	int numPoint = (int)(points.size());
	values_ = new Eigen::VectorXd(numPoint);
	values_->setZero();
	laplaceMat_ = new Eigen::SparseMatrix<double>(numPoint, numPoint);
	laplaceMat_->setZero();
	properties_.stencilSize = stencilSize;
	properties_.omega = omega;
	bcFlags_ = std::vector<int>(numPoint);
	properties_.iters = sorIters;
	source_ = source;
}
Grid::~Grid() {
	delete values_;
	delete laplaceMat_;
}

double distance(std::tuple<double, double, double> refPoint, std::tuple<double, double, double> queryPoint) {
	return std::sqrt(std::pow(std::get<0>(refPoint)-std::get<0>(queryPoint),2) + std::pow(std::get<1>(refPoint)-std::get<1>(queryPoint),2) 
		+ std::pow(std::get<2>(refPoint)*std::get<2>(queryPoint),2));
}

void Grid::setBCFlag(int bNum, std::string type, vector<double> boundValues) {
	if (type.compare("dirichlet") != 0 && !type.compare("neumann") != 0) {
		throw std::invalid_argument("Specified Boundary Condition on Boundary #" + std::to_string(bNum) + " is not dirichlet or neumann");
	}
	std::tuple<vector<int>, int, vector<double>> & bound = boundaries_.at(bNum);
	std::get<1>(bound) = type.compare("dirichlet") == 0 ? 1: 2;
	for (size_t i = 0; i < std::get<0>(bound).size(); i++) {
		bcFlags_[std::get<0>(bound)[i]] = std::get<1>(bound);
	}
	std::get<2>(bound) = boundValues;
}

void Grid::boundaryOp() {
	for (size_t i = 0; i < boundaries_.size(); i++) {
		//dirichlet
		if (std::get<1>(boundaries_[i]) == 1){
			for (size_t j = 0; j < std::get<0>(boundaries_[i]).size(); j++) {
				(*values_)(std::get<0>(boundaries_[i]).at(j)) = std::get<2>(boundaries_[i]).at(j);
			}
		}
		//implement neumann later.
	}
}

//Builds Sparse Matrix that approximates laplacian operator for SOR from the interpolation weights.
void Grid::build_laplacian() {
	vector<Eigen::Triplet<double>> tripletList;
	for (int i = 0; i < laplaceMatSize_; i++) {
		std::pair <Eigen::VectorXd, vector<int>> weights = laplaceWeights(i);
		for (size_t j = 0; j < weights.second.size(); j++) {
			
			tripletList.push_back(Eigen::Triplet<double>(i, weights.second[j], weights.first(j)));
		}
	}
	laplaceMat_->setFromTriplets(tripletList.begin(), tripletList.end());
}

void Grid::sor() {
	boundaryOp();
	
	double residual;
	double x_i;
	for (int it = 0; it < properties_.iters; it++) {		// set/enforce the boundary conditions
		//residual = 0;
		for (int i = 0; i < laplaceMatSize_; i++) {
			double x_i_old = values_->coeff(i);
			if (bcFlags_[i] != 0) {
				continue;
			}
			x_i = 0;
			double diagCoeff = laplaceMat_->coeff(i, i);
			x_i += (1 - properties_.omega)*values_->coeff(i) + properties_.omega / diagCoeff * source_.coeff(i);
			for (int j = 0; j < laplaceMatSize_; j++) {
				if (i != j) {
					x_i -= properties_.omega / diagCoeff * laplaceMat_->coeff(i, j)*values_->coeff(j);
				}
			}
			values_->coeffRef(i) = x_i;
			//residual += (x_i - x_i_old);
		}
		//residual = (residual / laplaceMatSize_);
		//std::cout << "residual: " << residual << std::endl;
	}
}
//Fix distance function if we want to extend code to 3D. Makes and sorts a distance array with int point IDs. Brute force, so O(n log n) time
vector<int> Grid::kNearestNeighbors(int pointID) {
	return kNearestNeighbors(points_[pointID]);
}

vector<int> Grid::kNearestNeighbors(std::tuple<double, double, double> refPoint) {
	vector<std::pair<double, int>> distances;
	for (size_t i = 0; i < points_.size(); i++) {
		std::tuple<double, double, double> queryPoint = points_[i];
		distances.push_back(std::pair<double, int>(distance(refPoint, queryPoint), i));
	}
	//std::sort sorts pairs by first (distance)
	std::sort(distances.begin(), distances.end());

	vector<int> nearestNeighbors;
	//Have to include point itself in the stencil since otherwise diag would be zeros.
	for (int i = 0; i < properties_.stencilSize; i++) {
		nearestNeighbors.push_back(distances[i].second);
	}
	return nearestNeighbors;
}
//Builds PHS coefficients matrix for a point
//returns pair with first = coeff matrix and second = list of neighbors, to save time.
std::pair<Eigen::MatrixXd, vector<int>> Grid::buildCoeffMatrix(std::tuple<double, double, double> point) {
	int polyTerms = (properties_.polyDeg + 1)*(properties_.polyDeg + 2) / 2;
	vector<int> neighbors = kNearestNeighbors(point);
	//build A-matrix
	Eigen::MatrixXd coeff_mat = Eigen::MatrixXd::Zero(properties_.stencilSize + polyTerms, properties_.stencilSize + polyTerms);
	for (int i = 0; i < properties_.stencilSize; i++) {
		for (int j = 0; j < properties_.stencilSize; j++) {
			double r = distance(points_[neighbors[i]], points_[neighbors[j]]);
			coeff_mat(i, j) = std::pow(r, properties_.rbfExp);
		}
	}
	//fill P-matrix sections.
	int colIndex;
	for (int row = 0; row < properties_.stencilSize; row++) {
		colIndex = properties_.stencilSize;
		for (int p = 0; p <= properties_.polyDeg; p++) {
			for (int q = 0; q <= p; q++) {
				double x = std::get<0>(points_.at(neighbors[row]));
				double y = std::get<1>(points_.at(neighbors[row]));
				coeff_mat(row, colIndex) = std::pow(x, p-q)*std::pow(y, q);
				coeff_mat(colIndex, row) = std::pow(x, p-q)*std::pow(y, q);
				colIndex++;
			}
		}
	}
	//std::cout << "Condition number: " << coeff_mat.norm()*coeff_mat.inverse().norm() << std::endl;
	//Rest of the matrix should automatically be zeros.
	return std::pair<Eigen::MatrixXd, vector<int>>(coeff_mat, neighbors);
}

std::pair<Eigen::MatrixXd, vector<int>> Grid::buildCoeffMatrix(int pointID) {
	return buildCoeffMatrix(points_[pointID]);
}

//Builds RHS and solves for the laplace weights for a point
std::pair<Eigen::VectorXd, vector<int>> Grid::laplaceWeights (int pointID) {

	std::pair<Eigen::MatrixXd, vector<int>> coeffs = buildCoeffMatrix(pointID);
	vector<int> neighbors = coeffs.second;
	int polyTerms = (properties_.polyDeg + 1)*(properties_.polyDeg + 2) / 2;
	Eigen::VectorXd rhs = Eigen::VectorXd::Zero(properties_.stencilSize + polyTerms);
	//Build RHS
	double xEval = std::get<0>(points_[pointID]);
	double yEval = std::get<1>(points_[pointID]);
	for (int i = 0; i < properties_.stencilSize; i++) {
		double delX = std::get<0>(points_[neighbors[i]]) - xEval;
		double delY = std::get<1>(points_[neighbors[i]]) - yEval;
		rhs(i) = properties_.rbfExp * properties_.rbfExp*std::pow(delX*delX + delY * delY, (double)(properties_.rbfExp) / 2.0 - 1);
	}
	//poly terms
	int rowIndex = properties_.stencilSize;
	for (int p = 0; p <= properties_.polyDeg; p++) {
		for (int q = 0; q <= p; q++) {
			double laplacePoly = 0;
			if (p - q - 2 >= 0) {
				laplacePoly += (p - q)*(p - q - 1)*std::pow(xEval, p - q - 2)*std::pow(yEval, q);
			}
			if (q - 2 >= 0) {
				laplacePoly += q * (q - 1)*std::pow(xEval, p - q)*std::pow(yEval, q - 2);
			}
			rhs(rowIndex) = laplacePoly;
			rowIndex++;
		}
	}
	//maybe add a condition number check here to use fullpiv or partialpiv
	Eigen::VectorXd weights = coeffs.first.fullPivLu().solve(rhs);
	return std::pair<Eigen::VectorXd, vector<int>>(weights, neighbors);
}
std::pair<Eigen::VectorXd, vector<int>> Grid::pointInterpWeights(std::tuple<double, double, double> point) {
	int polyTerms = (getPolyDeg() + 1)*(getPolyDeg() + 2);
	vector<int> neighbors = kNearestNeighbors(point);
	int matSize = getStencilSize() + polyTerms;
	Eigen::MatrixXd coeff_mat = Eigen::MatrixXd::Zero(matSize, matSize);
	//A-matrix
	for (int i = 0; i < properties_.stencilSize; i++) {
		for (int j = 0; j < properties_.stencilSize; j++) {
			double r = distance(points_[neighbors[i]], points_[neighbors[j]]);
			coeff_mat(i, j) = std::pow(r, properties_.rbfExp);
		}
	}
	//fill P-matrix sections.
	int colIndex;
	for (int row = 0; row < properties_.stencilSize; row++) {
		colIndex = properties_.stencilSize;
		for (int p = 0; p <= properties_.polyDeg; p++) {
			for (int q = 0; q <= p; q++) {
				double x = std::get<0>(points_.at(neighbors[row]));
				double y = std::get<1>(points_.at(neighbors[row]));
				coeff_mat(row, colIndex) = std::pow(x, p - q)*std::pow(y, q);
				coeff_mat(colIndex, row) = std::pow(x, p - q)*std::pow(y, q);
				colIndex++;
			}
		}
	}
	Eigen::VectorXd rhs = Eigen::VectorXd::Zero(properties_.stencilSize + polyTerms);
	//Build RHS rbf terms
	double xEval = std::get<0>(point);
	double yEval = std::get<1>(point);
	for (int i = 0; i < properties_.stencilSize; i++) {
		rhs(i) = std::pow(distance(point, points_[neighbors[i]]), properties_.rbfExp);
	}
	//RHS poly terms
	int rowIndex = properties_.stencilSize;
	for (int p = 0; p <= properties_.polyDeg; p++) {
		for (int q = 0; q <= p; q++) {
			rhs(rowIndex) = std::pow(xEval, p - q)*std::pow(yEval, q);
			rowIndex++;
		}
	}
	Eigen::VectorXd weights = coeff_mat.fullPivLu().solve(rhs);
	return std::pair<Eigen::VectorXd, vector<int>>(weights, neighbors);
}
int Grid::getSize() {
	return laplaceMatSize_;
}
int Grid::getStencilSize() {
	return properties_.stencilSize;
}
int Grid::getPolyDeg() {
	return properties_.polyDeg;
}
std::tuple<double, double, double> Grid::getPoint(int index) {
	return points_[index];
}
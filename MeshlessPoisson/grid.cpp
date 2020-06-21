#include "grid.h"
#include <stdexcept>
#include "math.h"
#include <iostream>
#include <algorithm>
Grid::Grid(vector<std::tuple<double, double, double>> points, vector<Boundary> boundaries ,
					GridProperties properties, Eigen::VectorXd source) {
	int numPoint = (int)(points.size());
	points_ = points;
	boundaries_ = boundaries;
	bcFlags_ = std::vector<int>(numPoint);
	properties_ = properties;
	source_ = source;



	laplaceMatSize_ = (int)(points.size());
	values_ = new Eigen::VectorXd(numPoint);
	values_->setZero();
	laplaceMat_ = new Eigen::SparseMatrix<double, Eigen::RowMajor>(numPoint, numPoint);
	laplaceMat_->setZero();
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
	Boundary & bound = boundaries_.at(bNum);
	(bound.type) = type.compare("dirichlet") == 0 ? 1: 2;
	for (size_t i = 0; i < (bound.bcPoints).size(); i++) {
		bcFlags_[(bound.bcPoints)[i]] = (bound.type);
	}
	(bound.values) = boundValues;
}

void Grid::boundaryOp() {
	for (size_t i = 0; i < boundaries_.size(); i++) {
		//dirichlet
		if ((boundaries_[i]).type == 1){
			for (size_t j = 0; j < (boundaries_[i].bcPoints).size(); j++) {
				(*values_)((boundaries_[i].bcPoints).at(j)) = (boundaries_[i].values).at(j);
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
	laplaceMat_->makeCompressed();
}

void Grid::sor() {
	boundaryOp();
	
	double x_i, diagCoeff;

	const double* laplaceValues = laplaceMat_->valuePtr();
	//column indices of values
	const int* innerValues = laplaceMat_->innerIndexPtr();
	//index in values of first nonzero entry of each row, has size of (laplaceMatSize + 1) with last value being the # of nonzeros/end flag.
	const int* outerValues = laplaceMat_->outerIndexPtr();

	int valueIdx, innerIdx, outerIdx, rowStartIdx, rowEndIdx;
	int numNonZeros = laplaceMat_->nonZeros();

	for (int it = 0; it < properties_.iters; it++) {
		valueIdx = 0;
		innerIdx = 0;
		outerIdx = 0;

		for (int i = 0; i < laplaceMatSize_; i++){
			
			if (bcFlags_[i] != 0) {
				outerIdx++;
				continue;
			}
		
			x_i = 0;
			diagCoeff = 0;
			rowStartIdx = outerValues[outerIdx];
			rowEndIdx = outerValues[outerIdx + 1];
			//sum coeffs*x_i_old on the row i
			for (int j = rowStartIdx; j < rowEndIdx; j++) {
				if (innerValues[j] == i) {
					diagCoeff = laplaceValues[j];
					continue;
				}
				x_i -= laplaceValues[j] * values_->coeff(innerValues[j]);
			}
			x_i += source_.coeff(i);
			x_i *= properties_.omega / diagCoeff;
			x_i += (1 - properties_.omega)*values_->coeff(i);
			values_->coeffRef(i) = x_i;
			outerIdx++;
		}
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
	//std::sort(distances.begin(), distances.end());

	
	//Partial Selection Sort algorithm, O(kn) time, somehow takes more time than std::sort on VS.
	for (int i = 0; i < properties_.stencilSize; i++) {
		int minIndex = i;
		std::pair<double, int> minDistPoint = distances[i];
		for (size_t j = i+1; j < points_.size(); j++) {
			if (distances[j] < minDistPoint) {
				minIndex = j;
				minDistPoint = distances[j];
				std::swap(distances[i], distances[minIndex]);
			}
		}
	}
	
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
	//Exploit Symmetry and Diagonals being zero.
	double a_coeff;
	for (int i = 0; i < properties_.stencilSize; i++) {
		for (int j = i; j < properties_.stencilSize; j++) {
			double r = distance(points_[neighbors[i]], points_[neighbors[j]]);
			a_coeff = std::pow(r, properties_.rbfExp);
			coeff_mat(i, j) = a_coeff;
			coeff_mat(j, i) = a_coeff;
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
	Eigen::VectorXd weights = coeffs.first.partialPivLu().solve(rhs);
	return std::pair<Eigen::VectorXd, vector<int>>(weights, neighbors);
}
std::pair<Eigen::VectorXd, vector<int>> Grid::pointInterpWeights(std::tuple<double, double, double> point) {
	int polyTerms = (getPolyDeg() + 1)*(getPolyDeg() + 2);
	vector<int> neighbors = kNearestNeighbors(point);
	int matSize = getStencilSize() + polyTerms;
	Eigen::MatrixXd coeff_mat = Eigen::MatrixXd::Zero(matSize, matSize);
	//A-matrix
	double a_coeff;
	for (int i = 0; i < properties_.stencilSize; i++) {
		for (int j = i; j < properties_.stencilSize; j++) {
			double r = distance(points_[neighbors[i]], points_[neighbors[j]]);
			a_coeff = std::pow(r, properties_.rbfExp);
			coeff_mat(i, j) = a_coeff;
			coeff_mat(j, i) = a_coeff;
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
	Eigen::VectorXd weights = coeff_mat.partialPivLu().solve(rhs);
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
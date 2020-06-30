#include "grid.h"
#include <stdexcept>
#include "math.h"
#include <iostream>
#include <algorithm>

using std::cout;
using std::endl;
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
	residuals_ = new Eigen::VectorXd(numPoint);
	residuals_->setZero();
	laplaceMat_ = new Eigen::SparseMatrix<double, Eigen::RowMajor>(numPoint, numPoint);
	laplaceMat_->setZero();
}
Grid::~Grid() {
	delete values_;
	delete laplaceMat_;
	delete residuals_;
}

double distance(std::tuple<double, double, double> refPoint, std::tuple<double, double, double> queryPoint) {
	return std::sqrt(std::pow(std::get<0>(refPoint)-std::get<0>(queryPoint),2) + std::pow(std::get<1>(refPoint)-std::get<1>(queryPoint),2) 
		+ std::pow(std::get<2>(refPoint)-std::get<2>(queryPoint),2));
}

void Grid::setBCFlag(int bNum, std::string type, vector<double> boundValues) {
	Boundary & bound = boundaries_.at(bNum);
	(bound.type) = type.compare("dirichlet") == 0 ? 1 : 2;
	for (size_t i = 0; i < (bound.bcPoints).size(); i++) {
		bcFlags_[(bound.bcPoints)[i]] = (bound.type);
	}
	(bound.values) = boundValues;
}

void Grid::boundaryOp(std::string coarse) {
	for (size_t i = 0; i < boundaries_.size(); i++) {
		//dirichlet
		if ((boundaries_[i]).type == 1){
			for (size_t j = 0; j < (boundaries_[i].bcPoints).size(); j++) {
				(*values_)((boundaries_[i].bcPoints).at(j)) = coarse.compare("coarse") == 0 ? 0 : (boundaries_[i].values).at(j);
			}
		}
		//implement neumann later.
		if (boundaries_[i].type == 2) {
			if (coarse.compare("coarse") == 0) {
				for (int i = 0; i < laplaceMat_->rows(); i++) {
					for (size_t j = 0; j < (boundaries_[i].bcPoints).size(); j++) {
						Eigen::VectorXd newCoeff = -(laplaceMat_->col(i))*(*values_);
					}
				}
			}
			else {

			}
		}
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
	laplaceMat_->makeCompressed();;
}

void Grid::modifyCoeffNeumann() {
}
void Grid::sor(Eigen::SparseMatrix<double, 1>* matrix, Eigen::VectorXd* values, Eigen::VectorXd* rhs) {
	double x_i, diagCoeff;

	const double* laplaceValues = matrix->valuePtr();
	//column indices of values
	const int* innerValues = matrix->innerIndexPtr();
	//index in values of first nonzero entry of each row, has size of (laplaceMatSize + 1) with last value being the # of nonzeros/end flag.
	const int* outerValues = matrix->outerIndexPtr();

	int valueIdx, innerIdx, outerIdx, rowStartIdx, rowEndIdx;
	int numNonZeros = matrix->nonZeros();

	for (int it = 0; it < properties_.iters; it++) {
		valueIdx = 0;
		innerIdx = 0;
		outerIdx = 0;

		for (int i = 0; i < matrix->rows(); i++){
			//dirichlet bc
			if (bcFlags_[i] == 1) {
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
				x_i -= laplaceValues[j] * values->coeff(innerValues[j]);
			}
			//residual here
			//take last residual and restrict
			//source - ap phi_nb
			x_i += rhs->coeff(i);

			x_i *= properties_.omega / diagCoeff;
			x_i += (1 - properties_.omega)*values->coeff(i);
			values->coeffRef(i) = x_i;
			outerIdx++;
		}
	}
}
Eigen::VectorXd Grid::residual() {
	return source_ - (*laplaceMat_) * (*values_);
}
//Fix distance function if we want to extend code to 3D. Makes and sorts a distance array with int point IDs. Brute force, so O(n log n) time
vector<int> Grid::kNearestNeighbors(int pointID) {
	return kNearestNeighbors(points_[pointID]);
}

vector<int> Grid::kNearestNeighbors(std::tuple<double, double, double> refPoint) {
	vector<std::pair<double, int>> distances;
	for (int i = 0; i < laplaceMatSize_; i++) {
		std::tuple<double, double, double> queryPoint = points_[i];
		distances.push_back(std::pair<double, int>(distance(refPoint, queryPoint), i));
	}

	
	//max-heap selection algorithm
	int k = properties_.stencilSize;
	vector<std::pair<double, int>> maxHeap;
	for (int i = 0; i < k; i++) {
		maxHeap.push_back(distances[i]);
	}
	std::make_heap(maxHeap.begin(), maxHeap.end());
	for (int i = k; i < laplaceMatSize_; i++) {
		if (distances[i] < maxHeap.front()) {
			std::pop_heap(maxHeap.begin(), maxHeap.end());
			maxHeap.pop_back();

			maxHeap.push_back(distances[i]);
			std::push_heap(maxHeap.begin(), maxHeap.end());
		}
	}
	std::sort_heap(maxHeap.begin(), maxHeap.end());
	vector<int> nearestNeighbors;
	//Have to include point itself in the stencil since otherwise diag would be zeros.
	for (int i = 0; i < k; i++) {
		nearestNeighbors.push_back(maxHeap[i].second);
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
	double yRef, xRef, D;
	double L = (double)(properties_.rbfExp);
	for (int i = 0; i < properties_.stencilSize; i++) {
		xRef = std::get<0>(points_[neighbors[i]]);
		yRef = std::get<1>(points_[neighbors[i]]);
		D = (xEval*xEval - 2 * xEval*xRef + xRef * xRef + yEval * yEval - 2 * yEval*yRef + yRef * yRef);
		if (D > 0) {
			rhs(i) = (std::pow(2 * xEval - 2 * xRef, 2) + std::pow(2 * yEval - 2 * yRef, 2))*(L / 2)*(L / 2 - 1)*std::pow(D, L / 2 - 2)
				+ 2 * L*std::pow(D, L / 2 - 1);
		}
		//cout << rhs(i) << endl;
		
		//double delX = std::get<0>(points_[neighbors[i]]) - xEval;
		//double delY = std::get<1>(points_[neighbors[i]]) - yEval;
		//rhs(i) = properties_.rbfExp * properties_.rbfExp*std::pow(delX*delX + delY * delY, (double)(properties_.rbfExp) / 2.0 - 1);
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
	int polyTerms = (getPolyDeg() + 1)*(getPolyDeg() + 2)/2;
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
	//std::cout << "interp condition number: " << coeff_mat.norm() * coeff_mat.inverse().norm() << std::endl;
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
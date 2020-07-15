#include "grid.h"
#include <stdexcept>
#include "math.h"
#include <iostream>
#include <algorithm>
#include <queue>
#define pi 3.141592653589793238462643383279
using std::cout;
using std::endl;
Grid::Grid(vector<Point> points, vector<Boundary> boundaries ,
					GridProperties properties, Eigen::VectorXd source) {
	int numPoint = (int)(points.size());
	points_ = points;
	boundaries_ = boundaries;
	properties_ = properties;
	source_ = source;
	cond_rbf = vector<double>();
	neumannFlag_ = false;
	setNeumannFlag();

	laplaceMatSize_ = numPoint;
	int A_size = neumannFlag_ ? numPoint + 1 : numPoint;
	bcFlags_ = std::vector<int>(numPoint);
	laplaceMat_ = new Eigen::SparseMatrix<double, Eigen::RowMajor>(A_size, A_size);
	laplaceMat_->setZero();
	values_ = new Eigen::VectorXd(A_size);
	values_->setZero();

}
Grid::~Grid() {
	delete values_;
	delete laplaceMat_;
}

double distance(Point refPoint, Point queryPoint) {
	return std::sqrt(std::pow(std::get<0>(refPoint)-std::get<0>(queryPoint),2) + std::pow(std::get<1>(refPoint)-std::get<1>(queryPoint),2));
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
	}
}
void Grid::setNeumannFlag() {
	for (size_t i = 0; i < boundaries_.size(); i++) {
		if (boundaries_[i].type == 2) {
			neumannFlag_ = true;
			return;
		}
	}
	neumannFlag_ = false;
}
//This is inefficient and just written so I don't have to worry about sparse matrix traversal for now.
//Fix once single grid converges.
void Grid::modify_coeff_neumann() {
	
	/*
	vector<int> neighbors;
	double A_ij, A_jj, A_ik, A_jk;
	for (int i = 0; i < laplaceMat_->rows() - 1; i++) {
		if (bcFlags_[i] != 0) {
			continue;
		}
		for (int j = 0; j < laplaceMat_->cols() - 1; j++) {
			//cout << j << endl;
			if (bcFlags_[j] == 2) {
				neighbors = kNearestNeighbors(j, true);
				A_jj = laplaceMat_->coeff(j, j);
				A_ij = laplaceMat_->coeff(i, j);
			
				for (int k = 0; k < neighbors.size(); k++) {
					A_ik = laplaceMat_->coeff(i, neighbors[k]);
					A_jk = laplaceMat_->coeff(j, neighbors[k]);
					//cout << A_jj << endl;
					laplaceMat_->coeffRef(i,neighbors[k]) = A_ik - A_ij * A_jk / A_jj;
				}
				laplaceMat_->coeffRef(i, j) = 0;
			}
		}
	}
	*/
	
	for (int i = 0; i < source_.rows() - 1; i++) {
		if (bcFlags_[i] == 2) {
			source_(i) = 0;
		}
	}
	laplaceMat_->makeCompressed();
}
void Grid::bound_eval_neumann() {
	const double* laplaceValues = laplaceMat_->valuePtr();
	//column indices of values
	const int* innerValues = laplaceMat_->innerIndexPtr();
	//index in values of first nonzero entry of each row, has size of (laplaceMatSize + 1) with last value being the # of nonzeros/end flag.
	const int* outerValues = laplaceMat_->outerIndexPtr();
	int rowStart, rowEnd, curr;
	double diag, boundValue;
	for (size_t i = 0; i < boundaries_.size(); i++) {
		if ((boundaries_[i]).type == 2) {
			for (size_t j = 0; j < (boundaries_[i].bcPoints).size(); j++) {
				curr = boundaries_[i].bcPoints[j];
				rowStart = outerValues[curr];
				rowEnd = outerValues[curr + 1];
				boundValue = 0;
				//cout << "initial boundvalue: " << boundValue << endl;
				for (int k = rowStart; k < rowEnd; k++) {
					if (innerValues[k] == curr) {
						diag = laplaceValues[k];
						continue;
					}
					boundValue -= values_->coeff(innerValues[k]) * laplaceValues[k];
				}
				boundValue /= diag;
				values_->coeffRef(curr) = boundValue;
				
			}
		}
	}
	//print_check_bc_normal_derivs();
}
void Grid::sor(Eigen::SparseMatrix<double, Eigen::RowMajor>* matrix, Eigen::VectorXd* values, Eigen::VectorXd* rhs) {
	double x_i, diagCoeff;
	const double* laplaceValues = matrix->valuePtr();
	//column indices of values
	const int* innerValues = matrix->innerIndexPtr();
	//index in values of first nonzero entry of each row, has size of (laplaceMatSize + 1) with last value being the # of nonzeros/end flag.
	const int* outerValues = matrix->outerIndexPtr();
	int valueIdx, innerIdx, outerIdx, rowStartIdx, rowEndIdx;
	int numNonZeros = matrix->nonZeros();
	int inner;
	for (int it = 0; it < properties_.iters; it++) {
		valueIdx = 0;
		innerIdx = 0;
		outerIdx = 0;

		for (int i = 0; i < matrix->rows(); i++){
			//dirichlet bc
			if (!(neumannFlag_ && i == matrix->rows()-1) && bcFlags_[i] != 0) {
				outerIdx++;
				continue;
			}
			x_i = 0;
			diagCoeff = 0;
			//diagCoeff = laplaceMat_->coeff(i, i);
			rowStartIdx = outerValues[outerIdx];
			rowEndIdx = outerValues[outerIdx + 1];
			//sum coeffs*x_i_old on the row i
			for (int j = rowStartIdx; j < rowEndIdx; j++) {
				if (innerValues[j] == i) {
					diagCoeff = laplaceValues[j];
					continue;
				}
				
				inner = innerValues[j];
				x_i -= laplaceValues[j] * values->coeff(inner);
			}
			//residual here
			//take last residual and restrict
			//source - ap phi_nb
			x_i += rhs->coeff(i);
	
			x_i *= properties_.omega / diagCoeff;
			x_i += (1 - properties_.omega)*values->coeff(i);
			//cout << "diag: " << diagCoeff << " bcflag: " << bcFlags_[i] << endl;
			//cout << x_i << endl;
			values->coeffRef(i) = x_i;
			outerIdx++;
		}
		bound_eval_neumann();
	}
}
void Grid::directSolve() {
	
	//modifyCoeffDirichlet();
	Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;
	solver.analyzePattern(*laplaceMat_);
	solver.factorize(*laplaceMat_);
	*values_ = solver.solve(source_);
	
}
Eigen::VectorXd Grid::residual() { 
	Eigen::VectorXd res = source_ - (*laplaceMat_) * (*values_);
	fix_vector_bound_coarse(&res);
	return res;
}
void Grid::print_bc_values(Eigen::VectorXd vec) {
	for (size_t i = 0; i < boundaries_.size(); i++) {
		if ((boundaries_[i]).type == 1) {
			for (size_t j = 0; j < (boundaries_[i].bcPoints).size(); j++) {
				cout << "bc value: " << vec.coeff(boundaries_[i].bcPoints.at(j)) << endl;
			}
		}
	}
}
void Grid::print_check_bc_normal_derivs() {
	double* laplaceValues = laplaceMat_->valuePtr();
	//column indices of values
	const int* innerValues = laplaceMat_->innerIndexPtr();
	//index in values of first nonzero entry of each row, has size of (laplaceMatSize + 1) with last value being the # of nonzeros/end flag.
	const int* outerValues = laplaceMat_->outerIndexPtr();
	int valueIdx, innerIdx, outerIdx, rowStartIdx, rowEndIdx;
	double sum;
	valueIdx = 0;
	innerIdx = 0;
	outerIdx = 0;

	for (int i = 0; i < laplaceMat_->rows() - 1; i++) {
		if (bcFlags_[i] != 2) {
			outerIdx++;
			continue;
		}
		sum = 0;
		rowStartIdx = outerValues[outerIdx];
		rowEndIdx = outerValues[outerIdx + 1];
		for (int j = rowStartIdx; j < rowEndIdx; j++) {
			sum += laplaceValues[j] * values_->coeff(innerValues[j]);
		}
		cout << "d/dn: " << sum << " rhs: " << source_.coeff(i) << endl;
		outerIdx++;
	}
}
void Grid::fix_vector_bound_coarse(Eigen::VectorXd* vec) {
	for (size_t i = 0; i < boundaries_.size(); i++) {
		//dirichlet
		if ((boundaries_[i]).type == 1) {
			for (size_t j = 0; j < (boundaries_[i].bcPoints).size(); j++) {
				vec->coeffRef(boundaries_[i].bcPoints.at(j)) = 0;
			}
		}
	}
}
//Fix distance function if we want to extend code to 3D. Makes and sorts a distance array with int point IDs. Brute force, so O(n log n) time
vector<int> Grid::kNearestNeighbors(int pointID, bool neumann) {
	return kNearestNeighbors(points_[pointID], neumann, (bcFlags_[pointID] != 0));
}
// THIS IS BUGGED AND RESULTS IN RCM GIVING N - 1 order, MEANING WE DONT VISIT ALL POINTS.
vector<int> Grid::kNearestNeighbors(Point refPoint, bool neumann, bool pointBCFlag) {
	vector<std::pair<double, int>> distances;
	Point queryPoint;
	int samePoint = -1;
	for (int i = 0; i < laplaceMatSize_; i++) {
		queryPoint = points_[i];
		distances.push_back(std::pair<double, int>(distance(refPoint, queryPoint), i));
		
		if (distances[i].first == 0) {
			samePoint = i;
		}
		
	}

	
	//max-heap selection algorithm
	int k = properties_.stencilSize;
	int lastInitIndex = k;
	vector<std::pair<double, int>> maxHeap;
	for (int i = 0; i < lastInitIndex; i++) {
		if (i != samePoint && (pointBCFlag && neumann && bcFlags_[i] != 0)) {
			lastInitIndex++;
			continue;
		}
		maxHeap.push_back(distances[i]);
	}
	std::make_heap(maxHeap.begin(), maxHeap.end());
	//cout << maxHeap.size() << " " << properties_.stencilSize << endl;
	for (int i = lastInitIndex; i < laplaceMatSize_; i++) {
		//cout << i << endl;
		
		if (i == samePoint || !(pointBCFlag && neumann && bcFlags_[i] != 0) && distances[i] < maxHeap.front()) {
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
	
	/*
	double x, y;
	int numBoundPoints = 0;
	x = std::get<0>(refPoint);
	y = std::get<1>(refPoint);
	if (true) {
		for (int i = 0; i < properties_.stencilSize; i++) {
			if (bcFlags_[nearestNeighbors[i]] != 0) {
				numBoundPoints++;
			}
		}
	}
	
	//cout << "maxheap size: " << maxHeap.size() << "nearest neighbors size: " << nearestNeighbors.size() << endl;

	cout << " num bc points: " << numBoundPoints << endl;
	*/
	return nearestNeighbors;
	
}
//Builds PHS coefficients matrix for a point
//returns pair with first = coeff matrix and second = list of neighbors, to save time.
std::tuple<Eigen::MatrixXd, vector<int>, vector<Point>> Grid::buildCoeffMatrix(Point point, bool neumann, bool pointBCFlag) {
	int polyTerms = (properties_.polyDeg + 1)*(properties_.polyDeg + 2) / 2;
	vector<int> neighbors = kNearestNeighbors(point, neumann, pointBCFlag);
	vector<Point> scaledPoints = shifting_scaling(neighbors, point);
	//build A-matrix
	Eigen::MatrixXd coeff_mat = Eigen::MatrixXd::Zero(properties_.stencilSize + polyTerms, properties_.stencilSize + polyTerms);
	
	double a_coeff, r;
	for (int i = 0; i < properties_.stencilSize; i++) {
		for (int j = i; j < properties_.stencilSize; j++) {
			r = distance(scaledPoints[i], scaledPoints[j]);
			a_coeff = std::pow(r, properties_.rbfExp);
			coeff_mat(i, j) = a_coeff;
			coeff_mat(j, i) = a_coeff;
		}
	}
	//fill P-matrix sections.
	int colIndex;
	double x, y, p_coeff;
	for (int row = 0; row < properties_.stencilSize; row++) {
		colIndex = properties_.stencilSize;
		for (int p = 0; p <= properties_.polyDeg; p++) {
			for (int q = 0; q <= p; q++) {
				//cout << "x " << p - q << "y " << q << endl;
				x = std::get<0>(scaledPoints[row]);
				y = std::get<1>(scaledPoints[row]);
				p_coeff = std::pow(x, p - q)*std::pow(y, q);
				coeff_mat(row, colIndex) = p_coeff;
				coeff_mat(colIndex, row) = p_coeff;
				colIndex++;
			}
		}
	}
	//cout << coeff_mat << "\n" << endl;
	//double cond = coeff_mat.lpNorm<2>()*coeff_mat.inverse().lpNorm<2>();
	//cout << cond << endl;
	//cond_rbf.push_back(cond);
	return std::tuple<Eigen::MatrixXd, vector<int>, vector<Point>>(coeff_mat, neighbors, scaledPoints);
}

std::tuple<Eigen::MatrixXd, vector<int>, vector<Point>> Grid::buildCoeffMatrix(int pointID, bool neumann) {
	return buildCoeffMatrix(points_[pointID], neumann, (bcFlags_[pointID] != 0));
}
std::pair<Eigen::VectorXd, vector<int>> Grid::derivx_weights(int pointID) {
	std::tuple<Eigen::MatrixXd, vector<int>, vector<Point>> coeffs = buildCoeffMatrix(pointID, neumannFlag_);
	vector<int> neighbors = std::get<1>(coeffs);
	vector<Point> scaledPoints = std::get<2>(coeffs);
	int polyTerms = (properties_.polyDeg + 1)*(properties_.polyDeg + 2) / 2;
	Eigen::VectorXd rhs = Eigen::VectorXd::Zero(properties_.stencilSize + polyTerms);

	//Build RHS. xRef, yRef, xEval, yEval are already shifted and scaled.
	Point evalPoint = scaledPoints.at(scaledPoints.size() - 1);
	double xEval = std::get<0>(evalPoint);
	double yEval = std::get<1>(evalPoint);
	double yRef, xRef;
	double M = (double)(properties_.rbfExp);
	double coeff;
	for (int i = 0; i < properties_.stencilSize; i++) {
		xRef = std::get<0>(scaledPoints[i]);
		yRef = std::get<1>(scaledPoints[i]);
		if (i > 0) {
			rhs(i) = M * std::pow(distance(scaledPoints[i], evalPoint), M - 2) * (xEval - xRef);
		}
	}
	//poly terms
	int rowIndex = properties_.stencilSize;
	for (int p = 0; p <= properties_.polyDeg; p++) {
		for (int q = 0; q <= p; q++) {
			double laplacePoly = 0;
			if (p - q - 1 >= 0) {
				laplacePoly += (p - q)*std::pow(xEval, p - q - 1)*std::pow(yEval, q);
			}
			rhs(rowIndex) = laplacePoly;
			rowIndex++;
		}
	}
	//maybe add a condition number check here to use fullpiv or partialpiv
	Eigen::VectorXd weights = std::get<0>(coeffs).fullPivLu().solve(rhs);
	//cout << std::get<0>(coeffs) << endl;
	//cout << rhs << endl;
	double scale = std::get<0>(scaledPoints[scaledPoints.size() - 2]);
	for (int i = 0; i < weights.rows(); i++) {
		weights(i) /= scale;
	}
	//cout << weights << endl;
	return std::pair<Eigen::VectorXd, vector<int>>(weights, neighbors);
}
std::pair<Eigen::VectorXd, vector<int>> Grid::derivy_weights(int pointID) {
	std::tuple<Eigen::MatrixXd, vector<int>, vector<Point>> coeffs = buildCoeffMatrix(pointID, neumannFlag_);
	vector<int> neighbors = std::get<1>(coeffs);
	vector<Point> scaledPoints = std::get<2>(coeffs);
	int polyTerms = (properties_.polyDeg + 1)*(properties_.polyDeg + 2) / 2;
	Eigen::VectorXd rhs = Eigen::VectorXd::Zero(properties_.stencilSize + polyTerms);

	//Build RHS. xRef, yRef, xEval, yEval are already shifted and scaled.
	Point evalPoint = scaledPoints.at(scaledPoints.size() - 1);
	double xEval = std::get<0>(evalPoint);
	double yEval = std::get<1>(evalPoint);
	double yRef, xRef;
	double M = (double)(properties_.rbfExp);
	double coeff;
	for (int i = 0; i < properties_.stencilSize; i++) {
		xRef = std::get<0>(scaledPoints[i]);
		yRef = std::get<1>(scaledPoints[i]);
		if (i > 0) {
			rhs(i) = M * std::pow(distance(scaledPoints[i], evalPoint), M - 2) * (yEval - yRef);
		}
	}
	//poly terms
	int rowIndex = properties_.stencilSize;
	for (int p = 0; p <= properties_.polyDeg; p++) {
		for (int q = 0; q <= p; q++) {
			double laplacePoly = 0;
			if (q - 1 >= 0) {
				laplacePoly += q*std::pow(xEval, p - q)*std::pow(yEval, q-1);
			}
			rhs(rowIndex) = laplacePoly;
			rowIndex++;
		}
	}
	//maybe add a condition number check here to use fullpiv or partialpiv
	Eigen::VectorXd weights = std::get<0>(coeffs).fullPivLu().solve(rhs);
	//cout << std::get<0>(coeffs) << endl;
	//cout << rhs << endl;
	double scale = std::get<0>(scaledPoints[scaledPoints.size() - 2]);
	for (int i = 0; i < weights.rows(); i++) {
		weights(i) /= scale;
	}
	//cout << weights << endl;
	return std::pair<Eigen::VectorXd, vector<int>>(weights, neighbors);
}
//Builds RHS and solves for the laplace weights for a point
std::pair<Eigen::VectorXd, vector<int>> Grid::laplaceWeights(int pointID) {
	std::tuple<Eigen::MatrixXd, vector<int>, vector<Point>> coeffs = buildCoeffMatrix(pointID, neumannFlag_);
	vector<int> neighbors = std::get<1>(coeffs);
	vector<Point> scaledPoints = std::get<2>(coeffs);
	int polyTerms = (properties_.polyDeg + 1)*(properties_.polyDeg + 2) / 2;
	Eigen::VectorXd rhs = Eigen::VectorXd::Zero(properties_.stencilSize + polyTerms);
	 
	//Build RHS. xRef, yRef, xEval, yEval are already shifted and scaled.
	Point evalPoint = scaledPoints.at(scaledPoints.size() - 1);
	double xEval = std::get<0>(evalPoint);
	double yEval = std::get<1>(evalPoint);
	double yRef, xRef, D;
	double M = (double)(properties_.rbfExp);
	double coeff;
	for (int i = 0; i < properties_.stencilSize; i++) {
		xRef = std::get<0>(scaledPoints[i]);
		yRef = std::get<1>(scaledPoints[i]);
		D = (xEval*xEval - 2 * xEval*xRef + xRef * xRef + yEval * yEval - 2 * yEval*yRef + yRef * yRef);
		if (D > 0) {
			rhs(i) = (std::pow(2 * xEval - 2 * xRef, 2) + std::pow(2 * yEval - 2 * yRef, 2))*(M / 2)*(M / 2 - 1)*std::pow(D, M / 2 - 2)
				+ 2 * M*std::pow(D, M / 2 - 1);
		}
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
	Eigen::VectorXd weights = std::get<0>(coeffs).fullPivLu().solve(rhs);
	//cout << std::get<0>(coeffs) << endl;
	//cout << rhs << endl;
	double scale = std::get<0>(scaledPoints[scaledPoints.size() - 2]);
	for (int i = 0; i < weights.rows(); i++) {
		weights(i) /= std::pow(scale, 2);
	}
	//cout << weights << endl;
	return std::pair<Eigen::VectorXd, vector<int>>(weights, neighbors);
}
//normal vector is hard coded for square for now. fix.
void Grid::build_deriv_normal_bound() {
	double x, y;
	std::pair<Eigen::VectorXd, vector<int>> coeffs;
	int currPoint;
	deriv_normal_coeffs_ = vector<deriv_normal_bc>();
	deriv_normal_bc bound;
	for (size_t i = 0; i < boundaries_.size(); i++) {
		//dirichlet
		if ((boundaries_[i]).type == 2) {
			for (size_t j = 0; j < (boundaries_[i].bcPoints).size(); j++) {
				currPoint = boundaries_[i].bcPoints[j];
				x = std::get<0>(points_[currPoint]);
				y = std::get<1>(points_[currPoint]);
				if (x == 0 || x == 1) {
					coeffs = derivx_weights(currPoint);
				}
				else if (y == 0 || y == 1) {
					coeffs = derivy_weights(currPoint);
				}
				bound.pointID = currPoint;
				bound.value = boundaries_[i].values[j];
				bound.weights = coeffs.first;
				bound.neighbors = coeffs.second;
				deriv_normal_coeffs_.push_back(bound);
			}
		}
	}
	
}
//Builds Sparse Matrix that approximates laplacian operator for SOR from the interpolation weights.
void Grid::build_laplacian() {
	vector<Eigen::Triplet<double>> tripletList;
	for (int i = 0; i < laplaceMatSize_; i++) {
		std::pair <Eigen::VectorXd, vector<int>> weights = laplaceWeights(i);
		//Handle interior points first.
		if (bcFlags_[i] != 2) {
			for (size_t j = 0; j < weights.second.size(); j++) {
				tripletList.push_back(Eigen::Triplet<double>(i, weights.second[j], weights.first(j)));
			}
		}
		//Add column of ones if Neumann
		if (neumannFlag_ && bcFlags_[i] != 2) {
			tripletList.push_back(Eigen::Triplet<double>(i, laplaceMatSize_, 1));
		}
	}
	//Add row of ones if Neumann
	if (neumannFlag_) {
		for (int i = 0; i < laplaceMatSize_ + 1; i++) {
			if ((i == laplaceMatSize_) || bcFlags_[i] != 2) {
				tripletList.push_back(Eigen::Triplet<double>(laplaceMatSize_, i, 1));
			}
		}
	}
	//Add neumann boundary rows: d/dn = 0
	deriv_normal_bc bound;
	if (neumannFlag_) {
		for (size_t i = 0; i < deriv_normal_coeffs_.size(); i++) {
			bound = deriv_normal_coeffs_[i];
			for (size_t j = 0; j < bound.neighbors.size(); j++) {
				tripletList.push_back(Eigen::Triplet<double>(bound.pointID, bound.neighbors[j], bound.weights.coeff(j)));
			}
		}
	}
	laplaceMat_->setFromTriplets(tripletList.begin(), tripletList.end());
	laplaceMat_->makeCompressed();
}
std::pair<double, double> Grid::minMaxCoord(vector<int> pointIDs, char coord) {
	double min;
	double max;
	double tempCoord;
	if (coord == 'x') {
		min = std::get<0>(points_[pointIDs[0]]);
		max = min;
	}
	else if (coord == 'y') {
		min = std::get<1>(points_[pointIDs[0]]);
		max = min;
	}
	for (size_t i = 0; i < pointIDs.size(); i++) {
		if (coord == 'x') {
			tempCoord = std::get<0>(points_[pointIDs[i]]);
		}
		else if (coord == 'y') {
			tempCoord = std::get<1>(points_[pointIDs[i]]);
		}
		if (tempCoord > max) {
			max = tempCoord;
		}
		else if (tempCoord < min) {
			min = tempCoord;
		}
	}
	return std::pair<double, double>(min, max);
}
//SECOND TO LAST element in the returned vector is a Point with coords (scale, scale, scale). This is messy but is here to prevent very complex return types.
//LAST ELEMENT returned here is the evaluation point. First k elements are those in the cloud.
//Therefore the returned vector has k+2 elements.
std::vector<Point> Grid::shifting_scaling(vector<int> pointIDs, Point evalPoint) {
	std::pair<double, double> minMaxX, minMaxY;
	minMaxX = minMaxCoord(pointIDs, 'x');
	minMaxY = minMaxCoord(pointIDs, 'y');
	double x, y, x_ss, y_ss, scale, minX, minY, maxX, maxY;
	minX = minMaxX.first;
	maxX = minMaxX.second;
	minY = minMaxY.first;
	maxY = minMaxY.second;
	scale = std::max(maxX - minX, maxY - minY);
	std::vector<Point> scaledPoints;
	for (size_t i = 0; i < pointIDs.size(); i++) {
		x = std::get<0>(points_[pointIDs[i]]);
		y = std::get<1>(points_[pointIDs[i]]);
		x_ss = (x - minX) / scale;
		y_ss = (y - minY) / scale;
		scaledPoints.push_back(Point(x_ss, y_ss, 0));
	}
	scaledPoints.push_back(Point(scale, scale, scale));
	x = std::get<0>(evalPoint);
	y = std::get<1>(evalPoint);
	x_ss = (x - minX) / scale;
	y_ss = (y - minY) / scale;
	scaledPoints.push_back(Point(x_ss, y_ss, 0));
	return scaledPoints;
}
void Grid::diagonal_scaling(Eigen::SparseMatrix<double, Eigen::RowMajor>* matrix, Eigen::VectorXd* rhs) {
	double diagCoeff;
	double* laplaceValues = matrix->valuePtr();
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

	for (int i = 0; i < matrix->rows(); i++) {
		diagCoeff = 0;
		rowStartIdx = outerValues[outerIdx];
		rowEndIdx = outerValues[outerIdx + 1];
		//Search for diagonal on every row
		for (int j = rowStartIdx; j < rowEndIdx; j++) {
			if (innerValues[j] == i) {
				diagCoeff = laplaceValues[j];
				break;
			}
		}
		for (int j = rowStartIdx; j < rowEndIdx; j++) {
			laplaceValues[j] /= diagCoeff;
		}
		//Modify Source
		rhs->coeffRef(i) /= diagCoeff;
		outerIdx++;
	}
}
//TODO: Fix multigrid stuff for Neumann after single grid works.
std::pair<Eigen::VectorXd, vector<int>> Grid::pointInterpWeights(Point point) {
	std::tuple<Eigen::MatrixXd, vector<int>, vector<Point>> coeffs = buildCoeffMatrix(point, false, false);
	int polyTerms = (properties_.polyDeg + 1)*(properties_.polyDeg + 2) / 2;
	Eigen::VectorXd rhs = Eigen::VectorXd::Zero(properties_.stencilSize + polyTerms);
	Eigen::MatrixXd coeff_mat = std::get<0>(coeffs);
	vector<int> neighbors = std::get<1>(coeffs);
	vector<Point> scaledPoints = std::get<2>(coeffs);
	//Build RHS rbf terms	
	Point evalPoint = scaledPoints.at(scaledPoints.size() - 1);
	double xEval = std::get<0>(evalPoint);
	double yEval = std::get<1>(evalPoint);
	for (int i = 0; i < properties_.stencilSize; i++) {
		//cout << distance(evalPoint, scaledPoints[i]) << endl;
		rhs(i) = std::pow(distance(evalPoint, scaledPoints[i]), properties_.rbfExp);
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
	double sum = 0;
	for (int i = 0; i < getStencilSize(); i++) {
		sum += weights.coeff(i);
	}
	return std::pair<Eigen::VectorXd, vector<int>>(weights, neighbors);
}
void Grid::cuthill_mckee_ordering(vector<vector<int>> &adjacency, vector<int> &order) {
	
	vector<bool> visited(laplaceMatSize_, false);
	//stores original values
	std::queue<int> queue;
	vector<int> new_order;
	int startPoint = 0;
	visited[startPoint] = true;
	queue.push(startPoint);
	int currPoint, adjPoint;
	while (!queue.empty()) {
		currPoint = queue.front();
		new_order.push_back(currPoint);
		queue.pop();
		for (size_t i = 0; i < adjacency[currPoint].size(); i++) {
			adjPoint = adjacency[currPoint].at(i);
			if (visited[adjPoint] == false) {
				visited[adjPoint] = true;
				queue.push(adjPoint);
			}
		}
	}
	order = new_order;
}
void Grid::reverse_cuthill_mckee_ordering(vector<vector<int>> &adjacency, vector<int> &order) {
	cuthill_mckee_ordering(adjacency, order);
	std::reverse(order.begin(), order.end());
}
void Grid::rcm_order_points() {
	vector<vector<int>> adjacency;
	vector<int> order; // = orderFromTxt("rcm/finer_grid_rcm.txt", laplaceMatSize_);
	
	vector<int> neighbors;
	for (size_t i = 0; i < points_.size(); i++) {
		neighbors = kNearestNeighbors(points_[i], neumannFlag_, (bcFlags_[i] != 0));
		adjacency.push_back(neighbors);
	}
	reverse_cuthill_mckee_ordering(adjacency, order);
	
	vector<Point> newPoints = points_;
	Eigen::VectorXd newSource = source_;
	vector<int> newBCFlag = bcFlags_;
	vector<int> oldToRCMPtrs(points_.size());

	for (size_t i = 0; i < points_.size(); i++) {
		newPoints[i] = points_[order[i]];
		newBCFlag[i] = bcFlags_[order[i]];
		newSource(i) = source_.coeff(order[i]);
		oldToRCMPtrs[order[i]] = i;
	}
	points_ = newPoints;
	source_ = newSource;
	bcFlags_ = newBCFlag;

	vector<int> newBCPoints;
	vector<double> newBCVals;
	
	for (size_t i = 0; i < boundaries_.size(); i++) {
		newBCPoints.clear();
		newBCVals.clear();
		for(size_t j = 0; j < boundaries_[i].bcPoints.size(); j++){
			newBCPoints.push_back(oldToRCMPtrs[boundaries_[i].bcPoints.at(j)]);
			newBCVals.push_back(boundaries_[i].values.at(j));
		}
		boundaries_[i].bcPoints = newBCPoints;
		boundaries_[i].values = newBCVals;
	}
	
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
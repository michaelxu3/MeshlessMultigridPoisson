#include "grid.h"
#include <stdexcept>
#include <Eigen/Dense>
#include "math.h"
#include <algorithm>
Grid::Grid(vector<std::tuple<double, double, double>> points, vector<std::tuple<vector<int>, int, double>> boundaries, int rbfExp, int polyDeg) {
	points_ = points;
	boundaries_ = boundaries;
	rbfExp_ = rbfExp;
	polyDeg_ = polyDeg;
	laplaceMatSize_ = (int)(points.size());
	int numPoint = (int)(points.size());
	values_ = new Eigen::VectorXd(numPoint);
	laplaceMat_ = new Eigen::MatrixXd(numPoint, numPoint);
}
Grid::~Grid() {
	delete values_;
	delete laplaceMat_;
}

double distance(std::tuple<double, double, double> refPoint, std::tuple<double, double, double> queryPoint) {
	return std::sqrt(std::get<0>(refPoint)*std::get<0>(queryPoint) + std::get<1>(refPoint)*std::get<1>(queryPoint) + std::get<2>(refPoint)*std::get<2>(queryPoint));
}

void Grid::setBC(int bNum, std::string type, double boundValue) {
	if (!type.compare("dirichlet") && !type.compare("neumann")) {
		throw std::invalid_argument("Specified Boundary Condition on Boundary #" + std::to_string(bNum) + " is not dirichlet or neumann");
	}
	std::tuple<vector<int>, int, double> bound = boundaries_.at(bNum);
	std::get<1>(bound) = type.compare("dirichlet") ? 0 : 1;
	std::get<2>(bound) = boundValue;
	if (type.compare("dirichlet")) {
		for (size_t i = 0; i < std::get<0>(bound).size(); i++) {
			(*values_)(std::get<0>(bound).at(i)) = boundValue;
		}
	}
}
void Grid::build_laplacian() {

}

//Fix distance function if we want to extend code to 3D. Makes and sorts a distance array with int point IDs. Brute force, so O(n log n) time
vector<int> Grid::kNearestNeighbors(int pointID, int k) {
	vector<std::pair<double, int>> distances;
	std::tuple<double, double, double> refPoint = points_[pointID];
	for (size_t i = 0; i < points_.size(); i++) {
		std::tuple<double, double, double> queryPoint = points_[i];
		distances.push_back(std::pair<double, int>(distance(refPoint, queryPoint), i));
	}
	//std::sort sorts pairs by first (distance)
	std::sort(distances.begin(), distances.end());

	vector<int> nearestNeighbors;
	for (int i = 1; i <= k; i++) {
		nearestNeighbors.push_back(distances[i].second);
	}
	return nearestNeighbors;
}
//Builds PHS coefficients matrix for a point
//returns pair with first = coeff matrix and second = list of neighbors, to save time.
std::pair<Eigen::MatrixXd, vector<int>> Grid::buildCoeffMatrix(int pointID, int k) {
	int polyTerms = (polyDeg_ + 1)*(polyDeg_ + 2) / 2;
	vector<int> neighbors = kNearestNeighbors(pointID, k);
	//build A-matrix
	Eigen::MatrixXd coeff_mat = Eigen::MatrixXd(k + polyTerms, k + polyTerms);
	for (int i = 0; i < k; i++) {
		for (int j = 0; j < k; j++) {
			double r = distance(points_[neighbors[i]], points_[neighbors[j]]);
			coeff_mat(i, j) = std::pow(r, rbfExp_);
			/*
			if (rbfExp_ % 2 == 0) {
				coeff_mat(i, j) *= std::log(r);
			}*/
		}
	}
	//fill P-matrix sections.
	int colIndex;
	for (int row = 0; row < k; row++) {
		colIndex = k;
		for (int p = 0; p <= polyDeg_; p++) {
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

//Builds RHS and solves for the laplace weights for a point
std::pair<Eigen::VectorXd, vector<int>> Grid::laplaceWeights (int pointID, int k) {
	std::pair<Eigen::MatrixXd, vector<int>> coeffs = buildCoeffMatrix(pointID, k);
	vector<int> neighbors = coeffs.second;
	int polyTerms = (polyDeg_ + 1)*(polyDeg_ + 2) / 2;
	Eigen::VectorXd rhs = Eigen::VectorXd(k + polyTerms);
	//Build RHS
	double xEval = std::get<0>(points_[pointID]);
	double yEval = std::get<1>(points_[pointID]);
	for (int i = 0; i < k; i++) {
		double delX = std::get<0>(points_[neighbors[i]]) - xEval;
		double delY = std::get<1>(points_[neighbors[i]]) - yEval;
		rhs(i) = rbfExp_ * rbfExp_*std::pow(delX*delX + delY * delY, (double)(rbfExp_) / 2.0 - 1);
	}
	//poly terms
	int rowIndex = k;
	for (int p = 0; p <= polyDeg_; p++) {
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
	Eigen::VectorXd weights = coeffs.first.partialPivLu().solve(rhs);
	return std::pair<Eigen::VectorXd, vector<int>>(weights, neighbors);
}

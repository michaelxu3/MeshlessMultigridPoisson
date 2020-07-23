#pragma once
#include "testing_functions.hpp"
int main() {
	//run_tests();
	testGmshSingleGrid();
	/*
	vector<Eigen::Triplet<double>> triples;
	Eigen::SparseMatrix<double, 1> mat (4,4);
	triples.push_back(Eigen::Triplet<double>(2, 2, 3));
	triples.push_back(Eigen::Triplet<double>(1, 1, 0.234));
	triples.push_back(Eigen::Triplet<double>(1, 1, -0.234));
	triples.push_back(Eigen::Triplet<double>(2, 2, 9));
	triples.push_back(Eigen::Triplet<double>(2, 2, -4));

	mat.setFromTriplets(triples.begin(), triples.end());

	cout << mat.toDense() << endl;
	*/
	return 0;
}
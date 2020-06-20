#include "multigrid.h"
#include <iostream>
Multigrid::Multigrid() {

	grids_ = vector<std::pair<int, Grid*>>(); //int = grid size
}
Multigrid::~Multigrid() {
	for (size_t i = 0; i < grids_.size(); i++) {
		delete grids_.at(i).second;
		delete prolongMatrices_.at(i);
		delete restrictionMatrices_.at(i);
	}
}
Eigen::SparseMatrix<double>* Multigrid::buildInterpMatrix(Grid* baseGrid, Grid* targetGrid) {
	//matrix is n x m where m is size of target grid, n is size of base grid
	Eigen::SparseMatrix<double>* interpMatrix = new Eigen::SparseMatrix<double>(targetGrid->getSize(), baseGrid->getSize());
	vector<Eigen::Triplet<double>> tripletList;

	for (int i = 0; i < targetGrid->getSize(); i++) {
		std::pair<Eigen::VectorXd, vector<int>> pointWeights = baseGrid->pointInterpWeights(targetGrid->getPoint(i));
		for (size_t j = 0; j < pointWeights.second.size(); j++) {

			tripletList.push_back(Eigen::Triplet<double>(i, pointWeights.second[j], pointWeights.first(j)));
		}
	}
	interpMatrix->setFromTriplets(tripletList.begin(), tripletList.end());
	return interpMatrix;
}

void Multigrid::buildProlongMatrices() {
	prolongMatrices_.resize(grids_.size());
	for (size_t i = 0; i < grids_.size() - 1; i++) {
		prolongMatrices_[i] = (buildInterpMatrix(grids_[i].second, grids_[i+1].second));
	}
	prolongMatrices_[prolongMatrices_.size() - 1] = NULL;
}
void Multigrid::buildRestrictionMatrices() {
	restrictionMatrices_.resize(grids_.size());
	restrictionMatrices_[0] = NULL;
	for (size_t i = 1; i < grids_.size(); i++) {
		restrictionMatrices_[i] = buildInterpMatrix(grids_[i].second, grids_[i - 1].second);
	}
}
void Multigrid::buildMatrices() {
	buildProlongMatrices();
	buildRestrictionMatrices();
}
void Multigrid::vCycle() {
	//Restriction
	Eigen::VectorXd  u_old = *(grids_[grids_.size() - 1].second->values_);
	for (size_t i = grids_.size() - 1; i > 0;  i--) {
		grids_[i].second->sor();
		*(grids_[i-1].second->values_) = (*(restrictionMatrices_[i])) * (*(grids_[i].second->values_));
	}
	//prolongation
	for (size_t i = 0; i < grids_.size() - 1; i++) {
		grids_[i].second->sor();
		*(grids_[i+1].second->values_) = (*(prolongMatrices_[i])) * (*grids_[i].second->values_);
	}
	double residual = (*(grids_[grids_.size() - 1].second->values_) - u_old).norm()/grids_[grids_.size() - 1].first;
	std::cout << "Residual: " << residual << std::endl;

}
void Multigrid::addGrid(Grid* grid) {
	grids_.push_back(std::pair<int, Grid*>(grid->getSize(), grid));
	sortGridsBySize();
}
void Multigrid::sortGridsBySize() {
	std::sort(grids_.begin(), grids_.end());
}
#include "multigrid.h"
#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;
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

	Eigen::SparseMatrix<double>* interpMatrix = new Eigen::SparseMatrix<double>(targetGrid->getSize(), baseGrid->getSize());
	vector<Eigen::Triplet<double>> tripletList;
	//tripletList: <row, col, value>
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
	//Eigen::VectorXd  u_old = *(grids_[grids_.size() - 1].second->values_);
	
	/*
	Grid* currGrid;
	currGrid = grids_[grids_.size() - 1].second;
	double resid_norm = currGrid->residual().lpNorm<1>() / currGrid->source_.lpNorm<1>();
	std::cout << "Residual: " << resid_norm << std::endl;
	Grid* coarsegrid = grids_[0].second;
	Grid* finegrid = grids_[1].second;
	finegrid->boundaryOp("fine");
	finegrid->sor(finegrid->laplaceMat_, finegrid->values_, &finegrid->source_);
	coarsegrid->source_ = *restrictionMatrices_[1] * finegrid->residual();
	//std::cout << finegrid->residual().norm() << std::endl;
	cout << grids_[0].second->source_.norm() << endl;
	coarsegrid-> boundaryOp("coarse");
	coarsegrid->values_->setZero();
	coarsegrid->sor(coarsegrid->laplaceMat_, coarsegrid->values_, &coarsegrid->source_);
	*(finegrid->values_) += (*prolongMatrices_[0]) * (*coarsegrid->values_);
	*/
	
	Grid* currGrid;
	currGrid = grids_[grids_.size() - 1].second;
	double resid_norm = residual();
	std::cout << std::setprecision(12) << "Residual: " << resid_norm << std::endl;
	//Restriction
	for (size_t i = grids_.size() - 1; i > 0; i--) {
		currGrid = grids_[i].second;
		std::string gridType = i == grids_.size() - 1 ? "fine" : "coarse";
		//cout << gridType << endl;
		if (i != grids_.size() - 1) {
			currGrid->values_->setZero();
		}
		currGrid->boundaryOp(gridType);
		currGrid->sor(currGrid->laplaceMat_, currGrid->values_, &(currGrid->source_));
		grids_[i - 1].second->source_ = (*(restrictionMatrices_[i])) * currGrid->residual();
	}
	//cout << grids_[0].second->source_.norm() << endl;
	//Iterate on coarsest grid
	currGrid->boundaryOp("coarse");
	currGrid = grids_[0].second;
	currGrid->values_->setZero();
	currGrid->sor(currGrid->laplaceMat_, currGrid->values_, &(currGrid->source_));
	currGrid->sor(currGrid->laplaceMat_, currGrid->values_, &(currGrid->source_));

	//correction and prolongation
	
	for (size_t i = 1; i < grids_.size(); i++) {
		currGrid = grids_[i].second;
		//correction
		*currGrid->values_ += (*prolongMatrices_[i-1]) * (*grids_[i-1].second->values_);
		//smoother
		currGrid->sor(currGrid->laplaceMat_, currGrid->values_, &(currGrid->source_));
	}
	
	
}

double Multigrid::residual() {
	Grid* finegrid = grids_[grids_.size() - 1].second;
	return finegrid->residual().lpNorm<1>() / finegrid->source_.lpNorm<1>();
}
void Multigrid::addGrid(Grid* grid) {
	grids_.push_back(std::pair<int, Grid*>(grid->getSize(), grid));
	sortGridsBySize();
}
void Multigrid::sortGridsBySize() {
	std::sort(grids_.begin(), grids_.end());
}
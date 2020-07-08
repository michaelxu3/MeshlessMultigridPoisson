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
		std::pair<Eigen::VectorXd, vector<int>> pointWeights = baseGrid->pointInterpWeights(targetGrid->points_[i]);
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
	/*
	Grid* currGrid;
	for (size_t i = 0; i < grids_.size() - 1; i++) {
		currGrid = grids_[i].second;
		currGrid->diagonal_scaling(currGrid->laplaceMat_, &currGrid->source_);
	}
	*/
}
void Multigrid::vCycle() {
	//Restriction
	/*
	Grid* currGrid;
	currGrid = grids_[grids_.size() - 1].second;
	double resid_norm = currGrid->residual().lpNorm<1>() / currGrid->source_.lpNorm<1>();
	std::cout << "Residual: " << resid_norm << std::endl;
	Grid* coarsegrid = grids_[0].second;
	Grid* finegrid = grids_[1].second;
	finegrid->boundaryOp("fine");
	finegrid->sor(finegrid->laplaceMat_, finegrid->values_, &finegrid->source_);
	//cout << "fine grid raw resid norm: " << finegrid->residual().lpNorm<1>() << endl;
	coarsegrid->source_ = *restrictionMatrices_[1] * finegrid->residual();
	//cout << "Restricted residual norm: " << (*restrictionMatrices_[1] * finegrid->residual()).lpNorm<1>() << endl;
	coarsegrid->boundaryOp("coarse");
	coarsegrid->values_->setZero();
	coarsegrid->sor(coarsegrid->laplaceMat_, coarsegrid->values_, &coarsegrid->source_);
	//cout << "Coarse grid residual: " << coarsegrid->residual().lpNorm<1>()/coarsegrid->source_.lpNorm<1>() << endl;
	*(finegrid->values_) += (*prolongMatrices_[0]) * (*coarsegrid->values_);
	finegrid->sor(finegrid->laplaceMat_, finegrid->values_, &finegrid->source_);
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
		grids_[i - 1].second->fix_vector_bound_coarse(&grids_[i-1].second->source_);

	}
	//cout << grids_[0].second->source_.norm() << endl;
	//Iterate on coarsest grid
	currGrid->boundaryOp("coarse");
	currGrid = grids_[0].second;
	currGrid->values_->setZero();
	currGrid->sor(currGrid->laplaceMat_, currGrid->values_, &(currGrid->source_));

	//correction and prolongation
	Eigen::VectorXd correction;
	for (size_t i = 1; i < grids_.size(); i++) {
		currGrid = grids_[i].second;
		//correction
		correction = (*prolongMatrices_[i - 1]) * (*grids_[i - 1].second->values_);
		currGrid->fix_vector_bound_coarse(&correction);
		*currGrid->values_ += correction;
		//smoother
		currGrid->sor(currGrid->laplaceMat_, currGrid->values_, &(currGrid->source_));
	}
	
	/*
	//code just to test prolong error
	Grid* coarsegrid = grids_[0].second;
	Grid* finegrid = grids_[1].second;
	coarsegrid->sor(coarsegrid->laplaceMat_, coarsegrid->values_, &coarsegrid->source_);
	*(finegrid->values_) = (*prolongMatrices_[0]) * (*coarsegrid->values_);
	finegrid->fix_vector_bound_coarse(finegrid->values_);
	cout << "Coarse grid residual: " << coarsegrid->residual().norm() / coarsegrid->source_.norm() << endl;
	cout << "Fine grid residual: " << finegrid->residual().norm() / finegrid->source_.norm() << endl;
	cout << "Ratio: " << (finegrid->values_->lpNorm<1>()) / (coarsegrid->values_->lpNorm<1>()) << endl;
	*/
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
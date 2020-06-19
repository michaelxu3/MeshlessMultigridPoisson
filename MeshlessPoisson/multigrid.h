#ifndef MULTIGRID_H
#define MULTIGRID_H
#include "grid.h"
class Multigrid {
private:
	vector<std::pair<int, Grid*>> grids_;
	vector<int> sorGridIters_;
	vector<Eigen::SparseMatrix<double>*> restrictionMatrices_;
	vector<Eigen::SparseMatrix<double>*> prolongMatrices_;
	void sortGridsBySize();
	Eigen::SparseMatrix<double>* buildInterpMatrix(Grid* baseGrid, Grid* targetGrid);
	void buildRestrictionMatrices();
	void buildProlongMatrices();
	void vCycle();

public:
	Multigrid();
	~Multigrid();
	void addGrid(Grid* grid);
};
#endif
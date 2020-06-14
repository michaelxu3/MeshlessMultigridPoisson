#ifndef MULTIGRID_H
#define MULTIGRID_H
#include "grid.h"
class Multigrid {
private:
	vector<Grid*> grids_;

public:
	Multigrid(vector<Grid*> grids);
	void interp_grid(Grid* baseGrid, Grid* targetGrid);
};
#endif
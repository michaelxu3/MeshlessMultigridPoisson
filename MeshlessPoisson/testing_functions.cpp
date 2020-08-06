#include "testing_functions.hpp"

double calc_l1_error(Grid* grid, bool neumannFlag, int k1, int k2) {
	vector<Point> points = grid->points_;
	Eigen::VectorXd actual(points.size());
	double x, y;
	for (size_t i = 0; i < points.size(); i++) {
		x = std::get<0>(points[i]);
		y = std::get<1>(points[i]);
		actual(i) = neumannFlag? std::cos(k1*pi*x)*std::cos(k2*pi*y) : std::sin(k1*pi*x)*std::sin(k2*pi*y);
	}
	double solMean = 0;
	double manufacturedMean = 0;
	if (!neumannFlag) {
		return (*grid->values_ - actual).lpNorm<1>() / grid->laplaceMatSize_;
	}

	for (int i = 0; i < grid->laplaceMatSize_; i++) {
		solMean += grid->values_->coeff(i) / grid->laplaceMatSize_;
		manufacturedMean += actual.coeff(i) / grid->laplaceMatSize_;
	}

	for (int i = 0; i < grid->laplaceMatSize_; i++) {
		grid->values_->coeffRef(i) += (manufacturedMean - solMean);
	}

	double error = 0;
	for (int i = 0; i < points.size(); i++) {
		error += std::abs(grid->values_->coeff(i) - actual.coeff(i));
	}
	error /= points.size();
	return error;
}
Grid* genGmshGridDirichlet(string geomtype, const char* filename, GridProperties props, std::string filetype, int k1, int k2) {
	vector<std::tuple<double, double, double>> points;
	if (filetype.compare("txt") == 0) {
		points = pointsFromTxts(filename);
	}
	else {
		points = pointsFromMshFile(filename);
	}
	double x, y;
	vector<int> bPts, bPts_inner;
	vector<double> bValues, bValues_inner;

	Eigen::VectorXd source(points.size());
	if (geomtype.compare("square") == 0) {
		for (size_t i = 0; i < points.size(); i++) {
			x = std::get<0>(points[i]);
			y = std::get<1>(points[i]);
			source(i) = -(k1*k1 + k2 * k2) * pi*pi*std::sin(k1*pi*x)*std::sin(k2*pi*y);
			if (x == 0 || x == 1 || y == 0 || y == 1) {
				bPts.push_back(i);
				bValues.push_back(0.0);
			}
		}
	}
	else if (geomtype.compare("square_with_circle") == 0) {
		for (size_t i = 0; i < points.size(); i++) {
			x = std::get<0>(points[i]);
			y = std::get<1>(points[i]);
			source(i) = -(k1*k1 + k2 * k2) * pi*pi*std::sin(k1*pi*x)*std::sin(k1*pi*y);
			if (x == 0 || x == 1 || y == 0 || y == 1) {
				bPts.push_back(i);
				bValues.push_back(0);
			}
			else if (std::abs(0.0625 - (x - 0.5)*(x - 0.5) - (y - 0.5)*(y - 0.5)) <= std::pow(10,-10)) {
				//cout << i << endl;
				bPts_inner.push_back(i);
				bValues_inner.push_back(std::sin(k1*pi*x)*std::sin(k1*pi*y));
			}
		}
	}
	else if (geomtype.compare("concentric_circles") == 0) {
		for (size_t i = 0; i < points.size(); i++) {
			x = std::get<0>(points[i]);
			y = std::get<1>(points[i]);
			source(i) = -(k1*k1 + k2 * k2) * pi*pi*std::sin(k1*pi*x)*std::sin(k1*pi*y);
			if (std::abs(0.25 - (x - 0.5)*(x - 0.5) - (y - 0.5)*(y - 0.5)) <= std::pow(10, -10)) {
				bPts.push_back(i);
				bValues.push_back(std::sin(k1*pi*x)*std::sin(k1*pi*y));

			}
			else if (std::abs(0.0625 - (x - 0.5)*(x - 0.5) - (y - 0.5)*(y - 0.5)) <= std::pow(10, -10)) {
				bPts_inner.push_back(i);
				bValues_inner.push_back(std::sin(k1*pi*x)*std::sin(k1*pi*y));

			}
		}
	}
	Boundary boundary;
	boundary.bcPoints = bPts;
	boundary.type = 1;
	boundary.values = bValues;
	vector<Boundary> bcs;
	bcs.push_back(boundary);
	if (geomtype.compare("square") != 0) {
		Boundary inner;
		inner.bcPoints = bPts_inner;
		inner.type = 1;
		inner.values = bValues_inner;
		bcs.push_back(inner);
	}

	Grid* grid = new Grid(points, bcs, props, source);
	grid->implicitFlag_ = false;
	grid->setBCFlag(0, std::string("dirichlet"), bValues);
	if (geomtype.compare("square") != 0) {
		//cout << grid->boundaries_[1].bcPoints.size() << endl;
		grid->setBCFlag(1, std::string("dirichlet"), bValues_inner);
	}
	//grid->print_bc_values();
	grid->rcm_order_points();
	grid->build_laplacian();
	return grid;
}
Grid* genGmshGridNeumann(string geomtype, const char* filename, GridProperties props, std::string filetype, int k1, int k2) {
	vector<std::tuple<double, double, double>> points;
	if (filetype.compare("txt") == 0) {
		points = pointsFromTxts(filename);
	}
	else {
		points = pointsFromMshFile(filename);
	}
	double x, y;
	vector<int> bPts, bPts_inner;
	vector<double> bValues, bValues_inner;
	Eigen::VectorXd actual(points.size() + 1);
	Eigen::VectorXd source(points.size() + 1);
	if (geomtype.compare("square") == 0) {
		for (size_t i = 0; i < points.size(); i++) {
			x = std::get<0>(points[i]);
			y = std::get<1>(points[i]);
			actual(i) = std::cos(k1*pi*x)*std::cos(k2*pi*y);
			source(i) = -(k1*k1 + k2 * k2)*pi*pi*std::cos(k1*pi*x)*std::cos(k2*pi*y);
			if (x == 0 || x == 1 || y == 0 || y == 1) {
				bPts.push_back(i);
				bValues.push_back(0.0);
			}
		}
	}
	else if (geomtype.compare("square_with_circle") == 0) {
		double normX, normY, norm;
		for (size_t i = 0; i < points.size(); i++) {
			x = std::get<0>(points[i]);
			y = std::get<1>(points[i]);
			source(i) = -(k1*k1 + k2 * k2)*pi*pi*std::cos(k1*pi*x)*std::cos(k2*pi*y);
			if (x == 0 || x == 1 || y == 0 || y == 1) {
				bPts.push_back(i);
				bValues.push_back(0);
			}
			else if (std::abs(0.0625 - (x - 0.5)*(x - 0.5) - (y - 0.5)*(y - 0.5)) <= std::pow(10, -10)) {
				normX = x - 0.5;
				normY = y - 0.5;
				norm = std::sqrt(normX*normX + normY * normY);
				normX /= norm;
				normY /= norm;
				bPts_inner.push_back(i);
				bValues_inner.push_back(-normX*k1*std::sin(k1*pi*x)*std::cos(k2*pi*y)
										-normY*k2*std::cos(k1*pi*x)*std::sin(k2*pi*y));
			}
		}
	}
	actual(actual.rows() - 1) = 0;
	source(source.rows() - 1) = 0;
	Boundary boundary;
	boundary.bcPoints = bPts;
	boundary.type = 2;
	boundary.values = bValues;
	vector<Boundary> bcs;
	bcs.push_back(boundary);
	if (geomtype.compare("square") != 0) {
		Boundary inner;
		inner.bcPoints = bPts_inner;
		inner.type = 2;
		inner.values = bValues_inner;
		bcs.push_back(inner);
	}
	Grid* grid = new Grid(points, bcs, props, source);
	grid->implicitFlag_ = true;
	grid->setBCFlag(0, std::string("neumann"), bValues);
	if (geomtype.compare("square") != 0) {
		grid->setBCFlag(1, std::string("neumann"), bValues_inner);

	}
	//grid->print_bc_values();
	grid->modify_coeff_neumann("fine");
	grid->build_normal_vecs(filename, geomtype);
	grid->rcm_order_points();
	grid->build_deriv_normal_bound();
	grid->build_laplacian();
	//grid->push_inhomog_to_rhs();
	return grid;
}
void write_temp_contour(Grid* testGrid, string directory, string extension) {
	//Write necessary things for temperature contour
	vector<Point> points = testGrid->points_;
	double x, y;
	vector<double> xv, yv;
	for (size_t i = 0; i < points.size(); i++) {
		x = std::get<0>(points[i]);
		y = std::get<1>(points[i]);
		xv.push_back(x);
		yv.push_back(y);
	}
	string x_str = "x_";
	string y_str = "y_";
	string temp_str = "temp_";
	string cond_str = "cond_rbf_";
	string txt = ".txt";

	writeVectorToTxt(xv, (directory + x_str + extension + txt).c_str());
	writeVectorToTxt(yv, (directory + y_str + extension + txt).c_str());

	vector<double> temp;
	for (int i = 0; i < testGrid->values_->rows() - 1; i++) {
		temp.push_back(testGrid->values_->coeff(i));
	}
	writeVectorToTxt(temp, (directory + temp_str + extension + txt).c_str());
}
void write_mg_resid(Multigrid& mg, string directory, string extension) {
	string res_str = "resid_";
	string txt = ".txt";
	writeVectorToTxt(mg.residuals_, (directory + res_str + extension + txt).c_str()); 
}
void write_l1error_cond(Multigrid& mg, MultigridParameters params, string directory, string extension) {
	string res_str = "cond_error_";
	string txt = ".txt";
	vector<double> vec;
	Grid* grid = mg.grids_.at(mg.grids_.size() - 1).second;
	vec.push_back(calc_l1_error(grid, grid->neumannFlag_, params.k1, params.k2));
	vec.push_back(grid->cond_L());
	writeVectorToTxt(mg.residuals_, (directory + res_str + extension + txt).c_str());
}
void run_mg_sim(MultigridParameters params){
	Multigrid mg;
	for (int i = 0; i < (int)(params.filenames.size()); i++) {
		if (params.neumann == true) {
			mg.addGrid(genGmshGridNeumann(params.geomtype, (params.directory + params.filenames[i]).c_str(), params.props[i], params.filetypes[i], params.k1, params.k2));
		}
		else {
			mg.addGrid(genGmshGridDirichlet(params.geomtype, (params.directory + params.filenames[i]).c_str(), params.props[i], params.filetypes[i], params.k1, params.k2));
		}
	}
	mg.buildMatrices();
	for (int i = 0; i < params.num_v_cycle; i++) {
		mg.vCycle();
	}
	write_mg_resid(mg, params.directory, params.extension);
	write_temp_contour(mg.grids_.back().second, params.directory, params.extension);

	//write cond number and l1 error

}
MultigridParameters gen_mg_param(string geom, int numGrids, int k, int poly_deg, int vcyc, bool neumann) {
	string dir;
	vector<string> msh_files;
	if (geom.compare("square") == 0) {
		msh_files = { "square_45.msh", "square_170.msh", "square_600.msh", "square_2.5k.msh", "square_10k.msh" };
		dir = "square_test_geometries/";
	}
	else if (geom.compare("square_with_circle") == 0) {
		msh_files = {"square_hole_44.msh",  "square_hole_176.msh", "square_hole_640.msh", "square_hole_2532.msh", "square_hole_10197.msh" };
		dir = "square_hole_geoms/";
	}
	else if (geom.compare("concentric_circles") == 0) {
		msh_files = { "concentric_circles_59.msh",  "concentric_circles_188.msh" , "concentric_circles_650.msh" , "concentric_circles_2581.msh" , "concentric_circles_10207.msh" };
		dir = "concentric_circle_geoms/";
	}
	vector<string> filenames, filetypes;
	for (int i = 0; i < numGrids; i++) {
		filenames.push_back(msh_files.at(i));
		filetypes.push_back("msh");
		//filetypes.push_back("txt");
	}
	vector<GridProperties> props(numGrids);
	for (int i = 0; i < numGrids; i++) {
		props[i].iters = 5;
		props[i].polyDeg = (i == numGrids - 1) ? poly_deg : 3;
		props[i].omega = 1;
		props[i].rbfExp = 3;
		props[i].stencilSize = (int)(1.5 * (props[i].polyDeg + 1) * (props[i].polyDeg + 2) / 2);
	}
	string bctype = (neumann) ? "neumann" : "dirichlet";
	MultigridParameters params;
	params.geomtype = geom;
	params.directory = dir;
	params.extension = bctype + "_L=" + std::to_string(poly_deg) + "_K=" + std::to_string(k);
	params.filenames = filenames;
	params.filetypes = filetypes;
	params.k1 = k;
	params.k2 = k;
	params.neumann = neumann;
	params.num_v_cycle = vcyc;
	params.props = props;
	return params;
}
void run_tests() {
	MultigridParameters param = gen_mg_param("square_with_circle", 5, 1, 5, 150, true);
	run_mg_sim(param);
}

void testGmshSingleGrid() {
	GridProperties props;
	
	props.iters = 5;
	props.polyDeg = 5;
	props.omega = 1.4;
	props.rbfExp = 3;
	props.stencilSize = (int)(1.5 * (props.polyDeg + 1) * (props.polyDeg + 2) / 2);
	
	//Grid* testGrid = genGmshGridNeumann("square_with_circle", "square_hole_geoms/square_hole_10197.msh", props, "msh", 1, 1);
	Grid* testGrid = genGmshGridNeumann("square", "square_test_geometries/square_10k.msh", props, "msh", 1, 1);
	//Grid* testGrid = genGmshGridDirichlet("concentric_circles", "concentric_circle_geoms/concentric_circles_10207.msh", props, "msh", 1, 1);
	
	vector<double> res;
	testGrid->boundaryOp("fine");
	for (int i = 0; i < 4000; i++) {
		cout << "residual: " << testGrid->residual().lpNorm<1>() / testGrid->source_.lpNorm<1>() << endl;
		res.push_back(testGrid->residual().lpNorm<1>() / testGrid->source_.lpNorm<1>());
		testGrid->sor(testGrid->laplaceMat_, testGrid->values_, &testGrid->source_);
	}
	//testGrid->bound_eval_neumann();
	writeVectorToTxt(res, "residual.txt");
	cout << "L1 error: " << calc_l1_error(testGrid, testGrid->neumannFlag_, 1, 1) << endl;
	/*
	//band matrix plotting.
	
	Eigen::SparseMatrix<double, 1> * matrix = testGrid->laplaceMat_;
	const double* laplaceValues = matrix->valuePtr();
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

	vector<double> iv, jv;
	for (int i = 0; i < matrix->rows(); i++) {
		//diagCoeff = laplaceMat_->coeff(i, i);
		rowStartIdx = outerValues[outerIdx];
		rowEndIdx = outerValues[outerIdx + 1];
		//sum coeffs*x_i_old on the row i
		for (int j = rowStartIdx; j < rowEndIdx; j++) {
			iv.push_back(i);
			jv.push_back(innerValues[j]);
		}
		outerIdx++;
	}
	writeVectorToTxt(iv, "i.txt");
	writeVectorToTxt(jv, "j.txt");
	*/
	cout << testGrid->values_->maxCoeff() << endl;
	cout << testGrid->values_->minCoeff() << endl;
}
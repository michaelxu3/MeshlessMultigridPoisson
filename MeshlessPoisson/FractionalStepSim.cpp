#include "FractionalStepSim.hpp"
#include "testing_functions.hpp"
FractionalStepGrid* genFractionalStepGrid(const char* filename, GridProperties props, double dt, double mu, double rho, double ppe_conv, std::string coarse) {
	vector<std::tuple<double, double, double>> points;
	points = pointsFromMshFile(filename);
	double x, y;
	vector<int> bPts, bPts_inner;
	vector<double> bValues, bValues_inner;
	Eigen::VectorXd actual(points.size() + 1);
	Eigen::VectorXd source(points.size() + 1);
	double re = rho / mu;
	double lambda = 0.5 * re - std::sqrt(0.25 * re * re + 4 * EIGEN_PI * EIGEN_PI);
	for (size_t i = 0; i < points.size(); i++) {
		x = std::get<0>(points[i]);
		y = std::get<1>(points[i]);
		if (x == 0 || x == 1 || y == 0 || y == 1) {
			bPts.push_back(i);
			bValues.push_back(0.5 * std::exp(2 * lambda * x));//(0.0);
		}
	}
	source.setZero();
	Boundary boundary;
	boundary.bcPoints = bPts;
	boundary.type = 2;//1;
	boundary.values = bValues;
	vector<Boundary> bcs;
	bcs.push_back(boundary);

	FractionalStepGrid* grid = new FractionalStepGrid(points, bcs, props, source);
	cout << grid->source_.rows() << endl;

	grid->mu = mu;
	grid->rho = rho;
	grid->ppe_conv_res = ppe_conv;
	grid->dt = dt;
	grid->implicitFlag_ = true;//false
	grid->flowType = "kovasznay";
	grid->setBCFlag(0, std::string("neumann"/*"dirichlet"*/), bValues);
	grid->build_normal_vecs(filename, "square");
	grid->rcm_order_points();
	grid->build_deriv_normal_bound();
	grid->build_laplacian();
	grid->modify_coeff_neumann(coarse);
	grid->build_derivX_mat();
	grid->build_derivY_mat(); 
	grid->build_uv_laplace_mat();
	grid->push_inhomog_to_rhs();
	return grid;
}
FractionalStepParams gen_fracstep_param(int numGrids, int poly_deg, double dt, double mu, double rho, double ppe_conv) {
	string dir;
	vector<string> msh_files;
	msh_files = {/*"square_98.msh",*/ "square_170.msh", "square_600.msh", "square_2.5k.msh", "square_10k.msh" };
	dir = "square_test_geometries/";
	vector<string> filenames, filetypes;
	for (int i = 0; i < numGrids; i++) {
		filenames.push_back(msh_files.at(i));
	}
	vector<GridProperties> props(numGrids);
	for (int i = 0; i < numGrids; i++) {
		props[i].iters = 5;
		props[i].polyDeg = (i == numGrids - 1) ? poly_deg : 3;
		props[i].omega = 1.4;
		props[i].rbfExp = 3;
		props[i].stencilSize = (i == 0) ? (int)(2.5 * (props[i].polyDeg + 1) * (props[i].polyDeg + 2) / 2)
			: (int)(2.5 * (props[i].polyDeg + 1) * (props[i].polyDeg + 2) / 2);
	}
	FractionalStepParams params;
	params.directory = dir;
	params.extension = std::to_string(numGrids) + "grid_" + "_L=" +
		std::to_string(poly_deg) + "_fracstep";
	params.filenames = filenames;
	params.dt = dt;
	params.mu = mu;
	params.ppe_conv_res = ppe_conv;
	params.rho = rho;
	params.props = props;
	return params;
}
void check_derivs(FractionalStepGrid* grid) {
	grid->prescribe_soln();
	Eigen::VectorXd dudx = *grid->derivXMat_ * *grid->u;
	Eigen::VectorXd dudy = *grid->derivYMat_ * *grid->u;
	Eigen::VectorXd del2u = *grid->uvLaplaceMat_ * *grid->u;

	Eigen::VectorXd exact_dudx(grid->laplaceMatSize_);
	Eigen::VectorXd exact_dudy(grid->laplaceMatSize_);
	Eigen::VectorXd exact_del2u(grid->laplaceMatSize_);
	Eigen::VectorXd exact_del2p(grid->laplaceMatSize_);
	double x, y;
	double re = grid->rho / grid->mu;
	double lambda = 0.5 * re - std::sqrt(0.25 * re * re + 4 * EIGEN_PI * EIGEN_PI);
	for (int i = 0; i < grid->laplaceMatSize_; i++) {
		x = std::get<0>(grid->points_[i]);
		y = std::get<1>(grid->points_[i]);
		exact_dudx.coeffRef(i) = -lambda * std::exp(lambda * x) * std::cos(2 * EIGEN_PI * y);
		exact_dudy.coeffRef(i) = std::exp(lambda * x) * 2 * EIGEN_PI * std::sin(2 * EIGEN_PI * y);
		exact_del2u.coeffRef(i) = std::cos(2*EIGEN_PI*y)*std::exp(lambda*x)*(4*EIGEN_PI*EIGEN_PI - lambda*lambda);
		exact_del2p.coeffRef(i) = 2*lambda*lambda*std::exp(2*lambda*x);
	}
	cout << "dudx error: " << (exact_dudx - dudx).lpNorm<1>() / grid->laplaceMatSize_ << endl;
	cout << "dudy error: " << (exact_dudy - dudy).lpNorm<1>() / grid->laplaceMatSize_ << endl;
	cout << "laplacian error: " << (exact_del2u - del2u).lpNorm<1>() / grid->laplaceMatSize_ << endl;
	grid->calc_u_hat();
	Eigen::VectorXd actualUHat = *grid->u + grid->dt / grid->rho * *grid->derivXMat_ * grid->values_->head(grid->laplaceMatSize_);
	grid->calc_v_hat();
	Eigen::VectorXd actualVHat = *grid->v + grid->dt / grid->rho * *grid->derivYMat_ * grid->values_->head(grid->laplaceMatSize_);
	cout << "uhat error: " << (actualUHat - *grid->u_hat).lpNorm<1>() / grid->laplaceMatSize_ << endl;
	cout << "vhat error: " << (actualVHat - *grid->v_hat).lpNorm<1>() / grid->laplaceMatSize_ << endl;
	cout << "continuity: " << (dudx + *grid->derivYMat_ * *grid->v).lpNorm<1>() / grid->laplaceMatSize_ << endl;
	cout << "correction diff: " << (*grid->u - (*grid->u_hat - grid->dt / grid->rho * (*grid->derivXMat_ * grid->values_->head(grid->laplaceMatSize_)))).lpNorm<1>()/grid->laplaceMatSize_ << endl;
	//u cout << "pressure laplace error: " << ((grid->laplaceMat_->block(0, 0, grid->laplaceMatSize_, grid->laplaceMatSize_) * grid->values_->head(grid->laplaceMatSize_)) - exact_del2p).lpNorm<1>() / grid->laplaceMatSize_ << endl;
}
void run_fracstep_param(FractionalStepParams params, double endtime) {
	FractionalStepMultigrid mg;
	for (int i = 0; i < (int)(params.filenames.size()); i++) {
		std::string coarse = (i == params.filenames.size() - 1) ? "fine" : "coarse";
		mg.addGrid(genFractionalStepGrid((params.directory + params.filenames[i]).c_str(), params.props[i], params.dt, params.mu, params.rho, params.ppe_conv_res, coarse));
	}
	mg.buildMatrices();
	std::clock_t start = std::clock();
	double time = 0;
	FractionalStepGrid* finestGrid = mg.grids_[mg.grids_.size() - 1].second;
	//finestGrid->prescribe_soln();
	double resid, oldresid;
	resid = 100;
	oldresid = 1000;
	//check_derivs(finestGrid);
	int timesteps = 0;
	while (true) {
		*(finestGrid->u_old) = *(finestGrid->u);
		*(finestGrid->v_old) = *(finestGrid->v);
		finestGrid->set_uv_bound();
		finestGrid->calc_u_hat();
		finestGrid->calc_v_hat();
		finestGrid->set_ppe_source();
		finestGrid->push_inhomog_to_rhs();
		
		while (mg.residual() >= params.ppe_conv_res) {
			mg.vCycle();
			finestGrid->bound_eval_neumann();
		}
		
		finestGrid->correct_u();
		finestGrid->correct_v();
		finestGrid->set_uv_bound();
		resid = finestGrid->fs_residual();
		cout << "residual: " << std::abs(resid - oldresid) << endl;
		
		time += finestGrid->dt;
		if (timesteps > 2000) {
			break;
		}
		timesteps++;
		oldresid = resid;
	}
	cout << time << endl;
	Eigen::VectorXd actual(finestGrid->laplaceMatSize_);
	double x, y;
	double re = finestGrid->rho / finestGrid->mu;
	double lambda = 0.5 * re - std::sqrt(0.25 * re * re + 4 * EIGEN_PI * EIGEN_PI);
	for (int i = 0; i < finestGrid->laplaceMatSize_; i++) {
		x = std::get<0>(finestGrid->points_[i]);
		y = std::get<1>(finestGrid->points_[i]);
		actual.coeffRef(i) = 1 - std::exp(lambda * x) * std::cos(2 * EIGEN_PI * y);
	}
	Eigen::VectorXd error = (actual - *finestGrid->u);
	cout << error.lpNorm<1>() / finestGrid->laplaceMatSize_ << endl;
	for (int i = 0; i < finestGrid->laplaceMatSize_; i++) {
		//cout << "error: " << error.coeff(i) << " x " << std::get<0>(finestGrid->points_[i]) << " y " << std::get<1>(finestGrid->points_[i]) << endl;
	}
	vector<Point> points = finestGrid->points_;
	vector<double> xv, yv;
	for (size_t i = 0; i < points.size(); i++) {
		x = std::get<0>(points[i]);
		y = std::get<1>(points[i]);
		xv.push_back(x);
		yv.push_back(y);
	}

	writeVectorToTxt(xv, "fsx.txt");
	writeVectorToTxt(yv, "fsy.txt");

	vector<double> temp, uvec, vvec, pvec, actualvec, sourcevec;
	for (int i = 0; i < error.rows(); i++) {
		actualvec.push_back(actual.coeff(i));
		temp.push_back(error.coeff(i));
		uvec.push_back(finestGrid->u->coeff(i));
		vvec.push_back(finestGrid->v->coeff(i));
		pvec.push_back(finestGrid->values_->coeff(i));
		sourcevec.push_back(finestGrid->source_.coeff(i));
	}
	writeVectorToTxt(temp, "fserror.txt");
	writeVectorToTxt(uvec, "fsuvec.txt");
	writeVectorToTxt(vvec, "fsvvec.txt");
	writeVectorToTxt(pvec, "fspvec.txt");
	writeVectorToTxt(actualvec, "fsactual.txt");
	writeVectorToTxt(sourcevec, "fssource.txt");
	double vCycTime = (std::clock() - start) / (double)(CLOCKS_PER_SEC);
}
void run_frac_step_test() {
	FractionalStepParams param = gen_fracstep_param(4, 6, 0.0002, 0.025, 1, 1e-10);
	run_fracstep_param(param, 5);
}
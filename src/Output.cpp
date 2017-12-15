#include "Output.h"

void OutputPressure(Matrix p, int n, std::list<Circle> iList, Param par){

	std::ofstream output;
	std::string filename = par.WorkDir + "solution_pressure" + std::to_string(n) + ".plt";

	output.open(filename);

	output << "title = " << '"' << "sample mesh" << '"' << std::endl;
	output << "Variables = x y p" << std::endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << par.N1 + 1 << ", j=" << par.N2 + 1 << ", f=point" << std::endl;
	output << "SolutionTime = " << n << std::endl;

	for (int j = 0; j <= par.N2; ++j){
		for (int i = 0; i <= par.N1; ++i){
			GeomVec xp = x_p(i, j, par);
			output << xp[1] << ' ' << xp[2] << ' ' << p[i][j] << std::endl;
		}
	}

	for (auto& solid : iList){
		output << "zone T = circle" << ",  i=" << solid.Nn << ", f=point" << std::endl;
		output << "SolutionTime = " << n << std::endl;
		for (int i = 0; i < solid.Nn; ++i){
			output << solid.Nodes[i].x[1] << " "
			       << solid.Nodes[i].x[2] << " "
			       << 0 << std::endl;
		}
	}

	output.close();
}


void OutputVelocity_U(Matrix u, int n, std::list<Circle> iList, Param par){

	std::ofstream output;
	std::string filename = par.WorkDir + "solution_velocity_u" + std::to_string(n) + ".plt";

	output.open(filename);

	output << "title = " << '"' << "sample mesh" << '"' << std::endl;
	output << "Variables = x y u v" << std::endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << par.N1 << ", j=" << par.N2 + 1 << ", f=point" << std::endl;
	output << "SolutionTime = " << n << std::endl;

	for (int j = 0; j <= par.N2; ++j){
		for (int i = 0; i < par.N1; ++i){
			GeomVec xu = x_u(i, j, par);
			output << xu[1] << ' ' << xu[2] << ' ' << u[i][j] << ' ' << 0 << std::endl;
		}
	}

	for (auto& solid : iList){
		output << "zone T = circle" << ",  i=" << solid.Nn << ", f=point" << std::endl;
		output << "SolutionTime = " << n << std::endl;
		for (int i = 0; i < solid.Nn; ++i) {
			output << solid.Nodes[i].x[1] << " "
			       << solid.Nodes[i].x[2] << " "
			       << solid.Nodes[i].f[1] << " "
			       << solid.Nodes[i].f[2] << " " << std::endl;
		}
	}

	output.close();
}


void OutputVelocity_V(Matrix v, int n, std::list<Circle> iList, Param par){

	std::ofstream output;

	// Solution v

	std::string filename = par.WorkDir + "solution_velocity_v" + std::to_string(n) + ".plt";

	output.open(filename);

	output << "title = " << '"' << "sample mesh" << '"' << std::endl;
	output << "Variables = x y u v" << std::endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << par.N1 + 1 << ", j=" << par.N2 << ", f=point" << std::endl;
	output << "SolutionTime = " << n << std::endl;

	for (int j = 0; j < par.N2; ++j){
		for (int i = 0; i <= par.N1; ++i){
			GeomVec xv = x_v(i, j, par);
			output << xv[1] << ' ' << xv[2] << ' ' << 0 << ' ' << v[i][j] << std::endl;
		}
	}

	for (auto& solid : iList){
		output << "zone T = circle" << ",  i=" << solid.Nn << ", f=point" << std::endl;
		output << "SolutionTime = " << n << std::endl;
		for (int i = 0; i < solid.Nn; ++i) {
			output << solid.Nodes[i].x[1] << " "
			       << solid.Nodes[i].x[2] << " "
			       << solid.Nodes[i].f[1] << " "
			       << solid.Nodes[i].f[2] << " " << std::endl;
		}
	}

	output.close();
}

void Output(Matrix p, Matrix u, Matrix v, Matrix Fx, Matrix Fy, int n, std::list<Circle> iList, Param par) {

	std::ofstream output;
	std::string filename = par.WorkDir + "step" + std::to_string(n) + ".plt";

	output.open(filename);

	output << "title = " << '"' << "sample mesh" << '"' << std::endl;
	output << "Variables = x y p u v fx fy" << std::endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << par.N1 + 1 << ", j=" << par.N2 + 1 << ", f=point" << std::endl;
	output << "SolutionTime = " << n << std::endl;

	for (int j = 0; j <= par.N2; ++j) {
		for (int i = 0; i <= par.N1; ++i) {
			GeomVec xp = x_p(i, j, par);
			double u_, v_, Fx_, Fy_;
			if      (i == 0)      {u_ =  u[0][j]                     ;	Fx_ = Fx[0][j];}
			else if (i == par.N1) {u_ =  u[par.N1 - 1][j]            ;	Fx_ = Fx[par.N1 - 1][j];}
			else                  {u_ = (u[i - 1][j] + u[i][j]) * 0.5;	Fx_ = (Fx[i - 1][j] + Fx[i][j]) * 0.5;}
			if      (j == 0)      {v_ =  v[i][0]                     ;	Fy_ = Fy[i][0];}
			else if (j == par.N2) {v_ =  v[i][par.N2 - 1];           ;	Fy_ = Fy[i][par.N2 - 1];}
			else                  {v_ = (v[i][j - 1] + v[i][j]) * 0.5;	Fy_ = (Fy[i][j - 1] + Fy[i][j]) * 0.5;}
			output << xp[1] << " " << xp[2] << " " << p[i][j] << " " << u_ << " " << v_ << " " << Fx_ << " " << Fy_ << std::endl;
		}
	}

	for (auto& solid : iList) {
		output << "zone T = circle" << ",  i=" << solid.Nn+2 << ", f=point" << std::endl;
		output << "SolutionTime = " << n << std::endl;

			output << solid.xc[1]          << " "
			       << solid.xc[2]          << " "
			       << 0                    << " "
			       << solid.uc[1]          << " "
			       << solid.uc[2]          << " "
			       << solid.F_hd[1]        << " "
			       << solid.F_hd[2]        << " "
			       << std::endl;
		for (int i = 0; i < solid.Nn; ++i) {
			output << solid.Nodes[i].x[1]  << " "
			       << solid.Nodes[i].x[2]  << " "
			       << solid.Nodes[i].p     << " "
			       << solid.Nodes[i].uf[1] << " "
			       << solid.Nodes[i].uf[2] << " "
			       << solid.Nodes[i].f [1] << " "
			       << solid.Nodes[i].f [2] << " "
			       << std::endl;
		}
			output << solid.Nodes[0].x[1]  << " "
			       << solid.Nodes[0].x[2]  << " "
			       << solid.Nodes[0].p     << " "
			       << solid.Nodes[0].uf[1] << " "
			       << solid.Nodes[0].uf[2] << " "
			       << solid.Nodes[0].f [1] << " "
			       << solid.Nodes[0].f [2] << " "
			       << std::endl;
	}

	output.close();
}

void Output_dp(Matrix dp, int n, Param par) {

	std::ofstream output;
	std::string filename = par.WorkDir + "dp" + std::to_string(n) + ".plt";

	output.open(filename);

	output << "title = " << '"' << "sample mesh" << '"' << std::endl;
	output << "Variables = x y dp" << std::endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << dp.size() << ", j=" << dp[0].size() << ", f=point" << std::endl;
	output << "SolutionTime = " << n << std::endl;

	for (int j = 0; j <= par.N2; ++j) {
		for (int i = 0; i <= par.N1; ++i) {
			GeomVec xp = x_p(i, j, par);
			output << xp[1] << " " << xp[2] << " " << dp[i][j] << std::endl;
		}
	}

	output.close();
}

void Output_c(Matrix c, int n, Param par) {

	std::ofstream output;
	std::string filename = par.WorkDir + "c" + std::to_string(n) + ".plt";

	output.open(filename);

	output << "title = " << '"' << "sample mesh" << '"' << std::endl;
	output << "Variables = x y c" << std::endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << c.size() << ", j=" << c[0].size() << ", f=point" << std::endl;
	output << "SolutionTime = " << n << std::endl;

	for (int j = 0; j < par.N2; ++j) {
		for (int i = 0; i < par.N1; ++i) {
			GeomVec xc = x_c(i, j, par);
			output << xc[1] << " " << xc[2] << " " << c[i][j] << std::endl;
		}
	}

	output.close();
}

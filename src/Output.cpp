#include "Output.h"

void OutputPressure(Matrix p, int n, std::list<Circle> iList, Param par, std::string WorkDir){

	std::ofstream output;
	std::string filename = WorkDir + "solution_pressure" + std::to_string(n) + ".plt";

	output.open(filename.c_str());

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


void OutputVelocity_U(Matrix u, int n, std::list<Circle> iList, Param par, std::string WorkDir){

	std::ofstream output;
	std::string filename = WorkDir + "solution_velocity_u" + std::to_string(n) + ".plt";

	output.open(filename.c_str());

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
			       << solid.Nodes[i].us[1] << " "
			       << solid.Nodes[i].us[2] << " " << std::endl;
		}
	}

	output.close();
}


void OutputVelocity_V(Matrix v, int n, std::list<Circle> iList, Param par, std::string WorkDir){

	std::ofstream output;

	// Solution v

	std::string filename = WorkDir + "solution_velocity_v" + std::to_string(n) + ".plt";

	output.open(filename.c_str());

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
			       << solid.Nodes[i].us[1] << " "
			       << solid.Nodes[i].us[2] << " " << std::endl;
		}
	}

	output.close();
}

void Output(Matrix p, Matrix u, Matrix v, int n, std::list<Circle> iList, Param par, std::string WorkDir) {

	std::ofstream output;
	std::string filename = WorkDir + "step" + std::to_string(n) + ".plt";

	output.open(filename.c_str());

	output << "title = " << '"' << "sample mesh" << '"' << std::endl;
	output << "Variables = x y p u v" << std::endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << par.N1 + 1 << ", j=" << par.N2 + 1 << ", f=point" << std::endl;
	output << "SolutionTime = " << n << std::endl;

	for (int j = 0; j <= par.N2; ++j) {
		for (int i = 0; i <= par.N1; ++i) {
			GeomVec xp = x_p(i, j, par);
			double u_, v_;
			if      (i == 0)      u_ =  u[0][j];
			else if (i == par.N1) u_ =  u[par.N1 - 1][j];
			else                  u_ = (u[i - 1][j] + u[i][j]) * 0.5;
			if      (j == 0)      v_ =  v[i][0];
			else if (j == par.N2) v_ =  v[i][par.N2 - 1];
			else                  v_ = (v[i][j - 1] + v[i][j]) * 0.5;
			output << xp[1] << " " << xp[2] << " " << p[i][j] << " " << u_ << " " << v_ << std::endl;
		}
	}

	for (auto& solid : iList) {
		output << "zone T = circle" << ",  i=" << solid.Nn << ", f=point" << std::endl;
		output << "SolutionTime = " << n << std::endl;
		for (int i = 0; i < solid.Nn; ++i) {
			output << solid.Nodes[i].x[1] << " "
			       << solid.Nodes[i].x[2] << " "
			       << solid.omega[3]      << " "
			       << solid.Nodes[i].uf[1] << " "
			       << solid.Nodes[i].uf[2] << " " << std::endl;
		}
	}

	output.close();
}

void Output_dp(Matrix dp, int n, Param par, std::string WorkDir) {

	std::ofstream output;
	std::string filename = WorkDir + "dp" + std::to_string(n) + ".plt";

	output.open(filename.c_str());

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

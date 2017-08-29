#include "Output.h"

void OutputPressure(Matrix data, int n, std::list<Circle> iList, Param par){

	std::ofstream output;
	int i = 0;
	int j = 0;

	// Solution p

	std::string filename = "Result/solution_pressure" + std::to_string(n) + ".plt";

	output.open(filename.c_str());

	output << "title = " << '"' << "sample mesh" << '"' << std::endl;
	output << "Variables = x y p" << std::endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << par.N1 + 1 << ", j=" << par.N2 + 1 << ", f=point" << std::endl;
	output << "SolutionTime = " << n << std::endl;

	output << 0.0 << ' ' << 0.0 << ' ' << data[0][0] << std::endl;

	for (i = 1; i < par.N1; ++i){

		output << (i - 0.5)*par.d_x << ' ' << 0.0 << ' ' << data[i][0] << std::endl;
	}
	i = par.N1;
	output << (i - 1)*par.d_x << ' ' << 0.0 << ' ' << data[i][0] << std::endl;


	for (j = 1; j < par.N2; ++j){

		output << 0.0 << ' ' << (j - 0.5)*par.d_y << ' ' << data[0][j] << std::endl;

		for (i = 1; i < par.N1; ++i){
			output << (i - 0.5)*par.d_x << ' ' << (j - 0.5)*par.d_y << ' ' << data[i][j] << std::endl;
		}
		i = par.N1;
		output << (i - 1)*par.d_x << ' ' << (j - 0.5)*par.d_y << ' ' << data[i][j] << std::endl;
	}


	j = par.N2;
	output << 0.0 << ' ' << (j - 1)*par.d_y << ' ' << data[0][j] << std::endl;
	for (i = 1; i < par.N1; ++i){

		output << (i - 0.5)*par.d_x << ' ' << (j - 1)*par.d_y << ' ' << data[i][j] << std::endl;
	}
	i = par.N1;
	output << (i - 1)*par.d_x << ' ' << (j - 1)*par.d_y << ' ' << data[i][j] << std::endl;

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


void OutputVelocity_U(Matrix data, int n, std::list<Circle> iList, Param par){

	std::ofstream output;

	// Solution u

	std::string filename = "Result/solution_velocity_u" + std::to_string(n) + ".plt";

	output.open(filename.c_str());

	output << "title = " << '"' << "sample mesh" << '"' << std::endl;
	output << "Variables = x y u v" << std::endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << par.N1 << ", j=" << par.N2 + 1 << ", f=point" << std::endl;
	output << "SolutionTime = " << n << std::endl;

	for (int i = 0; i < par.N1; ++i){

		output << i*par.d_x << ' ' << 0.0 << ' ' << data[i][0] << ' ' << 0 << std::endl;
	}

	for (int j = 1; j < par.N2; ++j){
		for (int i = 0; i < par.N1; ++i){
			output << (i)*par.d_x << ' ' << (j - 0.5)*par.d_y << ' ' << data[i][j] << ' ' << 0 << std::endl;
		}
	}

	for (int i = 0; i < par.N1; ++i){

		int j = par.N2;
		output << i*par.d_x << ' ' << (j - 1)*par.d_y << ' ' << data[i][j] << ' ' << 0 << std::endl;
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


void OutputVelocity_V(Matrix data, int n, std::list<Circle> iList, Param par){

	std::ofstream output;

	// Solution v

	std::string filename = "Result/solution_velocity_v" + std::to_string(n) + ".plt";

	output.open(filename.c_str());

	output << "title = " << '"' << "sample mesh" << '"' << std::endl;
	output << "Variables = x y u v" << std::endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << par.N1 + 1 << ", j=" << par.N2 << ", f=point" << std::endl;
	output << "SolutionTime = " << n << std::endl;

	for (int j = 0; j < par.N2; ++j){


		output << 0.0 << ' ' << (j)*par.d_y << ' ' << 0 << ' ' << data[0][j] << std::endl;

		for (int i = 1; i < par.N1; ++i){

			output << (i - 0.5)*par.d_x << ' ' << (j)*par.d_y << ' ' << 0 << ' ' << data[i][j] << std::endl;

		}

		int i = par.N1;

		output << (i - 1)*par.d_x << ' ' << (j)*par.d_y << ' ' << 0 << ' ' << data[i][j] << std::endl;
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

#include "stdafx.h"
#include "Output.h"

void OutputPressure(Matrix data, int n, double output_step, list<Circle> iList, Grid grid){

	ofstream output;
	int i = 0;
	int j = 0;

	// Solution p

	string filename = "Result/solution_pressure" + to_string(n) + ".plt";

	output.open(filename.c_str());

	output << "title = " << '"' << "sample mesh" << '"' << endl;
	output << "Variables=" << '"' << 'x' << '"' << ',' << '"' << 'y' << '"' << ',' << '"' << 'p' << '"' << endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << grid.N1 + 1 << ", j=" << grid.N2 + 1 << ", f=point" << endl;


	output << 0.0 << ' ' << 0.0 << ' ' << data[0][0] << endl;

	for (i = 1; i < grid.N1; ++i){

		output << (i - 0.5)*grid.d_x << ' ' << 0.0 << ' ' << data[i][0] << endl;
	}
	i = grid.N1;
	output << (i - 1)*grid.d_x << ' ' << 0.0 << ' ' << data[i][0] << endl;


	for (j = 1; j < grid.N2; ++j){

		output << 0.0 << ' ' << (j - 0.5)*grid.d_y << ' ' << data[0][j] << endl;

		for (i = 1; i < grid.N1; ++i){
			output << (i - 0.5)*grid.d_x << ' ' << (j - 0.5)*grid.d_y << ' ' << data[i][j] << endl;
		}
		i = grid.N1;
		output << (i - 1)*grid.d_x << ' ' << (j - 0.5)*grid.d_y << ' ' << data[i][j] << endl;
	}


	j = grid.N2;
	output << 0.0 << ' ' << (j - 1)*grid.d_y << ' ' << data[0][j] << endl;
	for (i = 1; i < grid.N1; ++i){

		output << (i - 0.5)*grid.d_x << ' ' << (j - 1)*grid.d_y << ' ' << data[i][j] << endl;
	}
	i = grid.N1;
	output << (i - 1)*grid.d_x << ' ' << (j - 1)*grid.d_y << ' ' << data[i][j] << endl;

	for (auto& solid : iList){

		output << "GEOMETRY ZN = " << n / output_step + 1 << ", X=" << solid.Bound[0][0] << ", Y=" << solid.Bound[1][0] << ", T=LINE, LT = 0.3, CS=GRID" << endl;
		output << "1" << endl << grid.NF + 1 << endl;
		output << 0 << " " << 0 << endl;
		for (int i = 1; i < grid.NF; ++i){

			output << solid.Bound[0][i] - solid.Bound[0][0] << " " << solid.Bound[1][i] - solid.Bound[1][0] << endl;

		}
		output << 0 << " " << 0 << endl;
		output << endl;
	}


	output.close();


}


void OutputVelocity_U(Matrix data, int n, int output_step, list<Circle> iList, Grid grid){

	ofstream output;

	// Solution u

	string filename = "Result/solution_velocity_u" + to_string(n) + ".plt";

	output.open(filename.c_str());

	output << "title = " << '"' << "sample mesh" << '"' << endl;
	output << "Variables=" << '"' << 'x' << '"' << ',' << '"' << 'y' << '"' << ',' << '"' << 'u' << '"' << endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << grid.N1 << ", j=" << grid.N2 + 1 << ", f=point" << endl;

	for (int i = 0; i < grid.N1; ++i){

		output << i*grid.d_x << ' ' << 0.0 << ' ' << data[i][0] << endl;
	}

	for (int j = 1; j < grid.N2; ++j){
		for (int i = 0; i < grid.N1; ++i){
			output << (i)*grid.d_x << ' ' << (j - 0.5)*grid.d_y << ' ' << data[i][j] << endl;
		}
	}

	for (int i = 0; i < grid.N1; ++i){

		int j = grid.N2;
		output << i*grid.d_x << ' ' << (j - 1)*grid.d_y << ' ' << data[i][j] << endl;
	}

	for (auto& solid : iList){

		output << "GEOMETRY ZN = " << n / output_step + 1 << ", X=" << solid.Bound[0][0] << ", Y=" << solid.Bound[1][0] << ", T=LINE, LT = 0.3, CS=GRID" << endl;
		output << "1" << endl << grid.NF + 1 << endl;
		output << 0 << " " << 0 << endl;
		for (int i = 1; i < grid.NF; ++i){

			output << solid.Bound[0][i] - solid.Bound[0][0] << " " << solid.Bound[1][i] - solid.Bound[1][0] << endl;

		}
		output << 0 << " " << 0 << endl;
		output << endl;
	}

	output.close();


}


void OutputVelocity_V(Matrix& data, int n, int output_step, list<Circle> iList, Grid grid){

	ofstream output;

	// Solution v

	string filename = "Result/solution_velocity_v" + to_string(n) + ".plt";

	output.open(filename.c_str());

	output << "title = " << '"' << "sample mesh" << '"' << endl;
	output << "Variables=" << '"' << 'x' << '"' << ',' << '"' << 'y' << '"' << ',' << '"' << 'v' << '"' << endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << grid.N1 + 1 << ", j=" << grid.N2 << ", f=point" << endl;


	for (int j = 0; j < grid.N2; ++j){


		output << 0.0 << ' ' << (j)*grid.d_y << ' ' << data[0][j] << endl;

		for (int i = 1; i < grid.N1; ++i){

			output << (i - 0.5)*grid.d_x << ' ' << (j)*grid.d_y << ' ' << data[i][j] << endl;

		}

		int i = grid.N1;

		output << (i - 1)*grid.d_x << ' ' << (j)*grid.d_y << ' ' << data[i][j] << endl;
	}

	for (auto& solid : iList){

		output << "GEOMETRY ZN = " << n / output_step + 1 << ", X=" << solid.Bound[0][0] << ", Y=" << solid.Bound[1][0] << ", T=LINE, LT = 0.3, CS=GRID" << endl;
		output << "1" << endl << grid.NF + 1 << endl;
		output << 0 << " " << 0 << endl;
		for (int i = 1; i < grid.NF; ++i){

			output << solid.Bound[0][i] - solid.Bound[0][0] << " " << solid.Bound[1][i] - solid.Bound[1][0] << endl;

		}
		output << 0 << " " << 0 << endl;
		output << endl;
	}

	output.close();


}

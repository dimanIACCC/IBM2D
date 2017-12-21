#include "stdafx.h"
#include "Parameters.h"



Param::Param() {
	// default parameters
	Re = 20;
	L = 30;
	H = 7;
	N1 = 101;
	N2 = 21;
	d_t = 0.00125;
	Nn = 50;
	rho = 10;
	r = 0.5;
	output_step = 50;
	N_max = 5000000;
	N_Zeidel = 500000;
	Zeidel_eps = 1e-5;
	InelasticCollision = false;
	k_dist = 1.1;
	AddSolids_N = 0;
	AddSolids_start = 0;
	AddSolids_interval = 200;
	BC = u_inflow;
	SolidName_max = 0;
	WorkDir = "";

	d_x = L / (N1 - 1);
	d_y = H / (N2 - 1);
}

Param::Param(std::string WorkDir, std::string filename) : Param() {
	std::ifstream input;
	std::string line;
	WorkDir = WorkDir;

	input.open(WorkDir + filename);
	if (input.is_open()) {
		while (getline(input, line)) { // read line from file to string $line$
			std::string PAR, VALUE;
			GetParValue(line, PAR, VALUE);
			if (VALUE.size() > 0) {
				if (PAR == "Re")           Re = stod(VALUE);
				else if (PAR == "L")            L = stod(VALUE);
				else if (PAR == "H")            H = stod(VALUE);
				else if (PAR == "N1")           N1 = stoi(VALUE);
				else if (PAR == "N2")           N2 = stoi(VALUE);
				else if (PAR == "d_t")          d_t = stod(VALUE);
				else if (PAR == "Nn")           Nn = stoi(VALUE);
				else if (PAR == "rho")          rho = stod(VALUE);
				else if (PAR == "r")            r = stod(VALUE);
				else if (PAR == "output_step")  output_step = stoi(VALUE);
				else if (PAR == "N_max")        N_max = stoi(VALUE);
				else if (PAR == "N_Zeidel")     N_Zeidel = stoi(VALUE);
				else if (PAR == "Zeidel_eps")   Zeidel_eps = stod(VALUE);
				else if (PAR == "InelasticCollision")   InelasticCollision = bool(stoi(VALUE));
				else if (PAR == "k_dist")       k_dist = (stod(VALUE));
				else if (PAR == "AddSolids_N")          AddSolids_N = stoi(VALUE);
				else if (PAR == "AddSolids_start")      AddSolids_start = stoi(VALUE);
				else if (PAR == "AddSolids_interval")   AddSolids_interval = stoi(VALUE);
				else if (PAR == "BC")                   BC = string_to_BC(VALUE);
				else    std::cout << "unknown parameter " << PAR << std::endl;
			}
			else {
				std::cout << filename << ": no value inputed" << std::endl;
			}
		}
	}
	else {
		std::cout << "File " << filename << " is not found" << std::endl;
		std::cout << "Using default parameters" << std::endl;
	}

	this->WorkDir = WorkDir;
	d_x = L / (N1 - 1);
	d_y = H / (N2 - 1);

}
Param::Param(std::string WorkDir) : Param() {
	this->WorkDir = WorkDir;
}


Boundary_Conditions string_to_BC(std::string s) {
	Boundary_Conditions BC;
	if      (s == "u_infinity") BC = u_infinity;
	else if (s == "u_inflow")   BC = u_inflow;
	else if (s == "periodical") BC = periodical;
	else std::cout << "string_to_BC: unknown BC" << std::endl;
	return BC;
}

double ux_Poiseuille(double y, double H) {
	double ux = (pow(H / 2.0, 2) - pow(y - H / 2.0, 2)) / pow(H / 2.0, 2);
	return ux;
}

double dpdx_Poiseuille(double H, double Re) {
	return 8.0 / H / H / Re;
}

double dux_dy_Poiseuille(double y, double H) {
	double dux_dy = - 2.0 * (y - H / 2.0) / pow(H / 2.0, 2);
	return dux_dy;
}

GeomVec x_p(int i, int j, Param par) {
	GeomVec result;
	result[0] = 0.0;
	result[1] = (i - 0.5) * par.d_x;
	result[2] = (j - 0.5) * par.d_y;
	result[3] = 0.0;
	//if (i == 0     ) result[1] = 0.0;
	//if (j == 0     ) result[2] = 0.0;
	//if (i == par.N1) result[1] = (i - 1) * par.d_x;
	//if (j == par.N2) result[2] = (j - 1) * par.d_y;
	return result;
}

GeomVec x_u(int i, int j, Param par) {
	GeomVec result;
	result[0] = 0.0;
	result[1] =  i        * par.d_x;
	result[2] = (j - 0.5) * par.d_y;
	result[3] = 0.0;
	//if (j == 0     ) result[2] = 0.0;
	//if (j == par.N2) result[2] = (j - 1) * par.d_y;
	return result;
}

GeomVec x_v(int i, int j, Param par) {
	GeomVec result;
	result[0] = 0.0;
	result[1] = (i - 0.5) * par.d_x;
	result[2] =  j        * par.d_y;
	result[3] = 0.0;
	//if (i == 0) result[1] = 0.0;
	//if (i == par.N1) result[1] = (i - 1) * par.d_x;
	return result;
}

GeomVec x_c(int i, int j, Param par) {
	GeomVec result;
	result[0] = 0.0;
	result[1] = i * par.d_x;
	result[2] = j * par.d_y;
	result[3] = 0.0;

	return result;
}

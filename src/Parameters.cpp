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
	alpha_f = 0;
	beta_f = -2000;
	Nn = 50;
	rho = 10 / (M_PI * 0.5 * 0.5); // corresponds to old formula for force
	r = 0.5;
	output_step = 50;
	N_max = 5000000;
	N_Zeidel = 500000;
	Zeidel_eps = 1e-5;
	InelasticCollision = false;
	k_dist = 1.1;

	d_x = L / (N1 - 1);
	d_y = H / (N2 - 1);
}

Param::Param(std::string filename): Param(){
	std::ifstream input;
	std::string line;


	input.open(filename.c_str());
	if (input.is_open()) {
		while (getline(input, line)) { // read line from file to string $line$
			std::string PAR, VALUE;
			GetParValue(line, PAR, VALUE);
			if (VALUE.size() > 0) {
				if      (PAR == "Re")           Re = stod(VALUE);
				else if (PAR == "L")            L = stod(VALUE);
				else if (PAR == "H")            H = stod(VALUE);
				else if (PAR == "N1")           N1 = stoi(VALUE);
				else if (PAR == "N2")           N2 = stoi(VALUE);
				else if (PAR == "d_t")          d_t = stod(VALUE);
				else if (PAR == "alpha_f")      alpha_f = stod(VALUE);
				else if (PAR == "beta_f")       beta_f = stod(VALUE);
				else if (PAR == "Nn")           Nn = stoi(VALUE);
				else if (PAR == "rho")          rho = stod(VALUE);
				else if (PAR == "r")            r = stod(VALUE);
				else if (PAR == "output_step")  output_step = stoi(VALUE);
				else if (PAR == "N_max")        N_max = stoi(VALUE);
				else if (PAR == "N_Zeidel")     N_Zeidel = stoi(VALUE);
				else if (PAR == "Zeidel_eps")   Zeidel_eps = stod(VALUE);
				else if (PAR == "InelasticCollision")   InelasticCollision = bool(stoi(VALUE));
				else if (PAR == "k_dist")       k_dist = (stod(VALUE));
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

	d_x = L / (N1 - 1);
	d_y = H / (N2 - 1);

}

double ux_Poiseuille(double y, double H) {
	double ux = (pow(H / 2.0, 2) - pow(y - H / 2.0, 2));
	return ux;
}

GeomVec x_p(int i, int j, Param par) {
	GeomVec result;
	result[0] = 0.0;
	result[1] = (i - 0.5) * par.d_x;
	result[2] = (j - 0.5) * par.d_y;
	result[3] = 0.0;
	if (i == 0     ) result[1] = 0.0;
	if (j == 0     ) result[2] = 0.0;
	if (i == par.N1) result[1] = (i - 1) * par.d_x;
	if (j == par.N2) result[2] = (j - 1) * par.d_y;
	return result;
}

GeomVec x_u(int i, int j, Param par) {
	GeomVec result;
	result[0] = 0.0;
	result[1] =  i        * par.d_x;
	result[2] = (j - 0.5) * par.d_y;
	result[3] = 0.0;
	if (j == 0     ) result[2] = 0.0;
	if (j == par.N2) result[2] = (j - 1) * par.d_y;
	return result;
}

GeomVec x_v(int i, int j, Param par) {
	GeomVec result;
	result[0] = 0.0;
	result[1] = (i - 0.5) * par.d_x;
	result[2] =  j        * par.d_y;
	result[3] = 0.0;
	if (i == 0) result[1] = 0.0;
	if (i == par.N1) result[1] = (i - 1) * par.d_x;
	return result;
}

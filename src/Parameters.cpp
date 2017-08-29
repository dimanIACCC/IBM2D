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
	NF = 50;
	R = 0.5;
	output_step = 50;
	N_max = 5000000;
	N_Zeidel = 500000;
	Zeidel_eps = 1e-5;

	d_x = L / (N1 - 1);
	d_y = H / (N2 - 1);
}

Param::Param(std::string filename): Param(){
	std::ifstream input;
	std::string line;


	input.open(filename.c_str());
	if (input.is_open()) {
		while (getline(input, line)) { // read line from file to string $line$
			int i = line.find('=');
			std::string PAR(line, 0, i - 1);
			PAR = boost::trim_copy(PAR);
			if (i > 0) {
				std::string VALUE(line, i + 1);
				VALUE = boost::trim_copy(VALUE);
				if (PAR == "Re")                Re = stod(VALUE);
				else if (PAR == "L")            L = stod(VALUE);
				else if (PAR == "H")            H = stod(VALUE);
				else if (PAR == "N1")           N1 = stoi(VALUE);
				else if (PAR == "N2")           N2 = stoi(VALUE);
				else if (PAR == "d_t")          d_t = stod(VALUE);
				else if (PAR == "alpha_f")      alpha_f = stod(VALUE);
				else if (PAR == "beta_f")       beta_f = stod(VALUE);
				else if (PAR == "NF")           NF = stoi(VALUE);
				else if (PAR == "R")            R = stod(VALUE);
				else if (PAR == "output_step")  output_step = stoi(VALUE);
				else if (PAR == "N_max")        N_max = stoi(VALUE);
				else if (PAR == "N_Zeidel")     N_Zeidel = stoi(VALUE);
				else if (PAR == "Zeidel_eps")   Zeidel_eps = stod(VALUE);
				else    std::cout << "unknown parameter " << PAR << std::endl;
			}
			else {
				std::cout << "no value inputed" << std::endl;
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

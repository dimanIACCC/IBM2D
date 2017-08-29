#pragma once

#include "Parameters.h"

void InputData(Param& par) {
	std::ifstream input;
	std::string filename = "input.txt";
	std::string line;

	// default parameters
	par.Re = 20;
	par.L = 30;
	par.H = 7;
	par.N1 = 101;
	par.N2 = 21;
	par.d_t = 0.00125;
	par.alpha_f = 0;
	par.beta_f = -2000;
	par.NF = 50;
	par.R = 0.5;
	par.output_step = 50;
	par.N_max = 5000000;
	par.N_Zeidel = 500000;
	par.Zeidel_eps = 1e-5;


	input.open(filename.c_str());
	while (getline(input, line)) { // read line from file to string $line$
		int i = line.find('=');
		std::string PAR(line, 0, i - 1);
		PAR = boost::trim_copy(PAR);
		if (i > 0) {
			std::string VALUE(line, i + 1);
			VALUE = boost::trim_copy(VALUE);
			if (PAR == "Re")           par.Re = stod(VALUE);
			else if (PAR == "L")            par.L = stod(VALUE);
			else if (PAR == "H")            par.H = stod(VALUE);
			else if (PAR == "N1")           par.N1 = stoi(VALUE);
			else if (PAR == "N2")           par.N2 = stoi(VALUE);
			else if (PAR == "d_t")          par.d_t = stod(VALUE);
			else if (PAR == "alpha_f")      par.alpha_f = stod(VALUE);
			else if (PAR == "beta_f")       par.beta_f = stod(VALUE);
			else if (PAR == "NF")           par.NF = stoi(VALUE);
			else if (PAR == "R")            par.R = stod(VALUE);
			else if (PAR == "output_step")  par.output_step = stoi(VALUE);
			else if (PAR == "N_max")        par.N_max = stoi(VALUE);
			else if (PAR == "N_Zeidel")     par.N_Zeidel = stoi(VALUE);
			else if (PAR == "Zeidel_eps")   par.Zeidel_eps = stod(VALUE);
			else    std::cout << "unknown parameter " << PAR << std::endl;
		}
		else {
			std::cout << "no value inputed" << std::endl;
		}
	}

	par.d_x = par.L / (par.N1 - 1);
	par.d_y = par.H / (par.N2 - 1);

}

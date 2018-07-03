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
	N_step = -1;
	Nn = 50;
	rho = 10;
	r = 0.5;
	output_step = 50;
	N_max = 5000000;
	N_Zeidel = 500000;
	Zeidel_eps = 1e-5;
	eps_P = 1.e-3;
	InelasticCollision = false;
	k_dist = 1.1;
	AddSolids_N = 0;
	AddSolids_start = 0;
	AddSolids_interval = 200;
	BC = u_inflow;
	u_wall = 0;
	SolidName_max = 0;
	WorkDir = "";

	this->init();
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
				if      (PAR == "Re")           Re = stod(VALUE);
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
				else if (PAR == "eps_P")        eps_P = stod(VALUE);
				else if (PAR == "InelasticCollision")   InelasticCollision = bool(stoi(VALUE));
				else if (PAR == "k_dist")               k_dist = (stod(VALUE));
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

	this->init();

}
Param::Param(std::string WorkDir) : Param() {
	this->WorkDir = WorkDir;
}

void Param::init() {
	N1_u = N1 + 2;
	N2_u = N2 + 1;
	N1_v = N1 + 1;
	N2_v = N2 + 2;

	d_x = L / (N1 - 1);
	d_y = H / (N2 - 1);

	ldxdx = 1 / (d_x*d_x);
	ldydy = 1 / (d_y*d_y);

	x0[1] = 0.5 * L;
	x0[2] = 0.5 * H;
	if (BC == Line_Vortex) {
		x0[1] = - 0.5 * L;
		x0[2] = - 0.5 * H;
	}

	k[1] = M_PI / L;
	k[2] = M_PI / H;
}


boundary_conditions string_to_BC(std::string s) {
	boundary_conditions BC;
	if      (s == "u_infinity"   || s == "0") BC = u_infinity;
	else if (s == "u_inflow"     || s == "1") BC = u_inflow;
	else if (s == "periodical"   || s == "2") BC = periodical;
	else if (s == "Taylor_Green" || s == "3") BC = Taylor_Green;
	else if (s == "Lamb_Oseen"   || s == "4") BC = Lamb_Oseen;
	else if (s == "Line_Vortex"  || s == "4") BC = Line_Vortex;
	else std::cout << "string_to_BC: unknown BC" << std::endl;
	return BC;
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
	result[1] = (i - 1.0) * par.d_x;
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
	result[2] = (j - 1.0) * par.d_y;
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

double DeltaFunction(double x, double y, Param par) {
	return FunctionD(x / par.d_x) * FunctionD(y / par.d_y) / (par.d_x*par.d_y);
}

double FunctionD(double r) {
	if ((0.0 <= fabs(r)) && (fabs(r) < 1.0)) {
		return 1.0 / 8.0*(3.0 - 2.0 * fabs(r) + sqrt(1.0 + 4.0 * fabs(r) - 4.0 * r * r));
	}
	if ((1.0 <= fabs(r)) && (fabs(r) < 2.0)) {
		return 1.0 / 8.0*(5.0 - 2.0 * fabs(r) - sqrt(-7.0 + 12.0 * fabs(r) - 4.0 * r * r));
	}
	if (2.0 <= fabs(r)) {
		return 0.0;
	}
	return 0;
}

void GetInfluenceArea(int &i_min, int &i_max, int &j_min, int &j_max, size_t Ni, size_t Nj, GeomVec x, int size, Param par) {
	i_max = (int)(x[1] / par.d_x) + size;
	i_min = (int)(x[1] / par.d_x) - size;

	j_max = (int)(x[2] / par.d_y) + size;
	j_min = (int)(x[2] / par.d_y) - size;

	if (i_min < 0 && par.BC != periodical) {
		i_min = 0;
	}
	if (j_min < 0) {
		j_min = 0;
	}
	if (i_max > Ni && par.BC != periodical) {
		i_max = Ni;
	}
	if (j_max > Nj) {
		j_max = Nj;
	}
}

double Volume_Frac(GeomVec xc, double r, GeomVec x, double dx, double dy) {
	GeomVec x_i[5];
	for (int i = 1; i <= 4; ++i) {
		x_i[i] = x;
	}
	x_i[1][1] += dx;	x_i[1][2] += dy;
	x_i[2][1] += dx;	x_i[2][2] -= dy;
	x_i[3][1] -= dx;	x_i[3][2] -= dy;
	x_i[4][1] -= dx;	x_i[4][2] += dy;

	double result = 0;
	double sum_phi = 0;
	for (int i = 1; i <= 4; ++i) {
		double phi = length(x_i[i] - xc) - r;
		result -= phi * Heaviside(-phi);
		sum_phi += std::abs(phi);
	}
	result /= sum_phi;
	return result;
}

double Heaviside(double x) {
	double result;
	result = (x >= 0) ? 1 : 0;
	return result;
}

int i_real_u(int i, Param par) {
	int i_real = i;
	if (i_real < 1     ) i_real += par.N1 - 1;
	if (i_real > par.N1) i_real -= par.N1 - 1;
	if (i_real == 1    ) i_real  = par.N1;
	return i_real;
}

int i_real_v(int i, Param par) {
	int i_real = i;
	if (i_real < 1         ) i_real += par.N1 - 1;
	if (i_real > par.N1 - 1) i_real -= par.N1 - 1;
	return i_real;
}

#include "stdafx.h"
#include "Parameters.h"



Param::Param() {  // default parameters
	this->init();
}

void Param::init() {
	N1_u = N1 + 3;
	N2_u = N2 + 2;
	N1_v = N1 + 2;
	N2_v = N2 + 3;
	N1_p = N1 + 2;
	N2_p = N2 + 2;

	d_x = L / N1;
	d_y = H / N2;

	ldxdx = 1 / (d_x*d_x);
	ldydy = 1 / (d_y*d_y);

	x0[1] = 0.5 * L;
	x0[2] = 0.5 * H;
	if (BC == Line_Vortex) {
		x0[1] = -0.5 * L;
		x0[2] = -0.5 * H;
	}

	k[1] = 2 * M_PI / L;
	k[2] = 2 * M_PI / H;

	Gravity_angle *= M_PI / 180.;

	Gravity = ZeroVec();
	Gravity[1] =  Gravity_module*sin(Gravity_angle);
	Gravity[2] = -Gravity_module*cos(Gravity_angle);
}

void Param::read(std::ifstream &input) {
	std::string line = "<par>";
	do { //while (line != "</par>")
		getline(input, line);
		if (line != "</par>") this->read_line(line);
	} while (line != "</par>");
	this->init();
}

void Param::read_line(std::string line) {
	std::string PAR, VALUE;
	GetParValue(line, PAR, VALUE);
	if (VALUE.size() > 0) {
		if      (PAR == "N_step")               N_step = stoi(VALUE);
		else if (PAR == "time")                 time = stod(VALUE);
		else if (PAR == "BC")                   BC = string_to_BC(VALUE);

		// physical parameters
		else if (PAR == "d_t")                  d_t = stod(VALUE);
		else if (PAR == "L")                    L = stod(VALUE);
		else if (PAR == "H")                    H = stod(VALUE);
		else if (PAR == "Re")                   Re = stod(VALUE);
		else if (PAR == "grad_p_x")             grad_p_x = (stod(VALUE));
		else if (PAR == "Gravity_angle")        Gravity_angle = (stod(VALUE));
		else if (PAR == "Gravity_module")       Gravity_module = (stod(VALUE));

		// numerical parameters
		else if (PAR == "N1")                   N1 = stoi(VALUE);
		else if (PAR == "N2")                   N2 = stoi(VALUE);
		else if (PAR == "Nn_")                  Nn_ = stoi(VALUE);
		else if (PAR == "output_step")          output_step = stoi(VALUE);
		else if (PAR == "IBM")                  IBM = stoi(VALUE);
		else if (PAR == "DeltaP_method")        DeltaP_method = stoi(VALUE);
		else if (PAR == "N_Zeidel")             N_Zeidel = stoi(VALUE);
		else if (PAR == "Zeidel_eps")           Zeidel_eps = stod(VALUE);
		else if (PAR == "s_max")                s_max = stoi(VALUE);
		else if (PAR == "eps_P")                eps_P = stod(VALUE);

		// parameters for many particles
		else if (PAR == "rho")                  rho = stod(VALUE);
		else if (PAR == "shape")                shape = stoi(VALUE);
		else if (PAR == "r0")                   r0 = stod(VALUE);
		else if (PAR == "r")                    r = stod(VALUE);
		else if (PAR == "e")                    e = stod(VALUE);
		else if (PAR == "AddSolids_N")          AddSolids_N = stoi(VALUE);
		else if (PAR == "AddSolids_start")      AddSolids_start = stoi(VALUE);
		else if (PAR == "AddSolids_interval")   AddSolids_interval = stoi(VALUE);
		else if (PAR == "SolidName_max")        SolidName_max = stoi(VALUE);
		else if (PAR == "k_u_dist")             k_u_dist = (stod(VALUE));
		else if (PAR == "k_r_dist")             k_r_dist = (stod(VALUE));
		else if (PAR == "k_u_collide")          k_u_collide = (stod(VALUE));
		else if (PAR == "k_r_collide")          k_r_collide = (stod(VALUE));

		// parameters for special problems
		else if (PAR == "Lamb_Oseen_r0")        Lamb_Oseen_r0 = (stod(VALUE));
		else if (PAR == "u_in")                 u_in = (stod(VALUE));
		else if (PAR == "u_down")               u_down = (stod(VALUE));
		else if (PAR == "u_up")                 u_up   = (stod(VALUE));
		else if (PAR == "omega_BC")             omega_BC = (stod(VALUE));

		else    std::cout << "unknown parameter " << PAR << std::endl;
	}
	else {
		std::cout << "Param::read_line(): unknown parameter: " << line << std::endl;
	}

}

boundary_conditions string_to_BC(std::string s) {
	boundary_conditions BC;
	if      (s == "u_infinity"   || s == "0") BC = u_infinity;
	else if (s == "u_inflow"     || s == "1") BC = u_inflow;
	else if (s == "periodical"   || s == "2") BC = periodical;
	else if (s == "Taylor_Green" || s == "3") BC = Taylor_Green;
	else if (s == "Lamb_Oseen"   || s == "4") BC = Lamb_Oseen;
	else if (s == "Line_Vortex"  || s == "5") BC = Line_Vortex;
	else if (s == "box"          || s == "6") BC = box;
	else std::cout << "string_to_BC: unknown BC" << std::endl;
	return BC;
}

GeomVec x_p(int i, int j, Param par) {
	GeomVec result = ZeroVec();
	result[1] = (i - 0.5) * par.d_x;
	result[2] = (j - 0.5) * par.d_y;
	return result;
}

GeomVec x_u(int i, int j, Param par) {
	GeomVec result = ZeroVec();
	result[1] = (i - 1.0) * par.d_x;
	result[2] = (j - 0.5) * par.d_y;
	return result;
}

GeomVec x_v(int i, int j, Param par) {
	GeomVec result = ZeroVec();
	result[1] = (i - 0.5) * par.d_x;
	result[2] = (j - 1.0) * par.d_y;
	return result;
}

GeomVec x_c(int i, int j, Param par) {
	GeomVec result = ZeroVec();
	result[1] = i * par.d_x;
	result[2] = j * par.d_y;
	return result;
}


double* x_p_(int i, int j, double d_x, double d_y) restrict(amp) {
	double result[4];
	result[0] = 0.0;
	result[1] = (i - 0.5) * d_x;
	result[2] = (j - 0.5) * d_y;
	result[3] = 0.0;
	return result;
}

double* x_u_(int i, int j, double d_x, double d_y) restrict(amp) {
	double result[4];
	result[0] = 0.0;
	result[1] = (i - 1.0) * d_x;
	result[2] = (j - 0.5) * d_y;
	result[3] = 0.0;
	return result;
}

double* x_v_(int i, int j, double d_x, double d_y) restrict(amp) {
	double result[4];
	result[0] = 0.0;
	result[1] = (i - 0.5) * d_x;
	result[2] = (j - 1.0) * d_y;
	result[3] = 0.0;
	return result;
}

double* x_c_(int i, int j, double d_x, double d_y) restrict(amp) {
	double result[4];
	result[0] = 0.0;
	result[1] = i * d_x;
	result[2] = j * d_y;
	result[3] = 0.0;

	return result;
}

double DeltaFunction(double x, double y, double d_x, double d_y) {
	return FunctionD(x / d_x) * FunctionD(y / d_y);
}

double DeltaFunction_(double x, double y, double d_x, double d_y) restrict(amp) {
	return FunctionD_(x / d_x) * FunctionD_(y / d_y);
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
	return 0.;
}

inline double FunctionD_(double r) restrict(amp) {
	if ((0.0 <= fabs(r)) && (fabs(r) < 1.0)) {
		return 1.0 / 8.0*(3.0 - 2.0 * fabs(r) + sqrt(1.0 + 4.0 * fabs(r) - 4.0 * r * r));
	}
	if ((1.0 <= fabs(r)) && (fabs(r) < 2.0)) {
		return 1.0 / 8.0*(5.0 - 2.0 * fabs(r) - sqrt(-7.0 + 12.0 * fabs(r) - 4.0 * r * r));
	}
	if (2.0 <= fabs(r)) {
		return 0.0;
	}
	return 0.;
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

InfluenceArea GetInfluenceArea_(int Ni, int Nj, double* x, int size, boundary_conditions BC, double d_x, double d_y) restrict(amp) {
	InfluenceArea IA;
	IA.i_max = (int)(x[1] / d_x) + size;
	IA.i_min = (int)(x[1] / d_x) - size;

	IA.j_max = (int)(x[2] / d_y) + size;
	IA.j_min = (int)(x[2] / d_y) - size;

	if (IA.i_min < 0 && BC != periodical) {
		IA.i_min = 0;
	}
	if (IA.j_min < 0) {
		IA.j_min = 0;
	}
	if (IA.i_max > Ni && BC != periodical) {
		IA.i_max = Ni;
	}
	if (IA.j_max > Nj) {
		IA.j_max = Nj;
	}

	return IA;
}


double ellipse_rho(double r, double e, double phi) {
	double k = sqrt(1 - e*e);
	double a = r / sqrt(k);
	double b = r * sqrt(k);
	//double ellipse_rho = b / sqrt(1 - e*e*cos(phi)*cos(phi));
	double ellipse_rho = a*b / sqrt(a*a*sin(phi)*sin(phi) + b*b*cos(phi)*cos(phi));
	return ellipse_rho;
}

double Volume_Frac(GeomVec xc, double r, double e, double alpha0, GeomVec x, double dx, double dy) {
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
		GeomVec r_fluid = x_i[i] - xc;
		double alpha = atan(r_fluid[2] / r_fluid[1]);
		double rho_solid = ellipse_rho(r,e,alpha-alpha0);
		double phi = length(r_fluid) - rho_solid;
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
	if (i_real < 1         ) i_real += par.N1;
	if (i_real > par.N1 + 1) i_real -= par.N1;
	if (i_real == 1        ) i_real  = par.N1 + 1;
	return i_real;
}

int i_real_u_(int i, int N1) restrict(amp) {
	int i_real = i;
	if (i_real < 1) i_real += N1;
	if (i_real > N1 + 1) i_real -= N1;
	if (i_real == 1) i_real = N1 + 1;
	return i_real;
}

int i_real_v(int i, Param par) {
	int i_real = i;
	if (i_real < 1         ) i_real += par.N1;
	if (i_real > par.N1    ) i_real -= par.N1;
	return i_real;
}

int i_real_v_(int i, int N1) restrict(amp) {
	int i_real = i;
	if (i_real < 1) i_real += N1;
	if (i_real > N1) i_real -= N1;
	return i_real;
}

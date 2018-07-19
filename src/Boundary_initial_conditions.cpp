#include "Boundary_initial_conditions.h"

void Boundary_Conditions(Matrix &u, Param par, Direction Dir, double time) {

	size_t Nx = u.size();
	size_t Ny = u[0].size();

	// Up-Down BC
	if (par.BC == u_inflow || par.BC == u_infinity || par.BC == periodical) {
		for (size_t i = 0; i < Nx; ++i) {
			if (Dir == Du) {
				u[i][0           ] = 2 * par.u_wall - u[i][1];
				u[i][par.N2_u - 1] = 2 * par.u_wall - u[i][par.N2_u - 2];

				// zero derivative of the horizontal velocity
				if (par.BC == u_infinity) {
					u[i][0           ] = u[i][1           ];
					u[i][par.N2_u - 1] = u[i][par.N2_u - 2];
				}
			}
			else if (Dir == Dv) {
				u[i][0           ] = 2 * u[i][1           ] - u[i][2           ];
				u[i][par.N2_v - 1] = 2 * u[i][par.N2_v - 2] - u[i][par.N2_v - 3];
			}
		}
	}

	// Left-Right BC
	if (par.BC == u_infinity || par.BC == u_inflow) {
		for (size_t j = 0; j < Ny; ++j) {
			if (Dir == Du) {
				u[par.N1_u - 1][j] = u[par.N1_u - 2][j];
			}
			else if (Dir == Dv) {
				u[par.N1_v - 1][j] = u[par.N1_v - 2][j];
			}
		}
	}

	// Periodical BC
	if (par.BC == periodical) {
		for (size_t j = 0; j < Ny; ++j) {
			if (Dir == Du) {
				u[0         ][j] = u[par.N1][j];
				u[par.N1 + 2][j] = u[2     ][j];
				u[par.N1 + 1][j] = u[1     ][j];
			}
			else if (Dir == Dv) {
				u[0         ][j] = u[par.N1][j];
				u[par.N1 + 1][j] = u[1     ][j];
			}
		}
	}

	if ((par.BC == Taylor_Green || par.BC == Lamb_Oseen || par.BC == Line_Vortex) && time > 0) {
		// Up-Down BC
		for (size_t i = 0; i < Nx; ++i) {
			GeomVec x_D;
			GeomVec x_U;
			if (Dir == Du) {
				x_D = x_u(i, 0     , par);
				x_U = x_u(i, Ny - 1, par);
				u[i][0     ] = exact_u(x_D, par, time);
				u[i][Ny - 1] = exact_u(x_U, par, time);
			}
			else if (Dir == Dv) {
				x_D = x_v(i, 0, par);
				x_U = x_v(i, Ny - 1, par);
				u[i][0     ] = exact_v(x_D, par, time);
				u[i][Ny - 1] = exact_v(x_U, par, time);
			}
		}
		// Left-Right BC
		for (size_t j = 0; j < Ny; ++j) {
			GeomVec x_L;
			GeomVec x_R;
			if (Dir == Du) {
				x_L = x_u(0     , j, par);
				x_R = x_u(Nx - 1, j, par);
				u[0     ][j] = exact_u(x_L, par, time);
				u[Nx - 1][j] = exact_u(x_R, par, time);
			}
			else if (Dir == Dv) {
				x_L = x_v(0     , j, par);
				x_R = x_v(Nx - 1, j, par);
				u[0     ][j] = exact_v(x_L, par, time);
				u[Nx - 1][j] = exact_v(x_R, par, time);
			}
		}
	}
}

void fill_exact(Matrix &u, Matrix &v, Matrix &p, Param par, double time) {
	fill_exact_u(u, par, time);
	fill_exact_v(v, par, time);
	fill_exact_p(p, par, time);
}

void fill_exact_u(Matrix &u, Param par, double time) {
	for (size_t i = 0; i < u.size(); ++i) {
		for (size_t j = 0; j < u[0].size(); ++j) {
			u[i][j] = exact_u(x_u(i, j, par), par, time);
		}
	}
}

void fill_exact_v(Matrix &v, Param par, double time) {
	for (size_t i = 0; i < v.size(); ++i) {
		for (size_t j = 0; j < v[0].size(); ++j) {
			v[i][j] = exact_v(x_v(i, j, par), par, time);
		}
	}
}

void fill_exact_p(Matrix &p, Param par, double time) {
	for (size_t i = 0; i < p.size(); ++i) {
		for (size_t j = 0; j < p[0].size(); ++j) {
			p[i][j] = exact_p(x_p(i, j, par), par, time);
		}
	}
}

void BC_exact_p(Matrix &p, Param par, double time) {

	size_t Nx = p.size();
	size_t Ny = p[0].size();

	// Up-Down BC
	for (size_t i = 0; i < Nx; ++i) {
		GeomVec x_D =  x_p(i, 0     , par);
		GeomVec x_U =  x_p(i, Ny - 1, par);
		p[i][0     ] = exact_p(x_D, par, time);
		p[i][Ny - 1] = exact_p(x_U, par, time);

	}
	// Left-Right BC
	for (size_t j = 0; j < Ny; ++j) {
		GeomVec x_L =  x_p(0     , j, par);
		GeomVec x_R =  x_p(Nx - 1, j, par);
		p[0     ][j] = exact_p(x_L, par, time);
		p[Nx - 1][j] = exact_p(x_R, par, time);
	}
}

void Lamb_Oseen_p_test(Param par, double time) {

	GeomVec x0 = par.x0;
	GeomVec x = par.x0;


	std::ofstream Ei_bug;
	std::string filename = par.WorkDir + "/Ei_bug.plt";
	Ei_bug.open(filename);
	Ei_bug << "title = Ei_bug" << std::endl;
	Ei_bug << "Variables = x P_integral P_Ei d_P" << std::endl;

	for (size_t i = 0; i < 500; ++i) {
		x[1] += 0.001;
		double P_integral = Lamb_Oseen_p(x, x0, par.Re, time, par.Lamb_Oseen_r0, false);
		double P_Ei       = Lamb_Oseen_p(x, x0, par.Re, time, par.Lamb_Oseen_r0, true);
		Ei_bug << x[1] << " " << P_integral << " " << P_Ei << " " << P_integral - P_Ei << std::endl;;
	}
}

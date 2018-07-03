#include "Boundary_initial_conditions.h"

void Boundary_Conditions(Matrix &u, Param par, Direction Dir, int N_step) {

	size_t Nx = u.size();
	size_t Ny = u[0].size();

	// Up-Down BC
	if (par.BC == u_inflow || par.BC == u_infinity || par.BC == periodical) {
		for (size_t i = 0; i < Nx; ++i) {
			if (Dir == Du) {
				u[i][0] = 2 * par.u_wall - u[i][1];
				u[i][par.N2] = 2 * par.u_wall - u[i][par.N2 - 1];

				// zero derivative of the horizontal velocity
				if (par.BC == u_infinity) {
					u[i][0] = u[i][1];
					u[i][par.N2] = u[i][par.N2 - 1];
				}
			}
			else if (Dir == Dv) {
				u[i][0] = 2 * u[i][1] - u[i][2];
				u[i][par.N2 + 1] = 2 * u[i][par.N2] - u[i][par.N2 - 1];
			}
		}
	}

	// Left-Right BC
	if (par.BC == u_infinity || par.BC == u_inflow) {
		for (size_t j = 0; j < Ny; ++j) {
			if (Dir == Du) {
				u[par.N1 + 1][j] = u[par.N1 - 1][j];
			}
			else if (Dir == Dv) {
				u[par.N1][j] = u[par.N1 - 1][j];
			}
		}
	}

	// Periodical BC
	if (par.BC == periodical) {
		for (size_t j = 0; j < Ny; ++j) {
			if (Dir == Du) {
				u[0][j] = u[par.N1 - 1][j];
				u[par.N1 + 1][j] = u[2][j];
				u[par.N1][j] = u[1][j];
			}
			else if (Dir == Dv) {
				u[0][j] = u[par.N1 - 1][j];
				u[par.N1][j] = u[1][j];
			}
		}
	}

	if ((par.BC == Taylor_Green || par.BC == Lamb_Oseen || par.BC == Line_Vortex) && N_step > 0) {
		double time = par.d_t * N_step;
		// Up-Down BC
		for (size_t i = 0; i < Nx; ++i) {
			GeomVec x_D;
			GeomVec x_U;
			if (Dir == Du) {
				x_D = (x_u(i, 0, par) + x_u(i, 1, par)) * 0.5;
				x_U = (x_u(i, Ny - 1, par) + x_u(i, Ny - 2, par)) * 0.5;
				double u_D = exact_u(x_D, par, time);
				double u_U = exact_u(x_U, par, time);
				u[i][0] = (2 * u_D - u[i][1]);
				u[i][Ny - 1] = (2 * u_U - u[i][Ny - 2]);
			}
			else if (Dir == Dv) {
				x_D = x_v(i, 1, par);
				x_U = x_v(i, Ny - 2, par);
				double v_D = exact_v(x_D, par, time);
				double v_U = exact_v(x_U, par, time);
				u[i][0] = (2 * v_D - u[i][2]);
				u[i][Ny - 1] = (2 * v_U - u[i][Ny - 3]);
			}
		}
		// Left-Right BC
		for (size_t j = 0; j < Ny; ++j) {
			GeomVec x_L;
			GeomVec x_R;
			if (Dir == Du) {
				x_L = x_u(1, j, par);
				x_R = x_u(Nx - 2, j, par);
				double u_L = exact_u(x_L, par, time);
				double u_R = exact_u(x_R, par, time);
				u[0][j] = (2 * u_L - u[2][j]);
				u[Nx - 1][j] = (2 * u_R - u[Nx - 3][j]);
			}
			else if (Dir == Dv) {
				x_L = (x_v(0, j, par) + x_v(1, j, par)) * 0.5;
				x_R = (x_v(Nx - 1, j, par) + x_v(Nx - 2, j, par)) * 0.5;
				double v_L = exact_v(x_L, par, time);
				double v_R = exact_v(x_R, par, time);
				u[0][j] = (2 * v_L - u[1][j]);
				u[Nx - 1][j] = (2 * v_R - u[Nx - 2][j]);
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
		GeomVec x_D = (x_p(i, 0     , par) + x_p(i, 1     , par)) * 0.5;
		GeomVec x_U = (x_p(i, Ny - 1, par) + x_p(i, Ny - 2, par)) * 0.5;
		double p_D = exact_p(x_D, par, time);
		double p_U = exact_p(x_U, par, time);
		p[i][0     ] = (2 * p_D - p[i][1]);
		p[i][Ny - 1] = (2 * p_U - p[i][Ny - 2]);

	}
	// Left-Right BC
	for (size_t j = 0; j < Ny; ++j) {
		GeomVec x_L = (x_p(1     , j, par) + x_p(2     , j, par)) * 0.5;
		GeomVec x_R = (x_p(Nx - 1, j, par) + x_p(Nx - 2, j, par)) * 0.5;
		double p_L = exact_p(x_L, par, time);
		double p_R = exact_p(x_R, par, time);
		p[0     ][j] = (2 * p_L - p[2     ][j]);
		p[Nx - 1][j] = (2 * p_R - p[Nx - 2][j]);
	}
}


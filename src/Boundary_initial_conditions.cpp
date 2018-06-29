#include "Boundary_initial_conditions.h"

// Apply initial data for velocity and pressure
void ApplyInitialData(Matrix &u, Matrix &v, Matrix &p, Param par) {

	for (size_t i = 0; i < u.size(); ++i) {
		for (size_t j = 0; j < u[0].size(); ++j) {
			GeomVec xu = x_u(i, j, par);
			switch (par.BC) {
				case u_infinity:   u[i][j] = 1.0; break;
				case u_inflow:     u[i][j] = 1.0 * ux_Poiseuille(xu[2], par.H); break;
				case periodical:   u[i][j] = 1.0 * ux_Poiseuille(xu[2], par.H); break;
				case Taylor_Green: break;
				case Lamb_Oseen:   break;
				default: std::cout << "ApplyInitialData: unknown BC" << std::endl;
			}
		}
	}

	if (par.BC == u_inflow || par.BC == u_infinity || par.BC == periodical) {
		for (size_t i = 0; i < p.size(); ++i) {
			for (size_t j = 0; j < p[0].size(); ++j) {
				GeomVec xp = x_p(i, j, par);
				p[i][j] = (par.L - xp[1]) * dpdx_Poiseuille(par.H, par.Re);
			}
		}
	}

	if (par.BC == Taylor_Green) {
		Taylor_Green_exact(u, v, p, par, 0.);
	}

	if (par.BC == Lamb_Oseen) {
		Lamb_Oseen_exact_uv(u, Du, par, 0., false);
		Lamb_Oseen_exact_uv(v, Dv, par, 0., false);
		Lamb_Oseen_exact_p(p, par, 0.);
	}

}

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

	// Taylor-Green BC
	if (par.BC == Taylor_Green && N_step > 0) {
		double k1 = M_PI / par.L;
		double k2 = M_PI / par.H;
		double time_exp = exp(-(k1*k1 + k2*k2) / par.Re * par.d_t * N_step);
		// Up-Down BC
		for (size_t i = 0; i < Nx; ++i) {
			GeomVec x_D;
			GeomVec x_U;
			if (Dir == Du) {
				x_D = (x_u(i, 0, par) + x_u(i, 1, par)) * 0.5;
				x_U = (x_u(i, Ny - 1, par) + x_u(i, Ny - 2, par)) * 0.5;
				double u_D = Taylor_Green_u(x_D, k1, k2, time_exp);
				double u_U = Taylor_Green_u(x_U, k1, k2, time_exp);
				u[i][0] = (2 * u_D - u[i][1]);
				u[i][Ny - 1] = (2 * u_U - u[i][Ny - 2]);
			}
			else if (Dir == Dv) {
				x_D = x_v(i, 1, par);
				x_U = x_v(i, Ny - 2, par);
				double v_D = Taylor_Green_v(x_D, k1, k2, time_exp);
				double v_U = Taylor_Green_v(x_U, k1, k2, time_exp);
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
				double u_L = Taylor_Green_u(x_L, k1, k2, time_exp);
				double u_R = Taylor_Green_u(x_R, k1, k2, time_exp);
				u[0][j] = (2 * u_L - u[2][j]);
				u[Nx - 1][j] = (2 * u_R - u[Nx - 3][j]);
			}
			else if (Dir == Dv) {
				x_L = (x_v(0, j, par) + x_v(1, j, par)) * 0.5;
				x_R = (x_v(Nx - 1, j, par) + x_v(Nx - 2, j, par)) * 0.5;
				double v_L = Taylor_Green_v(x_L, k1, k2, time_exp);
				double v_R = Taylor_Green_v(x_R, k1, k2, time_exp);
				u[0][j] = (2 * v_L - u[1][j]);
				u[Nx - 1][j] = (2 * v_R - u[Nx - 2][j]);
			}
		}
	}

	if (par.BC == Lamb_Oseen && N_step > 0) {
		Lamb_Oseen_exact_uv(u, Dir, par, par.d_t * N_step, true);
	}

}

void Taylor_Green_exact(Matrix &u, Matrix &v, Matrix &p, Param par, double time) {

	double k1 = M_PI / par.L;
	double k2 = M_PI / par.H;
	double time_exp = exp(-(k1*k1 + k2*k2) / par.Re * time);
	double time_exp2 = time_exp*time_exp;

	for (size_t i = 0; i < u.size(); ++i) {
		for (size_t j = 0; j < u[0].size(); ++j) {
			GeomVec xu = x_u(i, j, par);
			u[i][j] = Taylor_Green_u(xu, k1, k2, time_exp);
		}
	}

	for (size_t i = 0; i < v.size(); ++i) {
		for (size_t j = 0; j < v[0].size(); ++j) {
			GeomVec xv = x_v(i, j, par);
			v[i][j] = Taylor_Green_v(xv, k1, k2, time_exp);
		}
	}

	for (size_t i = 0; i < par.N1 + 1; ++i) {
		for (size_t j = 0; j < par.N2 + 1; ++j) {
			GeomVec xp = x_p(i, j, par);
			p[i][j] = Taylor_Green_p(xp, k1, k2, time_exp2);
		}
	}
}

void Lamb_Oseen_exact_uv(Matrix &uv, Direction Dir, Param par, double time, bool boundary) {
	// The key boundary stands for the domain where apply the Lamb-Oseen solution: true - only in the boundary, false - in the whole domain

	GeomVec omega, r, r0;
	r0[1] = 0.5 * par.L;
	r0[2] = 0.5 * par.H;
	r0[3] = 0;

	omega[1] = 0;
	omega[2] = 0;

	if (Dir == Du) {
		for (size_t i = 0; i < uv.size(); ++i) {
			for (size_t j = 0; j < uv[0].size(); ++j) {
				if (!boundary || i == 0 || j == 0 || i == uv.size() - 1 || j == uv[0].size() - 1) {
					r = x_u(i, j, par) - r0;
					omega[3] = Lamb_Oseen_velocity(length(r), par.Re, time);
					GeomVec U = x_product(omega, r / length(r));
					uv[i][j] = U[1];
				}
			}
		}
	}

	if (Dir == Dv) {
		for (size_t i = 0; i < uv.size(); ++i) {
			for (size_t j = 0; j < uv[0].size(); ++j) {
				if (!boundary || i == 0 || j == 0 || i == uv.size() - 1 || j == uv[0].size() - 1) {
					r = x_v(i, j, par) - r0;
					omega[3] = Lamb_Oseen_velocity(length(r), par.Re, time);
					GeomVec V = x_product(omega, r / length(r));
					uv[i][j] = V[2];
				}
			}
		}
	}

}

// Calculate exact pressure for Lamb-Oseen vortex
void Lamb_Oseen_exact_p(Matrix &p, Param par, double time) {
	// The key boundary stands for the domain where apply the Lamb-Oseen solution: true - only in the boundary, false - in the whole domain

	GeomVec r, r0;
	r0[1] = 0.5 * par.L;
	r0[2] = 0.5 * par.H;
	r0[3] = 0;

	for (size_t i = 0; i < p.size(); ++i) {
		for (size_t j = 0; j < p[0].size(); ++j) {
			r = x_p(i, j, par) - r0;
			p[i][j] = Lamb_Oseen_pressure(length(r), par.Re, time);
		}
	}

}

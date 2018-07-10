#include "Boundary_initial_conditions.h"

void Boundary_Conditions(Matrix &u_n, Matrix &u, Param par, Direction Dir, double time) {

	size_t Nx = u.size();
	size_t Ny = u[0].size();

	// Up-Down BC
	if (par.BC == u_inflow || par.BC == u_infinity || par.BC == periodical) {
		for (size_t i = 0; i < Nx; ++i) {
			if (Dir == Du) {
				u[i][0     ] = 2 * par.u_wall - u[i][1];
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

	if ((par.BC == Taylor_Green || par.BC == Lamb_Oseen || par.BC == Line_Vortex) && time > 0) {
		// Up-Down BC
		for (size_t i = 0; i < Nx; ++i) {
			GeomVec x_D;
			GeomVec x_U;
			if (Dir == Du) {
				//x_D = (x_u(i, 0     , par) + x_u(i, 1     , par)) * 0.5;
				//x_U = (x_u(i, Ny - 1, par) + x_u(i, Ny - 2, par)) * 0.5;
				x_D = x_u(i, 0     , par);
				x_U = x_u(i, Ny - 1, par);
				double un_D = u_n[i][0     ];
				double un_U = u_n[i][Ny - 1];
				double u_D = 2 * exact_u(x_D, par, time) - un_D;
				double u_U = 2 * exact_u(x_U, par, time) - un_U;
				//u[i][0     ] = (2 * u_D - u[i][1     ]);
				//u[i][Ny - 1] = (2 * u_U - u[i][Ny - 2]);
				u[i][0     ] = u_D;
				u[i][Ny - 1] = u_U;
			}
			else if (Dir == Dv) {
				//x_D = x_v(i, 1     , par);
				//x_U = x_v(i, Ny - 2, par);
				x_D = x_v(i, 0, par);
				x_U = x_v(i, Ny - 1, par);
				double vn_D = u_n[i][0     ];
				double vn_U = u_n[i][Ny - 1];
				double v_D = 2 * exact_v(x_D, par, time) - vn_D;
				double v_U = 2 * exact_v(x_U, par, time) - vn_U;
				//u[i][0     ] = (2 * v_D - u[i][2     ]);
				//u[i][Ny - 1] = (2 * v_U - u[i][Ny - 3]);
				u[i][0     ] = v_D;
				u[i][Ny - 1] = v_U;
			}
		}
		// Left-Right BC
		for (size_t j = 0; j < Ny; ++j) {
			GeomVec x_L;
			GeomVec x_R;
			if (Dir == Du) {
				//x_L = x_u(1     , j, par);
				//x_R = x_u(Nx - 2, j, par);
				x_L = x_u(0     , j, par);
				x_R = x_u(Nx - 1, j, par);
				double un_L = u_n[0     ][j];
				double un_R = u_n[Nx - 1][j];
				double u_L = 2 * exact_u(x_L, par, time) - un_L;
				double u_R = 2 * exact_u(x_R, par, time) - un_R;
				//u[0     ][j] = (2 * u_L - u[2     ][j]);
				//u[Nx - 1][j] = (2 * u_R - u[Nx - 3][j]);
				u[0     ][j] = u_L;
				u[Nx - 1][j] = u_R;
			}
			else if (Dir == Dv) {
				//x_L = (x_v(0     , j, par) + x_v(1     , j, par)) * 0.5;
				//x_R = (x_v(Nx - 1, j, par) + x_v(Nx - 2, j, par)) * 0.5;
				x_L = x_v(0     , j, par);
				x_R = x_v(Nx - 1, j, par);
				double vn_L = u_n[0     ][j];
				double vn_R = u_n[Nx - 1][j];
				double v_L = 2 * exact_v(x_L, par, time) - vn_L;
				double v_R = 2 * exact_v(x_R, par, time) - vn_R;
				//u[0     ][j] = (2 * v_L - u[1     ][j]);
				//u[Nx - 1][j] = (2 * v_R - u[Nx - 2][j]);
				u[0     ][j] = v_L;
				u[Nx - 1][j] = v_R;
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

void BC_exact_p(Matrix &p_n, Matrix &p_new, Param par, double time) {

	size_t Nx = p_n.size();
	size_t Ny = p_n[0].size();

	// Up-Down BC
	for (size_t i = 0; i < Nx; ++i) {
		//GeomVec x_D = (x_p(i, 0     , par) + x_p(i, 1     , par)) * 0.5;
		//GeomVec x_U = (x_p(i, Ny - 1, par) + x_p(i, Ny - 2, par)) * 0.5;
		GeomVec x_D =  x_p(i, 0     , par);
		GeomVec x_U =  x_p(i, Ny - 1, par);
		double pn_D = 0.5 * (p_n[i][0     ] + p_n[i][0     ]);
		double pn_U = 0.5 * (p_n[i][Ny - 1] + p_n[i][Ny - 1]);
		double p_D = 2 * exact_p(x_D, par, time) - pn_D;
		double p_U = 2 * exact_p(x_U, par, time) - pn_U;
		//p_new[i][0     ] = (2 * p_D - p_new[i][1     ]);
		//p_new[i][Ny - 1] = (2 * p_U - p_new[i][Ny - 2]);
		p_new[i][0     ] = p_D;
		p_new[i][Ny - 1] = p_U;

	}
	// Left-Right BC
	for (size_t j = 0; j < Ny; ++j) {
		//GeomVec x_L = (x_p(0     , j, par) + x_p(1     , j, par)) * 0.5;
		//GeomVec x_R = (x_p(Nx - 1, j, par) + x_p(Nx - 2, j, par)) * 0.5;
		GeomVec x_L =  x_p(0     , j, par);
		GeomVec x_R =  x_p(Nx - 1, j, par);
		double pn_L = 0.5 * (p_n[0     ][j] + p_n[0     ][j]);
		double pn_R = 0.5 * (p_n[Nx - 1][j] + p_n[Nx - 1][j]);
		double p_L = 2 * exact_p(x_L, par, time) - pn_L;
		double p_R = 2 * exact_p(x_R, par, time) - pn_R;
		//p_new[0     ][j] = (2 * p_L - p_new[2     ][j]);
		//p_new[Nx - 1][j] = (2 * p_R - p_new[Nx - 2][j]);
		p_new[0     ][j] = p_L;
		p_new[Nx - 1][j] = p_R;
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

#include "stdafx.h"
#include "PredictVel.h"

// LHS of Navier-Stokes equation operator
// LHS = (1 / d_t - 1/Re * 0.5 \Delta ) * U_new
void Calculate_A(Template &A, Param par, double Re, Direction Dir) {

	double d_u, d_v;

	if      (Dir == Du) { d_u = par.d_x;	d_v = par.d_y;	}
	else if (Dir == Dv) { d_u = par.d_y;	d_v = par.d_x;	}

	double d_uu = 1.0 / (d_u*d_u);
	double d_vv = 1.0 / (d_v*d_v);

	A.C =  1.0 /      Re *(d_uu + d_vv) + 1.0 / par.d_t;
	A.U = -1.0 / (2.0*Re)* d_vv;
	A.R = -1.0 / (2.0*Re)* d_uu;
	A.D = -1.0 / (2.0*Re)* d_vv;
	A.L = -1.0 / (2.0*Re)* d_uu;
}

// LHS of Navier-Stokes equation
Matrix Operator_Ax(Template &A, Matrix &u, Param par, Direction Dir) {

	size_t Nx = u.size();
	size_t Ny = u[0].size();

	CreateMatrix(result, Nx, Ny);

	Boundary_Conditions(u, par, Dir, -1);

	// For direction Du  (U, R, D, L)  are  (up, right, down, left)  neighbour elements in matrix u[i][j]
	// For direction Dv  (U, R, D, L)  are  (up, right, down, left)  neighbour elements in transpose matrix (v[i][j])^T
	// In periodical problem near the boundary the neighbour element is taken from the other side of the domain
	for (size_t j = 1; j < Ny - 1; ++j) {
		for (size_t i = 1; i < Nx - 1; ++i) {
			result[i][j] = A.C * u[i][j]
			             + A.U * U(u, i, j, Dir)
			             + A.R * R(u, i, j, Dir)
			             + A.D * D(u, i, j, Dir)
			             + A.L * L(u, i, j, Dir);
		}
	}

	return result;

}

// RHS of Navier-Stokes equation
Matrix CalculateB(Matrix &u_n, Matrix &v_n, Matrix &u_s, Matrix &v_s, Matrix &p, Matrix &p_new, Param par, Direction Dir) {

	size_t Nx = u_n.size();
	size_t Ny = u_n[0].size();

	double d_u, d_v;

	if      (Dir == Du) { d_u = par.d_x;	d_v = par.d_y;}
	else if (Dir == Dv) { d_u = par.d_y;	d_v = par.d_x;}

	double d_uu = 1.0 / (d_u*d_u);
	double d_vv = 1.0 / (d_v*d_v);

	double alpha = 0.5;

	CreateMatrix(result, Nx, Ny);

	for (size_t i = 1; i < (Nx - 1); ++i) {
		for (size_t j = 1; j < (Ny - 1); ++j) {
			double advective_term_n    = advective_term(u_n, v_n, i, j, d_u, d_v, Dir);
			double advective_term_s    = advective_term(u_s, v_s, i, j, d_u, d_v, Dir);
			double diffusion_term_n    = par.ldxdx * (u_n[i + 1][j] - 2.0 * u_n[i][j] + u_n[i - 1][j])
			                           + par.ldydy * (u_n[i][j + 1] - 2.0 * u_n[i][j] + u_n[i][j - 1]);
			double pressure_term = 0.5 * (
			                                (p    [i][j] - L(p    , i, j, Dir)) / d_u
			                              + (p_new[i][j] - L(p_new, i, j, Dir)) / d_u
			                             );

			result[i][j] = -(       alpha  * advective_term_n
			               + (1.0 - alpha) * advective_term_s)
			                               - pressure_term
			                               + diffusion_term_n / (2.0*par.Re)
			                               + u_n[i][j] / par.d_t;
		}
	}

	Boundary_Conditions(u_s, par, Dir, par.N_step + 1);

	return result;
}

void Boundary_Conditions(Matrix &u, Param par, Direction Dir, int N_step) {

	size_t Nx = u.size();
	size_t Ny = u[0].size();

	// Up-Down BC
	if (par.BC == u_inflow || par.BC == u_infinity || par.BC == periodical) {
		for (size_t i = 0; i < Nx; ++i) {
			if (Dir == Du) {
				u[i][0     ] = 2 * par.u_wall - u[i][1         ];
				u[i][par.N2] = 2 * par.u_wall - u[i][par.N2 - 1];

				// zero derivative of the horizontal velocity
				if (par.BC == u_infinity) {
					u[i][0     ] = u[i][1         ];
					u[i][par.N2] = u[i][par.N2 - 1];
				}
			}
			else if (Dir == Dv) {
				u[i][0         ] = 2 * u[i][1     ] - u[i][2         ];
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
				u[par.N1    ][j] = u[par.N1 - 1][j];
			}
		}
	}

	// Periodical BC
	if (par.BC == periodical) {
		for (size_t j = 0; j < Ny; ++j) {
			if (Dir == Du) {
				u[0         ][j] = u[par.N1 - 1][j];
				u[par.N1 + 1][j] = u[2         ][j];
				u[par.N1    ][j] = u[1         ][j];
			}
			else if (Dir == Dv) {
				u[0         ][j] = u[par.N1 - 1][j];
				u[par.N1    ][j] = u[1         ][j];
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
				x_D = (x_u(i, 0     , par) + x_u(i, 1     , par)) * 0.5;
				x_U = (x_u(i, Ny - 1, par) + x_u(i, Ny - 2, par)) * 0.5;
				double u_D = sin(k1 * x_D[1]) * cos(k2 * x_D[2]) * time_exp;
				double u_U = sin(k1 * x_U[1]) * cos(k2 * x_U[2]) * time_exp;
				u[i][0     ] = (2 * u_D - u[i][1     ]);
				u[i][Ny - 1] = (2 * u_U - u[i][Ny - 2]);
			}
			else if (Dir == Dv) {
				x_D = x_v(i, 1, par);
				x_U = x_v(i, Ny - 2, par);
				double v_D = -k1 / k2 * sin(k2 * x_D[2]) * cos(k1 * x_D[1]) * time_exp;
				double v_U = -k1 / k2 * sin(k2 * x_U[2]) * cos(k1 * x_U[1]) * time_exp;
				u[i][0     ] = (2 * v_D - u[i][2     ]);
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
				double u_L = sin(k1 * x_L[1]) * cos(k2 * x_L[2]) * time_exp;
				double u_R = sin(k1 * x_R[1]) * cos(k2 * x_R[2]) * time_exp;
				u[0     ][j] = (2 * u_L - u[2     ][j]);
				u[Nx - 1][j] = (2 * u_R - u[Nx - 3][j]);
			}
			else if (Dir == Dv) {
				x_L = (x_v(0     , j, par) + x_v(1     , j, par)) * 0.5;
				x_R = (x_v(Nx - 1, j, par) + x_v(Nx - 2, j, par)) * 0.5;
				double v_L = -k1 / k2 * sin(k2 * x_L[2]) * cos(k1 * x_L[1]) * time_exp;
				double v_R = -k1 / k2 * sin(k2 * x_R[2]) * cos(k1 * x_R[1]) * time_exp;
				u[0     ][j] = (2 * v_L - u[1     ][j]);
				u[Nx - 1][j] = (2 * v_R - u[Nx - 2][j]);
			}
		}
	}

}

double advective_term(Matrix &u, Matrix &v, size_t i, size_t j, double d_x, double d_y, Direction Dir) {
	double v_help = 0.25 * (v[i][j] + L(v, i, j, Dir) + U(v, i, j, Dir) + UL(v, i, j, Dir));
	double result = u[i][j] * (R(u, i, j, Dir) - L(u, i, j, Dir)) / (2.0*d_x)
	               + v_help * (U(u, i, j, Dir) - D(u, i, j, Dir)) / (2.0*d_y);
	return result;
}

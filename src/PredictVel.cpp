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
	CreateMatrix(u_n, Nx, Ny);

	Boundary_Conditions(u_n, u, par, Dir, -1.);

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

	Boundary_Conditions(u_n, u_s, par, Dir, par.d_t * (par.N_step + 0.5));

	return result;
}

double advective_term(Matrix &u, Matrix &v, size_t i, size_t j, double d_x, double d_y, Direction Dir) {
	double v_help = 0.25 * (v[i][j] + L(v, i, j, Dir) + U(v, i, j, Dir) + UL(v, i, j, Dir));
	double result = u[i][j] * (R(u, i, j, Dir) - L(u, i, j, Dir)) / (2.0*d_x)
	               + v_help * (U(u, i, j, Dir) - D(u, i, j, Dir)) / (2.0*d_y);
	return result;
}

#include "stdafx.h"
#include "PredictVel.h"

// LHS of Navier-Stokes equation operator
// LHS = (1 / d_t - 1/Re * 0.5 \Delta ) * U_new
void Calculate_A(Template &A, Param par, double Re) {
	A.C =  1.0 /      Re *(par.ldxdx + par.ldydy) + 1.0 / par.d_t;
	A.U = -1.0 / (2.0*Re)* par.ldydy;
	A.R = -1.0 / (2.0*Re)* par.ldxdx;
	A.D = -1.0 / (2.0*Re)* par.ldydy;
	A.L = -1.0 / (2.0*Re)* par.ldxdx;
}

// LHS of Navier-Stokes equation
Matrix Operator_Ax(Template &A, Matrix &u, Param par, Direction Dir) {

	size_t Nx = u.size();
	size_t Ny = u[0].size();

	CreateMatrix(result, Nx, Ny);

	Boundary_Conditions(u, par, Dir, -1.);

	for (size_t j = 1; j < Ny - 1; ++j) {
		for (size_t i = 1; i < Nx - 1; ++i) {
			result[i][j] = A.C * u[i][j]
			             + A.U * u[i][j+1]
			             + A.R * u[i+1][j]
			             + A.D * u[i][j-1]
			             + A.L * u[i-1][j];
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

	double alpha = 0.5;

	CreateMatrix(result, Nx, Ny);

	for (size_t i = 1; i < (Nx - 1); ++i) {
		for (size_t j = 1; j < (Ny - 1); ++j) {
			double advective_term_s = 0.25 * ( advective_term(u_n, v_n, u_n, v_n, i, j, d_u, d_v, Dir)
			                                 + advective_term(u_s, v_s, u_s, v_s, i, j, d_u, d_v, Dir)
			                                 + advective_term(u_n, v_n, u_s, v_s, i, j, d_u, d_v, Dir)
			                                 + advective_term(u_s, v_s, u_n, v_n, i, j, d_u, d_v, Dir) );
			double diffusion_term_n    = par.ldxdx * (u_n[i + 1][j] - 2.0 * u_n[i][j] + u_n[i - 1][j])
			                           + par.ldydy * (u_n[i][j + 1] - 2.0 * u_n[i][j] + u_n[i][j - 1]);
			double pressure_term = 0.5 * (
			                                (p    [i][j] - L(p    , i, j, Dir)) / d_u
			                              + (p_new[i][j] - L(p_new, i, j, Dir)) / d_u
			                             );

			result[i][j] = - advective_term_s
			               - pressure_term
			               + diffusion_term_n / (2.0*par.Re)
			               + u_n[i][j] / par.d_t;

			/*int iii = 10, jjj = 3;
			std::cout << std::setprecision(18);
			if ((Dir == Du && i == iii && j == jjj) || (Dir == Dv && i == jjj && j == Nx / 2 + 1 - iii)) {
				std::cout << "Direction = " << Dir << std::endl;
				std::cout << "(i,j)      = (" << i << ", " << j << ")" << std::endl;
				std::cout << "convection = " << std::setw(22) << std::fixed << - (alpha  * advective_term_n + (1.0 - alpha) * advective_term_s) << std::endl;
				std::cout << "pressure   = " << std::setw(22) << std::fixed << - pressure_term << std::endl;
				std::cout << "diffusion  = " << std::setw(22) << std::fixed << diffusion_term_n / (2.0*par.Re) << std::endl;
				std::cout << "u_n / dt   = " << std::setw(22) << std::fixed << u_n[i][j] / par.d_t << std::endl;
				std::cout << "result     = " << std::setw(22) << std::fixed << result[i][j] << std::endl;
				std::cout << std::endl;
			}*/
		}
	}

	Boundary_Conditions(u_s, par, Dir, par.d_t * (par.N_step + 1));

	return result;
}

double advective_term(Matrix &ul, Matrix &vl, Matrix &ur, Matrix &vr, size_t i, size_t j, double d_x, double d_y, Direction Dir) {
	double v_help = 0.25 * (vl[i][j] + L(vl, i, j, Dir) + U(vl, i, j, Dir) + UL(vl, i, j, Dir));
	double result = ul[i][j] * (R(ur, i, j, Dir) - L(ur, i, j, Dir)) / (2.0*d_x)
	                + v_help * (U(ur, i, j, Dir) - D(ur, i, j, Dir)) / (2.0*d_y);
	return result;
}

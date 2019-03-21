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
Matrix CalculateB(Matrix &u_n, Matrix &v_n, Matrix &u_s, Matrix &v_s, Matrix &p, Matrix &p_new, Matrix &F, Param &par, Direction Dir) {

	size_t Nx = u_n.size();
	size_t Ny = u_n[0].size();

	double d_u, d_v;

	if      (Dir == Du) { d_u = par.d_x;	d_v = par.d_y;}
	else if (Dir == Dv) { d_u = par.d_y;	d_v = par.d_x;}

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
			double Laplace_F = par.d_t / (2.0*par.Re) * (par.ldxdx * (F[i + 1][j] - 2.0 * F[i][j] + F[i - 1][j])
			                                           + par.ldydy * (F[i][j + 1] - 2.0 * F[i][j] + F[i][j - 1]) ) ;

			result[i][j] = - advective_term_s
			               - pressure_term
			               + diffusion_term_n / (2.0*par.Re)
			               + Laplace_F
			               + u_n[i][j] / par.d_t;

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

void Output_eq_terms(std::string filename, int n, Matrix &u_n, Matrix &v_n, Matrix &u_s, Matrix &v_s, Matrix &p, Matrix &p_new, Matrix &F, Param par, Direction Dir) {

	size_t Nx = u_n.size();
	size_t Ny = u_n[0].size();

	double d_u, d_v;

	if      (Dir == Du) { d_u = par.d_x;	d_v = par.d_y; }
	else if (Dir == Dv) { d_u = par.d_y;	d_v = par.d_x; }


	std::ofstream output;
	filename = par.WorkDir + filename + std::to_string(n) + ".plt";

	output.open(filename);

	output << "title = " << '"' << filename << '"' << std::endl;
	output << "Variables = i j advection diffusion pressure Laplace_F du_dt" << std::endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << Nx-2 << ", j=" << Ny-2 << ", f=point" << std::endl;
	output << "SolutionTime = " << n << std::endl;

	for (size_t j = 1; j < (Ny - 1); ++j) {
		for (size_t i = 1; i < (Nx - 1); ++i) {
			double advective_term_s = 0.25 * (advective_term(u_n, v_n, u_n, v_n, i, j, d_u, d_v, Dir)
			                                + advective_term(u_s, v_s, u_s, v_s, i, j, d_u, d_v, Dir)
			                                + advective_term(u_n, v_n, u_s, v_s, i, j, d_u, d_v, Dir)
			                                + advective_term(u_s, v_s, u_n, v_n, i, j, d_u, d_v, Dir));
			double diffusion_term_n = par.ldxdx * (u_n[i + 1][j] - 2.0 * u_n[i][j] + u_n[i - 1][j])
			                        + par.ldydy * (u_n[i][j + 1] - 2.0 * u_n[i][j] + u_n[i][j - 1]);
			double pressure_term = 0.5 * (
			                               (p    [i][j] - L(p    , i, j, Dir)) / d_u
			                             + (p_new[i][j] - L(p_new, i, j, Dir)) / d_u
			                             );
			double Laplace_F =  par.ldxdx * (F[i + 1][j] - 2.0 * F[i][j] + F[i - 1][j])
			                  + par.ldydy * (F[i][j + 1] - 2.0 * F[i][j] + F[i][j - 1]);

			output << i << ' '
			       << j << ' '
				   << advective_term_s << ' '
			       << diffusion_term_n / (2.0*par.Re) << ' '
			       << pressure_term << ' '
			       << Laplace_F * par.d_t / (2.0*par.Re) << ' '
			       << (u_s[i][j] - u_n[i][j]) / par.d_t << ' '
			       << std::endl;
		}
	}

}

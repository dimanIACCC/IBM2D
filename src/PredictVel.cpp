#include "stdafx.h"
#include "PredictVel.h"
#include "helmholtz.h"

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

			if (Dir == Du) result[i][j] += par.grad_p_x;

		}
	}

	Boundary_Conditions(u_s, par, Dir, par.d_t * (par.N_step + 1));
	Boundary_Conditions(result, par, Dir, par.d_t * (par.N_step + 1));

	return result;
}

double advective_term(Matrix &ul, Matrix &vl, Matrix &ur, Matrix &vr, size_t i, size_t j, double d_x, double d_y, Direction Dir) {
	double v_help = 0.25 * (vl[i][j] + L(vl, i, j, Dir) + U(vl, i, j, Dir) + UL(vl, i, j, Dir));
	double result = ul[i][j] * (R(ur, i, j, Dir) - L(ur, i, j, Dir)) / (2.0*d_x)
	                + v_help * (U(ur, i, j, Dir) - D(ur, i, j, Dir)) / (2.0*d_y);
	return result;
}

void make_uv_RHS(Matrix &rhsu, Matrix &rhsv, Matrix &u0, Matrix &v0, Matrix &u1, Matrix &v1, Matrix &p0, double fu, double fv,
                 int nx, int ny, double hx, double hy, double tau, double Re) {

	double re_x = Re / hx;
	double re_y = Re / hy;
	double e_xx = 1. / (hx*hx);
	double e_yy = 1. / (hy*hy);
	double re2 = 2.*Re;
	double re2_t = 2.*Re / tau;

	for (int i = 1; i <= nx + 1; i++) {
		for (int j = 1; j <= ny; j++) {
			rhsu[i][j] = re2_t*u0[i][j]
				+ e_xx *(u0[i + 1][j] - 2.*u0[i][j] + u0[i - 1][j])
				+ e_yy *(u0[i][j + 1] - 2.*u0[i][j] + u0[i][j - 1])
				- re_x *(p0[i][j] - p0[i - 1][j])
				+ re2  *fu;
		}
	}
	for (int i = 1; i <= nx; i++) {
		for (int j = 1; j <= ny + 1; j++) {
			rhsv[i][j] = re2_t*v0[i][j] 
				+ e_xx *(v0[i + 1][j] - 2.*v0[i][j] + v0[i - 1][j])
				+ e_yy *(v0[i][j + 1] - 2.*v0[i][j] + v0[i][j - 1])
				- re_y *(p0[i][j] - p0[i][j - 1])
				+ re2  *fv;
		}
	}

}

void predict_uv(Matrix &u0, Matrix &v0, Matrix &u, Matrix &v,
	Matrix &p, Matrix &rhsu, Matrix &rhsv, int &nx, int &ny, double &hx, double &hy, double &tau, double &Re, Param &par) {

	double re2_t = 2.0 * Re / tau;
	double re_x  = Re / hx;
	double re_y  = Re / hy;
	double re_2x = Re / (2.0*hx);
	double re_2y = Re / (2.0*hy);

	//Allocate arrayes
	CreateMatrix(fu, nx + 3, ny + 2);
	CreateMatrix(fv, nx + 2, ny + 3);

	for (int i = 1; i <= nx + 1; i++) {
		for (int j = 1; j <= ny; j++) {
			double a1 = 0.5*(u[i][j] + u0[i][j]);
			double a2 = 0.125*(v[i - 1][j] + v[i - 1][j + 1] + v[i][j] + v[i][j + 1] + v0[i - 1][j] + v0[i - 1][j + 1] + v0[i][j] + v0[i][j + 1]);
			fu[i][j] = rhsu[i][j]
			- re_2x*a1*(u[i + 1][j] - u[i - 1][j] + u0[i + 1][j] - u0[i - 1][j])
			- re_2y*a2*(u[i][j + 1] - u[i][j - 1] + u0[i][j + 1] - u0[i][j - 1])
			- re_x *(p[i][j] - p[i - 1][j]);
		}
	}
	if (par.BC == periodical)
		for (int j = 1; j <= ny; j++) {
			fu[0   ][j] = fu[nx  ][j];
			fu[1   ][j] = fu[nx+1][j];
			fu[nx+2][j] = fu[2   ][j];
		}

	for (int i = 1; i <= nx; i++) {
		for (int j = 1; j <= ny + 1; j++) {
			double a1 = 0.125*(u[i][j - 1] + u[i + 1][j - 1] + u[i][j] + u[i + 1][j] + u0[i][j - 1] + u0[i + 1][j - 1] + u0[i][j] + u0[i + 1][j]);
			double a2 = 0.5*(v[i][j] + v0[i][j]);
			fv[i][j] = rhsv[i][j]
			- re_2x*a1*(v[i + 1][j] - v[i - 1][j] + v0[i + 1][j] - v0[i - 1][j])
			- re_2y*a2*(v[i][j + 1] - v[i][j - 1] + v0[i][j + 1] - v0[i][j - 1])
			- re_y *(p[i][j] - p[i][j - 1]);
		}
	}
	if (par.BC == periodical)
		for (int j = 1; j <= ny + 1; j++) {
			fv[0][j] = fv[nx][j];
			fv[nx+1][j] = fv[1][j];
		}
	//Output_U(fu, "fu", -555, par);
	//Output_V(fv, "fv", -555, par);
	//std::cin.get();

	CreateMatrix(ax_u, ny + 2, 4);
	CreateMatrix(bx_u, ny + 2, 4);
	CreateMatrix(ay_u, nx + 3, 4);
	CreateMatrix(by_u, nx + 3, 4);

	//Set up boundary conditions
	for (int i = 0; i <= ny + 1; i++) {
		ax_u[i][1] =  0.;   //Inflow boundary condition : du / dn = 0 (soft)
		ax_u[i][2] = -1.;
		ax_u[i][3] =  0.;
		bx_u[i][1] =  0.;   //Outflow boundary condition : du / dn = 0 (soft)
		bx_u[i][2] = -1.;
		bx_u[i][3] =  0.;
	}
	for (int i = 0; i <= nx + 2; i++) {
		ay_u[i][1] =  1.; //Bottom boundary condition : du / dn = 0 (soft)
		ay_u[i][2] =  0.;
		ay_u[i][3] =  0.;
		by_u[i][1] =  1.; //Top boundary condition : du / dn = 0 (soft)
		by_u[i][2] =  0.;
		by_u[i][3] =  0.;
	}

	char* BCtype = "DDDD";
	if (par.BC == u_inflow  ) BCtype = "NDDD";
	if (par.BC == u_infinity) BCtype = "NDNN";
	if (par.BC == periodical) BCtype = "PPDD";
	Helmholtz_MKL(u, fu, re2_t, par.d_x, par.d_y, 1, nx + 1, 0, ny + 1, BCtype);

	//Solve Boundary value problem for Helmholtz equation by Gauss–Seidel method(it is faster here than SOR method)
	//Helmholtz_SOR(u, ax_u, bx_u, ay_u, by_u, fu, re2_t, nx + 1, ny, hx, hy, 200000, 1.e-11, 1.);
	Boundary_Conditions(u, par, Du, par.d_t * (par.N_step + 1));


	CreateMatrix(ax_v, ny + 3, 4);
	CreateMatrix(bx_v, ny + 3, 4);
	CreateMatrix(ay_v, nx + 2, 4);
	CreateMatrix(by_v, nx + 2, 4);

	for (int i = 0; i <= ny + 2; i++) {
		ax_v[i][1] = -1.;   //Inflow boundary condition : v = 0
		ax_v[i][2] =  0.;
		ax_v[i][3] =  0.;
		bx_v[i][1] = -1.;   //Outflow boundary condition : dv / dn = 0 (soft)
		bx_v[i][2] =  0.;
		bx_v[i][3] =  0.;
	}
	for (int i = 0; i <= nx + 1; i++) {
		ay_v[i][1] =  0.;   //Bottom boundary condition : v = 0 (no - slip)
		ay_v[i][2] =  1.;
		ay_v[i][3] =  0.;
		by_v[i][1] =  0.;   //Top    boundary condition : v = 0 (no - slip)
		by_v[i][2] =  1.;
		by_v[i][3] =  0.;
	}

	if (par.BC == u_inflow  ) BCtype = "NDDD";
	if (par.BC == u_infinity) BCtype = "NDDD";
	if (par.BC == periodical) BCtype = "PPDD";
	Helmholtz_MKL(v, fv, re2_t, par.d_x, par.d_y, 0, nx, 0, ny + 2, BCtype);

	//Solve Boundary value problem for Helmholtz equation by Gauss–Seidel method(it is faster here than SOR method)
	//Helmholtz_SOR(v, ax_v, bx_v, ay_v, by_v, fv, re2_t, nx, ny + 1, hx, hy, 200000, 1.e-11, 1.);
	Boundary_Conditions(v, par, Dv, par.d_t * (par.N_step + 1));
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

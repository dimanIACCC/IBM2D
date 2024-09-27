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
			double Laplace_F =  par.ldxdx * (F[i + 1][j] - 2.0 * F[i][j] + F[i - 1][j])
			                  + par.ldydy * (F[i][j + 1] - 2.0 * F[i][j] + F[i][j - 1]) ;

			result[i][j] = - advective_term_s
			               - pressure_term
			               + diffusion_term_n    / (2.0*par.Re)
			               + par.d_t * Laplace_F / (2.0*par.Re)
			               + u_n[i][j] / par.d_t;

			if (Dir == Du) result[i][j] += par.grad_p_x;

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

void make_uv_RHS(Matrix &rhsu, Matrix &rhsv, Matrix &u0, Matrix &v0, Matrix &u, Matrix &v, Matrix &p, Matrix &Fx, Matrix &Fy, Param par) {

	double re2_x = 2. * par.Re / par.d_x;
	double re2_y = 2. * par.Re / par.d_y;
	double re_2x = par.Re / (2.0*par.d_x);
	double re_2y = par.Re / (2.0*par.d_y);
	double e_xx = 1. / (par.d_x*par.d_x);
	double e_yy = 1. / (par.d_y*par.d_y);
	double re2 = 2.*par.Re;
	double re2_t = 2.*par.Re / par.d_t;

	for (int i = 1; i <= par.N1_u - 2; i++) {
		for (int j = 1; j <= par.N2_u - 2; j++) {
			double a1 = 0.5*(u[i][j] + u0[i][j]);
			double a2 = 0.125*(v[i - 1][j] + v[i - 1][j + 1] + v[i][j] + v[i][j + 1] + v0[i - 1][j] + v0[i - 1][j + 1] + v0[i][j] + v0[i][j + 1]);
			double Laplace_F = par.ldxdx * (Fx[i + 1][j] - 2.0 * Fx[i][j] + Fx[i - 1][j])
				             + par.ldydy * (Fx[i][j + 1] - 2.0 * Fx[i][j] + Fx[i][j - 1]);

			rhsu[i][j] = re2_t*u0[i][j]
				+ e_xx *(u0[i + 1][j] - 2.*u0[i][j] + u0[i - 1][j])
				+ e_yy *(u0[i][j + 1] - 2.*u0[i][j] + u0[i][j - 1])

				- re_2x*a1*(u[i + 1][j] - u[i - 1][j] + u0[i + 1][j] - u0[i - 1][j])
				- re_2y*a2*(u[i][j + 1] - u[i][j - 1] + u0[i][j + 1] - u0[i][j - 1])

				- re2_x *(p [i][j] - p [i - 1][j])
				
				+ Fx[i][j] * re2
				// + Laplace_F * par.d_t
				+ par.grad_p_x * re2;
		}
	}
	for (int i = 1; i <= par.N1_v - 2; i++) {
		for (int j = 1; j <= par.N2_v - 2; j++) {
			double a1 = 0.125*(u[i][j - 1] + u[i + 1][j - 1] + u[i][j] + u[i + 1][j] + u0[i][j - 1] + u0[i + 1][j - 1] + u0[i][j] + u0[i + 1][j]);
			double a2 = 0.5*(v[i][j] + v0[i][j]);
			double Laplace_F = par.ldxdx * (Fy[i + 1][j] - 2.0 * Fy[i][j] + Fy[i - 1][j])
				             + par.ldydy * (Fy[i][j + 1] - 2.0 * Fy[i][j] + Fy[i][j - 1]);

			rhsv[i][j] = re2_t*v0[i][j]
				+ e_xx *(v0[i + 1][j] - 2.*v0[i][j] + v0[i - 1][j])
				+ e_yy *(v0[i][j + 1] - 2.*v0[i][j] + v0[i][j - 1])
				
				- re_2x*a1*(v[i + 1][j] - v[i - 1][j] + v0[i + 1][j] - v0[i - 1][j])
				- re_2y*a2*(v[i][j + 1] - v[i][j - 1] + v0[i][j + 1] - v0[i][j - 1])

				- re2_y *(p [i][j] - p [i][j - 1])

				+ Fy[i][j] * re2;
			    //+ Laplace_F * par.d_t;
		}
	}
	Boundary_Conditions(rhsu, par, Du, par.d_t * (par.N_step + 1));
	Boundary_Conditions(rhsv, par, Dv, par.d_t * (par.N_step + 1));
}

void predict_uv(Matrix &u, Matrix &v, Matrix &rhsu, Matrix &rhsv, Param &par) {

	#pragma omp parallel sections num_threads(1)
	{
      #pragma omp section
	  {
		//BiCGStab(u, A_u, B_u, par, Du, N_BiCGStab_u);                   // solving A_u * U_new = B_u
		prepare_solve_helmholtz_velocity(u, rhsu, 2 * par.Re / par.d_t, par, Du);
		Boundary_Conditions(u, par, Du, par.d_t * (par.N_step + 1));
	  }
	  #pragma omp section
	  {
		//BiCGStab(v, A_v, B_v, par, Dv, N_BiCGStab_v);                   // solving A_v * V_new = B_v
		prepare_solve_helmholtz_velocity(v, rhsv, 2 * par.Re / par.d_t, par, Dv);
		Boundary_Conditions(v, par, Dv, par.d_t * (par.N_step + 1));
	  }
	}

}

void prepare_solve_helmholtz_velocity(Matrix &A, Matrix &RHS, double q, Param par, Direction Dir) {
	char* BCtype = "DDDD";
	if (par.BC == u_inflow) BCtype = "DNDD";
	if (par.BC == u_infinity) BCtype = "DNNN";
	if (par.BC == periodical) BCtype = "PPDD";

	MKL_INT nx, ny;
	if (Dir == Du) {
		nx = par.N1 + 2;
		ny = par.N2 + 1;
		if (par.BC == periodical) nx = par.N1;
	}
	if (Dir == Dv) {
		nx = par.N1 + 1;
		ny = par.N2 + 2;
		if (par.BC == periodical) nx = par.N1;
	}

	double ax = 0.;
	double bx = par.d_x*nx;
	double ay = 0.;
	double by = par.d_y*ny;

	double *f_mkl = NULL, *bd_ax = NULL, *bd_bx = NULL, *bd_ay = NULL, *bd_by = NULL;
	f_mkl = (double*)mkl_malloc((nx + 1)*(ny + 1) * sizeof(double), 64);
	bd_ax = (double*)mkl_malloc((ny + 1) * sizeof(double), 64);
	bd_bx = (double*)mkl_malloc((ny + 1) * sizeof(double), 64);
	bd_ay = (double*)mkl_malloc((nx + 1) * sizeof(double), 64);
	bd_by = (double*)mkl_malloc((nx + 1) * sizeof(double), 64);

	if (Dir == Du) MatrixU_to_DoubleArray(RHS, f_mkl, par.BC);
	if (Dir == Dv) MatrixV_to_DoubleArray(RHS, f_mkl, par.BC);

	for (MKL_INT iy = 0; iy <= ny; iy++) {
		bd_ax[iy] = 0.;
		bd_bx[iy] = 0.;
		if (par.BC == u_infinity || par.BC == u_inflow) bd_ax[iy] = A[0][iy];
		if (par.BC == Taylor_Green || par.BC == Lamb_Oseen || par.BC == Line_Vortex) {
			bd_ax[iy] = A[0][iy];
			bd_bx[iy] = A[nx][iy];
		}
	}
	for (MKL_INT ix = 0; ix <= nx; ix++) {
		bd_ay[ix] = 0.;
		bd_by[ix] = 0.;
		if (par.BC == periodical || par.BC == u_inflow || par.BC == u_infinity) {
			if (Dir == Du) bd_ay[ix] = par.u_down;
			if (Dir == Du) bd_by[ix] = par.u_up;
		}
		if (par.BC == Taylor_Green || par.BC == Lamb_Oseen || par.BC == Line_Vortex) {
			bd_ay[ix] = A[ix][0];
			bd_by[ix] = A[ix][ny];
		}
	}

	Helmholtz_MKL(f_mkl, ax, bx, ay, by, bd_ax, bd_bx, bd_ay, bd_by, nx, ny, BCtype, 2 * par.Re / par.d_t, par.d_x, par.d_y);

	if (Dir == Du) DoubleArray_to_MatrixU(f_mkl, A, par.BC);
	if (Dir == Dv) DoubleArray_to_MatrixV(f_mkl, A, par.BC);

	mkl_free(f_mkl);
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

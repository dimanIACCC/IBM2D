#include "PredictVel.h"



void Calculate_A(ublas::matrix<Template> &A, Param par, double Re, Direction Dir) {

	size_t n1 = A.size1();
	size_t n2 = A.size2();
	size_t N;

	double d_u, d_v;

	if      (Dir == Du) { d_u = par.d_x;	d_v = par.d_y;	N = n2; }
	else if (Dir == Dv) { d_u = par.d_y;	d_v = par.d_x;	N = n1; }

	double d_uu = 1.0 / (d_u*d_u);
	double d_vv = 1.0 / (d_v*d_v);

	for (size_t j = 0; j < n2; ++j) {
		for (size_t i = 0; i < n1; ++i) {

			A(i, j).C =  1.0 /      Re *(d_uu + d_vv) + 1.0 / par.d_t;
			A(i, j).U = -1.0 / (2.0*Re)* d_vv;
			A(i, j).R = -1.0 / (2.0*Re)* d_uu;
			A(i, j).D = -1.0 / (2.0*Re)* d_vv;
			A(i, j).L = -1.0 / (2.0*Re)* d_uu;

		}
	}
}

Matrix Operator_Ax(ublas::matrix<Template> &A, Matrix &u, Param par, Direction Dir) {

	size_t Nx = u.size();
	size_t Ny = u[0].size();

	CreateMatrix(result, Nx, Ny);

	for (size_t j = 0; j < Ny; ++j) {
		for (size_t i = 0; i < Nx; ++i) {
			result[i][j] = A(i, j).C * u[i][j]
			             + A(i, j).U * U(u, i, j, Dir, Nx, Ny)
			             + A(i, j).R * R(u, i, j, Dir, Nx, Ny)
			             + A(i, j).D * D(u, i, j, Dir, Nx, Ny)
			             + A(i, j).L * L(u, i, j, Dir, Nx, Ny);
		}
	}

	if (par.BC == periodical) {
		if (Dir == Du) {
			for (size_t j = 1; j < Ny; ++j) {
				result[Nx - 1][j] = result[0][j];
			}
		}
		else if (Dir == Dv) {
			for (size_t j = 0; j < Ny; ++j) {
				result[0][j] = result[Nx - 2][j];
				result[Nx - 1][j] = result[1][j];
			}
		}
	}

	// Up-Down BC
	for (size_t i = 0; i < Nx; ++i) {

		size_t j = 0;
		result[i][j] = u[i][j];
		if ((par.BC == u_infinity) && (Dir == Du)) result[i][j] = (u[i][j + 1] - u[i][j]) / par.d_y;

		j = Ny - 1;
		result[i][j] = u[i][j];
		if ((par.BC == u_infinity) && (Dir == Du)) result[i][j] = (u[i][j] - u[i][j - 1]) / par.d_y;

	}

	for (size_t j = 0; j < Ny; ++j) {
		switch (par.BC) {
			case u_infinity:	result[0][j] = u[0][j];		result[Nx - 1][j] = (3.0 * u[Nx - 1][j] - 4.0 * u[Nx - 2][j] + 1.0 * u[Nx - 3][j]) / (2.0*par.d_x);		break;
			case u_inflow  :	result[0][j] = u[0][j];		result[Nx - 1][j] = (3.0 * u[Nx - 1][j] - 4.0 * u[Nx - 2][j] + 1.0 * u[Nx - 3][j]) / (2.0*par.d_x);		break;
		}
	}
	
	return result;

}

Matrix CalculateB(Matrix &u_n, Matrix &v_n, Matrix &u_s, Matrix &v_s, Matrix &p, Matrix &force, Param par, Direction Dir) {

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
			double advective_term_n    = advective_term(u_n, v_n, i, j, d_u, d_v, Dir, Nx, Ny);
			double advective_term_s    = advective_term(u_s, v_s, i, j, d_u, d_v, Dir, Nx, Ny);
			double diffusion_term_n    = diffusion_term(u_n     , i, j, d_uu, d_vv, Dir, Nx, Ny);
			double pressure_term = (R(p, i, j, Dir, par.N1 + 1, par.N2 + 1) - p[i][j]) / d_u;

			result[i][j] = -(       alpha  * advective_term_n
			               + (1.0 - alpha) * advective_term_s)
			                               - pressure_term
			                               + diffusion_term_n / (2.0*par.Re)
			                               + u_n[i][j] / par.d_t;
			                               //+ force[i][j];
		}
	}

	if (par.BC == periodical) {
		if (Dir == Du) {
			for (size_t j = 1; j < Ny; ++j) {
				size_t i = 0;
				double advective_term_n = advective_term(u_n, v_n, i, j, d_u, d_v, Dir, par.N1, par.N2);
				double advective_term_s = advective_term(u_s, v_s, i, j, d_u, d_v, Dir, par.N1, par.N2);
				double diffusion_term_n = diffusion_term(u_n, i, j, d_uu, d_vv, Dir, par.N1, par.N2 + 1);
				double pressure_term = (p[1][j] - p[par.N1 - 1][j] - par.L * dpdx_Poiseuille(par.H, par.Re)) / d_u;
				result[i][j] = -(alpha * advective_term_n
				       + (1.0 - alpha) * advective_term_s)
				                       - pressure_term
				                       + diffusion_term_n / (2.0*par.Re)
				                       + u_n[i][j] / par.d_t;
				                       //+ force[i][j];
				result[Nx - 1][j] = result[0][j];
			}
		}
		else if (Dir == Dv) {
			for (size_t j = 0; j < Ny; ++j) {
				result[0]     [j] = result[Nx - 2][j];
				result[Nx - 1][j] = result[1][j];
			}
		}
	}


	// Up-Down BC
	for (size_t i = 0; i < Nx; ++i) {
		if      (Dir == Du) {
			result[i][0] = -u_n[i][1];
			result[i][Ny - 1] = -u_n[i][Ny - 2];
			if (par.BC == u_infinity) {
				result[i][0] = 0;
				result[i][Ny - 1] = 0;
			}
		}
		else if (Dir == Dv) {
			result[i][0] = 0;
			result[i][Ny - 1] = 0;
		}
	}

	for (size_t j = 0; j < Ny; ++j) {
		switch (par.BC) {
			case u_infinity:	result[0][j] = u_n[0][j]     ;		result[Nx - 1][j] = 0;		break;
			case u_inflow  :	result[0][j] = u_n[0][j]     ;		result[Nx - 1][j] = 0;		break;
		}
	}

	return result;
}

double advective_term(Matrix &u, Matrix &v, size_t i, size_t j, double d_x, double d_y, Direction Dir, size_t N1, size_t N2) {
	double v_help = 0.25 * (v[i][j] + R(v, i, j, Dir, N1, N2) + D(v, i, j, Dir, N1, N2) + RD(v, i, j, Dir, N1, N2));
	double result = u[i][j] * (R(u, i, j, Dir, N1, N2) - L(u, i, j, Dir, N1, N2)) / (2.0*d_x)
	               + v_help * (U(u, i, j, Dir, N1, N2) - D(u, i, j, Dir, N1, N2)) / (2.0*d_y);
	return result;
}

double diffusion_term(Matrix &u, size_t i, size_t j, double d_xx, double d_yy, Direction Dir, size_t N1, size_t N2) {
	double result = d_xx * (R(u, i, j, Dir, N1, N2) - 2.0 * u[i][j] + L(u, i, j, Dir, N1, N2))
	              + d_yy * (U(u, i, j, Dir, N1, N2) - 2.0 * u[i][j] + D(u, i, j, Dir, N1, N2));
	return result;
}

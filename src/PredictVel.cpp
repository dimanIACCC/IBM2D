#include "PredictVel.h"



void Calculate_A(ublas::matrix<Template> &A, Param par, double Re, Direction Dir) {

	size_t n1 = A.size1();
	size_t n2 = A.size2();
	size_t N;

	double d_u, d_v;
	int ij;

	if      (Dir == Du) { d_u = par.d_x;	d_v = par.d_y;	N = n2; }
	else if (Dir == Dv) { d_u = par.d_y;	d_v = par.d_x;	N = n1; }

	double d_uu = 1.0 / (d_u*d_u);
	double d_vv = 1.0 / (d_v*d_v);

	for (int j = 1; j < (n2 - 1); ++j) {
		for (int i = 1; i < (n1 - 1); ++i) {

			A(i, j).C =  1.0 /      Re *(d_uu + d_vv) + 1.0 / par.d_t;
			A(i, j).U = -1.0 / (2.0*Re)* d_vv;
			A(i, j).R = -1.0 / (2.0*Re)* d_uu;
			A(i, j).D = -1.0 / (2.0*Re)* d_vv;
			A(i, j).L = -1.0 / (2.0*Re)* d_uu;

			if      (Dir == Du)  ij = j;
			else if (Dir == Dv)  ij = i;

			if (ij == 1) {
				A(i, j).C = 1.0 / Re* (d_uu + 2.0 * d_vv) + 1.0 / par.d_t;
				A(i, j).U = -2.0 / (3.0*Re)* d_vv;
				A(i, j).R = -1.0 / (2.0*Re)* d_uu;
				A(i, j).D = -4.0 / (3.0*Re)* d_vv;
				A(i, j).L = -1.0 / (2.0*Re)* d_uu;
			}

			if (ij == N - 2) {
				A(i, j).C = 1.0 / Re* (d_uu + 2.0 * d_vv) + 1.0 / par.d_t;
				A(i, j).U = -4.0 / (3.0*Re)* d_vv;
				A(i, j).R = -1.0 / (2.0*Re)* d_uu;
				A(i, j).D = -2.0 / (3.0*Re)* d_vv;
				A(i, j).L = -1.0 / (2.0*Re)* d_uu;
			}

		}
	}
}

Matrix Operator_Ax(ublas::matrix<Template> &A, Matrix &u, Param par, Direction Dir) {

	size_t Nx = u.size();
	size_t Ny = u[0].size();

	CreateMatrix(result, Nx, Ny);

	for (size_t j = 1; j < (Ny - 1); ++j) {
		for (size_t i = 1; i < (Nx - 1); ++i) {
			result[i][j] = A(i, j).C * u[i][j]
			             + A(i, j).U * U(u, i, j, Dir)
			             + A(i, j).R * R(u, i, j, Dir)
			             + A(i, j).D * D(u, i, j, Dir)
			             + A(i, j).L * L(u, i, j, Dir);
		}
	}

	// Up-Down BC
	for (size_t i = 1; i < Nx - 1; ++i) {

		size_t j = 0;
		result[i][j] = u[i][j];
		if ((par.BC == u_infinity) && (Dir == Du)) result[i][j] = (u[i][j + 1] - u[i][j]) / (par.d_y*0.5);

		j = Ny - 1;
		result[i][j] = u[i][j];
		if ((par.BC == u_infinity) && (Dir == Du)) result[i][j] = (u[i][j] - u[i][j - 1]) / (par.d_y*0.5);

	}

	for (size_t j = 0; j < Ny; ++j) {
		result[0][j] = u[0][j];                                                                                 // inflow u = u0
		result[Nx - 1][j] = (3.0 * u[Nx - 1][j] - 4.0 * u[Nx - 2][j] + 1.0 * u[Nx - 3][j]) / (2.0*par.d_x); 	// outflow du/dx = 0
	}
	
	return result;

}

Matrix CalculateB(Matrix &u_n, Matrix &v_n, Matrix &u_prev, Matrix &v_prev, Matrix &p, Matrix &force, Param par, Direction Dir) {

	size_t Nx = u_n.size();
	size_t Ny = u_n[0].size();
	size_t N;

	double d_u, d_v;
	size_t ij;

	if      (Dir == Du) { d_u = par.d_x;	d_v = par.d_y;	N = Ny;}
	else if (Dir == Dv) { d_u = par.d_y;	d_v = par.d_x;	N = Nx;}

	double d_uu = 1.0 / (d_u*d_u);
	double d_vv = 1.0 / (d_v*d_v);


	double k, kp, k0, km;

	CreateMatrix(result, Nx, Ny);

	for (size_t j = 1; j < (Ny - 1); ++j) {
		for (size_t i = 1; i < (Nx - 1); ++i) {
			if      (Dir == Du)  ij = j;
			else if (Dir == Dv)  ij = i;
			B_coefficients(ij, N, k, kp, k0, km);
			double advective_term_n    = advective_term(u_n   , v_n   , i, j, d_u, d_v, k, Dir);
			double advective_term_prev = advective_term(u_prev, v_prev, i, j, d_u, d_v, k, Dir);
			double diffusion_term_n    = diffusion_term(u_n           , i, j, d_uu, d_vv, kp, k0, km, Dir);
			double pressure_term       = (R(p, i, j, Dir) - p[i][j]) / d_u;
			result[i][j] = -(3.0 / 2.0 * advective_term_n
			               - 1.0 / 2.0 * advective_term_prev)
			                           - pressure_term
			                           + diffusion_term_n / (2.0*par.Re)
			                           + u_n[i][j] / par.d_t
			                           + force[i][j];
		}
	}

	// outflow du/dx = 0
	for (size_t j = 0; j < Ny; ++j) {
		result[Nx - 1][j] = 0.0;
		result[0][j] = u_n[0][j];
		if (par.BC == periodical) result[0][j] = u_n[Nx-1][j];
	}

	return result;
}

void B_coefficients(size_t ij, size_t N, double &k, double &kp, double &k0, double &km) {
	if (ij > 1 && ij < N - 2) { k = 2. ;	kp = 1.      ;	k0 = -2.;	km = 1.; }
	else {
		                        k = 1.5;
		if      (ij == 1)                 { kp = 4. / 3.;	k0 = -4.;	km = 8. / 3.; }
		else if (ij == N - 2)             { kp = 8. / 3.;	k0 = -4.;	km = 4. / 3.; }
	}
}

double advective_term(Matrix &u, Matrix &v, size_t i, size_t j, double d_x, double d_y, double k, Direction Dir) {
	double v_help = 0.25 * (v[i][j] + R(v, i, j, Dir) + D(v, i, j, Dir) + RD(v, i, j, Dir));
	double result = u[i][j] * (R(u, i, j, Dir) - L(u, i, j, Dir)) / (2.0*d_x)
	               + v_help * (U(u, i, j, Dir) - D(u, i, j, Dir)) / (k  *d_y);
	return result;
}

double diffusion_term(Matrix &u, size_t i, size_t j, double d_xx, double d_yy, double kp, double k0, double km, Direction Dir) {
	double result = d_xx * (     R(u, i, j, Dir) - 2.0 * u[i][j] +      L(u, i, j, Dir))
	              + d_yy * (kp * U(u, i, j, Dir) + k0  * u[i][j] + km * D(u, i, j, Dir));
	return result;
}

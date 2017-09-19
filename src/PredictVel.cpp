#include "PredictVel.h"


void ExplicPredVel(Matrix& U_predict, Matrix& V_predict, Matrix& U_n, Matrix& V_n, Matrix& P, Matrix& Force_x, Matrix& Force_y, Param par) {
	
	for (int i = 0; i <(int)U_predict.size(); i++) {
		U_predict[i][0]                    = U_n[i][0]; 
		U_predict[i][U_predict[0].size() - 1] = U_n[i][U_predict[0].size() -1];
	}
	//for (int i = 0; i < U_predict[0].size(); i++) { U_predict[0][i] = U_n[0][i]; U_predict[U_predict[0].size() -1 ][i] = U_n[U_predict[0].size() -1][i];}
	for (int i = 0; i < (int)V_predict.size(); i++) {
		V_predict[i][0] = V_n[i][0];
		V_predict[i][V_predict[0].size() - 1] = V_n[i][V_predict[0].size() - 1];
	}
	//for (int i = 0; i < U_predict[0].size(); i++);
	

	for (int i = 1; i < (int)U_predict.size()-1; i++) {
		for (int j = 1; j < (int)U_predict[0].size()-1; j++)
		{
			double LaplasU = (U_n[i + 1][j] - 2 * U_n[i][j] + U_n[i - 1][j]) / pow(par.d_x, 2) + (U_n[i][j + 1] - 2 * U_n[i][j] + U_n[i][j - 1]) / pow(par.d_y, 2);
			double CentrDiffx_U = (U_n[i + 1][j] - U_n[i - 1][j]) / (2 * par.d_x);
			double CentrDiffy_U = (U_n[i][j + 1] - U_n[i][j - 1]) / (2 * par.d_y);
			double GradPressX   = (P[i + 1][j] - P[i - 1][j]) / (2 * par.d_x);
			U_predict[i][j] = par.d_t*(LaplasU - U_n[i][j] * CentrDiffx_U - V_n[i][j] * CentrDiffy_U) + U_n[i][j] - GradPressX +Force_x[i][j];
		}
	}
	for (int i = 1; i < (int)V_predict.size()-1; i++) {
		for (int j = 1; j < (int)V_predict[0].size()-1; j++)
		{
			double LaplasV = (V_n[i + 1][j] - 2 * V_n[i][j] + V_n[i - 1][j]) / pow(par.d_x, 2) + (V_n[i][j + 1] - 2 * V_n[i][j] + V_n[i][j - 1]) / pow(par.d_y, 2);
			double CentrDiffx_V = (V_n[i + 1][j] - V_n[i - 1][j]) / (2 * par.d_x);
			double CentrDiffy_V = (V_n[i][j + 1] - V_n[i][j - 1]) / (2 * par.d_y);
			double GradPressY   = (P[i][j + 1] - P[i][j - 1]) / (2 * par.d_y);
			V_predict[i][j] = par.d_t*(LaplasV - U_n[i][j] * CentrDiffx_V - V_n[i][j] * CentrDiffy_V) + V_n[i][j] - GradPressY +Force_y[i][j];
		}
	}



}


void Calculate_A_u(Matrix A[5], Param par, double Re) {

	double d_xx = 1.0 / (par.d_x*par.d_x);
	double d_yy = 1.0 / (par.d_y*par.d_y);
	int const n1 = par.N1;
	int const n2 = par.N2 + 1;

	for (int i = 0; i < 5; i++) {
		A[i].resize(n1);
		for (int j = 0; j < n1; j++) {
			A[i][j].resize(n2);
			fill(A[i][j].begin(), A[i][j].end(), 0);
		}
	}

	for (int j = 1; j < (n2 - 1); ++j) {
		for (int i = 1; i < (n1 - 1); ++i) {

			A[0][i][j] = 1.0 / par.d_t + (1.0 / Re) * (d_xx + d_yy);
			A[1][i][j] = -1.0 / (2.0*Re*par.d_y*par.d_y);
			A[2][i][j] = -1.0 / (2.0*Re*par.d_x*par.d_x);
			A[3][i][j] = -1.0 / (2.0*Re*par.d_y*par.d_y);
			A[4][i][j] = -1.0 / (2.0*Re*par.d_x*par.d_x);


			if (j == 1) {
				A[0][i][j] = 1.0 / par.d_t + (1.0 / Re) * (d_xx + 2.0*d_yy);
				A[1][i][j] = -2.0 / (3.0*Re*par.d_y*par.d_y);
				A[2][i][j] = -1.0 / (2.0*Re*par.d_x*par.d_x);
				A[3][i][j] = -4.0 / (3.0*Re*par.d_y*par.d_y);
				A[4][i][j] = -1.0 / (2.0*Re*par.d_x*par.d_x);
			}

			if (j == n2 - 2) {
				A[0][i][j] = 1.0 / par.d_t + (1.0 / Re) * (d_xx + 2.0*d_yy);
				A[1][i][j] = -4.0 / (3.0*Re*par.d_y*par.d_y);
				A[2][i][j] = -1.0 / (2.0*Re*par.d_x*par.d_x);
				A[3][i][j] = -2.0 / (3.0*Re*par.d_y*par.d_y);
				A[4][i][j] = -1.0 / (2.0*Re*par.d_x*par.d_x);
			}
		}
	}

	for (int i = 1; i < n1 - 1; ++i) {
		int j = 0;
		A[0][i][j] = 1.0;
		A[1][i][j] = 1.0;
		j = n2 - 1;
		A[0][i][j] = 1.0;
		A[1][i][j] = 1.0;


	}

	// outflow du/dx = 0
	for (int j = 0; j < n2; ++j) {
		A[0][n1 - 1][j] = 3.0 / (2.0*par.d_x);
		A[1][n1 - 1][j] = -4.0 / (2.0*par.d_x);
		A[2][n1 - 1][j] = 1.0 / (2.0*par.d_x);
		A[0][0][j] = 1.0;
	}



}

void Calculate_A_v(Matrix A[5], Param par, double Re) {

	double d_xx = 1.0 / (par.d_x*par.d_x);
	double d_yy = 1.0 / (par.d_y*par.d_y);
	int const n1 = par.N1 + 1;
	int const n2 = par.N2;

	for (int i = 0; i < 5; i++) {
		A[i].resize(n1);
		for (int j = 0; j < n1; j++) {
			A[i][j].resize(n2);
			fill(A[i][j].begin(), A[i][j].end(), 0);
		}
	}


	for (int j = 1; j < (n2 - 1); ++j) {
		for (int i = 1; i < (n1 - 1); ++i) {
			A[0][i][j] = 1.0 / par.d_t + (1.0 / Re) * (d_xx + d_yy);
			A[1][i][j] = -1.0 / (2.0*Re*par.d_y*par.d_y);
			A[2][i][j] = -1.0 / (2.0*Re*par.d_x*par.d_x);
			A[3][i][j] = -1.0 / (2.0*Re*par.d_y*par.d_y);
			A[4][i][j] = -1.0 / (2.0*Re*par.d_x*par.d_x);

			if (i == 1) {
				A[0][i][j] = 1.0 / par.d_t + (1.0 / Re) * (2.0*d_xx + d_yy);
				A[1][i][j] = -1.0 / (2.0*Re*par.d_y*par.d_y);
				A[2][i][j] = -2.0 / (3.0*Re*par.d_x*par.d_x);
				A[3][i][j] = -1.0 / (2.0*Re*par.d_y*par.d_y);
				A[4][i][j] = -4.0 / (3.0*Re*par.d_x*par.d_x);
			}

			if (i == n1 - 2) {
				A[0][i][j] = 1.0 / par.d_t + (1.0 / Re) * (2.0*d_xx + d_yy);
				A[1][i][j] = -1.0 / (2.0*Re*par.d_y*par.d_y);
				A[2][i][j] = -4.0 / (3.0*Re*par.d_x*par.d_x);
				A[3][i][j] = -1.0 / (2.0*Re*par.d_y*par.d_y);
				A[4][i][j] = -2.0 / (3.0*Re*par.d_x*par.d_x);
			}

		}
	}

	for (int i = 1; i < n1 - 1; ++i) {


		int j = 0;
		A[0][i][j] = -par.d_y * 0.5;
		A[1][i][j] = 0.0;
		j = n2 - 1;
		A[0][i][j] = par.d_y * 0.5;
		A[1][i][j] = 0.0;


	}

	// outflow du/dx = 0
	for (int j = 0; j < n2; ++j) {
		A[0][n1 - 1][j] = 3.0 / (2.0*par.d_x);
		A[1][n1 - 1][j] = -4.0 / (2.0*par.d_x);
		A[2][n1 - 1][j] = 1.0 / (2.0*par.d_x);
		A[0][0][j] = 1.0;
	}

}

Matrix Operator_Ax(Matrix A[5], Matrix &v, int const n1, int const n2, Param par,bool OverFlow) {

	CreateMatrix(result, n1, n2);


	for (int j = 1; j < (n2 - 1); ++j) {
		for (int i = 1; i < (n1 - 1); ++i) {

			result[i][j] = A[0][i][j] * v[i][j] + A[1][i][j] * v[i][j + 1] + A[2][i][j] * v[i + 1][j] + A[3][i][j] * v[i][j - 1] + A[4][i][j] * v[i - 1][j];
		}
	}
	if (OverFlow) {
		for (int i = 1; i < n1 - 1; ++i) {
			int j = 0;
			result[i][j] = (A[1][i][j] * v[i][j + 1] - A[0][i][j] * v[i][j]) / (par.d_y*0.5);
			j = n2 - 1;
			result[i][j] = (A[0][i][j] * v[i][j] - A[1][i][j] * v[i][j - 1]) / (par.d_y*0.5);
		}

	}
	else {
		for (int i = 1; i < n1 - 1; ++i) {
			int j = 0;
			result[i][j] = v[i][j];
			j = n2 - 1;
			result[i][j] = v[i][j];
		}
	}

	// outflow du/dx = 0
	for (int j = 0; j < n2; ++j) {

		result[n1 - 1][j] = (3.0 * v[n1 - 1][j] - 4.0 * v[n1 - 2][j] + 1.0 * v[n1 - 3][j]) / (2.0*par.d_x);
		result[0][j] = v[0][j];
	}

	return result;

}


Matrix CalculateB(Matrix &u_n, Matrix &v_n, Matrix &u_prev, Matrix &v_prev, Matrix &p, Matrix &force, Param par, Direction Dir) {

	size_t Nx = u_n.size();
	size_t Ny = u_n[0].size();

	double d_u, d_v, ij, N;
	if      (Dir == Du) { d_u = par.d_x;	d_v = par.d_y;	N = Ny;}
	else if (Dir == Dv) { d_u = par.d_y;	d_v = par.d_x;	N = Nx;}

	double d_uu = 1.0 / (d_u*d_u);
	double d_vv = 1.0 / (d_v*d_v);


	double k, kp, k0, km;

	CreateMatrix(result, Nx, Ny);

	for (int j = 1; j < (Ny - 1); ++j) {
		for (int i = 1; i < (Nx - 1); ++i) {
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
	for (int j = 0; j < Ny; ++j) {
		result[Nx - 1][j] = 0.0;
		result[0][j] = u_n[0][j];
	}

	return result;
}

void B_coefficients(int ij, int N, double &k, double &kp, double &k0, double &km) {
	if (ij > 1 && ij < N - 2) { k = 2. ;	kp = 1.      ;	k0 = -2.;	km = 1.; }
	else {
		                        k = 1.5;
		if      (ij == 1)                 { kp = 4. / 3.;	k0 = -4.;	km = 8. / 3.; }
		else if (ij == N - 2)             { kp = 8. / 3.;	k0 = -4.;	km = 4. / 3.; }
	}
}

double advective_term(Matrix &u, Matrix &v, int i, int j, double d_x, double d_y, double k, Direction Dir) {
	double v_help = 0.25 * (v[i][j] + R(v, i, j, Dir) + D(v, i, j, Dir) + RD(v, i, j, Dir));
	double result = u[i][j] * (R(u, i, j, Dir) - L(u, i, j, Dir)) / (2.0*d_x)
	               + v_help * (U(u, i, j, Dir) - D(u, i, j, Dir)) / (k  *d_y);
	return result;
}

double diffusion_term(Matrix &u, int i, int j, double d_xx, double d_yy, double kp, double k0, double km, Direction Dir) {
	double result = d_xx * (     R(u, i, j, Dir) - 2.0 * u[i][j] +      L(u, i, j, Dir))
	              + d_yy * (kp * U(u, i, j, Dir) + k0  * u[i][j] + km * D(u, i, j, Dir));
	return result;
}

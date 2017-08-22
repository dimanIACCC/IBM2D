#include "stdafx.h"
#include "PredictVel.h"


void ExplicPredVel(Matrix& U_predict, Matrix& V_predict, Matrix& U_n, Matrix& V_n, Matrix& P, Matrix& Force_x, Matrix& Force_y, Grid grid) {
	
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
			double LaplasU = (U_n[i + 1][j] - 2 * U_n[i][j] + U_n[i - 1][j]) / pow(grid.d_x, 2) + (U_n[i][j + 1] - 2 * U_n[i][j] + U_n[i][j - 1]) / pow(grid.d_y, 2);
			double CentrDiffx_U = (U_n[i + 1][j] - U_n[i - 1][j]) / (2 * grid.d_x);
			double CentrDiffy_U = (U_n[i][j + 1] - U_n[i][j - 1]) / (2 * grid.d_y);
			double GradPressX   = (P[i + 1][j] - P[i - 1][j]) / (2 * grid.d_x);
			U_predict[i][j] = grid.d_t*(LaplasU - U_n[i][j] * CentrDiffx_U - V_n[i][j] * CentrDiffy_U) + U_n[i][j] - GradPressX +Force_x[i][j];
		}
	}
	for (int i = 1; i < (int)V_predict.size()-1; i++) {
		for (int j = 1; j < (int)V_predict[0].size()-1; j++)
		{
			double LaplasV = (V_n[i + 1][j] - 2 * V_n[i][j] + V_n[i - 1][j]) / pow(grid.d_x, 2) + (V_n[i][j + 1] - 2 * V_n[i][j] + V_n[i][j - 1]) / pow(grid.d_y, 2);
			double CentrDiffx_V = (V_n[i + 1][j] - V_n[i - 1][j]) / (2 * grid.d_x);
			double CentrDiffy_V = (V_n[i][j + 1] - V_n[i][j - 1]) / (2 * grid.d_y);
			double GradPressY   = (P[i][j + 1] - P[i][j - 1]) / (2 * grid.d_y);
			V_predict[i][j] = grid.d_t*(LaplasV - U_n[i][j] * CentrDiffx_V - V_n[i][j] * CentrDiffy_V) + V_n[i][j] - GradPressY +Force_y[i][j];
		}
	}



}


void Calculate_A_u(Matrix A[5], Grid grid, double Re) {

	double d_xx = 1.0 / (grid.d_x*grid.d_x);
	double d_yy = 1.0 / (grid.d_y*grid.d_y);
	int const n1 = grid.N1;
	int const n2 = grid.N2 + 1;


	for (int j = 1; j < (n2 - 1); ++j) {
		for (int i = 1; i < (n1 - 1); ++i) {

			A[0][i][j] = 1.0 / grid.d_t + (1.0 / Re) * (d_xx + d_yy);
			A[1][i][j] = -1.0 / (2.0*Re*grid.d_y*grid.d_y);
			A[2][i][j] = -1.0 / (2.0*Re*grid.d_x*grid.d_x);
			A[3][i][j] = -1.0 / (2.0*Re*grid.d_y*grid.d_y);
			A[4][i][j] = -1.0 / (2.0*Re*grid.d_x*grid.d_x);


			if (j == 1) {
				A[0][i][j] = 1.0 / grid.d_t + (1.0 / Re) * (d_xx + 2.0*d_yy);
				A[1][i][j] = -2.0 / (3.0*Re*grid.d_y*grid.d_y);
				A[2][i][j] = -1.0 / (2.0*Re*grid.d_x*grid.d_x);
				A[3][i][j] = -4.0 / (3.0*Re*grid.d_y*grid.d_y);
				A[4][i][j] = -1.0 / (2.0*Re*grid.d_x*grid.d_x);
			}

			if (j == n2 - 2) {
				A[0][i][j] = 1.0 / grid.d_t + (1.0 / Re) * (d_xx + 2.0*d_yy);
				A[1][i][j] = -4.0 / (3.0*Re*grid.d_y*grid.d_y);
				A[2][i][j] = -1.0 / (2.0*Re*grid.d_x*grid.d_x);
				A[3][i][j] = -2.0 / (3.0*Re*grid.d_y*grid.d_y);
				A[4][i][j] = -1.0 / (2.0*Re*grid.d_x*grid.d_x);
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
		A[0][n1 - 1][j] = 3.0 / (2.0*grid.d_x);
		A[1][n1 - 1][j] = -4.0 / (2.0*grid.d_x);
		A[2][n1 - 1][j] = 1.0 / (2.0*grid.d_x);
		A[0][0][j] = 1.0;
	}



}

void Calculate_A_v(Matrix A[5], Grid grid, double Re) {

	double d_xx = 1.0 / (grid.d_x*grid.d_x);
	double d_yy = 1.0 / (grid.d_y*grid.d_y);
	int const n1 = grid.N1 + 1;
	int const n2 = grid.N2;


	for (int j = 1; j < (n2 - 1); ++j) {
		for (int i = 1; i < (n1 - 1); ++i) {
			A[0][i][j] = 1.0 / grid.d_t + (1.0 / Re) * (d_xx + d_yy);
			A[1][i][j] = -1.0 / (2.0*Re*grid.d_y*grid.d_y);
			A[2][i][j] = -1.0 / (2.0*Re*grid.d_x*grid.d_x);
			A[3][i][j] = -1.0 / (2.0*Re*grid.d_y*grid.d_y);
			A[4][i][j] = -1.0 / (2.0*Re*grid.d_x*grid.d_x);

			if (i == 1) {
				A[0][i][j] = 1.0 / grid.d_t + (1.0 / Re) * (2.0*d_xx + d_yy);
				A[1][i][j] = -1.0 / (2.0*Re*grid.d_y*grid.d_y);
				A[2][i][j] = -2.0 / (3.0*Re*grid.d_x*grid.d_x);
				A[3][i][j] = -1.0 / (2.0*Re*grid.d_y*grid.d_y);
				A[4][i][j] = -4.0 / (3.0*Re*grid.d_x*grid.d_x);
			}

			if (i == n1 - 2) {
				A[0][i][j] = 1.0 / grid.d_t + (1.0 / Re) * (2.0*d_xx + d_yy);
				A[1][i][j] = -1.0 / (2.0*Re*grid.d_y*grid.d_y);
				A[2][i][j] = -4.0 / (3.0*Re*grid.d_x*grid.d_x);
				A[3][i][j] = -1.0 / (2.0*Re*grid.d_y*grid.d_y);
				A[4][i][j] = -2.0 / (3.0*Re*grid.d_x*grid.d_x);
			}

		}
	}

	for (int i = 1; i < n1 - 1; ++i) {


		int j = 0;
		A[0][i][j] = -grid.d_y * 0.5;
		A[1][i][j] = 0.0;
		j = n2 - 1;
		A[0][i][j] = grid.d_y * 0.5;
		A[1][i][j] = 0.0;


	}

	// outflow du/dx = 0
	for (int j = 0; j < n2; ++j) {
		A[0][n1 - 1][j] = 3.0 / (2.0*grid.d_x);
		A[1][n1 - 1][j] = -4.0 / (2.0*grid.d_x);
		A[2][n1 - 1][j] = 1.0 / (2.0*grid.d_x);
		A[0][0][j] = 1.0;
	}

}

Matrix Operator_Ax(Matrix A[5], Matrix &v, int const n1, int const n2, Grid g,bool OverFlow) {

	CreateMatrix(result, n1, n2);


	for (int j = 1; j < (n2 - 1); ++j) {
		for (int i = 1; i < (n1 - 1); ++i) {

			result[i][j] = A[0][i][j] * v[i][j] + A[1][i][j] * v[i][j + 1] + A[2][i][j] * v[i + 1][j] + A[3][i][j] * v[i][j - 1] + A[4][i][j] * v[i - 1][j];
		}
	}
	if (OverFlow) {
		for (int i = 1; i < n1 - 1; ++i) {
			int j = 0;
			result[i][j] = (A[1][i][j] * v[i][j + 1] - A[0][i][j] * v[i][j]) / (g.d_y*0.5);
			j = n2 - 1;
			result[i][j] = (A[0][i][j] * v[i][j] - A[1][i][j] * v[i][j - 1]) / (g.d_y*0.5);
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

		result[n1 - 1][j] = (3.0 * v[n1 - 1][j] - 4.0 * v[n1 - 2][j] + 1.0 * v[n1 - 3][j]) / (2.0*g.d_x);
		result[0][j] = v[0][j];
	}

	return result;

}


Matrix CalculateB_u(Matrix &u_n, Matrix &v_n, Matrix &u_prev, Matrix &v_prev, Matrix &p, Matrix &force, Grid grid, double Re) {

	double d_xx = 1.0 / (grid.d_x*grid.d_x);
	double d_yy = 1.0 / (grid.d_y*grid.d_y);
	int const n1 = grid.N1;
	int const n2 = grid.N2 + 1;

	double advective_term_n = 0.0;
	double advective_term_prev = 0.0;
	double diffusion_term = 0.0;
	double pressure = 0.0;
	double v_help = 0.0;


	CreateMatrix(result, n1, n2);

	for (int j = 1; j < (n2 - 1); ++j) {
		for (int i = 1; i < (n1 - 1); ++i) {

			v_help              = 0.25 * (v_n[i][j] + v_n[i + 1][j] + v_n[i][j - 1] + v_n[i + 1][j - 1]);
			advective_term_n    = u_n[i][j] * (u_n[i + 1][j] - u_n[i - 1][j]) / (2.0*grid.d_x) + v_help * (u_n[i][j + 1] - u_n[i][j - 1]) / (2.0*grid.d_y);
			v_help              = 0.25 * (v_prev[i][j] + v_prev[i + 1][j] + v_prev[i][j - 1] + v_prev[i + 1][j - 1]);
			advective_term_prev = u_prev[i][j] * (u_prev[i + 1][j] - u_prev[i - 1][j]) / (2.0*grid.d_x) + v_help * (u_prev[i][j + 1] - u_prev[i][j - 1]) / (2.0*grid.d_y);
			diffusion_term      = d_xx * (u_n[i + 1][j] - 2.0 * u_n[i][j] + u_n[i - 1][j]) + d_yy * (u_n[i][j + 1] - 2.0 * u_n[i][j] + u_n[i][j - 1]);
			pressure            = (p[i + 1][j] - p[i][j]) / (grid.d_x);
			
			result[i][j]        = -(3.0 / 2.0 * advective_term_n - 1.0 / 2.0 * advective_term_prev) - pressure + 1.0 / (2.0*Re) * (diffusion_term)+u_n[i][j] / grid.d_t;

			if (j == 1 || j == n2 - 2) {
				v_help              = 0.25 * (v_n[i][j] + v_n[i + 1][j] + v_n[i][j - 1] + v_n[i + 1][j - 1]);
				advective_term_n    = u_n[i][j] * (u_n[i + 1][j] - u_n[i - 1][j]) / (2.0*grid.d_x) + v_help * (u_n[i][j + 1] - u_n[i][j - 1]) / (1.5*grid.d_y);
				v_help              = 0.25 * (v_prev[i][j] + v_prev[i + 1][j] + v_prev[i][j - 1] + v_prev[i + 1][j - 1]);
				advective_term_prev = u_prev[i][j] * (u_prev[i + 1][j] - u_prev[i - 1][j]) / (2.0*grid.d_x) + v_help * (u_prev[i][j + 1] - u_prev[i][j - 1]) / (1.5*grid.d_y);
				if (j == 1) {

					diffusion_term = d_xx * (u_n[i + 1][j] - 2.0 * u_n[i][j] + u_n[i - 1][j]) + d_yy * (4.0*u_n[i][j + 1] - 12.0 * u_n[i][j] + 8.0*u_n[i][j - 1]) / 3.0;
				}
				if (j == n2 - 2) {

					diffusion_term = d_xx * (u_n[i + 1][j] - 2.0 * u_n[i][j] + u_n[i - 1][j]) + d_yy * (8.0*u_n[i][j + 1] - 12.0 * u_n[i][j] + 4.0*u_n[i][j - 1]) / 3.0;
				}

				result[i][j] = -(3.0 / 2.0 * advective_term_n - 1.0 / 2.0 * advective_term_prev) - pressure + 1.0 / (2.0*Re) * (diffusion_term)+u_n[i][j] / grid.d_t;
			}
		}
	}

	// outflow du/dx = 0
	for (int j = 0; j < n2; ++j) {
		result[n1 - 1][j] = 0.0;
		result[0][j] = u_n[0][j];
	}

	for (int i = 1; i < n1 - 1; ++i) {
		for (int j = 1; j < n2 - 1; ++j) {

			result[i][j] += force[i][j];
		}
	}


	return result;
}

Matrix CalculateB_v(Matrix &u_n, Matrix &v_n, Matrix &u_prev, Matrix &v_prev, Matrix &p, Matrix &force, Grid grid, double Re) {

	double d_xx = 1.0 / (grid.d_x*grid.d_x);
	double d_yy = 1.0 / (grid.d_y*grid.d_y);
	int const n1 = grid.N1 + 1;
	int const n2 = grid.N2;

	double advective_term_n = 0.0;
	double advective_term_prev = 0.0;
	double diffusion_term = 0.0;
	double pressure = 0.0;
	double u_help = 0.0;

	CreateMatrix(result, n1, n2);


	for (int j = 1; j < (n2 - 1); ++j) {
		for (int i = 1; i < (n1 - 1); ++i) {
			u_help = 0.25 * (u_n[i][j] + u_n[i - 1][j] + u_n[i][j + 1] + u_n[i - 1][j + 1]);
			advective_term_n = u_help * (v_n[i + 1][j] - v_n[i - 1][j]) / (2.0*grid.d_x) + v_n[i][j] * (v_n[i][j + 1] - v_n[i][j - 1]) / (2.0*grid.d_y);
			u_help = 0.25 * (u_prev[i][j] + u_prev[i - 1][j] + u_prev[i][j + 1] + u_prev[i - 1][j + 1]);
			advective_term_prev = u_help * (v_prev[i + 1][j] - v_prev[i - 1][j]) / (2.0*grid.d_x) + v_prev[i][j] * (v_prev[i][j + 1] - v_prev[i][j - 1]) / (2.0*grid.d_y);
			diffusion_term = d_xx * (v_n[i + 1][j] - 2.0 * v_n[i][j] + v_n[i - 1][j]) + d_yy * (v_n[i][j + 1] - 2.0 * v_n[i][j] + v_n[i][j - 1]);
			pressure = (p[i][j + 1] - p[i][j]) / (grid.d_y);
			result[i][j] = -(3.0 / 2.0 * advective_term_n - 1.0 / 2.0 * advective_term_prev) - pressure + 1.0 / (2.0*Re) * (diffusion_term)+v_n[i][j] / grid.d_t;

			if (i == 1 || i == n1 - 2) {


				u_help = 0.25 * (u_n[i][j] + u_n[i - 1][j] + u_n[i][j + 1] + u_n[i - 1][j + 1]);
				advective_term_n = u_help * (v_n[i + 1][j] - v_n[i - 1][j]) / (1.5*grid.d_x) + v_n[i][j] * (v_n[i][j + 1] - v_n[i][j - 1]) / (2.0*grid.d_y);
				u_help = 0.25 * (u_prev[i][j] + u_prev[i - 1][j] + u_prev[i][j + 1] + u_prev[i - 1][j + 1]);
				advective_term_prev = u_help * (v_prev[i + 1][j] - v_prev[i - 1][j]) / (1.5*grid.d_x) + v_prev[i][j] * (v_prev[i][j + 1] - v_prev[i][j - 1]) / (2.0*grid.d_y);
				if (i == 1) {

					diffusion_term = d_xx * (4.0*v_n[i + 1][j] - 12.0 * v_n[i][j] + 8.0*v_n[i - 1][j]) / 3.0 + d_yy * (v_n[i][j + 1] - 2.0 * v_n[i][j] + v_n[i][j - 1]);
				}
				if (i == n1 - 2) {

					diffusion_term = d_xx * (8.0*v_n[i + 1][j] - 12.0 * v_n[i][j] + 4.0*v_n[i - 1][j]) / 3.0 + d_yy * (v_n[i][j + 1] - 2.0 * v_n[i][j] + v_n[i][j - 1]);
				}

				result[i][j] = -(3.0 / 2.0 * advective_term_n - 1.0 / 2.0 * advective_term_prev) - pressure + 1.0 / (2.0*Re) * (diffusion_term)+v_n[i][j] / grid.d_t;
			}
		}

	}


	for (int j = 0; j < n2; ++j) {
		result[n1 - 1][j] = 0.0;
		result[0][j] = v_n[0][j];
	}

	for (int i = 1; i < n1 - 1; ++i) {
		for (int j = 1; j < n2 - 1; ++j) {
			result[i][j] += force[i][j];
		}
	}

	return result;

}
#include "Calculate_B.h"


Matrix CalculateB_u(Matrix &u_n, Matrix &v_n, Matrix &u_prev, Matrix &v_prev, Matrix &p, Matrix &force, Grid grid,double Re){

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

	for (int j = 1; j < (n2 - 1); ++j){
		for (int i = 1; i < (n1 - 1); ++i){

			v_help              = 0.25 * (v_n[i][j] + v_n[i + 1][j] + v_n[i][j - 1] + v_n[i + 1][j - 1]);
			advective_term_n    = u_n[i][j] * (u_n[i + 1][j] - u_n[i - 1][j]) / (2.0*grid.d_x) + v_help * (u_n[i][j + 1] - u_n[i][j - 1]) / (2.0*grid.d_y);
			v_help			    = 0.25 * (v_prev[i][j] + v_prev[i + 1][j] + v_prev[i][j - 1] + v_prev[i + 1][j - 1]);
			advective_term_prev = u_prev[i][j] * (u_prev[i + 1][j] - u_prev[i - 1][j]) / (2.0*grid.d_x) + v_help * (u_prev[i][j + 1] - u_prev[i][j - 1]) / (2.0*grid.d_y);
			diffusion_term      = d_xx * (u_n[i + 1][j] - 2.0 * u_n[i][j] + u_n[i - 1][j]) + d_yy * (u_n[i][j + 1] - 2.0 * u_n[i][j] + u_n[i][j - 1]);
			pressure            = (p[i + 1][j] - p[i][j]) / (grid.d_x);
			result[i][j]        = -(3.0 / 2.0 * advective_term_n - 1.0 / 2.0 * advective_term_prev) - pressure + 1.0 / (2.0*Re) * (diffusion_term)+u_n[i][j] / grid.d_t;

			if (j == 1 || j == n2 - 2){
				v_help = 0.25 * (v_n[i][j] + v_n[i + 1][j] + v_n[i][j - 1] + v_n[i + 1][j - 1]);
				advective_term_n = u_n[i][j] * (u_n[i + 1][j] - u_n[i - 1][j]) / (2.0*grid.d_x) + v_help * (u_n[i][j + 1] - u_n[i][j - 1]) / (1.5*grid.d_y);
				v_help = 0.25 * (v_prev[i][j] + v_prev[i + 1][j] + v_prev[i][j - 1] + v_prev[i + 1][j - 1]);
				advective_term_prev = u_prev[i][j] * (u_prev[i + 1][j] - u_prev[i - 1][j]) / (2.0*grid.d_x) + v_help * (u_prev[i][j + 1] - u_prev[i][j - 1]) / (1.5*grid.d_y);
				if (j == 1){

					diffusion_term = d_xx * (u_n[i + 1][j] - 2.0 * u_n[i][j] + u_n[i - 1][j]) + d_yy * (4.0*u_n[i][j + 1] - 12.0 * u_n[i][j] + 8.0*u_n[i][j - 1]) / 3.0;
				}
				if (j == n2 - 2){

					diffusion_term = d_xx * (u_n[i + 1][j] - 2.0 * u_n[i][j] + u_n[i - 1][j]) + d_yy * (8.0*u_n[i][j + 1] - 12.0 * u_n[i][j] + 4.0*u_n[i][j - 1]) / 3.0;
				}

				result[i][j] = -(3.0 / 2.0 * advective_term_n - 1.0 / 2.0 * advective_term_prev) - pressure + 1.0 / (2.0*Re) * (diffusion_term)+u_n[i][j] / grid.d_t;
			}
		}
	}

	// outflow du/dx = 0
	for (int j = 0; j < n2; ++j){
		result[n1 - 1][j] = 0.0;
		result[0][j] = u_n[0][j];
	}
	
	for (int i = 1; i < n1 - 1; ++i){
		for (int j = 1; j < n2 - 1; ++j){

			result[i][j] += force[i][j];
		}
	}


	return result;
}

Matrix CalculateB_v(Matrix &u_n, Matrix &v_n, Matrix &u_prev, Matrix &v_prev, Matrix &p, Matrix &force, Grid grid, double Re){

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


	for (int j = 1; j < (n2 - 1); ++j){
		for (int i = 1; i < (n1 - 1); ++i){
			u_help = 0.25 * (u_n[i][j] + u_n[i - 1][j] + u_n[i][j + 1] + u_n[i - 1][j + 1]);
			advective_term_n = u_help * (v_n[i + 1][j] - v_n[i - 1][j]) / (2.0*grid.d_x) + v_n[i][j] * (v_n[i][j + 1] - v_n[i][j - 1]) / (2.0*grid.d_y);
			u_help = 0.25 * (u_prev[i][j] + u_prev[i - 1][j] + u_prev[i][j + 1] + u_prev[i - 1][j + 1]);
			advective_term_prev = u_help * (v_prev[i + 1][j] - v_prev[i - 1][j]) / (2.0*grid.d_x) + v_prev[i][j] * (v_prev[i][j + 1] - v_prev[i][j - 1]) / (2.0*grid.d_y);
			diffusion_term = d_xx * (v_n[i + 1][j] - 2.0 * v_n[i][j] + v_n[i - 1][j]) + d_yy * (v_n[i][j + 1] - 2.0 * v_n[i][j] + v_n[i][j - 1]);
			pressure = (p[i][j + 1] - p[i][j]) / (grid.d_y);
			result[i][j] = -(3.0 / 2.0 * advective_term_n - 1.0 / 2.0 * advective_term_prev) - pressure + 1.0 / (2.0*Re) * (diffusion_term)+v_n[i][j] / grid.d_t;

			if (i == 1 || i == n1 - 2){


				u_help = 0.25 * (u_n[i][j] + u_n[i - 1][j] + u_n[i][j + 1] + u_n[i - 1][j + 1]);
				advective_term_n = u_help * (v_n[i + 1][j] - v_n[i - 1][j]) / (1.5*grid.d_x) + v_n[i][j] * (v_n[i][j + 1] - v_n[i][j - 1]) / (2.0*grid.d_y);
				u_help = 0.25 * (u_prev[i][j] + u_prev[i - 1][j] + u_prev[i][j + 1] + u_prev[i - 1][j + 1]);
				advective_term_prev = u_help * (v_prev[i + 1][j] - v_prev[i - 1][j]) / (1.5*grid.d_x) + v_prev[i][j] * (v_prev[i][j + 1] - v_prev[i][j - 1]) / (2.0*grid.d_y);
				if (i == 1){

					diffusion_term = d_xx * (4.0*v_n[i + 1][j] - 12.0 * v_n[i][j] + 8.0*v_n[i - 1][j]) / 3.0 + d_yy * (v_n[i][j + 1] - 2.0 * v_n[i][j] + v_n[i][j - 1]);
				}
				if (i == n1 - 2){

					diffusion_term = d_xx * (8.0*v_n[i + 1][j] - 12.0 * v_n[i][j] + 4.0*v_n[i - 1][j]) / 3.0 + d_yy * (v_n[i][j + 1] - 2.0 * v_n[i][j] + v_n[i][j - 1]);
				}

				result[i][j] = -(3.0 / 2.0 * advective_term_n - 1.0 / 2.0 * advective_term_prev) - pressure + 1.0 / (2.0*Re) * (diffusion_term)+v_n[i][j] / grid.d_t;
			}
		}

	}


	for (int j = 0; j < n2; ++j){
		result[n1 - 1][j] = 0.0;
		result[0][j] = v_n[0][j];
	}

	for (int i = 1; i < n1 - 1; ++i){
		for (int j = 1; j < n2 - 1; ++j){
			result[i][j] += force[i][j];
		}
	}

	return result;

}

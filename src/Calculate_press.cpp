#include "stdafx.h"
#include "Calculate_press.h"



double Calculate_Press_correction(Matrix &delta_p, Matrix &b_p, Param par, int &N_out){

	double eps;
	double delta_p_max;

	double dx2 = 1.0 / (par.d_x*par.d_x);
	double dy2 = 1.0 / (par.d_y*par.d_y);
	double A   = 1.0 / (2.0 * (dx2 + dy2));

	double help;

	size_t n1 = delta_p.size();
	size_t n2 = delta_p[0].size();

	for (int n = 0; n < par.N_Zeidel; n++) {
		eps = 0.0;
		delta_p_max = 0.0;

		for (size_t i = 1; i < n1-1; ++i){
			for (size_t j = 1; j < n2-1; ++j){

				help = A * (dx2 * (delta_p[i + 1][j] + delta_p[i - 1][j])
				          + dy2 * (delta_p[i][j + 1] + delta_p[i][j - 1]) - b_p[i][j]);

				// Fixed pressure
				if (i == n1 - 2) { // Right boundary
					if (par.BC == periodical) {
						if (j == 1)
							help = 0.0; // down corner
					}
					else if (par.BC == u_inflow || par.BC == u_infinity) {
						help = 0.0;
					}
				}

				if (par.BC == Taylor_Green && i == par.N1 / 2 && j == par.N2 / 2) {
					help = 0;
				}

				if (fabs(help) > delta_p_max) {
					delta_p_max = fabs(help);
				}
				if (fabs(help - delta_p[i][j]) > eps){
					eps = fabs(help - delta_p[i][j]);
				}

				delta_p[i][j] = help;

			}
		}

		if (par.BC == u_infinity || par.BC == u_inflow || par.BC == periodical || par.BC == Taylor_Green) {

			// Up-Down BC
			for (size_t i = 0; i < n1; ++i) {
				delta_p[i][0]      = delta_p[i][1];            // D
				delta_p[i][n2 - 1] = delta_p[i][n2 - 2];       // U
			}

			// Left-Right BC
			for (size_t j = 0; j < n2; ++j) {
				delta_p[0][j]      = delta_p[1][j];            // L
				delta_p[n1 - 1][j] = delta_p[n1 - 2][j];       // R
			}

			// Corners
			delta_p[0][0]           = delta_p[1][1];            // LD
			delta_p[0][n2 - 1]      = delta_p[1][n2 - 2];       // LU
			delta_p[n1 - 1][0]      = delta_p[n1 - 2][1];       // RD
			delta_p[n1 - 1][n2 - 1] = delta_p[n1 - 2][n2 - 2];  // RU
		}

		// Periodical Left-Right BC
		if (par.BC == periodical) {
			for (size_t j = 0; j < n2; ++j) {
				delta_p[0][j] = delta_p[n1 - 2][j];        // L
				delta_p[n1 - 1][j] = delta_p[1][j];        // R
			}
		}

		/*if (n % 100 == 0) {
			std::cout << eps << "  " << delta_p_max << std::endl;
			OutputPressure(delta_p, n, solidList, par);
		}*/

		if (eps/delta_p_max < par.Zeidel_eps){
			N_out = n;
			break;
		}

	}



	if (eps / delta_p_max > par.Zeidel_eps) {
		std::cout << "Zeidel has not converged, eps_p = " << eps << ",   delta_p_max = " << delta_p_max << std::endl;
		
	}

	return delta_p_max;

}


Matrix Calculate_Press_Right(Matrix &u, Matrix &v, Param par){
	double d = 0.0;

	size_t n1 = par.N1 + 1;
	size_t n2 = par.N2 + 1;
	CreateMatrix(result, n1, n2);


	for (size_t i = 1; i < n1 - 1; ++i){
		for (size_t j = 1; j < n2 - 1; ++j){

			double d = (1.0 / par.d_x) * (u[i + 1][j] - u[i][j])
			         + (1.0 / par.d_y) * (v[i][j + 1] - v[i][j]);

			result[i][j] = d / par.d_t;

		}

	}
	result[0][0] = 0.0;
	result[n1 - 1][0] = 0.0;
	result[0][n2 - 1] = 0.0;
	result[n1 - 1][n2 - 1] = 0.0;


	for (size_t i = 1; i < n1 - 1; ++i){
		result[i][0] = 0.0;
		result[i][n2 - 1] = 0.0;
	}
	for (size_t j = 1; j < n2 - 1; ++j){
		result[0][j] = 0.0;
		result[n1 - 1][j] = 0.0;
	}

	return result;

}

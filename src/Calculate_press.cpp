#include "stdafx.h"
#include "Calculate_press.h"



double Calculate_Press_correction(Matrix &delta_p, Matrix &b_p, Param par, int &N_out){

	double eps;
	double delta_p_max;

	double dx2 = 1.0 / (par.d_x*par.d_x);
	double dy2 = 1.0 / (par.d_y*par.d_y);
	double A   = 1.0 / (2.0 * (dx2 + dy2));

	double p_tmp, p_fix = 0;

	size_t n1 = delta_p.size();
	size_t n2 = delta_p[0].size();

	for (int n = 0; n < par.N_Zeidel; n++) {
		eps = 0.0;
		delta_p_max = 0.0;

		for (size_t i = 1; i < n1-1; ++i){
			for (size_t j = 1; j < n2-1; ++j){

				p_tmp = A * (dx2 * (delta_p[i + 1][j] + delta_p[i - 1][j])
				           + dy2 * (delta_p[i][j + 1] + delta_p[i][j - 1]) - b_p[i][j]);

				if (fabs(p_tmp) > delta_p_max) {
					delta_p_max = fabs(p_tmp);
				}
				if (fabs(p_tmp - delta_p[i][j]) > eps){
					eps = fabs(p_tmp - delta_p[i][j] - p_fix);
				}

				delta_p[i][j] = p_tmp;
			}
		}

		// subtract constant pressure to make delta_p = 0 in the centre of the domain
		//if (par.BC == periodical || par.BC == Taylor_Green) {
		if (par.BC == periodical) {
			p_fix = delta_p[par.N1 / 2][par.N2 / 2];
			for (size_t i = 1; i < n1 - 1; ++i)
			for (size_t j = 1; j < n2 - 1; ++j)
				delta_p[i][j] -= p_fix;
		}

		if (par.BC == u_infinity || par.BC == u_inflow || par.BC == periodical) {

			// Up-Down BC
			for (size_t i = 0; i < n1; ++i) {
				delta_p[i][0]      = delta_p[i][1];            // D
				delta_p[i][n2 - 1] = delta_p[i][n2 - 2];       // U
			}

			// Left-Right BC
			for (size_t j = 0; j < n2; ++j) {
				delta_p[0     ][j] = delta_p[1     ][j];       // L
				delta_p[n1 - 1][j] = delta_p[n1 - 2][j];       // R
				if (par.BC == u_inflow || par.BC == u_infinity) {
					delta_p[n1 - 1][j] = - delta_p[n1 - 2][j];
				}
			}

		}

		// Periodical Left-Right BC
		if (par.BC == periodical) {
			for (size_t j = 0; j < n2; ++j) {
				delta_p[0     ][j] = delta_p[n1 - 2][j];        // L
				delta_p[n1 - 1][j] = delta_p[1     ][j];        // R
			}
		}

		/*if (n % 100 == 0) {
			std::cout << eps << "  " << delta_p_max << std::endl;
			Output_P(delta_p, "delta_p", n, par);
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

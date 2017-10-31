#include "Calculate_press.h"



double Calculate_Press_correction(Matrix &delta_p, Matrix &b_p, Param par){

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
		for (size_t i = 0; i < n1; ++i){
			for (size_t j = 0; j < n2; ++j){

				help = A * (dx2 * (R(delta_p, i, j, Du, n1, n2) + L(delta_p, i, j, Du, n1, n2))
				          + dy2 * (U(delta_p, i, j, Du, n1, n2) + D(delta_p, i, j, Du, n1, n2)) - b_p[i][j]);

				if (i == 0     )                 help = delta_p[i + 1][j];       // L
				if (i == n1 - 1)                 help = delta_p[i - 1][j];       // R
				if (j == 0     )                 help = delta_p[i][j + 1];       // D
				if (j == n2 - 1)                 help = delta_p[i][j - 1];       // U

				if (i == 0      && j == 0     )  help = delta_p[i + 1][j + 1];   // LD
				if (i == 0      && j == n2 - 1)  help = delta_p[i + 1][j - 1];   // LU
				if (i == n1 - 1 && j == 0     )  help = delta_p[i - 1][j + 1];   // RD
				if (i == n1 - 1 && j == n2 - 1)  help = delta_p[i - 1][j - 1];   // RU

				if (par.BC == periodical && i == 0     ) help = delta_p[n1 - 2][j];
				if (par.BC == periodical && i == n1 - 1) help = delta_p[1][j];

				if (j == 0 || j == 1 || j == n2 - 2 || j == n2 - 1) {
					if ( i == n1 - 1 || i == n1 - 2)
						help = 0.0; // Right boundary condition
					if ((i == 0      || i == 1     ) && par.BC == periodical)
						help = 0.0; // Left boundary condition for periodical problem
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


		if (eps/delta_p_max < par.Zeidel_eps){
			break;
		}

	}

	if (eps / delta_p_max > par.Zeidel_eps) {
		std::cout << "Zeidel has not converged, eps_p = " << eps << ",   delta_p_max = " << delta_p_max << std::endl;
	}

	return delta_p_max;

}


Matrix Calculate_Press_Right(Matrix &u, Matrix &v, Matrix &Fx, Matrix &Fy, Param par){
	double d = 0.0;

	size_t n1 = par.N1 + 1;
	size_t n2 = par.N2 + 1;
	CreateMatrix(result, n1, n2);


	for (size_t i = 1; i < n1 - 1; ++i){
		for (size_t j = 1; j < n2 - 1; ++j){

			double d = (1.0 / par.d_x) * (u[i][j] - u[i - 1][j])
			         + (1.0 / par.d_y) * (v[i][j] - v[i][j - 1]);

			double dF = (1.0 / par.d_x) * (Fx[i][j] - Fx[i - 1][j])
				      + (1.0 / par.d_y) * (Fy[i][j] - Fy[i][j - 1]);

			result[i][j] = d / par.d_t - dF;

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

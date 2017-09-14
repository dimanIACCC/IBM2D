#include "Calculate_press.h"



double Calculate_Press_correction(Matrix& delta_p, Matrix &b_p, Param par, bool OverFlow){

	for (int i = 0; i < (int)delta_p.size(); ++i) {
		std::fill(delta_p[i].begin(), delta_p[i].end(), 0);
	}

	int n = 0;
	double eps = 0.0;

	double dx2 = 1.0 / (par.d_x*par.d_x);
	double dy2 = 1.0 / (par.d_y*par.d_y);
	double A   = 1.0 / (2.0 * (dx2 + dy2));

	double kxp, kxm, kyp, kym;

	double help;

	int const n1 = par.N1 + 1;
	int const n2 = par.N2 + 1;

	while (n < par.N_Zeidel){
		eps = 0.0;
		for (int i = 0; i < n1; ++i){
			for (int j = 0; j < n2; ++j){

				if (i == 0     )                 help = delta_p[i + 1][j];       // L
				if (i == n1 - 1)                 help = delta_p[i - 1][j];       // R
				if (j == 0     )                 help = delta_p[i][j + 1];       // D
				if (j == n2 - 1)                 help = delta_p[i][j - 1];       // U

				if (i == 0      && j == 0     )  help = delta_p[i + 1][j + 1];   // LD
				if (i == 0      && j == n2 - 1)  help = delta_p[i + 1][j - 1];   // LU
				if (i == n1 - 1 && j == 0     )  help = delta_p[i - 1][j + 1];   // RD
				if (i == n1 - 1 && j == n2 - 1)  help = delta_p[i - 1][j - 1];   // RU

				if (i == n1 - 1 && (!OverFlow))  help = 0.0; // Right boundary condition


				if     (i > 0 && i < n1 - 1 && j > 0 && j < n2 - 1){
					if (i > 1 && i < n1 - 2 && j > 1 && j < n2 - 2)
						help = A * (dx2 * (delta_p[i + 1][j] + delta_p[i - 1][j])
						          + dy2 * (delta_p[i][j + 1] + delta_p[i][j - 1]) - b_p[i][j]);
					else {
						if (i == 1     )                 { kxp = 4.;	kxm = 8.;	kyp = 1.;	kym = 1.; }   // L
						if (i == n1 - 2)                 { kxp = 8.;	kxm = 4.;	kyp = 1.;	kym = 1.; }   // R
						if (j == 1     )                 { kxp = 1.;	kxm = 1.;	kyp = 4.;	kym = 8.; }   // D
						if (j == n2 - 2)                 { kxp = 1.;	kxm = 1.;	kyp = 8.;	kym = 4.; }   // U

						if (1 == i      && j == 1     )  { kxp = 4.;	kxm = 8.;	kyp = 4.;	kym = 8.; }   // LD
						if (1 == i      && j == n2 - 2)  { kxp = 4.;	kxm = 8.;	kyp = 8.;	kym = 4.; }   // LU
						if (n1 - 2 == i && j == 1     )  { kxp = 8.;	kxm = 4.;	kyp = 4.;	kym = 8.; }   // RD
						if (n1 - 2 == i && j == n2 - 2)  { kxp = 8.;	kxm = 4.;	kyp = 8.;	kym = 4.; }   // RU

						help = (dx2 * (kxp * delta_p[i + 1][j] + kxm * delta_p[i - 1][j])
						      + dy2 * (kyp * delta_p[i][j + 1] + kym * delta_p[i][j - 1]) - b_p[i][j])
						     / (dx2 * (kxp + kxm)
						      + dy2 * (kyp + kym));
					}
				}

				if (fabs(help - delta_p[i][j]) > eps){
					eps = fabs(help - delta_p[i][j]);
				}

				delta_p[i][j] = help;

			}
		}


		if (eps < par.Zeidel_eps){
			break;
		}
		n++;
	}

	if (eps > par.Zeidel_eps) {
		std::cout << "Zeidel has not converged, eps_p = " << eps << std::endl;
	}

	return eps;

}


Matrix Calculate_Press_Right(Matrix &u, Matrix &v, Param par){
	double d = 0.0;

	int const n1 = par.N1 + 1;
	int const n2 = par.N2 + 1;
	CreateMatrix(result, n1, n2);


	for (int i = 1; i < n1 - 1; ++i){
		for (int j = 1; j < n2 - 1; ++j){

			d = (1.0 / par.d_x) * (u[i][j] - u[i - 1][j])
			  + (1.0 / par.d_y) * (v[i][j] - v[i][j - 1]);

			result[i][j] = d / par.d_t;

		}

	}
	result[0][0] = 0.0;
	result[n1 - 1][0] = 0.0;
	result[0][n2 - 1] = 0.0;
	result[n1 - 1][n2 - 1] = 0.0;


	for (int i = 1; i < n1 - 1; ++i){
		result[i][0] = 0.0;
		result[i][n2 - 1] = 0.0;
	}
	for (int j = 1; j < n2 - 1; ++j){
		result[0][j] = 0.0;
		result[n1 - 1][j] = 0.0;
	}

	return result;

}

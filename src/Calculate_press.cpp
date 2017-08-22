#include "stdafx.h"
#include "Calculate_press.h"



double Calculate_Press_correction(Matrix& delta_p, Matrix &b_p, int const N_Zeidel, double const Zeidel_eps, Grid grid,bool OverFlow){

	int n = 0;
	double eps = 0.0;

	double a = (1.0 / (grid.d_x*grid.d_x) + 1.0 / (grid.d_y*grid.d_y));
	double b = 1.0 / (grid.d_x*grid.d_x);
	double c = 1.0 / (grid.d_y*grid.d_y);
	double d = 0.0;
	double help;

	int const n1 = grid.N1 + 1;
	int const n2 = grid.N2 + 1;

	while (n < N_Zeidel){
		eps = 0.0;
		for (int i = 0; i < n1; ++i){
			for (int j = 0; j < n2; ++j){
				if (0 == i && 0 == j){
					help = delta_p[i + 1][j + 1];
				}

				if (0 == i && 0 != j && n2 - 1 != j){
					help = delta_p[i + 1][j];
				}

				if (0 == i && n2 - 1 == j){
					help = delta_p[i + 1][j - 1];
				}

				if (n1 - 1 == i && 0 != j && n2 - 1 != j) {
					if (OverFlow) help = delta_p[i - 1][j];
					else help = 0.0;
				}
				if (n1 - 1 == i && 0 == j) {
					if (OverFlow) help = delta_p[i - 1][j + 1];
					else help = 0.0;
				}
				if (n1 - 1 == i && n2 - 1 == j) {
					if (OverFlow) help = delta_p[i - 1][j - 1];
					else help = 0.0;
				}


				if (0 != i && n1 - 1 != i && 0 == j){
					help = delta_p[i][j + 1];
				}

				if (0 != i && n1 - 1 != i && n2 - 1 == j){
					help = delta_p[i][j - 1];
				}

				if (0 != i && n1 - 1 != i && 0 != j && n2 - 1 != j){
					help = (1.0 / (2.0*a)) * (b * (delta_p[i + 1][j] + delta_p[i - 1][j]) + c * (delta_p[i][j + 1] + delta_p[i][j - 1]) - b_p[i][j]);

					if (1 == i){
						help = (1.0 / (12.0*b + 2.0*c)) * (b * (4.0*delta_p[i + 1][j] + 8.0*delta_p[i - 1][j]) + c * (delta_p[i][j + 1] + delta_p[i][j - 1]) - b_p[i][j]);
					}
					if (n1 - 2 == i){
						help = (1.0 / (12.0*b + 2.0*c)) * (b * (8.0*delta_p[i + 1][j] + 4.0*delta_p[i - 1][j]) + c * (delta_p[i][j + 1] + delta_p[i][j - 1]) - b_p[i][j]);
					}
					if (1 == j){
						help = (1.0 / (2.0*b + 12.0*c)) * (b * (delta_p[i + 1][j] + delta_p[i - 1][j]) + c * (4.0*delta_p[i][j + 1] + 8.0*delta_p[i][j - 1]) - b_p[i][j]);
					}
					if (n2 - 2 == j){
						help = (1.0 / (2.0*b + 12.0*c)) * (b * (delta_p[i + 1][j] + delta_p[i - 1][j]) + c * (8.0*delta_p[i][j + 1] + 4.0*delta_p[i][j - 1]) - b_p[i][j]);
					}

					if (1 == i && 1 == j){
						help = (1.0 / (12.0*a)) * (b * (4.0*delta_p[i + 1][j] + 8.0*delta_p[i - 1][j]) + c * (4.0*delta_p[i][j + 1] + 8.0*delta_p[i][j - 1]) - b_p[i][j]);
					}
					if (1 == i && n2 - 2 == j){
						help = (1.0 / (12.0*a)) * (b * (4.0*delta_p[i + 1][j] + 8.0*delta_p[i - 1][j]) + c * (8.0*delta_p[i][j + 1] + 4.0*delta_p[i][j - 1]) - b_p[i][j]);
					}
					if (n1 - 2 == i && 1 == j){
						help = (1.0 / (12.0*a)) * (b * (8.0*delta_p[i + 1][j] + 4.0*delta_p[i - 1][j]) + c * (4.0*delta_p[i][j + 1] + 8.0*delta_p[i][j - 1]) - b_p[i][j]);
					}
					if (n1 - 2 == i && n2 - 2 == j){
						help = (1.0 / (12.0*a)) * (b * (8.0*delta_p[i + 1][j] + 4.0*delta_p[i - 1][j]) + c * (8.0*delta_p[i][j + 1] + 4.0*delta_p[i][j - 1]) - b_p[i][j]);
					}
				}

				if (fabs(help - delta_p[i][j]) > eps){
					eps = fabs(help - delta_p[i][j]);
				}

				delta_p[i][j] = help;

			}
		}


		if (eps < Zeidel_eps){
			break;
		}
		n++;
	}



	return eps;

}


Matrix Calculate_Press_Right(Matrix &u, Matrix &v, Grid grid){
	double d = 0.0;

	int const n1 = grid.N1 + 1;
	int const n2 = grid.N2 + 1;
	CreateMatrix(result, n1, n2);


	for (int i = 1; i < n1 - 1; ++i){
		for (int j = 1; j < n2 - 1; ++j){

			d = (1.0 / grid.d_x) * (u[i][j] - u[i - 1][j]) + (1.0 / grid.d_y) * (v[i][j] - v[i][j - 1]);

			result[i][j] = d / grid.d_t;

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
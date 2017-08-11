#define _USE_MATH_DEFINES
#include <math.h>
#include "CalculateForce.h"

double DeltaFunction(double x, double y, Grid grid){
	return 1.0 / (grid.d_x*grid.d_y) * FunctionD(x / grid.d_x) * FunctionD(y / grid.d_y);
}

double FunctionD(double r){
	if ((0.0 <= fabs(r)) && (fabs(r) < 1.0)){
		return 1.0 / 8.0*(3.0 - 2.0 * fabs(r) + sqrt(1.0 + 4.0 * fabs(r) - 4.0 * r * r));
	}
	if ((1.0 <= fabs(r)) && (fabs(r) < 2.0)){
		return 1.0 / 8.0*(5.0 - 2.0 * fabs(r) - sqrt(-7.0 + 12.0 * fabs(r) - 4.0 * r * r));
	}
	if (2.0 <= fabs(r)){
		return 0.0;
	}
	return 0;
}
void GetInfluenceArea(int& i_min, int& i_max, int& j_min, int& j_max, double x, double y, int size, Grid grid){

	i_max = (int)((x / grid.d_x) + size);
	i_min = (int)((x / grid.d_x) - size);

	j_max = (int)(y / grid.d_y) + size;
	j_min = (int)(y / grid.d_y) - size;

	if (i_min < 0){
		i_min = 0;
	}
	if (j_min < 0){
		j_min = 0;
	}
}

double CalculateForce_X(Matrix& force_x, list<Circle> &iList, Matrix& u, double r, double & Coeff, Grid grid, double alpha_f, double beta_f, double M){

	int const n1 = grid.N1;
	int	const n2 = grid.N2+1;

	vector<double> bound_Force_x;
	double bound_norm = 0.0;
	double new_x = 0.0;
	double new_u_bound = 0.0;
	double f1 = 0.0;

	bound_Force_x.resize(grid.NF);

	for (int i = 0; i < n1; ++i){
		for (int j = 0; j < n2; ++j){
			force_x[i][j] = 0.0;
		}
	}

	/*
	calculating force F for Lagrange
	in new_x and new_y calculated value of velocity in Lagrangian point. It calculates by using near points and discrete delta function
	*/

	for (auto& solid : iList){
		CreateMatrix(force_x_temp, n1, n2);
		for (int k = 0; k < grid.NF; ++k){

			new_x = 0.0;
			int i_max = 0;
			int i_min = 0;
			int j_max = 0;
			int j_min = 0;

			GetInfluenceArea(i_min, i_max, j_min, j_max, solid.Bound[0][k], solid.Bound[1][k], 3, grid);

			if (i_max >= n1){
				i_max = n1 - 1;
			}
			if (j_max >= n2){
				j_max = n2 - 1;
			}

			for (int i = i_min; i <= i_max; ++i){
				for (int j = j_min; j <= j_max; ++j){

					new_x += u[i][j] * DeltaFunction(i*grid.d_x - solid.Bound[0][k], (j - 0.5)*grid.d_y - solid.Bound[1][k],grid) * grid.d_x * grid.d_y;

				}
			}


			solid.Integral_x[k] += (new_x - solid.Uc[1]) * grid.d_t;
			bound_Force_x[k] = alpha_f * solid.Integral_x[k] + beta_f * (new_x - solid.Uc[1]);
		}

		//calculating force f for Euler points
		// spreading force by delta function
		for (int k = 0; k < grid.NF; ++k){

			int i_max = 0;
			int i_min = 0;

			int j_max = 0;
			int j_min = 0;

			GetInfluenceArea(i_min, i_max, j_min, j_max, solid.Bound[0][k], solid.Bound[1][k], 3,grid);

			if (i_max >= n1){
				i_max = n1 - 1;
			}
			if (j_max >= n2){
				j_max = n2 - 1;
			}

			for (int i = i_min; i <= i_max; ++i){
				for (int j = j_min; j <= j_max; ++j){

					force_x_temp[i][j] += bound_Force_x[k] * DeltaFunction(i*grid.d_x - solid.Bound[0][k], (j - 0.5)*grid.d_y - solid.Bound[1][k], grid) * solid.d_s * solid.d_s;

				}
			}
		}


		int i_max = 0;
		int i_min = n1;

		int j_max = 0;
		int j_min = n2;

		for (int k = 0; k < grid.NF; ++k){

			int i_max_temp = 0;
			int i_min_temp = 0;

			int j_max_temp = 0;
			int j_min_temp = 0;

			GetInfluenceArea(i_min_temp, i_max_temp, j_min_temp, j_max_temp, solid.Bound[0][k], solid.Bound[1][k], 3, grid);

			if (i_max_temp > i_max){
				i_max = i_max_temp;
			}
			if (i_min_temp < i_min){
				i_min = i_min_temp;
			}
			if (j_max_temp > j_max){
				j_max = j_max_temp;
			}
			if (j_min_temp < j_min){
				j_min = j_min_temp;
			}

		}

		if (i_max >= n1){
			i_max = n1 - 1;
		}
		if (j_max >= n2){
			j_max = n2 - 1;
		}

		Coeff = 0.0;

		for (int i = i_min; i <= i_max; ++i){
			for (int j = j_min; j <= j_max; ++j){


				force_x[i][j] += force_x_temp[i][j];
				//sum += force_x[i][j];
				Coeff += force_x_temp[i][j] * grid.d_x * grid.d_y;

			}

		}

		if (solid.moveSolid){

			solid.Uc[1] = solid.Uc[1] + (-Coeff * grid.d_t) / (M - M_PI * solid. r * solid.r);
		}
	}

	return 0;

}

double CalculateForce_Y(Matrix& force_y, list<Circle> &iList, Matrix& v, double r, double & Coeff, Grid grid, double alpha_f, double beta_f, double M){

	int const n1 = grid.N1 + 1;
	int	const n2 = grid.N2 ;
	vector<double> bound_Force_y;

	double bound_norm = 0.0;

	double new_y = 0.0;

	double new_v_bound = 0.0;
	CreateMatrix(force_y_temp, n1, n2);
	bound_Force_y.resize(grid.NF);

	for (int i = 0; i < n1; ++i){
		for (int j = 0; j < n2; ++j){
			force_y[i][j] = 0.0;
		}
	}
	/*
	for (int i = 0; i < n1; ++i){
		for (int j = 0; j < n2; ++j){
			force_y[i][j] = 0.0;
		}
	}
	*/

	/*calculating force F for Lagrange

	in new_x and new_y calculated value of velocity in Lagrangian point. It calculates by using near points and discrete delta function

	*/

	for (auto& solid : iList){
		/*
		for (int i = 0; i < n1; ++i){
			for (int j = 0; j < n2; ++j){
				force_y_temp[i][j] = 0.0;
			}
		}
		*/
		for (int k = 0; k < grid.NF; ++k){

			new_y = 0.0;

			int i_max = 0;
			int i_min = 0;

			int j_max = 0;
			int j_min = 0;

			GetInfluenceArea(i_min, i_max, j_min, j_max, solid.Bound[0][k], solid.Bound[1][k], 3,grid);

			if (i_max >= n1){
				i_max = n1 - 1;
			}
			if (j_max >= n2){
				j_max = n2 - 1;
			}

			for (int i = i_min; i <= i_max; ++i){
				for (int j = j_min; j <= j_max; ++j){


					new_y += v[i][j] * DeltaFunction((i - 0.5)*grid.d_x - solid.Bound[0][k], j*grid.d_y - solid.Bound[1][k],grid) * grid.d_x * grid.d_y;

				}
			}


			solid.Integral_y[k] += (new_y - solid.Uc[2]) * grid.d_t;

			bound_Force_y[k] = alpha_f * solid.Integral_y[k] + beta_f * (new_y - solid.Uc[2]);

		}



		//calculating force f for Euler points
		// spreading force by delta function
		for (int k = 0; k < grid.NF; ++k){

			int i_max = 0;
			int i_min = 0;

			int j_max = 0;
			int j_min = 0;

			GetInfluenceArea(i_min, i_max, j_min, j_max, solid.Bound[0][k], solid.Bound[1][k], 3,grid);

			if (i_max >= n1){
				i_max = n1 - 1;
			}
			if (j_max >= n2){
				j_max = n2 - 1;
			}

			for (int i = i_min; i <= i_max; ++i){
				for (int j = j_min; j <= j_max; ++j){

					force_y_temp[i][j] += bound_Force_y[k] * DeltaFunction((i - 0.5)*grid.d_x - solid.Bound[0][k], j*grid.d_y - solid.Bound[1][k], grid) * solid.d_s * solid.d_s;
				}
			}

		}

		int i_max = 0;
		int i_min = n1;

		int j_max = 0;
		int j_min = n2;

		for (int k = 0; k < grid.NF; ++k){

			int i_max_temp = 0;
			int i_min_temp = 0;

			int j_max_temp = 0;
			int j_min_temp = 0;

			GetInfluenceArea(i_min_temp, i_max_temp, j_min_temp, j_max_temp, solid.Bound[0][k], solid.Bound[1][k], 3,grid);

			if (i_max_temp > i_max){
				i_max = i_max_temp;
			}
			if (i_min_temp < i_min){
				i_min = i_min_temp;
			}
			if (j_max_temp > j_max){
				j_max = j_max_temp;
			}
			if (j_min_temp < j_min){
				j_min = j_min_temp;
			}
		}

		if (i_max >= n1){
			i_max = n1 - 1;
		}
		if (j_max >= n2){
			j_max = n2 - 1;
		}

		Coeff = 0.0;
		for (int i = i_min; i <= i_max; ++i){
			for (int j = j_min; j <= j_max; ++j){


				force_y[i][j] += force_y_temp[i][j];

				Coeff += force_y_temp[i][j] * grid.d_x * grid.d_y;

			}

		}

		if (solid.moveSolid){

			solid.Uc[2] = solid.Uc[2] + (-Coeff * grid.d_t) / (M - M_PI * solid.r * solid.r);
		}
	}


	return 0;

}
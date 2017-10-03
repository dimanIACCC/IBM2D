#include "BiCGStab.h"

void BiCGStab(Matrix &x, ublas::matrix<Template> &A, Matrix &b, Param par, Direction Dir){

	size_t n1 = x.size();
	size_t n2 = x[0].size();

	double eps = 0.0;
	double help_value = 0.0;

	int n = 0;

	CreateMatrix(r, n1, n2); // variables for BiCGStab ( Biconjugate gradient stabilized method )
	CreateMatrix(r_wave, n1, n2);
	CreateMatrix(v, n1, n2); // this is not velocity v.
	CreateMatrix(p, n1, n2); // this is not velocity p.
	CreateMatrix(s, n1, n2);
	CreateMatrix(t, n1, n2);
	


	double rho_prev = 1.0;
	double rho_current = 1.0;
	double alpha = 1.0;
	double beta = 0.0;
	double omega = 1.0;
	double r_norm = 0.0;
	double b_norm = 0.0;

	/* We want to solve Poisson equation for velocity.
	Au = b, where A Matrix of coefficents of leap-frog scheme applyied to poisson equation. and b right side
	preparation for CG
	Suppose in u first approximation ( in fact in u - velocity fromprevious step)
	r(0) = b - Au
	z(0) = r(0)
	in b_norm calculate Euclid norm of vector b*/

	// for( int i = 0; i < n1; ++i){
	// 	for( int j = 0; j < n2; ++j){

	// 		res(i,j) = 0.0;
	// 	}
	// }


	r = Operator_Ax(A, x, par, Dir);


	for (size_t j = 0; j < n2; ++j){
		for (size_t i = 0; i < n1; ++i){
			help_value = b[i][j] - r[i][j];

			r[i][j] = help_value;
			r_wave[i][j] = help_value;

			b_norm += pow(b[i][j], 2);
		}
	}

	b_norm = sqrt(b_norm);


	// iteration of method
	while (n < 10000 && b_norm != 0.0){
		r_norm = 0.0;

		// 1.

		rho_current = ScalarOperator(r_wave, r, n1, n2);


		// 2.

		beta = (rho_current / rho_prev) * (alpha / omega);

		// 3.

		for (size_t i = 0; i < n1; ++i){
			for (size_t j = 0; j < n2; ++j){

				p[i][j] = r[i][j] + beta * (p[i][j] - omega * v[i][j]);

			}
		}

		// 4.

		v = Operator_Ax(A, p, par, Dir);


		// 5.
		//cout << v;
		alpha = rho_current / ScalarOperator(r_wave, v, n1, n2);


		// 6.

		for (size_t j = 0; j < n2; ++j){
			for (size_t i = 0; i < n1; ++i){
				s[i][j] = r[i][j] - alpha * v[i][j];


			}
		}


		// 7. 

		t = Operator_Ax(A, s, par, Dir);

		// 8.

		omega = ScalarOperator(t, s, n1, n2) / ScalarOperator(t, t, n1, n2);


		// 9.
		// 10.

		for (size_t j = 0; j < n2; ++j){
			for (size_t i = 0; i < n1; ++i){

				// if ( i < M1 && j < M2 ){
				// 	continue;
				// }
				x[i][j] = x[i][j] + omega * s[i][j] + alpha * p[i][j];

				r[i][j] = s[i][j] - omega * t[i][j];

				r_norm += pow(r[i][j], 2);

			}
		}

		rho_prev = rho_current;


		//ApplyVelocityBoundary(u, v);

		// calculate error = ||r|| / ||b||
		eps = sqrt(r_norm) / b_norm;


		if (eps < 5e-10){

			//cout<<"vel ended at"<<n<<"  eps ="<<eps<<endl;
			break;

		}

		++n;


	}

	if (10000 == n){

		std::cout << "---------------- Ops. iteration exit in BiCGStab method (eps = " << eps << ") ----------------" << std::endl;
	}

}

double ScalarOperator(Matrix &a, Matrix &b, size_t n1, size_t n2){

	double result = 0;

	for (size_t i = 0; i < n1; ++i){
		for (size_t j = 0; j < n2; ++j){
			result += a[i][j] * b[i][j];
		}
	}

	return result;
}

#include "Calculate_press.h"

double Pressure_correction_solve_SOR(Matrix &delta_p, Matrix &rhs, Param par, int &N_out){

	double eps;
	double delta_p_max;

	double rxx1 = par.ldxdx;
	double ryy1 = par.ldydy;
	double rp   = 1.0 / (2.0 * (rxx1 + ryy1));

	double p_tmp, p_fix = 0;

	size_t n1 = delta_p.size();
	size_t n2 = delta_p[0].size();

	//Close to optimal value
	//See page 135 in "Khakimzyanov G.S. Cherny S.G. Method of computations.
	//Part 3. Numerical methods to solve problems for parabolic and elliptic
	//problems.Novosibirsk: NSU, 2007. 160 p.".
	double omega = 2. / (1. + sin(M_PI / std::min(n1, n2)));
	double omega1 = 1. - omega;

	for (int n = 0; n < par.N_Zeidel; n++) {
		eps = 0.0;
		delta_p_max = 0.0;

		for (size_t i = 1; i < n1-1; ++i){
			for (size_t j = 1; j < n2-1; ++j){

				p_tmp = rp * (rxx1 * (delta_p[i + 1][j] + delta_p[i - 1][j])
				            + ryy1 * (delta_p[i][j + 1] + delta_p[i][j - 1]) + rhs[i][j]);

				if (par.BC == periodical && i == n1-2 && j == 1) p_tmp = 0.0; // down corner
				
				if (fabs(p_tmp) > delta_p_max) {
					delta_p_max = fabs(p_tmp);
				}
				if (fabs(p_tmp - delta_p[i][j]) > eps){
					eps = fabs(p_tmp - delta_p[i][j]);
				}

				delta_p[i][j] = omega1*delta_p[i][j] + omega*p_tmp;

			}
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

		if (eps/delta_p_max < par.Zeidel_eps && n > 2){
			N_out = n;
			break;
		}

	}

	if (eps / delta_p_max > par.Zeidel_eps) {
		std::cout << "Zeidel has not converged, eps_p = " << eps << ",   delta_p_max = " << delta_p_max << std::endl;
	}

	return delta_p_max;

}


Matrix Pressure_RHS(Matrix &u, Matrix &v, Param par){
	double d = 0.0;

	CreateMatrix(result, par.N1_p, par.N2_p);


	for (size_t i = 1; i < par.N1_p - 1; ++i){
		for (size_t j = 1; j < par.N2_p - 1; ++j){

			double d = (1.0 / par.d_x) * (u[i + 1][j] - u[i][j])
			         + (1.0 / par.d_y) * (v[i][j + 1] - v[i][j]);

			result[i][j] = - 2. * d / par.d_t;

		}

	}
	result[0           ][0           ] = 0.0;
	result[par.N1_p - 1][0           ] = 0.0;
	result[0           ][par.N2_p - 1] = 0.0;
	result[par.N1_p - 1][par.N2_p - 1] = 0.0;


	for (size_t i = 1; i < par.N1_p - 1; ++i){
		result[i][0           ] = 0.0;
		result[i][par.N2_p - 1] = 0.0;
	}
	for (size_t j = 0; j < par.N2_p; ++j){
		result[0           ][j] = 0.0;
		result[par.N1_p - 1][j] = 0.0;
		if (par.BC == periodical) {
			result[0][j] = result[par.N1_p - 2][j];
			result[par.N1_p - 1][j] = result[1][j];
		}
	}

	return result;

}

// subroutine to solve Pressure correction equation
double Pressure_correction_solve(Matrix &delta_p, Matrix &rhs, Param par, double q, MKL_INT nx, MKL_INT ny, char* BCtype, int &N_DeltaP) {

	double dp_max;

	if (par.DeltaP_method == 0) { //SOR method
		dp_max = Pressure_correction_solve_SOR(delta_p, rhs, par, N_DeltaP);       // SOR method for solving Poisson equation
	}
	else if (par.DeltaP_method >= 1) { //MKL solver with first order approximation of boundary conditions

		Helmholtz_MKL(delta_p, rhs, q, par.d_x, par.d_y, 0, nx, 0, ny, BCtype);

		if (par.DeltaP_method == 2) { //Hybrid variant MKL + SOR with second order approximation of boundary conditions
			dp_max = Pressure_correction_solve_SOR(delta_p, rhs, par, N_DeltaP);       // SOR method for solving Poisson equation
		}
		dp_max = max(delta_p);
	}
	return dp_max;

} //Pressure_correction_solve



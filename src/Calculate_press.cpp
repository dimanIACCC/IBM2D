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
	for (size_t j = 1; j < par.N2_p - 1; ++j){
		result[0           ][j] = 0.0;
		result[par.N1_p - 1][j] = 0.0;
	}

	return result;

}

// subroutine to solve Pressure correction equation
double Pressure_correction_solve(Matrix &delta_p, Matrix &rhs, Param par, int &N_DeltaP) {

	if (par.DeltaP_method == 0) { //SOR method
		return Pressure_correction_solve_SOR(delta_p, rhs, par, N_DeltaP);       // SOR method for solving Poisson equation
	}
	else if (par.DeltaP_method >= 1) { //MKL solver with first order approximation of boundary conditions

		char* BCtype = "DDDD";

		if (par.BC == u_inflow || par.BC == u_infinity) BCtype = "NDNN";
		if (par.BC == periodical)                       BCtype = "PPNN";

		double q = 0.;

		MKL_INT nx = par.N1 + 1;
		MKL_INT ny = par.N2 + 1;

		//Boundaries
		double ax = 0.;
		double bx = par.d_x*nx;
		double ay = 0.;
		double by = par.d_y*ny;

		double *f_mkl = NULL, *bd_ax = NULL, *bd_bx = NULL, *bd_ay = NULL, *bd_by = NULL;
		f_mkl = (double*)mkl_malloc((nx + 1)*(ny + 1) * sizeof(double), 64);
		bd_ax = (double*)mkl_malloc((ny + 1) * sizeof(double), 64);
		bd_bx = (double*)mkl_malloc((ny + 1) * sizeof(double), 64);
		bd_ay = (double*)mkl_malloc((nx + 1) * sizeof(double), 64);
		bd_by = (double*)mkl_malloc((nx + 1) * sizeof(double), 64);

		Matrix_to_DoubleArray(rhs, f_mkl);

		for (MKL_INT iy = 0; iy <= ny; iy++) {
			bd_ax[iy] = 0.;
			bd_bx[iy] = 0.;
		}
		for (MKL_INT ix = 0; ix <= nx; ix++) {
			bd_ay[ix] = 0.;
			bd_by[ix] = 0.;
		}

		Helmholtz_MKL(f_mkl, ax, bx, ay, by, bd_ax, bd_bx, bd_ay, bd_by, nx, ny, BCtype, q, par.d_x, par.d_y);

		DoubleArray_to_Matrix(f_mkl, delta_p);

		mkl_free(bd_ax);
		mkl_free(bd_bx);
		mkl_free(bd_ay);
		mkl_free(bd_by);
		mkl_free(f_mkl);
		MKL_Free_Buffers();

		if (par.DeltaP_method == 2) { //Hybrid variant MKL + SOR with second order approximation of boundary conditions
			return Pressure_correction_solve_SOR(delta_p, rhs, par, N_DeltaP);       // SOR method for solving Poisson equation
		}
		return max(delta_p);
	}

} //Pressure_correction_solve

void Matrix_to_DoubleArray(Matrix &M, double* D) {

	int Nx = M.size();
	int Ny = M[0].size();

	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			D[i + j*Nx] = M[i][j];
		}
	}
}

void DoubleArray_to_Matrix(double* D, Matrix &M) {

	int Nx = M.size();
	int Ny = M[0].size();

	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			M[i][j] = D[i + j*Nx];
		}
	}
}

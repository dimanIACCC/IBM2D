
// Module to implement subroutines to solve boundary value problem for Helmholtz
// equation in the rectangular domain:
// - div(grad p) + q*p = rhs
// The equation is approximated by the second order difference schema:
// (q*I - Lambda_xx - Lambda_yy)p = rhs
// and obtained SLAE is solved by the Successive over-relaxation (SOR) method
// To start main subroutine 6 arrayes should be set:
//   p(0:nx+1,0:ny+1) - initial guess
//   ax(0:ny+1,1:3), bx(0:ny+1,1:3) - boundary conditions (left and right)
//   ay(0:nx+1,1:3), by(0:nx+1,1:3) - boundary conditions (bottom and top)
//   rhs(0:nx+1,0:ny+1) - right hand side function
// The boundary conditions realize in the following manner:
//   p(0,i)    + ax(i,1)*p(1,i)  + ax(i,2)*p(2,i)    = ax(i,3)
//   p(nx+1,i) + bx(i,1)*p(nx,i) + bx(i,2)*p(nx-1,i) = bx(i,3)
//   p(i,0)    + ay(i,1)*p(i,1)  + ay(i,2)*p(i,2)    = ay(i,3)
//   p(i,ny+1) + by(i,1)*p(i,ny) + by(i,2)*p(i,ny-1) = by(i,3)

///////////////////////////////////////////////////////////////////////////////
// Program: ibm2d_cyl -- SIMPLE + IBM method for incomressible viscous Navier -
// Stockes equations in 2D statement for flow past a cylinder
// Copyright (C): Denis V. Esipov (2019)
///////////////////////////////////////////////////////////////////////////////

#include "helmholtz.h"

DFTI_DESCRIPTOR_HANDLE xhandle = 0;

void Helmholtz_MKL(double*f, double &ax, double &bx, double &ay, double &by, 
                   double*bd_ax, double*bd_bx, double*bd_ay, double*bd_by, MKL_INT &nx, MKL_INT &ny, char *BCtype, double q, double hx, double hy)
{
	MKL_INT stat, ipar[128];
	double *dpar = NULL;

	dpar = (double*)mkl_malloc((13 * nx / 2 + 7) * sizeof(double), 64);

	d_init_Helmholtz_2D(&ax, &bx, &ay, &by, &nx, &ny, BCtype, &q, ipar, dpar, &stat);
	ipar[2] = 0;      // Disable warning for Neumann BC
	d_commit_Helmholtz_2D(f, bd_ax, bd_bx, bd_ay, bd_by, &xhandle, ipar, dpar, &stat);
	d_Helmholtz_2D(f, bd_ax, bd_bx, bd_ay, bd_by, &xhandle, ipar, dpar, &stat);
	free_Helmholtz_2D(&xhandle, ipar, &stat);

	mkl_free(dpar);
	//mkl_free(f);

	mkl_free(bd_ax);
	mkl_free(bd_bx);
	mkl_free(bd_ay);
	mkl_free(bd_by);
	MKL_Free_Buffers();
}

//subroutine to solve Helmholtz problem by the Successive over - relaxation(SOR) method
void Helmholtz_SOR(Matrix &p, Matrix &ax, Matrix &bx, Matrix &ay, Matrix &by, Matrix &rhs,
                         double q, int nx, int ny, double hx, double hy, int nit, double eps, double omega) {

	double omega1 = 1. - omega;

	//Compute ratios
	double rxx1 = 1. / (hx*hx);
	double ryy1 = 1. / (hy*hy);
	double rp = 1. / (2.*(rxx1 + ryy1) + q);

	double st;
	int k;
	//Main loop of SOR iterations
	for (k = 0; k <= nit; k++) {
		//Apply boundary conditions
		for (int i = 0; i <= ny + 1; i++) {
			p[0][i] = ax[i][3] - ax[i][1] * p[1][i] - ax[i][2] * p[2][i];
			p[nx + 1][i] = bx[i][3] - bx[i][1] * p[nx][i] - bx[i][2] * p[nx - 1][i];
		}
		for (int i = 0; i <= nx + 1; i++) {
			p[i][0] = ay[i][3] - ay[i][1] * p[i][1] - ay[i][2] * p[i][2];
			p[i][ny + 1] = by[i][3] - by[i][1] * p[i][ny] - by[i][2] * p[i][ny - 1];
		}
		//One SOR iteration
		st = 0.;
		double pmax = 0.;

		for (int i = 1; i <= nx; i++) {
			for (int j = 1; j <= ny; j++) {
				double c = (rxx1*(p[i - 1][j] + p[i + 1][j]) + ryy1*(p[i][j - 1] + p[i][j + 1]) + rhs[i][j])*rp;
				st = std::max(st, abs(p[i][j] - c));
				p[i][j] = omega1*p[i][j] + omega*c;
				pmax = std::max(pmax, abs(p[i][j]));
			}
		}
		if (pmax < eps) pmax = 1.;   //To avoid divide by zero

		//Check stop condition
		if (st*omega / pmax < eps)  return;

		//Write history of convergence
		//std::cout << "SOR: k = " << k << " st*omega/pmax = " << st*omega/pmax << std::endl;

	}

	// Write message on the screen if at the same time iteration process is stopped and stop condition is not fullfilled
	if (k == nit + 1) {
		std::cout << "Stop condition in the Successive over-relaxation (SOR) method is not fulfilled!" << std::endl;
		std::cout << "||p^" << nit << " - p^" << "nit-1" << st << std::endl;
	}
}

void Helmholtz_MKL(Matrix &f, Matrix &rhs, double q, double hx, double hy, MKL_INT nx_l, MKL_INT nx_r, MKL_INT ny_l, MKL_INT ny_r, char* BCtype) {

	MKL_INT nx = nx_r - nx_l;
	MKL_INT ny = ny_r - ny_l;

	//Boundaries
	double ax = 0.;
	double bx = hx*nx;
	double ay = 0.;
	double by = hy*ny;

	double *f_mkl = NULL, *bd_ax = NULL, *bd_bx = NULL, *bd_ay = NULL, *bd_by = NULL;
	f_mkl = (double*)mkl_malloc((nx + 1)*(ny + 1) * sizeof(double), 64);
	bd_ax = (double*)mkl_malloc((ny + 1) * sizeof(double), 64);
	bd_bx = (double*)mkl_malloc((ny + 1) * sizeof(double), 64);
	bd_ay = (double*)mkl_malloc((nx + 1) * sizeof(double), 64);
	bd_by = (double*)mkl_malloc((nx + 1) * sizeof(double), 64);

	Matrix_to_DoubleArray(rhs, f_mkl, 0, nx, 0, ny);

	for (MKL_INT iy = 0; iy <= ny; iy++) {
		bd_ax[iy] = 0.;
		bd_bx[iy] = 0.;
	}
	for (MKL_INT ix = 0; ix <= nx; ix++) {
		bd_ay[ix] = 0.;
		bd_by[ix] = 0.;
	}

	Helmholtz_MKL(f_mkl, ax, bx, ay, by, bd_ax, bd_bx, bd_ay, bd_by, nx, ny, BCtype, q, hx, hy);

	DoubleArray_to_Matrix(f_mkl, f, 0, nx, 0, ny);

	//Output_P(f, "dp", 555, par);

	mkl_free(f_mkl);
	MKL_Free_Buffers();

}

void Matrix_to_DoubleArray(Matrix &M, double* D, int Nx_l, int Nx_r, int Ny_l, int Ny_r) {
	for (int i = Nx_l; i <= Nx_r; i++) {
		for (int j = Ny_l; j <= Ny_r; j++) {
			D[(i - Nx_l) + (j - Ny_l)*(Nx_r - Nx_l + 1)] = M[i][j];
		}
	}

	//Output_2DArray(D, Nx - 1, Ny, "Result/", "Array", 555);
}

void DoubleArray_to_Matrix(double* D, Matrix &M, int Nx_l, int Nx_r, int Ny_l, int Ny_r) {
	for (int i = Nx_l; i <= Nx_r; i++) {
		for (int j = Ny_l; j <= Ny_r; j++) {
			M[i][j] = D[(i - Nx_l) + (j - Ny_l)*(Nx_r - Nx_l + 1)];
		}
	}
}
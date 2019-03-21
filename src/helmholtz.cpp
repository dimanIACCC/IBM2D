
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

void Helmholtz_MKL(double*f, double &ax, double &bx, double &ay, double &by, 
                   double*bd_ax, double*bd_bx, double*bd_ay, double*bd_by, MKL_INT &nx, MKL_INT &ny, char *BCtype, double q, double hx, double hy)
{
	MKL_INT stat, ipar[128];
	double *dpar = NULL;

	DFTI_DESCRIPTOR_HANDLE xhandle = 0;

	dpar = (double*)mkl_malloc((13 * nx / 2 + 7) * sizeof(double), 64);

	for (MKL_INT i = 0; i<128; i++) ipar[i] = 0;

	/* Initializing simple data structures of Poisson Library for 2D Helmholtz Solver */
	d_init_Helmholtz_2D(&ax, &bx, &ay, &by, &nx, &ny, BCtype, &q, ipar, dpar, &stat);
	ipar[2] = 0;      // Disable warning for Neumann BC
	if (stat != 0) goto end;
	/* Initializing complex data structures of Poisson Library for 2D Helmholtz Solver */
	d_commit_Helmholtz_2D(f, bd_ax, bd_bx, bd_ay, bd_by, &xhandle, ipar, dpar, &stat);
	if (stat != 0) goto end;
	/* Computing the approximate solution of 2D Helmholtz problem */
	d_Helmholtz_2D(f, bd_ax, bd_bx, bd_ay, bd_by, &xhandle, ipar, dpar, &stat);
	if (stat != 0) goto end;
	/* Cleaning the memory used by xhandle */
	free_Helmholtz_2D(&xhandle, ipar, &stat);
	if (stat != 0) goto end;
	/* Now we can use xhandle to solve another 2D Helmholtz problem */

end:
	/* Free Intel(R) MKL memory if any was allocated */
	mkl_free(dpar);
	//mkl_free(f);

	mkl_free(bd_ax);
	mkl_free(bd_bx);
	mkl_free(bd_ay);
	mkl_free(bd_by);
	MKL_Free_Buffers();
}

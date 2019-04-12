#pragma once

#include "mkl_service.h"
#include "mkl_poisson.h"
#include "GeomVec.h"
#include "Matrix.h"

void Helmholtz_MKL(double *f, double &ax, double &bx, double &ay, double &by, 
                   double*bd_ax, double*bd_bx, double*bd_ay, double*bd_by, MKL_INT &nx, MKL_INT &ny, char *BCtype, double q, double hx, double hy);
void Helmholtz_SOR(Matrix &p, Matrix &ax, Matrix &bx, Matrix &ay, Matrix &by, Matrix &rhs,
	                     double q, int nx, int ny, double hx, double hy, int nit, double eps, double omega = 1.);
void Helmholtz_MKL(Matrix &f, Matrix &rhs, double q, double hx, double hy, MKL_INT nx_l, MKL_INT nx_r, MKL_INT ny_l, MKL_INT ny_r, char* BCtype);
void Matrix_to_DoubleArray(Matrix &M, double* D, int nx_l, int nx_r, int ny_l, int ny_r);
void DoubleArray_to_Matrix(double* D, Matrix &M, int nx_l, int nx_r, int ny_l, int ny_r);

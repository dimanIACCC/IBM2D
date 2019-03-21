#pragma once

#include "mkl_service.h"
#include "mkl_poisson.h"


void Helmholtz_MKL(double *f, double &ax, double &bx, double &ay, double &by, 
                   double*bd_ax, double*bd_bx, double*bd_ay, double*bd_by, MKL_INT &nx, MKL_INT &ny, char *BCtype, double q, double hx, double hy);

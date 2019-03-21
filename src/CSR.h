#pragma once

#include "stdafx.h"
#include <mkl_pardiso.h>
#include <mkl_types.h>
#include <mkl_spblas.h>
class TCSRMatrix {
	int nnz, n;                                                   //Number of nonzeros in M and dimension of M
	double val[99999];                                            //Value of elements
	int col_ind[99999], row_ptr[99999];                           //Column indices and row pointers
};

MKL_INT solve_pardiso();
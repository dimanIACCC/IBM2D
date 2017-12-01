#pragma once

#include "GeomVec.h"
#include "Parameters.h"
#include "Matrix.h"
#include "Output.h"


double Calculate_Press_correction(Matrix& delta_p, Matrix &b_p, Param par);               // solve the Poisson equation:  Laplace delta_p = b_p

Matrix Calculate_Press_Right(Matrix& u, Matrix& v, Param par);                            // right-hand part of the Poisson equation


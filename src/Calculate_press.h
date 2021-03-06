#pragma once

#include "GeomVec.h"
#include "Parameters.h"
#include "Matrix.h"
#include "Output.h"
#include "helmholtz.h"

double Pressure_correction_solve    (Matrix &p, Matrix &rhs, Param par);
double Pressure_correction_solve_SOR(Matrix &delta_p, Matrix &rhs, Param par, int &N_DeltaP);    // solve the Poisson equation:  -Laplace delta_p = rhs

Matrix Pressure_RHS(Matrix& u, Matrix& v, Param par);                            // right-hand part of the Poisson equation

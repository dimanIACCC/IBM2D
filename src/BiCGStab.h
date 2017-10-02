#pragma once

#include "PredictVel.h"

/* We want to solve Poisson equation for velocity.
Au = b, where A Matrix of coefficents of leap-frog scheme applyied to poisson equation. and b right side
preparation for CG
Suppose in u first approximation ( in fact in u - velocity fromprevious step)
r(0) = b - Au
z(0) = r(0)
in b_norm calculate Euclid norm of vector b*/
void BiCGStab(Matrix& res, ublas::matrix<Template> &A, Matrix &b, Param par, Direction Dir);
double ScalarOperator(Matrix &a, Matrix &b, int const n1, int const n2);

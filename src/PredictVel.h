#pragma once

#include "Parameters.h"
#include "GeomVec.h"
#include "Matrix.h"

void Calculate_A(Template &A, Param par, double Re, Direction Dir);
Matrix Operator_Ax(Template &A, Matrix &x, Param par, Direction Dir);

Matrix CalculateB(Matrix &u_n, Matrix &v_n, Matrix &u_prev, Matrix &v_prev, Matrix &p, Matrix &p_new, Param par, Direction Dir);
double advective_term(Matrix &u, Matrix &v, size_t i, size_t j, double d_x, double d_y, Direction Dir, size_t N1, size_t N2);
double diffusion_term(Matrix &u, size_t i, size_t j, double d_xx, double d_yy, Direction Dir, size_t N1, size_t N2);

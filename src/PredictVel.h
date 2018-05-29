#pragma once

#include "Parameters.h"
#include "GeomVec.h"
#include "Matrix.h"
#include "Output.h"

void Calculate_A(Template &A, Param par, double Re, Direction Dir);
Matrix Operator_Ax(Template &A, Matrix &x, Param par, Direction Dir);

Matrix CalculateB(Matrix &u_n, Matrix &v_n, Matrix &u_prev, Matrix &v_prev, Matrix &p, Matrix &p_new, Param par, Direction Dir);
void Boundary_Conditions(Matrix &u, Param par, Direction Dir, int N_step);
double advective_term(Matrix &u, Matrix &v, size_t i, size_t j, double d_x, double d_y, Direction Dir);

void TaylorGreen_exact(Matrix &u, Matrix &v, Matrix &p, Param par, double time);
double Taylor_Green_u(GeomVec x, double k1, double k2, double time_exp);
double Taylor_Green_v(GeomVec x, double k1, double k2, double time_exp);
double Taylor_Green_p(GeomVec x, double k1, double k2, double time_exp2);

#pragma once

#include "Parameters.h"
#include "GeomVec.h"
#include "Matrix.h"
#include "Output.h"
#include "Exact_solutions.h"

void Calculate_A(Template &A, Param par, double Re, Direction Dir);
Matrix Operator_Ax(Template &A, Matrix &x, Param par, Direction Dir);

Matrix CalculateB(Matrix &u_n, Matrix &v_n, Matrix &u_prev, Matrix &v_prev, Matrix &p, Matrix &p_new, Param par, Direction Dir);
void Boundary_Conditions(Matrix &u, Param par, Direction Dir, int N_step);
double advective_term(Matrix &u, Matrix &v, size_t i, size_t j, double d_x, double d_y, Direction Dir);

void Taylor_Green_exact(Matrix &u, Matrix &v, Matrix &p, Param par, double time);
void Lamb_Oseen_exact(Matrix &uv, Direction Dir, Param par, double time, bool boundary);
void Lamb_Oseen_exact_p(Matrix &p, Param par, double time);

#pragma once

#include "Parameters.h"
#include "GeomVec.h"
#include "Matrix.h"
#include "Output.h"
#include "Boundary_initial_conditions.h"

void Calculate_A(Template &A, Param par, double Re);
Matrix Operator_Ax(Template &A, Matrix &x, Param par, Direction Dir);

Matrix CalculateB(Matrix &u_n, Matrix &v_n, Matrix &u_prev, Matrix &v_prev, Matrix &p, Matrix &p_new, Matrix &F, Param &par, Direction Dir);
double advective_term(Matrix &ul, Matrix &vl, Matrix &ur, Matrix &vr, size_t i, size_t j, double d_x, double d_y, Direction Dir);
void Output_eq_terms(std::string filename, int n, Matrix &u_n, Matrix &v_n, Matrix &u_s, Matrix &v_s, Matrix &p, Matrix &p_new, Matrix &F, Param par, Direction Dir);

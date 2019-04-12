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
void make_uv_RHS(Matrix &rhsu, Matrix &rhsv, Matrix &u0, Matrix &v0, Matrix &u1, Matrix &v1, Matrix &p0, double fu, double fv,
                 int nx, int ny, double hx, double hy, double tau, double Re);
void predict_uv(Matrix &u0, Matrix &v0, Matrix &u, Matrix &v,
                Matrix &p, Matrix &rhsu, Matrix &rhsv, int &nx, int &ny, double &hx, double &hy, double &tau, double &Re, Param &par);

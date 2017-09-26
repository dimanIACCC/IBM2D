#pragma once

#include "Parameters.h"
#include "GeomVec.h"
#include "Matrix.h"

void Calculate_A(ublas::matrix<Template> &A, Param par, double Re, Direction Dir);
Matrix Operator_Ax(ublas::matrix<Template> &A, Matrix &x, Param par, Direction Dir);

Matrix CalculateB(Matrix &u_n, Matrix &v_n, Matrix &u_prev, Matrix &v_prev, Matrix &p, Matrix &force, Param par, Direction Dir);
void B_coefficients(size_t ij, size_t N, double &k, double &kp, double &k0, double &km);
double advective_term(Matrix &u, Matrix &v, size_t i, size_t j, double d_x, double d_y, double k, Direction Dir);
double diffusion_term(Matrix &u, size_t i, size_t j, double d_xx, double d_yy, double kp, double k0, double km, Direction Dir);

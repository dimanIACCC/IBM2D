#pragma once

#include "Parameters.h"
#include "GeomVec.h"

void ExplicPredVel(Matrix& U_predict, Matrix& V_predict, Matrix& U_n, Matrix& V_n, Matrix& P, Matrix& Force_x, Matrix& Force_y, Param par);

void Calculate_A_u(Matrix A[5], Param par, double Re);
void Calculate_A_v(Matrix A[5], Param par, double Re);
Matrix Operator_Ax(Matrix A[5], Matrix &x, int n1, int n2, Param g,bool OverFlow);

Matrix CalculateB_u(Matrix &u_n, Matrix &v_n, Matrix &u_prev, Matrix &v_prev, Matrix &p, Matrix &force, Param par);
Matrix CalculateB_v(Matrix &u_n, Matrix &v_n, Matrix &u_prev, Matrix &v_prev, Matrix &p, Matrix &force, Param par);
void B_coefficients(int ij, int N, double &k, double &kp, double &k0, double &km);
double advective_term_u(Matrix &u, Matrix &v, int i, int j, double d_x, double d_y, double k);
double advective_term_v(Matrix &v, Matrix &u, int j, int i, double d_y, double d_x, double k);
double diffusion_term_u(Matrix &u, int i, int j, double d_xx, double d_yy, double kp, double k0, double km);
double diffusion_term_v(Matrix &v, int j, int i, double d_yy, double d_xx, double kp, double k0, double km);

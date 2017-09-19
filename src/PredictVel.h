#pragma once

#include "Parameters.h"
#include "GeomVec.h"

void ExplicPredVel(Matrix& U_predict, Matrix& V_predict, Matrix& U_n, Matrix& V_n, Matrix& P, Matrix& Force_x, Matrix& Force_y, Param par);

void Calculate_A_u(Matrix A[5], Param par, double Re);
void Calculate_A_v(Matrix A[5], Param par, double Re);
Matrix Operator_Ax(Matrix A[5], Matrix &x, int n1, int n2, Param g,bool OverFlow);

Matrix CalculateB_u(Matrix &u_n, Matrix &v_n, Matrix &u_prev, Matrix &v_prev, Matrix &p, Matrix &force, Param par);
Matrix CalculateB_v(Matrix &u_n, Matrix &v_n, Matrix &u_prev, Matrix &v_prev, Matrix &p, Matrix &force, Param par);

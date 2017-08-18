#pragma once

#include "stdafx.h"
#include "Grid.h"
#include "SolidBody.h"

using namespace std;

double DeltaFunction(double x, double y, Grid grid);
double FunctionD(double r);
void GetInfluenceArea(int& i_min, int& i_max, int& j_min, int& j_max, double x, double y, int size, Grid grid);

double CalculateForce_X(Matrix& force_x, list<Circle> &iList, Matrix& u, double r, double & Coeff, Grid grid, double alpha_f, double beta_f, double M);
double CalculateForce_Y(Matrix& force_y, list<Circle> &iList, Matrix& v, double r, double & Coeff, Grid grid, double alpha_f, double beta_f, double M);
Matrix Calculate_F_real_x(Matrix& u_n, Matrix& v_n, Matrix& u_prev, Matrix& p, Grid g, double Re);
Matrix Calculate_F_real_y(Matrix& u_n, Matrix& v_n, Matrix& v_prev, Matrix& p, Grid g, double Re);
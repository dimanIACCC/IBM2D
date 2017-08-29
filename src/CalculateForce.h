#pragma once

#include "stdafx.h"
#include "Parameters.h"
#include "SolidBody.h"

using namespace std;

double DeltaFunction(double x, double y, Grid grid);
double FunctionD(double r);
void GetInfluenceArea(int& i_min, int& i_max, int& j_min, int& j_max, double x, double y, int size, Grid grid);

double CalculateForce(Matrix& force_x, Matrix& force_y, list<Circle> &iList, Matrix& u, Matrix& v, Grid grid);
Matrix Calculate_F_real_x(Matrix& u_n, Matrix& v_n, Matrix& u_prev, Matrix& p, Grid g, double Re);
Matrix Calculate_F_real_y(Matrix& u_n, Matrix& v_n, Matrix& v_prev, Matrix& p, Grid g, double Re);

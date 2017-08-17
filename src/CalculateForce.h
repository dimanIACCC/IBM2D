#pragma once

#include "stdafx.h"
#include "Grid.h"
#include "SolidBody.h"

using namespace std;

double DeltaFunction(double x, double y, Grid grid);
double FunctionD(double r);
void GetInfluenceArea(int& i_min, int& i_max, int& j_min, int& j_max, double x, double y, int size, Grid grid);

double CalculateForce(Matrix& force_x, Matrix& force_y, list<Circle> &iList, Matrix& u, Matrix& v, double r, double & CoeffX, double & CoeffY, Grid grid, double alpha_f, double beta_f, double M);

#pragma once

#include "Parameters.h"
#include "SolidBody.h"
#include "Matrix.h"

double DeltaFunction(double x, double y, Param par);
double FunctionD(double r);
void GetInfluenceArea(int& i_min, int& i_max, int& j_min, int& j_max, double x, double y, int size, Param par);

void CalculateForce(Matrix& force_x, Matrix& force_y, std::list<Circle> &iList, Matrix& u, Matrix& v, Param par);

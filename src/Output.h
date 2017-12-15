#pragma once

#include "Parameters.h"
#include "SolidBody.h"
#include "Matrix.h"


void OutputPressure  (Matrix data, int n, std::list<Circle> iList, Param par);
void OutputVelocity_U(Matrix data, int n, std::list<Circle> iList, Param par);
void OutputVelocity_V(Matrix data, int n, std::list<Circle> iList, Param par);
void Output(Matrix p, Matrix u, Matrix v, Matrix Fx, Matrix Fy, int n, std::list<Circle> iList, Param par);
void Output_dp(Matrix dp, int n, Param par);
void Output_c (Matrix c,  int n, Param par);

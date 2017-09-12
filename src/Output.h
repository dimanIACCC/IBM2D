#pragma once

#include "Parameters.h"
#include "SolidBody.h"


void OutputPressure  (Matrix data, int n, std::list<Circle> iList, Param par);
void OutputVelocity_U(Matrix data, int n, std::list<Circle> iList, Param par);
void OutputVelocity_V(Matrix data, int n, std::list<Circle> iList, Param par);
void Output(Matrix p, Matrix u, Matrix v, int n, std::list<Circle> iList, Param par);

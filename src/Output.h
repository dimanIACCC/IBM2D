#pragma once

#include "stdafx.h"
#include "Parameters.h"
#include "SolidBody.h"
using namespace std;

void OutputPressure  (Matrix data, int n, list<Circle> iList, Param par);
void OutputVelocity_U(Matrix data, int n, list<Circle> iList, Param par);
void OutputVelocity_V(Matrix data, int n, list<Circle> iList, Param par);
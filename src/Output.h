#pragma once

#include "stdafx.h"
#include "Grid.h"
#include "SolidBody.h"
using namespace std;

void OutputPressure  (Matrix data, int n, list<Circle> iList, Grid grid);
void OutputVelocity_U(Matrix data, int n, list<Circle> iList, Grid grid);
void OutputVelocity_V(Matrix data, int n, list<Circle> iList, Grid grid);
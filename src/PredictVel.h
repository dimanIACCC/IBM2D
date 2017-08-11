#pragma once

#include "stdafx.h"
#include "Grid.h"

using namespace std;

void ExplicPredVel(Matrix& U_predict, Matrix& V_predict, Matrix& U_n, Matrix& V_n, int const n1, int const n2, Matrix operator_A[5], Matrix &b, Grid grid);
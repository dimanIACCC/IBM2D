#pragma once

#include "stdafx.h"
#include "Grid.h"
using namespace std;






void Calculate_A_u(Matrix A[5], Grid grid, double Re);
void Calculate_A_v(Matrix A[5], Grid grid, double Re);
Matrix Operator_Ax(Matrix A[5], Matrix &x, int n1, int n2, Grid grid);
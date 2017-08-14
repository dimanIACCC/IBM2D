#pragma once

#include "stdafx.h"
#include "Grid.h"

using namespace std;

void ExplicPredVel(Matrix& U_predict, Matrix& V_predict, Matrix& U_n, Matrix& V_n, Matrix& P, Matrix& Force_x, Matrix& Force_y, Grid grid);
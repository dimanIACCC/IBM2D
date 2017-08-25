#pragma once

#include "stdafx.h"
#include "Grid.h"
using namespace std;


double Calculate_Press_correction(Matrix& delta_p, Matrix &b_p, Grid grid,bool OverFlow);
Matrix Calculate_Press_Right(Matrix& u, Matrix& v, Grid grid);
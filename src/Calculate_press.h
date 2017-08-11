#pragma once

#include "stdafx.h"
#include "Grid.h"
using namespace std;


double Calculate_Press_correction(Matrix& delta_p, Matrix &b_p, int const N_Zeidel, double const Zeidel_eps, Grid grid);
Matrix Calculate_Press_Right(Matrix& u, Matrix& v, Grid grid);
#pragma once

#include "stdafx.h"
#include "Parameters.h"
using namespace std;


double Calculate_Press_correction(Matrix& delta_p, Matrix &b_p, Param par,bool OverFlow);
Matrix Calculate_Press_Right(Matrix& u, Matrix& v, Param par);
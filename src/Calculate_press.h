#pragma once

#include "GeomVec.h"
#include "Parameters.h"

double Calculate_Press_correction(Matrix& delta_p, Matrix &b_p, Param par, std::ostream& log, bool OverFlow);
Matrix Calculate_Press_Right(Matrix& u, Matrix& v, Param par);

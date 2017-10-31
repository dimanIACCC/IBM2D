#pragma once

#include "GeomVec.h"
#include "Parameters.h"
#include "Matrix.h"
#include "Output.h"

<<<<<<< HEAD
double Calculate_Press_correction(Matrix& delta_p, Matrix &b_p, Param par, std::ostream& log, bool OverFlow);
Matrix Calculate_Press_Right(Matrix& u, Matrix& v, Param par);
=======
double Calculate_Press_correction(Matrix& delta_p, Matrix &b_p, Param par);               // solve the Poisson equation:  Laplace delta_p = b_p
Matrix Calculate_Press_Right(Matrix& u, Matrix& v, Matrix& Fx, Matrix& Fy, Param par);    // right-hand part of the Poisson equation
>>>>>>> refs/remotes/origin/Kuranakov

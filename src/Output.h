#pragma once

#include "Parameters.h"
#include "SolidBody.h"
#include <boost/filesystem.hpp>
#include "Matrix.h"

namespace fs = boost::filesystem;


void CreateDirectory(fs::path directory);
void SetLog(std::ostream &log, Param par);
void PushLog(std::ostream &log, int n, double eps_u, double eps_v);

void OutputPressure  (Matrix data, int n, std::list<Circle> iList, Param par);
void OutputVelocity_U(Matrix data, int n, std::list<Circle> iList, Param par);
void OutputVelocity_V(Matrix data, int n, std::list<Circle> iList, Param par);
void Output(Matrix p, Matrix u, Matrix v, Matrix Fx, Matrix Fy, int n, std::list<Circle> iList, Param par);
void Output_dp(Matrix dp, int n, Param par);
void Output_c (Matrix c,  int n, Param par);


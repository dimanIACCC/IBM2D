#pragma once

#include "Parameters.h"
#include "SolidBody.h"
#include <boost/filesystem.hpp>
#include "Matrix.h"

namespace fs = boost::filesystem;

void MakeResultDir(fs::path dir_Result);
void SetLog(std::ostream &log, Param par);
void PushLog(std::ostream &log, int n, double eps_u, double eps_v);
void OutputPressure  (Matrix data, int n, std::list<Circle> iList, Param par, std::string WorkDir);
void OutputVelocity_U(Matrix data, int n, std::list<Circle> iList, Param par, std::string WorkDir);
void OutputVelocity_V(Matrix data, int n, std::list<Circle> iList, Param par, std::string WorkDir);
void Output(Matrix p, Matrix u, Matrix v, int n, std::list<Circle> iList, Param par, std::string WorkDir);
void Output_dp(Matrix dp, int n, Param par, std::string WorkDir);

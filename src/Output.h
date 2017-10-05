#pragma once

#include "Parameters.h"
#include "SolidBody.h"
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

void OutputPressure  (Matrix data, int n, std::list<Circle> iList, Param par);
void OutputVelocity_U(Matrix data, int n, std::list<Circle> iList, Param par);
void OutputVelocity_V(Matrix data, int n, std::list<Circle> iList, Param par);
void Output(Matrix p, Matrix u, Matrix v, int n, std::list<Circle> iList, Param par, fs::path ResultFolder);
void MakeResultDir(fs::path dir_Result);
void SetLog(std::ostream &log, Param par);
void PushLog(std::ostream &log, int n, double eps_u, double eps_v);
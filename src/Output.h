#pragma once

#include "Parameters.h"
#include "SolidBody.h"
#include <boost/filesystem.hpp>
#include "Matrix.h"

namespace fs = boost::filesystem;


void CreateDirectory(fs::path directory);
void SetLog(std::ostream &log, Param par);
void PushLog(std::ostream &log, int n, double eps_u, double eps_v);

void Output_U(Matrix data, std::string filename, int n, Param par);
void Output_V(Matrix data, std::string filename, int n, Param par);
void Output(Matrix p, Matrix u, Matrix v, Matrix Fx, Matrix Fy, int n, std::list<Circle> iList, Param par);
void Output_P(Matrix dp, std::string filename, int n, Param par);
void Output_c(Matrix c , std::string filename, int n, Param par);
void MakeHibernationFile(int n, Param& par, std::list<Circle>& solidList, Matrix& U_n, Matrix& V_n, Matrix& P);

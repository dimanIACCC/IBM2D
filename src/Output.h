#pragma once

#include "Parameters.h"
#include "SolidBody.h"
#include <boost/filesystem.hpp>
#include "Matrix.h"

namespace fs = boost::filesystem;


void CreateDir(fs::path directory);

void history_init(std::string WorkDir, std::string file, boundary_conditions BC);
void history_log(std::string WorkDir, std::string file, double t, double var1, double var2, double var3);

void Output_U(Matrix data, std::string filename, int n, Param par);
void Output_V(Matrix data, std::string filename, int n, Param par);
void Output(Matrix p, Matrix u, Matrix v, Matrix Fx, Matrix Fy, int n, std::list<Circle> iList, Param par);
void Output_P(Matrix dp, std::string filename, int n, Param par);
bool Read_plt(std::string filename, Param &par, std::list<Circle>& solidList);
void Output_c(Matrix c , std::string filename, int n, Param par);
void Output_Matrix(Matrix A, std::string WorkDir, std::string Variable, int n);
void Output_Matrix_mid(Matrix A, std::string WorkDir, std::string Variable, int n);
void Output_2DArray(double* A, int Nx, int Ny, std::string WorkDir, std::string Variable, int n);

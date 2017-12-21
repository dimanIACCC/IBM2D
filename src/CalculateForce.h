#pragma once

#include "Parameters.h"
#include "SolidBody.h"
#include "Matrix.h"

double DeltaFunction(double x, double y, Param par);
double FunctionD(double r);
void GetInfluenceArea(int &i_min, int &i_max, int &j_min, int &j_max, size_t Ni, size_t Nj, GeomVec x, int size, Param par);

void CalculateForce(Matrix& force_x, Matrix& force_y, std::list<Circle> &iList, Matrix& u, Matrix& v, Param par);
void deformation_velocity(Matrix &u, Matrix &v, Matrix &Exx, Matrix &Eyy, Matrix &Exy, Param par);
void Solids_deformation_velocity_pressure(std::list<Circle> &Solids, Matrix &Exx, Matrix &Eyy, Matrix &Exy, Matrix &p, Param par);
void Solids_Force(std::list<Circle> &Solids, double Re);

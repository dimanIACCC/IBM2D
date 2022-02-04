#pragma once

#include "Parameters.h"
#include "SolidBody.h"
#include "Matrix.h"

void Multidirect_Forcing_Method(Matrix &Fx, Matrix &Fy, Matrix &u, Matrix &v, std::list<Circle> &solidList, Param par);
void CalculateForce(Matrix& force_x, Matrix& force_y, std::list<Circle> &iList, Matrix& u, Matrix& v, Param par);
void deformation_velocity(Matrix &u, Matrix &v, Matrix &Exx, Matrix &Eyy, Matrix &Exy, Param par);
void Solids_deformation_velocity_pressure(std::list<Circle> &Solids, Matrix &Exx, Matrix &Eyy, Matrix &Exy, Matrix &p, Param par);
void Solids_Force(std::list<Circle> &Solids, double Re);
void uf_in_Nodes(std::vector<Node>& Nodes, Matrix u, Matrix v, Param par, int Nn);
void uf_in_Nodes_old(std::vector<Node>& Nodes, Matrix u, Matrix v, Param par, int Nn);

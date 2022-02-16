#pragma once

#include "Parameters.h"
#include "SolidBody.h"
#include "Matrix.h"

void CalculateForce(Matrix& force_x, Matrix& force_y, std::vector<Circle> &iList, Matrix& u, Matrix& v, Param par);
void deformation_velocity(Matrix &u, Matrix &v, Matrix &Exx, Matrix &Eyy, Matrix &Exy, Param par);
void Solids_deformation_velocity_pressure(std::vector<Circle> &Solids, Matrix &Exx, Matrix &Eyy, Matrix &Exy, Matrix &p, Param par);
void uf_in_Nodes    (std::vector<Node>& Nodes, Matrix &u, Matrix &v, Param par, int Nn);
void uf_in_Nodes_old(std::vector<Node>& Nodes, Matrix &u, Matrix &v, Param par, int Nn);
void F_to_Euler_grid    (std::vector<Node>& Nodes, Matrix &Fx_temp, Matrix &Fy_temp, Param par, int Nn);
void F_to_Euler_grid_old(std::vector<Node>& Nodes, Matrix &Fx_temp, Matrix &Fy_temp, Param par, int Nn);

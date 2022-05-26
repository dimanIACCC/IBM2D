#pragma once

#include "Parameters.h"
#include "SolidBody.h"
#include "Matrix.h"

void CalculateForce(Matrix& dFx, Matrix& dFy, float* Ax_beg, float* Ax_end, float* Ay_beg, float* Ay_end, std::vector<Circle> &iList, std::vector<Node> &Nodes, Matrix& u, Matrix& v, Param par);
void deformation_velocity(Matrix &u, Matrix &v, Matrix &Exx, Matrix &Eyy, Matrix &Exy, Param par);
void Solids_deformation_velocity_pressure(std::vector<Circle> &Solids, std::vector<Node> &Nodes, Matrix &Exx, Matrix &Eyy, Matrix &Exy, Matrix &p, Param par);
void uf_in_Nodes    (std::vector<Node>& Nodes, Matrix &u, Matrix &v, Param par, int Nn);
void uf_in_Nodes_old(std::vector<Node>& Nodes, Matrix &u, Matrix &v, Param par, int Nn);
void Make_interaction_Matrix(float* A_beg, float* A_end, int N1, int N2, double d_x, double d_y, std::vector<Node>& Nodes, int Nn_max, Direction Dir);
void F_to_Euler_grid    (std::vector<Node>& Nodes, Matrix &Fx_temp, Matrix &Fy_temp, float* Ax_beg, float* Ax_end, float* Ay_beg, float* Ay_end, Param par, int Nn);
void F_to_Euler_grid_old(std::vector<Node>& Nodes, Matrix &Fx_temp, Matrix &Fy_temp, Param par, int Nn);

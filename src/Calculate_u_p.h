#pragma once

#include "Calculate_press.h"
#include "PredictVel.h"
#include "BiCGStab.h"
#include "CalculateForce.h"

// calculate velocity $U_new, V_new$ and pressure $P$ at the new time step
void Calculate_u_p(Matrix &U_n , Matrix &U_new,
                   Matrix &V_n , Matrix &V_new,
                   Matrix &P_n , Matrix &P_new,
                   Matrix &Fx_n, Matrix &Fx_new,
                   Matrix &Fy_n, Matrix &Fy_new,
                   Template A_u, Template A_v,
                   std::list<Circle> &solidList, Param par);

void ApplyInitialData(Matrix &u, Matrix &v, Matrix &p, Param par);
void Zero_velocity_in_Solids(Matrix &u, Param par, std::list<Circle> iList);

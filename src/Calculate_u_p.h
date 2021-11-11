#pragma once

#include "Calculate_press.h"
#include "PredictVel.h"
#include "BiCGStab.h"
#include "CalculateForce.h"

// calculate velocity $U_new, V_new$ and pressure $P$ at the new time step
void Calculate_u_p(Matrix &U_n , Matrix &U_new,
                   Matrix &V_n , Matrix &V_new,
                   Matrix &P,
                   Matrix &Fx,
                   Matrix &Fy,
                   std::list<Circle> &solidList, Param par);

void Zero_velocity_in_Solids(Matrix &u, Param par, std::list<Circle> iList);

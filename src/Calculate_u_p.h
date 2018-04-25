#pragma once

#include "Calculate_press.h"
#include "PredictVel.h"
#include "BiCGStab.h"
#include "CalculateForce.h"

// calculate velocity $U_new, V_new$ and pressure $P$ at the new time step
void Calculate_u_p(Matrix &U_n, Matrix &V_n,
	Matrix &U_new, Matrix &V_new,
	Matrix &P, Matrix &P_new,
	Matrix &Fx, Matrix &Fy,
    Template A_u, Template A_v,
    std::list<Circle> &solidList, Param par, int N_step);

void ApplyInitialData(Matrix &u, Matrix &p, Param par);
void Zero_velocity_in_Solids(Matrix &u, Param par, std::list<Circle> iList);

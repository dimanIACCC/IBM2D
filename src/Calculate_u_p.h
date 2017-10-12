#pragma once

#include "Calculate_press.h"
#include "PredictVel.h"
#include "BiCGStab.h"
#include "CalculateForce.h"

void Calculate_u_p(Matrix &U_n, Matrix &V_n,
	Matrix &U_new, Matrix &V_new,
	Matrix &P,
	Matrix &Fx, Matrix &Fy,
	ublas::matrix<Template> A_u,
	ublas::matrix<Template> A_v, std::list<Circle> solidList, Param par);
void ApplyInitialData(Matrix &u, Matrix &p, Param par);
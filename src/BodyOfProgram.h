#pragma once
#include "stdafx.h"
#include "CalculateForce.h"
#include "Calculate_u_p.h"
#include "Output.h"
#include "PredictVel.h"


void BodyOfProgram(Param par, std::list<Circle> solidList, Matrix U_n, Matrix V_n, Matrix P, int n0 = 0, bool NeedNewLog = true);
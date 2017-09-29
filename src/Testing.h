#pragma once
#include "stdafx.h"
#include "GeomVec.h"
#include <cmath>
#include "CalculateForce.h"
#include "Output.h"
#include "BiCGStab.h"
#include "Calculate_press.h"
#include "PredictVel.h"
void DoTestForce(int Re);
void DoTesting();
double Sum(Matrix& f);
void ApplyInitialVelocity(Matrix &u, Param par);
void CalcForceDrugLift(Matrix& f, int n, std::ostream &stream);
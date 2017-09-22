#pragma once
#include "stdafx.h"
#include "GeomVec.h"
#include <cmath>
#include "CalculateForce.h"
#include "Output.h"
#include "BiCGStab.h"
#include "Calculate_press.h"
#include "PredictVel.h"
void DoTestForce();
void DoSomeTest();
double diff(Matrix A, Matrix B);
void ApplyInitialVelocity(Matrix &u, Param par);
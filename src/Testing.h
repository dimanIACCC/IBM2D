#pragma once
#include "stdafx.h"
#include "GeomVec.h"
#include <cmath>
#include "BodyOfProgram.h"


void DoTesting();
double Sum(Matrix& f);
void ApplyInitialVelocity(Matrix &u, Param par);
void CalcForceDrugLift(Matrix& f, int n, std::ostream &stream);
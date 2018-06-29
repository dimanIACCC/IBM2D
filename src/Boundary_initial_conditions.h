#pragma once

#include "Parameters.h"
#include "GeomVec.h"
#include "Matrix.h"
#include "Output.h"
#include "Exact_solutions.h"

void ApplyInitialData(Matrix &u, Matrix &v, Matrix &p, Param par);
void Boundary_Conditions(Matrix &u, Param par, Direction Dir, int N_step);

void Taylor_Green_exact(Matrix &u, Matrix &v, Matrix &p, Param par, double time);
void Lamb_Oseen_exact_uv(Matrix &uv, Direction Dir, Param par, double time, bool boundary);
void Lamb_Oseen_exact_p (Matrix &p, Param par, double time);

#pragma once

#include "Parameters.h"
#include "GeomVec.h"
#include "Matrix.h"
#include "Output.h"
#include "Exact_solutions.h"

void Boundary_Conditions(Matrix &u, Param par, Direction Dir, double time);

void fill_exact  (Matrix &u, Matrix &v, Matrix &p, Param &par, double time_uv, double time_p);
void fill_exact_u(Matrix &u, Param par, double time);
void fill_exact_v(Matrix &v, Param par, double time);
void fill_exact_p(Matrix &p, Param par, double time);
void BC_exact_p  (Matrix &p, Param par, double time);

void Lamb_Oseen_p_test(Param par, double time);

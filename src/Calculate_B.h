#pragma once
#include "stdafx.h"
#include "Grid.h"
using namespace std;


Matrix CalculateB_u(Matrix &u_n, Matrix &v_n, Matrix &u_prev, Matrix &v_prev, Matrix &p, Matrix &force, Grid grid, double Re);
Matrix CalculateB_v(Matrix &u_n, Matrix &v_n, Matrix &u_prev, Matrix &v_prev, Matrix &p, Matrix &force, Grid grid, double Re);

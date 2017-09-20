#pragma once

#include "stdafx.h"
#include <vector>

//Type of matrix with double
typedef  std::vector<std::vector<double>> Matrix;

//Macros which make matrix type of double size of n*m
#define CreateMatrix(name, n, m) Matrix name(n,std::vector<double>(m, 0))

enum Direction {
	Du, Dv
};

class Template {
public:
	double C, L, R, U, D;
	Template();
};

double L(Matrix &A, int i, int j, Direction dir);
double R(Matrix &A, int i, int j, Direction dir);
double D(Matrix &A, int i, int j, Direction dir);
double U(Matrix &A, int i, int j, Direction dir);

double RD(Matrix &A, int i, int j, Direction dir);

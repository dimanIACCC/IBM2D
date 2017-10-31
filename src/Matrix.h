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

double L(Matrix &A, size_t i, size_t j, Direction dir, size_t Ni, size_t Nj);
double R(Matrix &A, size_t i, size_t j, Direction dir, size_t Ni, size_t Nj);
double D(Matrix &A, size_t i, size_t j, Direction dir, size_t Ni, size_t Nj);
double U(Matrix &A, size_t i, size_t j, Direction dir, size_t Ni, size_t Nj);

double RD(Matrix &A, size_t i, size_t j, Direction dir, size_t Ni, size_t Nj);

double max(Matrix &A);
double diff(Matrix &A, Matrix &B);
double Summ(Matrix& A);

Matrix &operator+=(Matrix &A, const Matrix &B);
Matrix operator*(const Matrix &A, const double &b);

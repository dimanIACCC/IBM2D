#pragma once

#include "stdafx.h"
#include "Parameters.h"
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

double L(Matrix &A, size_t i, size_t j, Direction dir);
double R(Matrix &A, size_t i, size_t j, Direction dir);
double D(Matrix &A, size_t i, size_t j, Direction dir);
double U(Matrix &A, size_t i, size_t j, Direction dir);

double UL(Matrix &A, size_t i, size_t j, Direction dir);

double Matrix_max(const Matrix &A);
double diff(Matrix &A, Matrix &B);
double Summ(Matrix& A);

Matrix &operator+=(Matrix &A, const Matrix &B);
Matrix operator+(const Matrix& A, const Matrix& B);
Matrix operator-(const Matrix& A, const Matrix& B);
Matrix operator*(const Matrix &A, const double &b);

void Matrix_to_DoubleArray(Matrix &M, double* D, boundary_conditions BC);
void DoubleArray_to_Matrix(double* D, Matrix &M, boundary_conditions BC);
void MatrixU_to_DoubleArray(Matrix &M, double* D, boundary_conditions BC);
void DoubleArray_to_MatrixU(double* D, Matrix &M, boundary_conditions BC);

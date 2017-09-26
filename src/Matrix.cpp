﻿#include "Matrix.h"


Template::Template() {
	C = L = R = D = U = 0.0;
}

double L(Matrix &A, size_t i, size_t j, Direction dir) {
	double result;
	if      (dir == Du)
		result = A[i - 1][j];
	else if (dir == Dv)
		result = A[i][j - 1];
	else
		std::cout << "L: wrong direction" << std::endl;
	return result;
}
double R(Matrix &A, size_t i, size_t j, Direction dir) {
	double result;
	if (dir == Du)
		result = A[i + 1][j];
	else if (dir == Dv)
		result = A[i][j + 1];
	else
		std::cout << "R: wrong direction" << std::endl;
	return result;
}
double D(Matrix &A, size_t i, size_t j, Direction dir) {
	double result;
	if (dir == Du)
		result = A[i][j - 1];
	else if (dir == Dv)
		result = A[i - 1][j];
	else
		std::cout << "D: wrong direction" << std::endl;
	return result;
}
double U(Matrix &A, size_t i, size_t j, Direction dir) {
	double result;
	if (dir == Du)
		result = A[i][j + 1]; 
	else if (dir == Dv)
		result = A[i + 1][j];
	else
		std::cout << "U: wrong direction" << std::endl;
	return result;
}

double RD(Matrix &A, size_t i, size_t j, Direction dir) {
	double result;
	if (dir == Du)
		result = A[i + 1][j - 1];
	else if (dir == Dv)
		result = A[i - 1][j + 1];
	else
		std::cout << "RD: wrong direction" << std::endl;
	return result;
}

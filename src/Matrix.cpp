#include "Matrix.h"


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

double max(Matrix &A) {
	double eps = 0.0;
	for (size_t i = 0; i < A.size(); ++i) {
		for (size_t j = 0; j < A[0].size(); ++j) {
			const double eps_tmp = fabs(A[i][j]);
			if (eps_tmp > eps) {
				eps = eps_tmp;
			}
		}
	}
	return eps;
}

double diff(Matrix &A, Matrix &B) {
	double eps = 0.0;
	for (size_t i = 0; i < A.size(); ++i) {
		for (size_t j = 0; j < A[0].size(); ++j) {
			const double eps_tmp = fabs(A[i][j] - B[i][j]);
			if (eps_tmp > eps) {
				eps = eps_tmp;
			}
		}
	}
	return eps;
}

double Summ(Matrix& A) {
	double sum = 0;
	for (size_t i = 0; i < A.size(); i++) {
		for (size_t j = 0; j < A[0].size(); j++) {
			sum += A[i][j];
		}
	}
	return sum;
}

Matrix &operator+=(Matrix& A, const Matrix& B) {
	for (size_t i = 0; i < A.size(); i++) {
		for (size_t j = 0; j < A[0].size(); j++) {
			A[i][j] += B[i][j];
		}
	}
	return A;
}

Matrix operator*(const Matrix &A, const double &b) {
	CreateMatrix(C, A.size(), A[0].size());
	for (size_t i = 0; i < C.size(); i++) {
		for (size_t j = 0; j < C[0].size(); j++) {
			C[i][j] = A[i][j] * b;
		}
	}
	return C;
}

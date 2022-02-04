#include "Matrix.h"


Template::Template() {
	C = L = R = D = U = 0.0;
}

double L(Matrix &A, size_t i, size_t j, Direction dir) {
	double result;
	if (dir == Du) {
		result = A[i - 1][j];
	}
	else if (dir == Dv) {
		result = A[i][j - 1];
	}
	else
		std::cout << "L: wrong direction" << std::endl;
	return result;
}
double R(Matrix &A, size_t i, size_t j, Direction dir) {
	double result;
	if (dir == Du) {
		result = A[i + 1][j];
	}
	else if (dir == Dv) {
		result = A[i][j + 1];
	}
	else
		std::cout << "R: wrong direction" << std::endl;
	return result;
}
double D(Matrix &A, size_t i, size_t j, Direction dir) {
	double result;
	if (dir == Du) {
		result = A[i][j - 1];
	}
	else if (dir == Dv) {
		if (i > 0) result = A[i - 1][j];
	}
	else
		std::cout << "D: wrong direction" << std::endl;
	return result;
}
double U(Matrix &A, size_t i, size_t j, Direction dir) {
	double result;
	if (dir == Du) {
		result = A[i][j + 1];
	}
	else if (dir == Dv) {
		result = A[i + 1][j];
	}
	else
		std::cout << "U: wrong direction" << std::endl;
	return result;
}

double UL(Matrix &A, size_t i, size_t j, Direction dir) {
	double result;
	if (dir == Du) {
		result = A[i - 1][j + 1];
	}
	else if (dir == Dv) {
		result = A[i + 1][j - 1];
	}
	else
		std::cout << "RD: wrong direction" << std::endl;
	return result;
}

double Matrix_max(const Matrix &A) {
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

double Summ(Matrix &A) {
	double sum = 0;
	for (size_t i = 0; i < A.size(); i++) {
		for (size_t j = 0; j < A[0].size(); j++) {
			sum += A[i][j];
		}
	}
	return sum;
}

Matrix &operator+=(Matrix &A, const Matrix &B) {
	for (size_t i = 0; i < A.size(); i++) {
		for (size_t j = 0; j < A[0].size(); j++) {
			A[i][j] += B[i][j];
		}
	}
	return A;
}

Matrix operator+(const Matrix &A, const Matrix &B) {
	CreateMatrix(C, A.size(), A[1].size());
	for (size_t i = 0; i < C.size(); i++) {
		for (size_t j = 0; j < C[0].size(); j++) {
			C[i][j] = A[i][j] + B[i][j];
		}
	}
	return C;
}

Matrix operator-(const Matrix &A, const Matrix &B) {
	CreateMatrix(C, A.size(), A[0].size());
	for (size_t i = 0; i < C.size(); i++) {
		for (size_t j = 0; j < C[0].size(); j++) {
			C[i][j] = A[i][j] - B[i][j];
		}
	}
	return C;
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

void Matrix_to_DoubleArray(Matrix &M, double* D, boundary_conditions BC) {

	int Nx = M.size();
	int Ny = M[0].size();

	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			if (BC == periodical) {
				if (i == Nx - 1); else  D[i + j*(Nx - 1)] = M[i][j];
			}
			else D[i + j*Nx] = M[i][j];
		}
	}

	//Output_2DArray(D, Nx - 1, Ny, "Result/", "Array", 555);
}

void DoubleArray_to_Matrix(double* D, Matrix &M, boundary_conditions BC) {

	int Nx = M.size();
	int Ny = M[0].size();

	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			if (BC == periodical) {
				if (i == Nx - 1) M[i][j] = D[1 + j*(Nx - 1)];
				else           M[i][j] = D[i + j*(Nx - 1)];
			}
			else M[i][j] = D[i + j*Nx];
		}
	}
}


void MatrixU_to_DoubleArray(Matrix &M, double* D, boundary_conditions BC) {

	int Nx = M.size();
	int Ny = M[0].size();

	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			if (BC == periodical) {
				if (i == 0);
				else if (i == Nx - 1);
				else  D[(i - 1) + j*(Nx - 2)] = M[i][j];
			}
			else D[i + j*Nx] = M[i][j];
		}
	}

}

void DoubleArray_to_MatrixU(double* D, Matrix &M, boundary_conditions BC) {

	int Nx = M.size();
	int Ny = M[0].size();

	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			if (BC == periodical) {
				if (i == 0)      M[i][j] = D[(Nx - 2) + j*(Nx - 2)];
				else if (i == Nx - 1) M[i][j] = D[0 + j*(Nx - 2)];
				else                  M[i][j] = D[(i - 1) + j*(Nx - 2)];
			}
			else M[i][j] = D[i + j*Nx];
		}
	}
}


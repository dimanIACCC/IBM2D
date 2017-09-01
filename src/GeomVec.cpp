#include "GeomVec.h"

double length(GeomVec x) {
	double result = 0.0;
	for (int i = 1; i <= 3; i++) {
		result += pow(x[i], 2);
	}
	result = sqrt(result);
	return result;
}

double dot_product(GeomVec v1, GeomVec v2) {
	double result;
	result = v1[1]*v2[1] + v1[2]*v2[2] + v1[3]*v2[3];
	return result;
}


GeomVec x_product(GeomVec v1, GeomVec v2){
	GeomVec result;
	result[0] = 0;
	result[1] = v1[2] * v2[3] - v1[3] * v2[2];
	result[2] = v1[3] * v2[1] - v1[1] * v2[3];
	result[3] = v1[1] * v2[2] - v1[2] * v2[1];

	return result;
}

GeomMat x_product_Vec_Mat(GeomVec v, GeomMat A) {
	GeomMat result;
	GeomVec column;

	for (int i = 0; i <= 3; i++) {
		column = x_product(v, get_column(A, i)); // vector product of Vector $v$ and $i$-th column of Tensor $A$
		put_column(result, i, column);           // puts column to result Tensor
	}
	return result;
}

GeomVec get_column(GeomMat A, int i) {
	GeomVec v;
	for (int j = 0; j <= 3; j++) {
		v[j] = A(i, j);
	}
	return v;
}

void put_column(GeomMat& A, int i, GeomVec v) {
	for (int j = 0; j <= 3; j++) {
		A(i, j) = v[j];
	}
}

GeomMat diada(GeomVec v1, GeomVec v2) {
	GeomMat A;
	for (int i = 0; i <= 3; i++) {
		for (int j = 0; j <= 3; j++){
			A(i,j) = v1[i] * v2[j];
		}
	}
	return A;
}


GeomMat E = InitE(); // initialize unit matrix

GeomMat M_rotate(GeomVec& o){
	double angle;
	angle = length(o);
	if (length(o)>1.e-14) o = o / angle;
	GeomMat O;
	O = (1 - cos(angle)) * diada(o, o)
		   + cos(angle)  * E
		   + sin(angle)  * x_product_Vec_Mat(o, E);
	return O;
}

GeomVec rotate_Vector_around_vector(GeomVec v, GeomVec o) {
	GeomVec v_out;
	GeomMat Q, Qinv;

	Q = M_rotate(o);
	Qinv = InvertMat(Q);
	v_out = ublas::prod (Qinv, v);

	return v_out;
}

ublas::matrix<double> InvertMat(ublas::matrix<double> A) {
	typedef ublas::permutation_matrix<std::size_t> pmatrix;
	ublas::matrix<double> inverse(ublas::identity_matrix<double>(A.size1()));

	pmatrix pm(A.size1());                  	// create a permutation matrix for the LU-factorization
	int res = ublas::lu_factorize(A, pm);       // perform LU-factorization
	if (res != 0) std::cout << "InvertGeomMat: cannot inverse Matrix" << std::endl;
	lu_substitute(A, pm, inverse);             // backsubstitute to get the inverse matrix

	return inverse;
}

GeomMat InitE() {
	GeomMat E;
	for (int i = 0; i <= 3; i++) {
		for (int j = 0; j <= 3; j++) {
			if (i == j) {
				E(i, j) = 1;
			}
			else {
				E(i, j) = 0;
			}
		}
	}
	return E;
}

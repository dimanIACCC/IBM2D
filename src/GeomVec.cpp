#include "GeomVec.h"

double length(GeomVec x) {
	double result = 0.0;
	for (int i = 1; i <= DIM; i++) {
		result += pow(x[i], 2);
	}
	result = sqrt(result);
	return result;
}
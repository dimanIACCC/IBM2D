#include "Testing.h"

void DoSomeTest() {
	DoTestForce();
	//
	//TO DO: realize some tests here!
}
void DoTestForce() {
	std::list<Circle> solidList;
	Param par;
	par.Nn = 50;
	par.r = 0.5;
	par.alpha_f = -8000;
	par.beta_f = -2000;
	Circle c(1.5, 5.1, par);
	solidList.push_back(c);
	CreateMatrix(Force_x, par.N1, par.N2 + 1);
	CreateMatrix(Force_y, par.N1 + 1, par.N2);
	CreateMatrix(U, par.N1, par.N2 + 1);
	CreateMatrix(V, par.N1 + 1, par.N2);
	CreateMatrix(Force_x_prev, par.N1, par.N2 + 1);
	CreateMatrix(Force_y_prev, par.N1 + 1, par.N2);
	for (int i = 0; i < 4; i++)
	{
		double coef=0;
		par.N1 = 10*pow(2,i);
		par.N2 = 10* pow(2, i);

		CalculateForce(Force_x, Force_y, solidList, U, V, par);
		if (i > 0) {
			double dif = diff(Force_x, Force_x_prev);
			std::cout << dif / (pow(2, 2) - 1);
		}
		
		Force_x_prev = Force_x;
		Force_y_prev = Force_y;

	}
}
double diff(Matrix A, Matrix B) {
	double dif = 0;
	for (int i = 0; i < (int)A.size(); i++) {
		for (int j = 0; j < (int)A[0].size(); j++)
		{
			if(abs(A[i][j] - B[i][j]) > dif) dif = abs(A[i][j] - B[i][j]);
		}
	}
	return dif;
}

double Summ(Matrix& force) {
	double sum = 0;
	for (int i = 0; i < (int)force.size(); i++)
		for (int j = 0; j < (int)force[0].size(); j++)
		{
			sum += force[i][j];
		}


	return abs(sum);
}

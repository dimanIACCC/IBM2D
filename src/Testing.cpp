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
	Circle c(1.5, 5.1, par);
	solidList.push_back(c);
	CreateMatrix(Force_x, par.N1, par.N2 + 1);
	CreateMatrix(Force_y, par.N1 + 1, par.N2);
	CreateMatrix(U, par.N1, par.N2 + 1);
	CreateMatrix(V, par.N1 + 1, par.N2);
	CreateMatrix(Force_x_prev, par.N1, par.N2 + 1);
	CreateMatrix(Force_y_prev, par.N1 + 1, par.N2);
	CreateMatrix(P, par.N1 + 1, par.N2 + 1);

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


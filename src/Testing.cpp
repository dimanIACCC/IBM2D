#include "stdafx.h"
#include "Testing.h"
#include "CalculateForce.h"
using namespace std;
void DoSomeTest() {
	DoTestForce();
	//
	//TO DO: realize some tests here!
}
void DoTestForce() {
	list<Circle> solidList;
	Grid grid;
	grid.NF = 50;
	Circle c(1.5, 5.1, 0.5, 0, grid);
	solidList.push_back(c);
	CreateMatrix(Force_x, grid.N1, grid.N2 + 1);
	CreateMatrix(Force_y, grid.N1 + 1, grid.N2);
	CreateMatrix(U, grid.N1, grid.N2 + 1);
	CreateMatrix(V, grid.N1 + 1, grid.N2);
	CreateMatrix(Force_x_prev, grid.N1, grid.N2 + 1);
	CreateMatrix(Force_y_prev, grid.N1 + 1, grid.N2);
	for (int i = 0; i < 4; i++)
	{
		double coef=0;
		grid.N1 = 10*pow(2,i);
		grid.N2 = 10* pow(2, i);


		CalculateForce_X(Force_x, solidList, U, coef, grid, -8000, -2000, 10);
		CalculateForce_Y(Force_y, solidList, V, coef, grid, -8000, -2000, 10);
		if (i > 0) {
			double dif = diff(Force_x, Force_x_prev);
			cout << dif / (pow(2, 2) - 1);
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
#pragma once


#include "stdafx.h"
#include "Grid.h"

using namespace std;
class SolidBody
{
public:
	int start_n; // number of iteration, when solid added
	bool moveSolid;
	bool eraseSolid;
	double U, V, x,y;


	vector<double> Bound[2];
	vector<double> Integral_x; // value of integral when calculating force
	vector<double> Integral_y;

	
	SolidBody(double x, double y, int n);
	~SolidBody();

 	
};
class Circle : public SolidBody{
public:
	double r;
	double d_s;
	Circle(double x, double y,  double r, int n, Grid grid);
	//void AddSolid(list<Circle> &iList);
	~Circle();
};

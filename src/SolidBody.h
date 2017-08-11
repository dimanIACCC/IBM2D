#pragma once


#include "stdafx.h"
#include "Grid.h"

using namespace std;

class Node
{
public:
	double x[3];        // coordinates
	double u[3];        // velocity
	double integral[3]; // value of integral when calculating force
};

class SolidBody
{
public:
	int start_n; // number of iteration, when solid added
	bool moveSolid;
	bool eraseSolid;
	double xc[3];     // coordinates of the mass center
	double Uc[3];     // velocity of the mass center
	double omega;     // angular velocity


	vector<double> Bound[2];
	vector<double> Integral_x; // value of integral when calculating force
	vector<double> Integral_y;

	vector<Node> Nodes;
	
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


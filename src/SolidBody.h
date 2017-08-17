#pragma once


#include "stdafx.h"
#include "Grid.h"

using namespace std;


class Node
{
public:
	GeomVec x;        // coordinates
	GeomVec U;        // velocity of the fluid in the Node

	GeomVec f;        // force
	GeomVec Integral; // value of integral when calculating force
};

class SolidBody
{

public:
	int start_n; // number of iteration, when solid added
	bool moveSolid;
	bool eraseSolid;
	GeomVec xc;     // coordinates of the mass center
	GeomVec uc;     // velocity of the mass center
	GeomVec omega;  // angular velocity
	vector<Node> Nodes;   // Nodes of the SolidBody mesh
	
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


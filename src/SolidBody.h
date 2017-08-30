#pragma once

#include "GeomVec.h"
#include "Parameters.h"
#include <boost/algorithm/string.hpp>

class Node
{
public:
	GeomVec x;        // coordinates
	GeomVec uf;       // velocity of the fluid in the Node
	GeomVec us;       // velocity of the SolidBody in the Node
	GeomVec f;        // force
	GeomVec Integral; // value of integral when calculating force
};

class SolidBody
{
public:
	bool moveSolid;
	GeomVec xc;     // coordinates of the mass center
	GeomVec uc;     // velocity of the mass center
	GeomVec omega;  // angular velocity
	GeomVec f;      // force applied to the whole SolidBody
	GeomVec tau;    // torque, moment of force applied to the whole SolidBody
	double I;       // angular momentum
	double rho;     // density
	double V;       // volume
	std::vector<Node> Nodes;   // Nodes of the SolidBody mesh
	int Nn;                    // Number of Nodes
	
	SolidBody(double x, double y, double ux, double uy, double omega, double rho, int Nn, bool moveSolid);
	~SolidBody();
	void velocities();      // calculates the velocities in all Nodes of the SolidBody
};

class Circle : public SolidBody{
public:
	double r;
	double d_s;
	Circle(double x, double y, double ux, double uy, double omega, double rho, int Nn, bool moveSolid, double r);
	Circle::Circle(double x, double y, Param par);
	~Circle();
};

void Read_Solids(std::string filename, std::list<Circle>& Solids, Param par);


#include "stdafx.h"
#include "SolidBody.h"


SolidBody::SolidBody(double x, double y)
{
	 moveSolid = false;
	 std::fill(uc.begin()   , uc.end()   , 0.0); // fill vector uc    with zeros
	 std::fill(omega.begin(), omega.end(), 0.0); // fill vector omega with zeros
	 xc[1] = x;
	 xc[2] = y;
}


SolidBody::~SolidBody()
{
}

Circle::Circle(double x, double y, double r, int NF) : SolidBody(x, y){
	this->r = r;
	this->Nn = NF;
	this->d_s = (2.0*M_PI*r) / Nn;
	Nodes.resize(Nn);

	for (int i = 0; i < Nn; ++i){
		Nodes[i].x[1] = x + cos(i * 2.0 * M_PI / Nn) * r;
		Nodes[i].x[2] = y + sin(i * 2.0 * M_PI / Nn) * r;
	}
	V = M_PI * r * r;
	rho = 10.0 / V; // corresponds to old formula for force
	I =  V * r * r / 2.0; // angular momentum for unit density
}
 
Circle::~Circle()
{
}

void SolidBody::velocities() {
	for (int i = 0; i < Nn; ++i) {
		GeomVec r;
		r = Nodes[i].x - xc;
		Nodes[i].us = uc + x_product(omega, r);
	}
}
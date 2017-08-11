
#include "stdafx.h"
#include "SolidBody.h"


SolidBody::SolidBody(double x, double y, int n)
{
	 this->start_n = n; // number of iteration, when solid added
	 moveSolid = false;
	 moveSolid = false;
	 Uc[1] = 0.0;
	 Uc[2] = 0.0;
	 this->xc[1] = x;
	 this->xc[2] = y;
}


SolidBody::~SolidBody()
{
}

Circle::Circle(double x, double y, double r, int n, Grid grid) : SolidBody(x, y, n){
	this->r = r;
	this->d_s = (2.0*M_PI*r) / grid.NF;
	Bound[0].resize(grid.NF);
	Bound[1].resize(grid.NF);
	Integral_x.resize(grid.NF);
	Integral_y.resize(grid.NF);
	for (int i = 0; i < grid.NF; ++i){
		Bound[0][i] = x + cos(i * 2.0 * M_PI / grid.NF) * r;
		Bound[1][i] = y + sin(i * 2.0 * M_PI / grid.NF) * r;
	}
}
 
Circle::~Circle()
{
}

// void Circle::AddSolid(list<Circle> &iList){
//	 if (!iList.empty()) {
//		 int a = prev(iList.end())->first;
//		 iList.insert(make_pair(a + 1, this));
//	 }
//	 else
//		iList.insert(make_pair(0, this));
//}

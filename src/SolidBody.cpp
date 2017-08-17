
#include "stdafx.h"
#include "SolidBody.h"


SolidBody::SolidBody(double x, double y, int n)
{
	 start_n = n; // number of iteration, when solid added
	 moveSolid = false;
	 moveSolid = false;
	 std::fill(uc.begin()   , uc.end()   , 0.0); // fill vector uc    with zeros
	 std::fill(omega.begin(), omega.end(), 0.0); // fill vector omega with zeros
	 xc[1] = x;
	 xc[2] = y;
}


SolidBody::~SolidBody()
{
}

Circle::Circle(double x, double y, double r, int n, Grid grid) : SolidBody(x, y, n){
	this->r = r;
	this->d_s = (2.0*M_PI*r) / grid.NF;
	Nodes.resize(grid.NF);

	for (int i = 0; i < grid.NF; ++i){
		Nodes[i].x[1] = x + cos(i * 2.0 * M_PI / grid.NF) * r;
		Nodes[i].x[2] = y + sin(i * 2.0 * M_PI / grid.NF) * r;
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

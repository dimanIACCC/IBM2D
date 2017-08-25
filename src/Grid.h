#pragma once

struct Grid{
	double L = 0; //length
	double H = 0; //height
	int N1 = 0; // number of points in x-direction
	int N2 = 0; // number of points in y-direction
	double d_x = 0.0;
	double d_y = 0.0;
	double d_t = 0.0;
	int NF = 0; // number of points of Immersed Boundary
};
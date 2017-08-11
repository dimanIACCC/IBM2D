#pragma once
//Type of matrix with double
typedef  std::vector<std::vector<double>> Matrix;

//Macros which make matrix type of double size of n*m
#define CreateMatrix(name, n, m) Matrix name(n,std::vector<double>(m, 0))


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
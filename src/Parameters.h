#pragma once

#include "GeomVec.h"
#include "String.h"

class Param{
public:
	double Re;                // Reynolds number
	double L;                 // length
	double H;                 // height
	int N1;                   // number of points in x-direction
	int N2;                   // number of points in y-direction
	double d_x;               // mesh step in x-direction
	double d_y;               // mesh step in y-direction
	double d_t;               // time step
	double alpha_f;           // Force parameter
	double beta_f;            // Force parameter
	int Nn;                   // number of nodes of Solid
	double rho;               // default density of Solid
	double r;                 // default radius of Circle
	int output_step = 0;      // frequency of output
	int N_max;                // number of total iterations
	int N_Zeidel;             // number of iterations in Zeidel method
	double Zeidel_eps;        // tolerance for Zeidel method
	bool InelasticCollision;  // Inelastic - true, elastic - false
	double k_dist;            // coefficient for minimal distance between Solids
	Param();
	Param(std::string filename);
};

double ux_Poiseuille(double y, double H);
GeomVec x_p(int i, int j, Param par); // coordinates of (i,j)-th node for pressure p mesh
GeomVec x_u(int i, int j, Param par); // coordinates of (i,j)-th node for velocity u mesh
GeomVec x_v(int i, int j, Param par); // coordinates of (i,j)-th node for velocity v mesh

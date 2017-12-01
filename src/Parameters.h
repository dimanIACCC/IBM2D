#pragma once

#include "GeomVec.h"
#include "String.h"
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

enum Boundary_Conditions {
	u_infinity, u_inflow, periodical
};

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
	int AddSolids_N;          // number of added Solids
	int AddSolids_start;      // step when Solids start to add
	int AddSolids_interval;   // interval for Solids adding
	Boundary_Conditions BC;   // Boundary Conditions
	int SolidName_max;        // Maximal Name of Solids
	std::string WorkDir;      // WorkDir
	Param();

	Param(std::string WorkDir, std::string filename);

};

Boundary_Conditions string_to_BC(std::string s);
double ux_Poiseuille(double y, double H);
double dpdx_Poiseuille(double H, double Re);
GeomVec x_p(int i, int j, Param par); // coordinates of (i,j)-th node for pressure p mesh
GeomVec x_u(int i, int j, Param par); // coordinates of (i,j)-th node for velocity u mesh
GeomVec x_v(int i, int j, Param par); // coordinates of (i,j)-th node for velocity v mesh

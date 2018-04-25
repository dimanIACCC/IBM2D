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
	int N1_u, N2_u;           // sizes for u-direction arrays
	double d_x;               // mesh step in x-direction
	double d_y;               // mesh step in y-direction
	double d_t;               // time step
	int Nn;                   // number of nodes of Solid
	double rho;               // default density of Solid
	double r;                 // default radius of Circle
	int output_step = 0;      // frequency of output
	int N_max;                // number of total iterations
	int N_Zeidel;             // number of iterations in Zeidel method
	double Zeidel_eps;        // tolerance for Zeidel method
	double eps_P;             // tolerance for Pressure correction
	bool InelasticCollision;  // Inelastic - true, elastic - false
	double k_dist;            // coefficient for minimal distance between Solids
	int AddSolids_N;          // number of added Solids
	int AddSolids_start;      // step when Solids start to add
	int AddSolids_interval;   // interval for Solids adding
	Boundary_Conditions BC;   // Boundary Conditions
	double u_wall;            // velocity of the channel walls
	int SolidName_max;        // Maximal Name of Solids
	std::string WorkDir;      // WorkDir
	Param();
	Param(std::string WorkDir);
	Param(std::string WorkDir, std::string filename);

};

Boundary_Conditions string_to_BC(std::string s);
double ux_Poiseuille(double y, double H);
double dux_dy_Poiseuille(double y, double H);
double dpdx_Poiseuille(double H, double Re);
GeomVec x_p(int i, int j, Param par); // coordinates of (i,j)-th node for pressure p mesh
GeomVec x_u(int i, int j, Param par); // coordinates of (i,j)-th node for velocity u mesh
GeomVec x_v(int i, int j, Param par); // coordinates of (i,j)-th node for velocity v mesh
GeomVec x_c(int i, int j, Param par); // coordinates of (i,j)-th node for corners    mesh

double DeltaFunction(double x, double y, Param par);
double FunctionD(double r);
void GetInfluenceArea(int &i_min, int &i_max, int &j_min, int &j_max, size_t Ni, size_t Nj, GeomVec x, int size, Param par);
double Volume_Frac(GeomVec xc, double r, GeomVec x, double dx, double dy);
double Heaviside(double x);

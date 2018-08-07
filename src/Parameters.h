#pragma once

#include "GeomVec.h"
#include "String.h"
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

enum boundary_conditions {
	u_infinity, u_inflow, periodical, Taylor_Green, Lamb_Oseen, Line_Vortex
};

class Param{
public:
	double Re;                // Reynolds number
	double L;                 // length
	double H;                 // height
	int N1;                   // number of points in x-direction
	int N2;                   // number of points in y-direction
	int N1_u, N2_u;           // sizes for u-direction arrays
	int N1_v, N2_v;           // sizes for v-direction arrays
	int N1_p, N2_p;           // sizes for pressure (cell centers) arrays
	double d_x;               // mesh step in x-direction
	double d_y;               // mesh step in y-direction
	double ldxdx, ldydy;      // coefficients in Laplace approximation
	double d_t;               // time step
	int N_step;               // step number
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
	boundary_conditions BC;   // Boundary Conditions
	GeomVec x0;               // Special point in test problems Lamb_Oseen and Line_Vortex
	GeomVec k;                // Spatial frequencies for Taylor_Green vortices
	double Lamb_Oseen_r0;     // initial radius of Lamb_Oseen vortex
	double u_wall;            // velocity of the channel walls
	int SolidName_max;        // Maximal Name of Solids
	std::string WorkDir;      // WorkDir
	Param();
	Param(std::string WorkDir);
	Param(std::string WorkDir, std::string filename);
	void init();
};

boundary_conditions string_to_BC(std::string s);
GeomVec x_p(int i, int j, Param par); // coordinates of (i,j)-th node for pressure p mesh
GeomVec x_u(int i, int j, Param par); // coordinates of (i,j)-th node for velocity u mesh
GeomVec x_v(int i, int j, Param par); // coordinates of (i,j)-th node for velocity v mesh
GeomVec x_c(int i, int j, Param par); // coordinates of (i,j)-th node for corners    mesh

double DeltaFunction(double x, double y, Param &par);
double FunctionD(double r);
void GetInfluenceArea(int &i_min, int &i_max, int &j_min, int &j_max, size_t Ni, size_t Nj, GeomVec x, int size, Param par);
double Volume_Frac(GeomVec xc, double r, GeomVec x, double dx, double dy);
double Heaviside(double x);
int i_real_u(int i, Param par);
int i_real_v(int i, Param par);

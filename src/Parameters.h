#pragma once

#include "GeomVec.h"
#include "String.h"
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

enum boundary_conditions {
	u_infinity, u_inflow, periodical, Taylor_Green, Lamb_Oseen, Line_Vortex, box
};

class Param{
public:
	std::string WorkDir;      // WorkDir
	int N_step;               // step number
	boundary_conditions BC;   // Boundary Conditions

	// physical parameters
	double d_t;               // time step
	double L;                 // length
	double H;                 // height
	double Re;                // Reynolds number
	double grad_p_x;          // pressure gradient in x direction
	double Gravity_angle;     // Angle for rotation of Gravity vector
	double Gravity_module;    // Gravity vector module

	// numerical parameters
	int N1;                   // number of points in x-direction
	int N2;                   // number of points in y-direction
	int Nn;                   // number of nodes of Solid
	int output_step;          // frequency of output
	int IBM;                  // type of Immersed Boundary Method
	int DeltaP_method;        // method for pressure correction equation
	int N_Zeidel;             // number of iterations in Zeidel method
	double Zeidel_eps;        // tolerance for Zeidel method
	int s_max;                // max number of iterations for pressure and IBM force
	double eps_P;             // tolerance for Pressure correction
	int N_Force;              // number of iterations for force

	// parameters for many particles
	double rho;               // default density of Solid
	double r;                 // default radius of Solid
	int AddSolids_N;          // number of added Solids
	int AddSolids_start;      // step when Solids start to add
	int AddSolids_interval;   // interval for Solids adding
	int SolidName_max;        // Maximal Name of Solids
	double k_dist;            // coefficient for minimal distance between Solids

	// parameters for special problems
	double Lamb_Oseen_r0;     // initial radius of Lamb_Oseen vortex
	double u_wall;            // velocity of the channel walls
	double omega_BC;          // omega of the circle driven by external force

	// calculated parameters (dependent on basic ones)
	int N1_u, N2_u;           // sizes for u-direction arrays
	int N1_v, N2_v;           // sizes for v-direction arrays
	int N1_p, N2_p;           // sizes for pressure (cell centers) arrays
	double d_x;               // mesh step in x-direction
	double d_y;               // mesh step in y-direction
	double ldxdx, ldydy;      // coefficients in Laplace approximation
	GeomVec x0;               // Special point in test problems Lamb_Oseen and Line_Vortex
	GeomVec k;                // Spatial frequencies for Taylor_Green vortices
	GeomVec Gravity;          // Gravity vector

	//constructors
	Param();
	Param(std::string WorkDir);
	Param(std::string WorkDir, std::string filename);

	//methods
	void init();
	void read_line(std::string);
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

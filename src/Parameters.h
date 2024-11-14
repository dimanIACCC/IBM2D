#pragma once

#include "GeomVec.h"
#include "String.h"
#include <boost/filesystem.hpp>
#include <amp.h>
#include <amp_math.h>
using namespace concurrency;
using namespace precise_math;

namespace fs = boost::filesystem;

enum boundary_conditions {
	u_infinity, u_inflow, periodical, Taylor_Green, Lamb_Oseen, Line_Vortex, box
};

class Param{
public:
	std::string WorkDir = "";              // WorkDir
	int N_step = 0 ;                       // step number
	double time = 0.;                      // time
	boundary_conditions BC = Taylor_Green; // Boundary Conditions

	// physical parameters
	double d_t = 0.0001;               // time step
	double L = 1.    ;                 // length
	double H = 1.;                     // height
	double Re = 1.;                    // Reynolds number
	double grad_p_x = 0.;              // pressure gradient in x direction
	double Gravity_angle = 00.0;       // Angle for rotation of Gravity vector
	double Gravity_module = 0.0;       // Gravity vector module
	
	// numerical parameters
	int N1 = 100;                    // number of points in x-direction
	int N2 = 100;                    // number of points in y-direction
	int output_step = 10;            // frequency of output
	int IBM = 1;                     // type of Immersed Boundary Method
	int DeltaP_method = 1;           // method for pressure correction equation
	int N_Zeidel = 500000;           // number of iterations in Zeidel method
	double Zeidel_eps = 1e-5;        // tolerance for Zeidel method
	int s_max = 20;                  // max number of iterations for pressure and IBM force
	double eps_P = 1.e-5;            // tolerance for Pressure correction
	bool AMP = true;                 // using of AMP parallel velocity and force calculation

	// parameters for many particles
	double rho = 3.;                 // default density of Solid
	int shape = 0;                   // shape
	double r = 0.05;                 // default radius of Solid
	double r0 = 0.025;               // default little radius of Solid
	double e = 0.;                   // default eccentricity of Solid
	int Nn_max = 0;                  // number of nodes
	int Nn_ = 12;                    // number of nodes for basic segment of added Solid
	int AddSolids_N = 0;             // number of added Solids
	int AddSolids_start = 0;         // step when Solids start to add
	int AddSolids_interval = 5000000;// interval for Solids adding
	int SolidName_max = 0;           // Maximal Name of Solids
	double k_u_dist = 4.0;           // coefficient for minimal distance between Solids for collision force $au_collide$
	double k_r_dist = 1.0;           // coefficient for minimal distance between Solids for collision force $ar_collide$
	double k_u_collide = 1.;         // coefficient in collision force $au_collide$
	double k_r_collide = 1.;         // coefficient in collision force $ar_collide$

	// parameters for special problems
	double Lamb_Oseen_r0 = 0.1;      // initial radius of Lamb_Oseen vortex
	double u_in = 0.;                // velocity of down channel wall
	double u_down = 0.;              // velocity of down channel wall
	double u_up = 0.;                // velocity of up   channel wall
	double omega_BC = 0.;            // omega of the circle driven by external force

	// calculated parameters (dependent on basic ones)
	int N1_u, N2_u;           // sizes for u-direction arrays
	int N1_v, N2_v;           // sizes for v-direction arrays
	int N1_p, N2_p;           // sizes for pressure (cell centers) arrays
	double d_x;               // mesh step in x-direction
	double d_y;               // mesh step in y-direction
	double ldxdx, ldydy;      // coefficients in Laplace approximation
	GeomVec x0 = ZeroVec();   // Special point in test problems Lamb_Oseen and Line_Vortex
	GeomVec k;                // Spatial frequencies for Taylor_Green vortices
	GeomVec Gravity;          // Gravity vector

	//constructors
	Param();

	//methods
	void read(std::ifstream &input);
	void init();
	void read_line(std::string);
};

class InfluenceArea
{
public:
	int i_min, i_max, j_min, j_max;
};

boundary_conditions string_to_BC(std::string s);
GeomVec x_p(int i, int j, Param par); // coordinates of (i,j)-th node for pressure p mesh
GeomVec x_u(int i, int j, Param par); // coordinates of (i,j)-th node for velocity u mesh
GeomVec x_v(int i, int j, Param par); // coordinates of (i,j)-th node for velocity v mesh
GeomVec x_c(int i, int j, Param par); // coordinates of (i,j)-th node for corners    mesh

double* x_p_(int i, int j, double d_x, double d_y) restrict(amp);
double* x_u_(int i, int j, double d_x, double d_y) restrict(amp);
double* x_v_(int i, int j, double d_x, double d_y) restrict(amp);
double* x_c_(int i, int j, double d_x, double d_y) restrict(amp);

double DeltaFunction(double x, double y, double d_x, double d_y);
double DeltaFunction_(double x, double y, double d_x, double d_y) restrict(amp);
double FunctionD(double r);
double FunctionD_(double r) restrict(amp);
void GetInfluenceArea(int &i_min, int &i_max, int &j_min, int &j_max, size_t Ni, size_t Nj, GeomVec x, int size, Param par);

InfluenceArea GetInfluenceArea_(int Ni, int Nj, double* x, int size, boundary_conditions BC, double d_x, double d_y)  restrict(amp);
double ellipse_rho(double r, double e, double phi);
double Volume_Frac(GeomVec xc, double r, double e, double alpha0, GeomVec x, double dx, double dy);
double Heaviside(double x);
int i_real_u(int i, Param par);
int i_real_u_(int i, int N1) restrict(amp);
int i_real_v(int i, Param par);
int i_real_v_(int i, int N1) restrict(amp);

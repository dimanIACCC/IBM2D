#pragma once

#include "GeomVec.h"
#include "Parameters.h"
#include "Matrix.h"
#include "Exact_solutions.h"


enum SolidShape {
	circle, line
};

class Node
{
public:
	GeomVec x, x_n, x_s; // coordinates
	GeomVec uf;       // velocity of the fluid in the Node
	GeomVec us;       // velocity of the Solid in the Node
	GeomVec f;        // force
	double ds;        // size of Lagrange element
	double p;         // pressure
};

class Solid
{
public:
	int moving;              // 0 - Solid does not move, 1 - Solid is driven by the fluid, 2 - Solid moves according to the given law
	GeomVec x_n, x, x_n_plt, x_plt;          // coordinates of the mass center
	GeomVec u_n, u, u_s;          // velocity of the mass center
	GeomVec omega, omega_n, omega_s;  // angular velocity
	GeomVec alpha;                    // angular orientation
	GeomVec d_uv_collide, d_ur_collide, d_omega_collide;   // velocity and andgular velocity corrections due to the collision
	GeomVec f_new, f;      // force applied to the whole Solid
	GeomVec tau_new, tau;    // torque, moment of force applied to the whole Solid
	double Fr;      // average radial Force applied to Solid
	double S;       // length of contour for radial Force averaging
	double I;       // moment of inertia
	double rho;     // density
	double V;       // volume
	GeomVec integralV_du_dt;      // integral of the du/dt       over the Volume of the Solid
	GeomVec integralV_dur_dt;     // integral of the d(u x r)/dt over the Volume of the Solid
	//std::vector<Node> Nodes;   // Nodes of the Solid mesh
	std::vector<int> IndNodes; // Indices of the Nodes
	int Nn;                    // Number of Nodes
	int name;                     // integer name of the Solid
	int shape;                    // shape of the Solid
	double r;                     // radius
	double e;                     // eccentricity 
	Solid(double x, double y, double ux, double uy, double alpha, double omega, double rho, int Nn, int moving, int name, int shape, double r, double e);
	Solid(double x, double y, Param &par);
	~Solid();
	void add_Nodes(std::vector<Node> &Nodes, const int Nn_max);
	void log_init(std::string WorkDir);
	void log(std::string WorkDir, int n);
	void integrals(Matrix U_n, Matrix V_n, Matrix U_new, Matrix V_new, Param par);
};

void fill_solid_coordinates(std::vector<Node> &Nodes, const int Nn_max, const int Nn, const int shape, const double r, const double e, const double alpha, const double dxy);
void line_segment(std::vector<Node> &Nodes, const int N_start, const int Nn, GeomVec Xbeg, GeomVec Xend);
void circular_segment(std::vector<Node> &Nodes, const int N_start, const int Nn, GeomVec Xc, double r, double alpha_beg, double alpha_end);
void copy_solid_mesh(std::vector<Node> &Nodes, const int N_beg_from, const int N_beg_to, const int Nn);
void fill_solid_ds(std::vector<Node> &Nodes, const int Nn_max, const int Nn, const int shape, const double dxy);
void Read_Solids(std::string filename, std::vector<Solid>& Solids, std::vector<Node>& Nodes, Param &par);
void Add_Solids(std::vector<Solid>& Solids, std::vector<Node>& Nodes, Param &par);
bool Collide(Solid& s1, Solid& s2, std::vector<Node> &Nodes, Param par, double alpha, double beta, double friction, double kr, double F_collide);
void Solids_move(std::vector<Solid> &solidList, std::vector<Node> &Nodes, Param par);
void Solids_collide(std::vector<Solid> &solidList, std::vector<Node> &Nodes, Param par);
void h_average_of_Solids_Layer(std::vector<Solid> &solidList, Param par, double& h_average);
void Solids_zero_force(std::vector<Solid>& Solids, std::vector<Node> &Nodes, int N_max);
void Solids_velocity_new(std::vector<Solid>& Solids, Param par);
void Solids_position_new(std::vector<Solid>& Solids, std::vector<Node> &Nodes, Param par);

void velocities(std::vector<Solid>::iterator &Solid, std::vector<Node> &Nodes);      // calculates the velocities  in all Nodes of the Solid
void coordinates(std::vector<Solid>::iterator &Solid, std::vector<Node> &Nodes);     // calculates the coordinates of all Nodes of the Solid in global coordinate system


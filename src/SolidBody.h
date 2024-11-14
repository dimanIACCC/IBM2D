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
	GeomVec f_r_collide; // force of particle collision
	double ds;        // size of Lagrange element
	double p;         // pressure
};

class Solid
{
public:
	int moving = 1;              // 0 - Solid does not move, 1 - Solid is driven by the fluid, 2 - Solid moves according to the given law
	GeomVec x_n = ZeroVec();     // mass center at $n$ step
	GeomVec x   = ZeroVec();     // mass center at $n+1$ step
	GeomVec x_s = ZeroVec();     // mass center at $n+1$ step $s$ iteration
	GeomVec x_n_plt = ZeroVec(); // mass center at $n$ step
	GeomVec x_plt   = ZeroVec(); // mass center at $n+1$ step
	GeomVec u_n = ZeroVec();     // velocity of the mass center at $n$ step
	GeomVec u   = ZeroVec();     // velocity of the mass center at $n+1$ step
	GeomVec u_s = ZeroVec();     // velocity of the mass center at $n+1$ step $s$ iteration
	GeomVec omega   = ZeroVec(); // angular velocity at $n$ step
	GeomVec omega_n = ZeroVec(); // angular velocity at $n+1$ step
	GeomVec omega_s = ZeroVec(); // angular velocity at $n+1$ step $s$ iteration
	GeomVec alpha_n = ZeroVec();   // angular orientation
	GeomVec alpha   = ZeroVec();   // angular orientation
	GeomVec alpha_s = ZeroVec();   // angular orientation
	GeomVec f_collide = ZeroVec();         // acceleration of the collision
	GeomVec tau_collide = ZeroVec();       // angular acceleration of the collision
	GeomVec f_new = ZeroVec();             // force applied to the whole Solid
	GeomVec f     = ZeroVec();             // force applied to the whole Solid
	GeomVec tau_new = ZeroVec();           // torque, moment of force applied to the whole Solid
	GeomVec tau     = ZeroVec();           // torque, moment of force applied to the whole Solid
	double Fr = 0.;      // average radial Force applied to Solid
	double S;       // length of contour for radial Force averaging
	double I;       // moment of inertia
	double rho;     // density
	double V;       // volume
	GeomVec integralV_du_dt;      // integral of the du/dt       over the Volume of the Solid
	GeomVec integralV_dur_dt;     // integral of the d(u x r)/dt over the Volume of the Solid
	//std::vector<Node> Nodes;   // Nodes of the Solid mesh
	std::vector<int> IndNodes; // Indices of the Nodes
	int Nn_r0, Nn_r;           // Number of Nodes for $r0$ and $r$ segments of particle mesh
	int Nn;                    // Number of Nodes
	int name;                     // integer name of the Solid
	int shape;                    // shape of the Solid
	double r;                     // radius
	double r0;                    // little radius
	double e;                     // eccentricity
	bool Poiseuille = false;      // key for initial ux, uy and omega_new corresponding to Poiseuille flow
	//constructor
	Solid(Param &par);
	//Methods
	void Init();
	void read_line(std::string);
	void add_Nodes(std::vector<Node> &Nodes, const int Nn_max);
	void log_init(std::string WorkDir);
	void log(std::string WorkDir, int n);
	void integrals(Matrix U_n, Matrix V_n, Matrix U_new, Matrix V_new, Param par);
};

void fill_solid_coordinates(std::vector<Node> &Nodes, const int Nn_max, const int Nn, const int Nn_e, const int Nn_r, const int shape, const double r0, const double r, const double e, const double alpha);
void line_segment(std::vector<Node> &Nodes, const int N_start, const int Nn, GeomVec Xbeg, GeomVec Xend);
void circular_segment(std::vector<Node> &Nodes, const int N_start, const int Nn, GeomVec Xc, double r, double alpha_beg, double alpha_end);
void copy_solid_mesh(std::vector<Node> &Nodes, const int N_beg_from, const int N_beg_to, const int Nn);
void fill_solid_ds(std::vector<Node> &Nodes, const int Nn_max, const int Nn, const int shape, const double dxy);
void Read_Solids(std::ifstream &input, std::vector<Solid>& Solids, std::vector<Node> &Nodes, Param &par);
void Add_Solids(std::vector<Solid>& Solids, std::vector<Node>& Nodes, Param &par);
void Collide_2Solids(Solid& s1, Solid& s2, std::vector<Node> &Nodes, Param par, double dist_u, double dist_r, double alpha, double beta, double friction);
void Collide_Walls(Solid& s1, std::vector<Node> &Nodes, Param& par, double dist_u, double dist_r, double alpha, double beta, double friction);
void Solids_move(std::vector<Solid> &solidList, std::vector<Node> &Nodes, Param par);
void Solids_collide(std::vector<Solid> &solidList, std::vector<Node> &Nodes, Param par);
void h_average_of_Solids_Layer(std::vector<Solid> &solidList, Param par, double& h_average);
void Solids_zero_force(std::vector<Solid>& Solids, std::vector<Node> &Nodes, int N_max);
void Solids_velocity_new(std::vector<Solid>& Solids, Param par);
void Solids_position_new(std::vector<Solid>& Solids, std::vector<Node> &Nodes, Param par);

void velocities(std::vector<Solid>::iterator &Solid, std::vector<Node> &Nodes);      // calculates the velocities  in all Nodes of the Solid
void coordinates(std::vector<Solid>::iterator &Solid, std::vector<Node> &Nodes);     // calculates the coordinates of all Nodes of the Solid in global coordinate system


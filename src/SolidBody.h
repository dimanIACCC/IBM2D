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
	GeomVec us;       // velocity of the SolidBody in the Node
	GeomVec f;        // force
	GeomVec n;        // norm
	GeomMat Eps;      // deformation velocity
	double ds;        // size of Lagrange element
	double p;         // pressure
};

class SolidBody
{
public:
	int moving;              // 0 - Solid does not move, 1 - Solid is driven by the fluid, 2 - Solid moves according to the given law
	GeomVec x_n, x;          // coordinates of the mass center
	GeomVec u_n, u, u_s;          // velocity of the mass center
	GeomVec omega, omega_n, omega_s;  // angular velocity
	GeomVec d_uv_collide, d_ur_collide, d_omega_collide;   // velocity and andgular velocity corrections due to the collision
	GeomVec f_new, f;      // force applied to the whole SolidBody
	GeomVec tau_new, tau;    // torque, moment of force applied to the whole SolidBody
	double Fr;      // average radial Force applied to SolidBody
	double S;       // length of contour for radial Force averaging
	double I;       // moment of inertia
	double rho;     // density
	double V;       // volume
	GeomVec integralV_du_dt;      // integral of the du/dt       over the Volume of the Solid
	GeomVec integralV_dur_dt;     // integral of the d(u x r)/dt over the Volume of the Solid
	//std::vector<Node> Nodes;   // Nodes of the SolidBody mesh
	std::vector<int> IndNodes; // Indices of the Nodes
	int Nn;                    // Number of Nodes
	int name;                     // integer name of the Solid
	//int shape;                    // shape of the Solid
	double r;
	SolidBody(double x, double y, double ux, double uy, double omega, double rho, int Nn, int moving, int name);
	~SolidBody();
	void add_Nodes(std::vector<Node> &Nodes, const int Nn_max);
	void log_init(std::string WorkDir);
	void log(std::string WorkDir, int n);
};

class Circle : public SolidBody{
public:
	Circle(double x, double y, double ux, double uy, double omega, double rho, int Nn, int moving, int name, double r);
	Circle(double x, double y, Param &par);
	~Circle();
	void integrals(Matrix U_n, Matrix V_n, Matrix U_new, Matrix V_new, Param par);
};

void fill_circle_coordinates(std::vector<Node> &Nodes, const int Nn_max, const int Nn, double r, double e);
void fill_ds_x(std::vector<Node> &Nodes, const int Nn_max, const int Nn, const double r, const double d_x, const double d_y);

void Read_Solids(std::string filename, std::vector<Circle>& Solids, std::vector<Node>& Nodes, Param &par);
void Add_Solids(std::vector<Circle>& Solids, std::vector<Node>& Nodes, Param &par);
bool Collide(Circle& s1, Circle& s2, std::vector<Node> &Nodes, Param par, double alpha, double beta, double friction, double kr);
void Solids_move(std::vector<Circle> &solidList, std::vector<Node> &Nodes, Param par);
void Solids_collide(std::vector<Circle> &solidList, std::vector<Node> &Nodes, Param par);
void h_average_of_Solids_Layer(std::vector<Circle> &solidList, Param par, double& h_average);
void Solids_zero_force(std::vector<Circle>& Solids, std::vector<Node> &Nodes, int N_max);
void Solids_velocity_new(std::vector<Circle>& Solids, Param par);
void Solids_position_new(std::vector<Circle>& Solids, std::vector<Node> &Nodes, Param par);

void velocities(std::vector<Circle>::iterator &Solid, std::vector<Node> &Nodes);      // calculates the velocities  in all Nodes of the SolidBody
void coordinates(std::vector<Circle>::iterator &Solid, std::vector<Node> &Nodes);     // calculates the coordinates of all Nodes of the SolidBody in global coordinate system

GeomVec Circle_Equation(GeomVec xc, double r, double theta);

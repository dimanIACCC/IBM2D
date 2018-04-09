#pragma once

#include "GeomVec.h"
#include "Parameters.h"
#include "Matrix.h"

class Node
{
public:
	GeomVec x, xn;    // coordinates
	GeomVec uf;       // velocity of the fluid in the Node
	GeomVec us;       // velocity of the SolidBody in the Node
	GeomVec f, f_tmp; // force and temporary force in iterations
	GeomVec n;        // norm
	GeomVec t;        // traction vector
	GeomMat Eps;      // deformation velocity
	double p;         // pressure
};

class SolidBody
{
public:
	bool moving;
	GeomVec xc, xc_n;        // coordinates of the mass center
	GeomVec uc, uc_n, uc_s;  // velocity of the mass center
	GeomVec omega, omega_n, omega_s;  // angular velocity
	GeomVec f;      // force applied to the whole SolidBody
	double Fr, Fr_all;      // average radial Force applied to SolidBody
	GeomVec F_hd;   // Force calculated from hydrodynamics
	GeomVec tau_hd; // torque, moment of force calculated from hydrodynamics
	double S;       // length of contour for radial Force averaging
	GeomVec tau;    // torque, moment of force applied to the whole SolidBody
	double I;       // moment of inertia
	double rho;     // density
	double V;       // volume
	GeomVec integralV_du_dt;      // integral of the du/dt       over the Volume of the Solid
	GeomVec integralV_dur_dt;     // integral of the d(u x r)/dt over the Volume of the Solid
	std::vector<Node> Nodes;   // Nodes of the SolidBody mesh
	size_t Nn;                    // Number of Nodes
	int name;                     // integer name of the solid
	
	SolidBody(double x, double y, double ux, double uy, double omega, double rho, int Nn, bool moving, int &name);
	~SolidBody();
	void velocities();      // calculates the velocities in all Nodes of the SolidBody
	double ds(size_t i);
	void log_init(std::string WorkDir);
	void log(std::string WorkDir, int n);
};

class Circle : public SolidBody{
public:
	double r;
	Circle(double x, double y, double ux, double uy, double omega, double rho, int Nn, bool moving, int name, double r);
	Circle(double x, double y, Param &par);
	~Circle();
	void integrals(Matrix U_n, Matrix V_n, Matrix U_new, Matrix V_new, Param par);
};

void Read_Solids(std::string filename, std::list<Circle>& Solids, Param &par);
void Add_Solids(std::list<Circle>& Solids, int n, Param &par);
bool Collide(Circle& s1, Circle& s2, Param par);
void Solids_move(std::list<Circle> &solidList, Param par, int n);
void Solids_zero_force(std::list<Circle>& Solids);
void Solids_velocity_new(std::list<Circle>& Solids, Param par);
GeomVec Circle_Equation(GeomVec xc, double r, double theta);

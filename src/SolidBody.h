#pragma once

#include "GeomVec.h"
#include "Parameters.h"

class Node
{
public:
	GeomVec x;        // coordinates
	GeomVec uf;       // velocity of the fluid in the Node
	GeomVec us;       // velocity of the SolidBody in the Node
	GeomVec f, f_tmp; // force and temporary force in iterations
	GeomVec Integral; // value of integral when calculating force
};

class SolidBody
{
public:
	bool moving;
	GeomVec xc;     // coordinates of the mass center
	GeomVec uc, uc_n;        // velocity of the mass center
	GeomVec omega, omega_n;  // angular velocity
	GeomVec f;      // force applied to the whole SolidBody
	double Fr, Fr_all;      // average radial Force applied to SolidBody
	double S;       // length of contour for radial Force averaging
	GeomVec tau;    // torque, moment of force applied to the whole SolidBody
	double I;       // angular momentum
	double rho;     // density
	double V;       // volume
	std::vector<Node> Nodes;   // Nodes of the SolidBody mesh
	size_t Nn;                    // Number of Nodes
	int name;                     // integer name of the solid
	
	SolidBody(double x, double y, double ux, double uy, double omega, double rho, int Nn, bool moving, int &name);
	~SolidBody();
	void velocities();      // calculates the velocities in all Nodes of the SolidBody
	void move(double d_t); // move Solid using $uc$ and $omega$
	double ds(size_t i);
	void log_init(std::string WorkDir);
	void log(std::string WorkDir, int n);
};

class Circle : public SolidBody{
public:
	double r;
	Circle(double x, double y, double ux, double uy, double omega, double rho, int Nn, bool moving, int &name, double r);
	Circle(double x, double y, Param &par);
	~Circle();
};

void Read_Solids(std::string filename, std::list<Circle>& Solids, Param &par);
void Add_Solids(std::list<Circle>& Solids, int n, Param &par);
bool Collide(Circle& s1, Circle& s2, Param par);
void Solids_move(std::list<Circle> &solidList, Param par, int n);
void Solids_zero_force(std::list<Circle>& Solids);
void Solids_velocity_new(std::list<Circle>& Solids, Param par);

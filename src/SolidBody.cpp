#include "SolidBody.h"


SolidBody::SolidBody(double x, double y, double ux, double uy, double omega, double rho, int Nn, bool moveSolid)
{
	std::fill(this->xc.begin()   , this->xc.end(),    0.0); // fill vector xc    with zeros
	std::fill(this->uc.begin()   , this->uc.end()   , 0.0); // fill vector uc    with zeros
	std::fill(this->omega.begin(), this->omega.end(), 0.0); // fill vector omega with zeros
	this->xc[1] = x;
	this->xc[2] = y;
	this->uc[1] = ux;
	this->uc[2] = uy;
	this->omega[3] = omega;
	this->rho = rho;
	this->Nn = Nn;
	Nodes.resize(Nn);
	moveSolid = false;
}


SolidBody::~SolidBody()
{

}

Circle::Circle(double x, double y, double ux, double uy, double omega, double rho, int Nn, bool moveSolid, double r) :
     SolidBody(       x,        y,        ux,        uy,        omega,        rho,     Nn,      moveSolid) {
	this->r = r;
	this->d_s = (2.0*M_PI*r) / Nn;

	for (int i = 0; i < Nn; ++i){
		Nodes[i].x[1] = x + cos(i * 2.0 * M_PI / Nn) * r;
		Nodes[i].x[2] = y + sin(i * 2.0 * M_PI / Nn) * r;
	}
	V = M_PI * r * r;
	//rho = 10.0 / V; // corresponds to old formula for force
	I =  V * r * r / 2.0; // angular momentum for unit density
}

Circle::Circle(double x, double y, Param par): Circle(x, y, 0.0, 0.0, 0.0, par.rho, par.Nn, false, par.r) {
	this->uc[1] = ux_Poiseuille(y, par.H);
}
 
Circle::~Circle()
{
}

void SolidBody::velocities() {
	for (int i = 0; i < Nn; ++i) {
		GeomVec r;
		r = Nodes[i].x - xc;
		Nodes[i].us = uc + x_product(omega, r);
	}
}

void Read_Solids(std::string filename, std::list<Circle>& Solids, Param par) {
	std::ifstream input;
	std::string line;

	input.open(filename.c_str());
	if (input.is_open()) {
		while (getline(input, line)) { // read line from file to string $line$
			if (line == "circle{") {
				double x = par.L*0.1;
				double y = par.H*0.5;
				double ux = ux_Poiseuille(y, par.H);
				double uy = 0;
				double omega = 0;
				double rho = par.rho;
				int Nn = par.Nn;
				bool move = false;
				double r = par.r;

				while (line != "}") {
					getline(input, line);
					if (line == "}") break;
					std::string PAR, VALUE;
					GetParValue(line, PAR, VALUE);
					if (VALUE.size() > 0) {
						if      (PAR == "x")          x            = stod(VALUE);
						else if (PAR == "y")          y            = stod(VALUE);
						else if (PAR == "ux")         ux           = stod(VALUE);
						else if (PAR == "uy")         uy           = stod(VALUE);
						else if (PAR == "omega")      omega        = stod(VALUE);
						else if (PAR == "rho")        rho          = stod(VALUE);
						else if (PAR == "Nn")         Nn           = stoi(VALUE);
						else if (PAR == "move")       move         = stoi(VALUE);
						else if (PAR == "r")          r            = stod(VALUE);
						else    std::cout << "Read_Solids: unknown parameter " << PAR << std::endl;
					}
					else {
						std::cout << "Read_Solids: no value inputed" << std::endl;
					}
				}
				ux = ux_Poiseuille(y, par.H);  // recalculate $ux$ for new $y$
				Circle c(x, y, ux, uy, omega, rho, Nn, move, r);
				Solids.push_back(c);
			}
		}
	}
	else {
		std::cout << "Read_Solids: File " << filename << " is not found" << std::endl;
	}

}

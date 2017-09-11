#include "SolidBody.h"


SolidBody::SolidBody(double x, double y, double ux, double uy, double omega, double rho, int Nn, bool moving)
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
	this->moving = moving;
}


SolidBody::~SolidBody()
{

}

Circle::Circle(double x, double y, double ux, double uy, double omega, double rho, int Nn, bool moving, double r) :
     SolidBody(       x,        y,        ux,        uy,        omega,        rho,     Nn,      moving) {
	this->r = r;
	this->d_s = (2.0*M_PI*r) / Nn;

	for (int i = 0; i < Nn; ++i){
		Nodes[i].x[1] = x + cos(i * 2.0 * M_PI / Nn) * r;
		Nodes[i].x[2] = y + sin(i * 2.0 * M_PI / Nn) * r;
	}
	V = M_PI * r * r;
	I =  V * r * r / 2.0; // angular momentum for unit density
}

Circle::Circle(double x, double y, Param par):
	    Circle(       x,        y, 0.0, 0.0, 0.0, par.rho, par.Nn, true, par.r) {
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

void SolidBody::move(double d_t) {
	if (moving) {
		//update position
		for (int k = 0; k < Nn; ++k) {
			//rotate
			GeomVec r = Nodes[k].x - xc;
			GeomVec x_temp = rotate_Vector_around_vector(r, omega  * length(r) * d_t); //
			Nodes[k].x = xc + x_temp; // rotate solid by angle $omega$ * $dt$
			Nodes[k].x += uc * d_t; // move
		}
		xc += uc * d_t;
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
				bool moving = true;
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
						else if (PAR == "moving")     moving       = bool(stoi(VALUE));
						else if (PAR == "r")          r            = stod(VALUE);
						else    std::cout << "Read_Solids: unknown parameter " << PAR << std::endl;
					}
					else {
						std::cout << "Read_Solids: no value inputed" << std::endl;
					}
				}
				if (moving == false) {
					ux = 0;
					uy = 0;
				}
				Circle c(x, y, ux, uy, omega, rho, Nn, moving, r);
				Solids.push_back(c);
			}
		}
	}
	else {
		std::cout << "Read_Solids: File " << filename << " is not found" << std::endl;
	}

}

void Add_Solids(std::list<Circle>& Solids, int nSolids, int n, int n_start, int n_interval, Param par) {
	if ((n % n_interval) == n_start) { //create new solids starting from $n_start$ iteration with interval of $n_interval$ iterations
		for (int i = 0; i < nSolids; i++) { // add $nSolids$ solids
			GeomVec x;
			x[0] = 0;
			x[1] = par.L / 10 + par.L / 10 * 0.95 * (double(rand()) - RAND_MAX / 2) / RAND_MAX;
			x[2] = par.H / 2  + par.H / 2  * 0.95 * (double(rand()) - RAND_MAX / 2) / RAND_MAX;
			x[3] = 0;
			Circle c(x[1], x[2], par);
			// check if new Solid does not cross other Solids
			bool add = true;
			for (auto solid = Solids.begin(); solid != Solids.end(); solid++) {
				if (length(x - solid->xc) < par.k_dist * (par.r + solid->r)) {
					add = false;
					break;
				}
			}
			if (add) {
				Solids.push_back(c);
			}
		}
	}

}

bool Collide(Circle& s1, Circle& s2, Param par) {
	bool result = false;
	GeomVec r = s1.xc - s2.xc;
	double distance = length(r); //<----distance between two particles
	if (distance <= par.k_dist*(s1.r + s2.r)) {
		result = true;
		if (par.InelasticCollision) { //Perfectly inelastic collision
			s1.uc = (s1.uc + s2.uc) / 2;
			s2.uc = s1.uc;
		}
		else { //Perfectly elastic collision
			r = r / distance;
			double u1_before = dot_product(s1.uc, r);
			double u2_before = dot_product(s2.uc, r);
			if (u1_before - u2_before < 0.0) { // if Solids move to each other
				double m1 = s1.rho * s1.V;
				double m2 = s2.rho * s2.V;
				double u1_after = (u1_before * (m1 - m2) + 2 * m2 * u2_before) / (m1 + m2);
				double u2_after = (u2_before * (m2 - m1) + 2 * m1 * u1_before) / (m1 + m2);
				s1.uc += (u1_after - u1_before) * r;
				s2.uc += (u2_after - u2_before) * r;
			}
		}
	}
	return result;
}

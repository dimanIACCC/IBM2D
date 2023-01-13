#include "stdafx.h"
#include "SolidBody.h"

SolidBody::SolidBody(double x, double y, double ux, double uy, double omega, double rho, int Nn, int moving, int name)
{
	std::fill(this->x_n.begin()   , this->x_n.end(),    0.0); // fill vector x_n    with zeros
	std::fill(this->u_n.begin()   , this->u_n.end()   , 0.0); // fill vector u_n    with zeros
	std::fill(this->omega_n.begin(), this->omega_n.end(), 0.0); // fill vector omega_n with zeros
	this->x_n[1] = x;
	this->x_n[2] = y;
	this->u_n[1] = ux;
	this->u_n[2] = uy;
	this->omega_n[3] = omega;
	this->rho = rho;
	this->Nn = Nn;
	this->moving = moving;
	this->name = name;
	this->x     = this->x_n;
	this->u     = this->u_n;
	this->omega = this->omega_n;
	this->Fr = 0.;
}


SolidBody::~SolidBody()
{

}

void SolidBody::add_Nodes(std::vector<Node> &Nodes, const int Nn_max) {
	this->IndNodes.resize(this->Nn);
	Nodes.resize(Nn_max + this->Nn);
	for (size_t i = 0; i < Nn; ++i) {
		int Ind = Nn_max + i;
		this->IndNodes[i] = Ind;
	}
}

GeomVec Circle_Equation(GeomVec xc, double r, double theta) {
	GeomVec xy;
	xy[1] = xc[1] + cos(theta) * r;
	xy[2] = xc[2] + sin(theta) * r;
	xy[3] = 0;
	return xy;
}

bool operator <(const SolidBody& a, const SolidBody& b) {
	if (length(a.x_n) < length(b.x_n)) return true;
	else return false;
}

bool operator >(const SolidBody& a, const SolidBody& b) {
	if (length(a.x_n) > length(b.x_n)) return true;
	else return false;
}

Circle::Circle(double x, double y, double ux, double uy, double omega, double rho,  int Nn, int moving, int name, double r) :
     SolidBody(       x,        y,        ux,        uy,        omega,        rho,      Nn,     moving,     name) {
	this->r = r;
	V = M_PI * r * r;
	I =  V * r * r / 2.0; // angular momentum for unit density
}

Circle::Circle(double x, double y, Param &par):
	    Circle(       x,        y, 0.0, 0.0, 0.0, par.rho, par.Nn, true, par.SolidName_max+1, par.r) {
	this->u_n[1] = 0.0;  //  ux_Poiseuille(y, par.H);
	this->omega_n[3] = 0.0; // -dux_dy_Poiseuille(y, par.H);
	this->omega     = this->omega_n;
	this->u    = this->u_n;
}
 
Circle::~Circle()
{
}

void fill_solid_coordinates(std::vector<Node> &Nodes, const int Nn_max, const int Nn, const double r, const double e, const double alpha, const double dxy) {
	GeomVec o;
	o[3] = M_PI / 180 * alpha;

	for (size_t i = 0; i < Nn; ++i) {
		int Ind = Nn_max + i;
		if (e < 0.999) {    // ellipse
			double phi = i * 2.0 * M_PI / Nn;
			Nodes[Ind].x_n[1] = cos(phi) * r / sqrt(1 - e*e*cos(phi)*cos(phi)) * sqrt(1 - e*e);
			Nodes[Ind].x_n[2] = sin(phi) * r / sqrt(1 - e*e*cos(phi)*cos(phi)) * sqrt(1 - e*e);
			//Nodes[Ind].ds = 2.0 * M_PI * r / Nn * dxy;
		}
		else{  // line
			double x = 2 * r * (i - 0.5*Nn) / Nn;
			Nodes[Ind].x_n[1] = x;
			Nodes[Ind].x_n[2] = 0.;
			//Nodes[Ind].ds = 2. * r / Nn * dxy;
		}
		Nodes[Ind].x_n = rotate_Vector_around_vector(Nodes[Ind].x_n, o);
		Nodes[Ind].x = Nodes[Ind].x_n;
	}

	for (size_t i = 1; i < Nn-1; ++i) {
		int Ind = Nn_max + i;
		Nodes[Ind].ds = 0.5 * length(Nodes[Ind + 1].x_n - Nodes[Ind - 1].x_n) * dxy;
	}
	if (e < 0.999) {    // ellipse
		Nodes[Nn_max + 0     ].ds = 0.5 * length(Nodes[Nn_max + 1].x_n - Nodes[Nn_max + Nn - 1].x_n) * dxy;
		Nodes[Nn_max + Nn - 1].ds = 0.5 * length(Nodes[Nn_max + 0].x_n - Nodes[Nn_max + Nn - 2].x_n) * dxy;
	}
	else {  // line
		Nodes[Nn_max + 0     ].ds = length(Nodes[Nn_max      + 1].x_n - Nodes[Nn_max      + 0].x_n) * dxy;
		Nodes[Nn_max + Nn - 1].ds = length(Nodes[Nn_max + Nn - 1].x_n - Nodes[Nn_max + Nn - 2].x_n) * dxy;
	}
}

void fill_solid_ds(std::vector<Node> &Nodes, const int Nn_max, const int Nn, const double e, const double dxy) {

	for (size_t i = 1; i < Nn - 1; ++i) {
		int Ind = Nn_max + i;
		Nodes[Ind].ds = 0.5 * length(Nodes[Ind + 1].x_n - Nodes[Ind - 1].x_n) * dxy;
	}
	if (e < 0.999) {    // ellipse
		Nodes[Nn_max + 0].ds = 0.5 * length(Nodes[Nn_max + 1].x_n - Nodes[Nn_max + Nn - 1].x_n) * dxy;
		Nodes[Nn_max + Nn - 1].ds = 0.5 * length(Nodes[Nn_max + 0].x_n - Nodes[Nn_max + Nn - 2].x_n) * dxy;
	}
	else {  // line
		Nodes[Nn_max + 0].ds = length(Nodes[Nn_max + 1].x_n - Nodes[Nn_max + 0].x_n) * dxy;
		Nodes[Nn_max + Nn - 1].ds = length(Nodes[Nn_max + Nn - 1].x_n - Nodes[Nn_max + Nn - 2].x_n) * dxy;
	}

}

void velocities(std::vector<Circle>::iterator &Solid, std::vector<Node> &Nodes) {
	for (size_t k = 0; k < Solid->Nn; ++k) {
		int Ind = Solid->IndNodes[k];
		Nodes[Ind].us = 0.5 * (Solid->u + Solid->u + x_product(Solid->omega+ Solid->omega, Nodes[Ind].x));
	}
}

void coordinates(std::vector<Circle>::iterator &Solid, std::vector<Node> &Nodes) {
	for (size_t k = 0; k < Solid->Nn; ++k) {
		int Ind = Solid->IndNodes[k];
		Nodes[Ind].x_s = Solid->x + Nodes[Ind].x;
	}
}

void SolidBody::log_init(std::string WorkDir) {
	std::ofstream output;
	std::string filename = WorkDir + "Solids/" + std::to_string(name) + ".plt";
	output.open(filename);

	output << "title = " << '"' << "Solid" << name << '"' << std::endl;
	output << "Variables = n x y u v fx fy omega tau" << std::endl;
}

void SolidBody::log(std::string WorkDir, int n) {
	std::ofstream output;
	std::string filename = WorkDir + "Solids/" + std::to_string(name) + ".plt";
	output.open(filename, std::ios::app);

	output << std::setprecision(8);
	output << n << "   " 
		   << x_n[1] << "   "
	       << x_n[2] << "   "
	       << u_n[1] << "   "
	       << u_n[2] << "   "
	       << f_new[1] << "   "
	       << f_new[2] << "   "
	       << omega_n[3] << "   "
	       << tau_new[3]   << "   "
	       << std::endl;
}

void Read_Solids(std::string filename, std::vector<Circle>& Solids, std::vector<Node> &Nodes, Param &par) {
	std::ifstream input;
	std::string line;

	input.open(filename.c_str());
	if (input.is_open()) {
		while (getline(input, line)) { // read line from file to string $line$
			if (line == "circle{") {
				double x = par.L*0.1;
				double y = par.H*0.5;
				double ux = 0;
				double uy = 0;
				double omega = 0;
				double rho = par.rho;
				int Nn = par.Nn;
				int moving = 1;
				double r = par.r;
				double e = par.e;
				double alpha = 0.;
				bool Poiseuille = false;   //key for initial ux, uy and omega_new corresponding to Poiseuille flow

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
						else if (PAR == "moving")     moving       = stoi(VALUE);
						else if (PAR == "Poiseuille") Poiseuille   = bool(stoi(VALUE));
						else if (PAR == "r")          r            = stod(VALUE);
						else if (PAR == "e")          e            = stod(VALUE);
						else if (PAR == "alpha")      alpha        = stod(VALUE);
						else    std::cout << "Read_Solids: unknown parameter " << PAR << std::endl;
					}
					else {
						std::cout << "Read_Solids: no value inputed" << std::endl;
					}
				}
				if (Poiseuille) {
					ux = ux_Poiseuille(y, par.H);
					uy = 0;
					omega = - dux_dy_Poiseuille(y, par.H);
				}
				Circle c(x, y, ux, uy, omega, rho, Nn, moving, par.SolidName_max+1, r);
				c.add_Nodes(Nodes, par.Nn_max);
				fill_solid_coordinates(Nodes, par.Nn_max, c.Nn, c.r, e, alpha, 0.5*(par.d_x+par.d_y));
				par.Nn_max += c.Nn;
				if (c.name > par.SolidName_max) par.SolidName_max = c.name;
				Solids.push_back(c);
				c.log_init(par.WorkDir);
			}
		}
	}
	else {
		std::cout << "Read_Solids: File " << filename << " is not found" << std::endl;
	}

}

void Add_Solids(std::vector<Circle>& Solids, std::vector<Node>& Nodes, Param &par) {
	if (   (par.N_step - par.AddSolids_start) % par.AddSolids_interval == 0   &&   (par.N_step >= par.AddSolids_start)) { //create new solids starting from $AddSolids_start$ iteration with interval of $AddSolids_interval$ iterations
		for (int i = 0; i < par.AddSolids_N; i++) { // add $AddSolids_N$ solids
			GeomVec x;
			x[0] = 0;
			x[1] = (par.L)  * (0.5 + 0.5*(1 - 4 * par.r / par.L) * (double(rand()) - RAND_MAX / 2) / RAND_MAX);
			x[2] = (par.H)  * (0.5 + 0.5*(1 - 6 * par.r / par.H) * (double(rand()) - RAND_MAX / 2) / RAND_MAX);
			x[3] = 0;

			// check if new Solid does not cross other Solids
			bool add = true;
			for (auto solid = Solids.begin(); solid != Solids.end(); solid++) {
				//if (length(x - solid->x) < 1.1 * (par.r + solid->r)) {
				if (length(x - solid->x) < 1.25 * (par.r + solid->r) && solid->moving == 1) { // flow inside Solid
					add = false;
					i--;
					break;
				}
			}
			if (add) {
				Circle c(x[1], x[2], par);
				c.add_Nodes(Nodes, par.Nn_max);
				fill_solid_coordinates(Nodes, par.Nn_max, c.Nn, par.r, par.e, 0., 0.5*(par.d_x + par.d_y));
				par.Nn_max += c.Nn;
				if (c.name > par.SolidName_max) par.SolidName_max = c.name;
				Solids.push_back(c);
				c.log_init(par.WorkDir);
			}
		}
	}

}

double Distance_2Solids(SolidBody& s1, SolidBody& s2, Param& par, GeomVec& r) {
	double dist = 0;
	r = s1.x - s2.x;
	double rrr = length(r);      // distance between particle centers
	if (par.BC == periodical) {
		GeomVec r_plus  = r;
		GeomVec r_minus = r;
		r_plus[1]  += par.L;
		r_minus[1] -= par.L;
		if (length(r_plus ) < rrr) { r = r_plus ; rrr = length(r_plus ); };
		if (length(r_minus) < rrr) { r = r_minus; rrr = length(r_minus); };
	}

	//if (s1.shape == circle && s2.shape == circle) {

		double r1 = std::fmin(s1.r, s2.r);
		double r2 = std::fmax(s1.r, s2.r);

		if      (rrr >= r2) dist = rrr - (r2 + r1);
		else if (rrr <= r2) { dist = (r2 - r1) - rrr; r = -r; } // circle r1 is inside circle r2
	//}
	r = r / rrr;

	return dist;
}

double Distance_2Solids_(SolidBody& s1, SolidBody& s2, std::vector<Node> Nodes, Param& par, GeomVec& r_out, GeomVec& x1, GeomVec& x2) {
	double dist = +1.e99;

	for (size_t k1 = 0; k1 < s1.Nn; ++k1) {
		int Ind1 = s1.IndNodes[k1];
		for (size_t k2 = 0; k2 < s2.Nn; ++k2) {
			int Ind2 = s2.IndNodes[k2];
			GeomVec r = Nodes[Ind1].x_s - Nodes[Ind2].x_s;
			double rrr = length(r);      // distance between particle centers
			if (par.BC == periodical) {
				GeomVec r_plus = r;
				GeomVec r_minus = r;
				r_plus[1] += par.L;
				r_minus[1] -= par.L;
				if (length(r_plus) < rrr) { r = r_plus; rrr = length(r_plus); };
				if (length(r_minus) < rrr) { r = r_minus; rrr = length(r_minus); };
			}
			if (rrr < dist) {
				dist = rrr;
				x1 = Nodes[Ind1].x_s - s1.x_n;
				x2 = Nodes[Ind2].x_s - s2.x_n;
				r_out = r / rrr;
			}
		}
	}

	return dist;
}

bool Collide(Circle& s1, Circle& s2, std::vector<Node> &Nodes, Param par, double alpha, double beta, double friction, double kr) {
	bool result = false;

	GeomVec r = s1.x - s2.x;
	GeomVec x1 = s1.x;
	GeomVec x2 = s2.x;
	GeomVec omega0;
	omega0[3] = 1.;
	double dist = Distance_2Solids(s1, s2, par, r);
	if (dist < 10 * par.d_x) dist = Distance_2Solids_(s1, s2, Nodes, par, r, x1, x2);

	double dist_u = par.k_dist*par.d_x;
	double dist_r = par.k_dist*par.d_x * kr;
	if (dist <= dist_u) {
		result = true;
		double u1_before = dot_product(s1.u_n, r);
		double u2_before = dot_product(s2.u_n, r);
		double omega1_before = s1.omega_n[3];
		double omega2_before = s2.omega_n[3];
		double m1 = s1.rho * s1.V;
		double m2 = s2.rho * s2.V;
		if (u1_before - u2_before < 0.0) { // if Solids move to each other
			double u1_after = (u1_before * (m1 - m2) + 2 * m2 * u2_before) / (m1 + m2);
			double u2_after = (u2_before * (m2 - m1) + 2 * m1 * u1_before) / (m1 + m2);
			double omega1_after = (omega1_before * (m1 - m2) + 2 * m2 * omega2_before) / (m1 + m2);
			double omega2_after = (omega2_before * (m1 - m2) + 2 * m2 * omega1_before) / (m1 + m2);
			if (s1.moving == 1) s1.d_uv_collide += alpha*(u1_after - u1_before)*r;
			if (s2.moving == 1) s2.d_uv_collide += alpha*(u2_after - u2_before)*r;
			if (s1.moving == 1) s1.d_omega_collide[3] += friction*(omega1_after - omega1_before);
			if (s2.moving == 1) s2.d_omega_collide[3] += friction*(omega2_after - omega2_before);
			//if (Debug) std::cout << "Collision u" << std::endl;
		}
		if (dist <= dist_r) {
			double Delta_u = (dist_r - dist)*(dist_r - dist) / dist_r / dist_r * par.Gravity_module * par.d_t;  // distance force
			if (s1.moving == 1) s1.d_ur_collide += beta*m2 / (m1 + m2)*Delta_u*r;
			if (s2.moving == 1) s2.d_ur_collide -= beta*m1 / (m1 + m2)*Delta_u*r;

			if (s1.moving == 1) s1.d_omega_collide[3] += beta*m2 / (m1 + m2)*Delta_u*dot_product(x_product(x1, r), omega0) * s1.V / s1.I;
			if (s2.moving == 1) s2.d_omega_collide[3] -= beta*m1 / (m1 + m2)*Delta_u*dot_product(x_product(x2, r), omega0) * s2.V / s2.I;

			//if (Debug) std::cout << "Collision r" << std::endl;
		}
	}
	return result;
}

void Solids_collide(std::vector<Circle> &solidList, std::vector<Node> &Nodes, Param par) {

	double kr = 0.25;    // fraction of distance where distance force switches on
	double dist_u = par.k_dist*par.d_x;
	double dist_r = par.k_dist*par.d_x * kr;

	double alpha = 1500 * par.d_t; // coefficient for the collision force based on velocity value
	double beta  = 50.; // coefficient for the collision force based on distance between particles value
	double friction = 0.1; // coefficient for the friction force based on velocity value

	///-------------collision with walls--------------------------
	double dist;
	for (auto one = solidList.begin(); one != solidList.end(); one++) {
		if (one->moving == 1){
			dist = 2*(par.H - one->x[2] - one->r); //<-------distance to upper wall
			if (dist < dist_u) {
				if (one->u[2] > 0) one->d_uv_collide[2] -= alpha * one->u_n[2];   // velocity force
				one->d_omega_collide[3] -= friction * one->omega_n[3];
				if (Debug) std::cout << "Collision with upper wall,   Delta_u_v =" << alpha * one->u[2] << std::endl;
			}
			if (dist < dist_r) {
				double Delta_u = (dist_r - dist)*(dist_r - dist) / dist_r / dist_r * par.Gravity_module * par.d_t;  // distance force
				one->d_ur_collide[2] -= beta*Delta_u;
				if (Debug) std::cout << "Collision with upper wall,   Delta_u_r =" << beta*Delta_u << std::endl;
			}
			dist = 2*(one->x[2] - one->r);//<-------distance to lower wall
			if (dist < dist_u) {
				if (one->u[2] < 0) one->d_uv_collide[2] -= alpha * one->u_n[2];   // velocity force
				one->d_omega_collide[3] -= friction * one->omega_n[3];
				if (Debug) std::cout << "Collision with lower wall,   Delta_u_v =" << alpha * one->u[2] << std::endl;
			}
			if (dist < dist_r) {
				double Delta_u = (dist_r - dist)*(dist_r - dist) / dist_r / dist_r * par.Gravity_module * par.d_t;  // distance force
				one->d_ur_collide[2] += beta*Delta_u;
				if (Debug) std::cout << "Collision with lower wall,   Delta_u_r =" << beta*Delta_u << std::endl;
			}

			if (par.BC == box) {
				dist = 2*(par.L - one->x[1] - one->r);//<-------distance to right wall
				if (dist < dist_u) {
					if (one->u[1] > 0) one->d_uv_collide[1] -= alpha * one->u_n[1];  // velocity force
					one->d_omega_collide[3] -= friction * one->omega_n[3];
					if (Debug) std::cout << "Collision with right wall,   Delta_u_v =" << alpha * one->u[2] << std::endl;
				}
				if (dist < dist_r) {
					double Delta_u = (dist_r - dist)*(dist_r - dist) / dist_r / dist_r * par.Gravity_module * par.d_t;  // distance force
					one->d_ur_collide[1] -= beta*Delta_u;
					if (Debug) std::cout << "Collision with right wall,   Delta_u_r =" << alpha * one->u[2] << std::endl;
				}
				dist = 2*(one->x[1] - one->r);//<-------distance to left  wall
				if (dist < dist_u) {
					if (one->u[1] < 0) one->d_uv_collide[1] -= alpha *one->u_n[1];   // velocity force
					one->d_omega_collide[3] -= friction * one->omega_n[3];
					if (Debug) std::cout << "Collision with left wall,   Delta_u_v =" << alpha * one->u[2] << std::endl;
				}
				if (dist < dist_r) {
					double Delta_u = (dist_r - dist)*(dist_r - dist) / dist_r / dist_r * par.Gravity_module * par.d_t;  // distance force
					one->d_ur_collide[1] += beta*Delta_u;
					if (Debug) std::cout << "Collision with left wall,   Delta_u_r =" << alpha * one->u[2] << std::endl;
				}
			}
		}
	}

	///--------------collisions between particles-----------------
	for (auto one = solidList.begin(); one != solidList.end(); one++) {
		for (auto two = next(one); two != solidList.end(); two++) {
			if (Collide(*one, *two, Nodes, par, alpha, beta, friction, kr)) if (Debug) std::cout << "Collision of particles" << std::endl;
		}
	}

}

void Solids_move(std::vector<Circle> &solidList, std::vector<Node> &Nodes, Param par) {
	for (auto it = solidList.begin(); it != solidList.end();) {
			if (it->moving > 0) {
				it->u_n    = it->u;
				it->omega_n = it->omega;
				it->x_n    = it->x;
				for (size_t k = 0; k < it->Nn; ++k) {
					int Ind = it->IndNodes[k];
					Nodes[Ind].x_n = Nodes[Ind].x;
				}
			}
			it->log(par.WorkDir, par.N_step);

			//Right boundary conditions for Solids
			if (it->x_n[1] < par.L) {
				it++;
			}
			else {
				if (par.BC == periodical) {
					it->x_n[1] -= par.L;
					it++;
				}
				else
					it = solidList.erase(it);
			}
		}
	}


void h_average_of_Solids_Layer(std::vector<Circle> &solidList, Param par, double& h_average) {
	//solidList.sort(std::greater<Circle>());
	sort(solidList.begin(), solidList.end(), std::greater<Circle>());
	h_average = 0.;
	int i = 0;
	int i_max = 10;

	for (auto it = solidList.begin(); it != solidList.end();) {
		if (i < i_max &&  par.L - it->x_n[1] > it->r * 5.0 ) {
			i++;
			h_average += length(it->x_n);
		}
		it++;
	}
	h_average /= i_max;
}


void Solids_zero_force(std::vector<Circle>& Solids, std::vector<Node>& Nodes, int Nn_max) {
	for (auto& it : Solids) {
		std::fill(it.f.begin(), it.f.end(), 0.0);
		std::fill(it.tau.begin(), it.tau.end(), 0.0);
		std::fill(it.d_uv_collide.begin(), it.d_uv_collide.end(), 0.0);
		std::fill(it.d_ur_collide.begin(), it.d_ur_collide.end(), 0.0);
		std::fill(it.d_omega_collide.begin(), it.d_omega_collide.end(), 0.0);
		
	}
	for (size_t k = 0; k < Nn_max; ++k) {
		std::fill(Nodes[k].f.begin(), Nodes[k].f.end(), 0.0);
	}
}

void Solids_velocity_new(std::vector<Circle>& Solids, Param par) {

	for (auto& it : Solids) {
		if (it.moving == 0) {
			it.u[1] = 0.;
			it.u[2] = 0.;
			it.omega[3] = 0.;
		}
		else if (it.moving == 1) {
			it.u_s     = it.u;
			it.omega_s = it.omega;

			GeomVec d_uc    = par.d_t * (it.integralV_du_dt  - it.f_new  ) / it.V / it.rho * 1 + par.Gravity * (it.rho - 1.) / it.rho * par.d_t + it.d_uv_collide + it.d_ur_collide;  // fluid density equals 1
			GeomVec d_omega = par.d_t * (it.integralV_dur_dt - it.tau_new) / it.I / it.rho * 1 + it.d_omega_collide;  // angular moment I is normalized with density

			it.u     = it.u_n     + d_uc;
			it.omega = it.omega_n + d_omega;
		}
		else if (it.moving == 2) {
			it.u[1] = 0.;
			it.u[2] = 0.;
			double t = 100*par.d_t*(par.N_step + 1);
			if (t < 0.25) {
				it.omega[3] = sin(2.*M_PI*t)
					        * sin(2.*M_PI*t) * par.omega_BC;
			}
			else {
				it.omega[3] = par.omega_BC;
			}
		}
	}
}

void Solids_position_new(std::vector<Circle>& Solids, std::vector<Node>& Nodes, Param par) {
	for (auto& it : Solids) {
		if (it.moving>0) {
			it.x = it.x_n + 0.5 * (it.u_n + it.u) * par.d_t;
			for (size_t k = 0; k < it.Nn; ++k) {
				int Ind = it.IndNodes[k];
				Nodes[Ind].x = rotate_Vector_around_vector(Nodes[Ind].x_n, 0.5 * (it.omega_n + it.omega) * par.d_t); //rotate
			}
		}
	}
}

void Circle::integrals(Matrix U_n, Matrix V_n, Matrix U_new, Matrix V_new, Param par) {
	GeomVec integralV_un, integralV_unew, integralV_un_r, integralV_unew_r;
	std::fill(integralV_un.begin()    , integralV_un.end()    , 0.0);
	std::fill(integralV_unew.begin()  , integralV_unew.end()  , 0.0);
	std::fill(integralV_un_r.begin()  , integralV_un_r.end()  , 0.0);
	std::fill(integralV_unew_r.begin(), integralV_unew_r.end(), 0.0);

	int nx1 = U_n.size();
	int nx2 = U_n[0].size();
	int ny1 = V_n.size();
	int ny2 = V_n[0].size();
	int i_max, i_min;
	int j_max, j_min;

	GetInfluenceArea(i_min, i_max, j_min, j_max, nx1 - 1, nx2 - 1, x_n, int(r / par.d_x) + 4, par);
	for (int i = i_min; i <= i_max; ++i) {
		for (int j = j_min; j <= j_max; ++j) {
			int i_real = i_real_u(i, par);
			GeomVec xu = x_u(i, j, par);
			double Frac_n = par.d_x * par.d_y * Volume_Frac(x_n, r, xu, par.d_x, par.d_y);
			integralV_un[1] += Frac_n * U_n[i_real][j];
			GeomVec un;
			un[1] = U_n[i_real][j];
			un[2] = 0.0;
			un[3] = 0.0;
			integralV_un_r   += Frac_n * x_product(xu-x_n, un);
		}
	}

	GetInfluenceArea(i_min, i_max, j_min, j_max, nx1 - 1, nx2 - 1, x, int(r / par.d_x) + 4, par);
	for (int i = i_min; i <= i_max; ++i) {
		for (int j = j_min; j <= j_max; ++j) {
			int i_real = i_real_u(i, par);
			GeomVec xu = x_u(i, j, par);
			double Frac = par.d_x * par.d_y * Volume_Frac(x, r, xu, par.d_x, par.d_y);
			integralV_unew[1] += Frac * U_new[i_real][j];
			GeomVec unew;
			unew[1] = U_new[i_real][j];
			unew[2] = 0.0;
			unew[3] = 0.0;
			integralV_unew_r += Frac * x_product(xu - x, unew);
		}
	}

	GetInfluenceArea(i_min, i_max, j_min, j_max, ny1 - 1, ny2 - 1, x_n, int(r / par.d_y) + 4, par);
	for (int i = i_min; i <= i_max; ++i) {
		for (int j = j_min; j <= j_max; ++j) {
			int i_real = i_real_v(i, par);
			GeomVec xv = x_v(i, j, par);
			double Frac_n = par.d_x * par.d_y * Volume_Frac(x_n, r, xv, par.d_x, par.d_y);
			integralV_un[2] += Frac_n * V_n[i_real][j];
			GeomVec vn;
			vn[1] = 0.0;
			vn[2] = V_n[i_real][j];
			vn[3] = 0.0;
			integralV_un_r   += Frac_n * x_product(xv - x_n, vn);
		}
	}

	GetInfluenceArea(i_min, i_max, j_min, j_max, ny1 - 1, ny2 - 1, x, int(r / par.d_y) + 4, par);
	for (int i = i_min; i <= i_max; ++i) {
		for (int j = j_min; j <= j_max; ++j) {
			int i_real = i_real_v(i, par);
			GeomVec xv = x_v(i, j, par);
			double Frac = par.d_x * par.d_y * Volume_Frac(x, r, xv, par.d_x, par.d_y);
			integralV_unew[2] += Frac * V_new[i_real][j];
			GeomVec vnew;
			vnew[1] = 0.0;
			vnew[2] = V_new[i_real][j];
			vnew[3] = 0.0;
			integralV_unew_r += Frac * x_product(xv - x, vnew);
		}
	}

	//std::cout << integralV_unew_r[3] << "   " << integralV_un_r[3] << std::endl;

	integralV_du_dt  = (integralV_unew   - integralV_un  ) / par.d_t;
	integralV_dur_dt = (integralV_unew_r - integralV_un_r) / par.d_t;

}


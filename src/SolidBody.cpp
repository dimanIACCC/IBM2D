#include "stdafx.h"
#include "SolidBody.h"

Solid::Solid(double x, double y, double ux, double uy, double alpha, double omega, double rho, int Nn_r0, int moving, int name, int shape, double r0, double r, double e)
{
	this->x_n = ZeroVec();
	this->u_n = ZeroVec();
	this->omega_n = ZeroVec();
	this->alpha = ZeroVec();
	this->x_n[1] = x;
	this->x_n[2] = y;
	this->u_n[1] = ux;
	this->u_n[2] = uy;
	this->omega_n[3] = omega;
	this->alpha[3] = alpha;
	this->rho = rho;
	this->Nn_r0 = Nn_r0;
	this->Nn_r = 0;
	this->Nn = 0;
	this->moving = moving;
	this->name = name;
	this->x_n_plt = this->x_n;
	this->x_plt   = this->x_n;
	this->x     = this->x_n;
	this->u     = this->u_n;
	this->omega = this->omega_n;
	this->Fr = 0.;
	this->shape = shape;
	this->r0 = r0;
	this->r = r;
	this->e = e;

	if (shape == 0) {
		this->Nn_r = 0;
		this->Nn = 4 * this->Nn_r0;
	}
	else if (shape == 2){
		this->Nn_r = int(this->Nn_r0 * 0.5 * r / r0);
		this->Nn = 2 * this->Nn_r + 4 * this->Nn_r0;
	}
	else if (shape == 4) {
		this->Nn_r = int(this->Nn_r0 * 0.5 * r / r0);
		this->Nn = 8 * this->Nn_r + 8 * this->Nn_r0;
	}
	else if (shape == 8) {
		this->Nn_r = int(this->Nn_r0 * 0.5 * r / r0);
		this->Nn = 16 * this->Nn_r + 12 * this->Nn_r0;
	}
	else if (shape == 16) {
		this->Nn_r = int(this->Nn_r0 * 0.5 * r / r0);
		this->Nn = 8 * this->Nn_r + 16 * this->Nn_r0;
	}
	V = M_PI * r * r;
	I = V * r * r / 4.0 * (2.0 - e*e) / sqrt(1.0 - e*e); // angular momentum for unit density
}

Solid::Solid(double x, double y, Param &par) :
	Solid(x, y, 0.0, 0.0, 0.0, 0.0, par.rho, par.Nn_, true, par.SolidName_max + 1, par.shape, par.r0, par.r, par.e) {
	this->u_n[1] = 0.0;  //  par.u_in * ux_Poiseuille(y, par.H);
	this->omega_n[3] = 0.0; // -dux_dy_Poiseuille(y, par.H);
	this->omega = this->omega_n;
	this->u = this->u_n;
}


Solid::~Solid()
{

}

void Solid::add_Nodes(std::vector<Node> &Nodes, const int Nn_max) {
	this->IndNodes.resize(this->Nn);
	Nodes.resize(Nn_max + this->Nn);
	for (size_t i = 0; i < Nn; ++i) {
		int Ind = Nn_max + i;
		this->IndNodes[i] = Ind;
	}
}

bool operator <(const Solid& a, const Solid& b) {
	if (length(a.x_n) < length(b.x_n)) return true;
	else return false;
}

bool operator >(const Solid& a, const Solid& b) {
	if (length(a.x_n) > length(b.x_n)) return true;
	else return false;
}


 

void fill_solid_coordinates(std::vector<Node> &Nodes, const int Nn_max, const int Nn, const int Nn_r0, const int Nn_r,
                            const int shape, const double r0, const double r, const double e, const double alpha, const double dxy) {
	int Nn_ = 0;
	GeomVec Xbeg = ZeroVec();
	GeomVec Xend = ZeroVec();
	GeomVec Xc = ZeroVec();
	double alpha_beg, alpha_end;

	if (shape == 0) {    // ellipse
		for (size_t i = 0; i < Nn; ++i) {
			int Ind = Nn_max + i;
			double phi = i * 2.0 * M_PI / Nn;
			Nodes[Ind].x_n[1] = cos(phi) * r / sqrt(1 - e*e*cos(phi)*cos(phi)) * pow(1 - e*e, 0.25);
			Nodes[Ind].x_n[2] = sin(phi) * r / sqrt(1 - e*e*cos(phi)*cos(phi)) * pow(1 - e*e, 0.25);
		}
	}
	else if (shape == 2) {  // line
		// down line
		Xbeg[1] = -r;  Xbeg[2] = -r0;
		Xend[1] =  r;  Xend[2] = -r0;
		line_segment(Nodes, Nn_max + Nn_, Nn_r, Xbeg, Xend);
		Nn_ += Nn_r;

		// right half-circle
		Xc[1] = r;  Xc[2] = 0;
		alpha_beg = -M_PI / 2;   alpha_end = M_PI / 2;
		circular_segment(Nodes, Nn_max + Nn_, 2 * Nn_r0, Xc, r0, alpha_beg, alpha_end);
		Nn_ += 2 * Nn_r0;

		// up line
		Xbeg[1] =  r;  Xbeg[2] = r0;
		Xend[1] = -r;  Xend[2] = r0;
		line_segment(Nodes, Nn_max + Nn_, Nn_r, Xbeg, Xend);
		Nn_ += Nn_r;

		// left half-circle
		Xc[1] = -r;  Xc[2] = 0;
		alpha_beg = M_PI / 2;   alpha_end = 3 * M_PI / 2;
		circular_segment(Nodes, Nn_max + Nn_, 2 * Nn_r0, Xc, r0, alpha_beg, alpha_end);
		Nn_ += 2 * Nn_r0;
	}
	else if (shape == 4) {  // cross
		// R down line
		Xbeg[1] =  r0;  Xbeg[2] = -r0;
		Xend[1] =  r ;  Xend[2] = -r0;
		line_segment(Nodes, Nn_max + Nn_, Nn_r, Xbeg, Xend);
		Nn_ += Nn_r;

		// R  half-circle
		Xc[1] = r;  Xc[2] = 0;
		alpha_beg = -M_PI / 2;   alpha_end =  M_PI / 2;
		circular_segment(Nodes, Nn_max + Nn_, 2 * Nn_r0, Xc, r0, alpha_beg, alpha_end);
		Nn_ += 2 * Nn_r0;

		// R up line
		Xbeg[1] = r ;  Xbeg[2] = r0;
		Xend[1] = r0;  Xend[2] = r0;
		line_segment(Nodes, Nn_max + Nn_, Nn_r, Xbeg, Xend);
		Nn_ += Nn_r;

		// U
		copy_solid_mesh(Nodes, Nn_max + 0 * Nn_, Nn_max + 1 * Nn_, Nn_);

		// L
		copy_solid_mesh(Nodes, Nn_max + 1 * Nn_, Nn_max + 2 * Nn_, Nn_);

		// D
		copy_solid_mesh(Nodes, Nn_max + 2 * Nn_, Nn_max + 3 * Nn_, Nn_);

	}
	else if (shape == 8) {  // svastic
		// R
		Xbeg[1] = r0;  Xbeg[2] = -r0;
		Xend[1] = r ;  Xend[2] = -r0;
		line_segment(Nodes, Nn_max + Nn_, Nn_r, Xbeg, Xend);
		Nn_ += Nn_r;

		Xc[1] = r;  Xc[2] = 0;
		alpha_beg = -M_PI/2;   alpha_end = 0;
		circular_segment(Nodes, Nn_max + Nn_, Nn_r0, Xc, r0, alpha_beg, alpha_end);
		Nn_ += Nn_r0;

		Xbeg[1] = r + r0;  Xbeg[2] = 0;
		Xend[1] = r + r0;  Xend[2] = r;
		line_segment(Nodes, Nn_max + Nn_, Nn_r, Xbeg, Xend);
		Nn_ += Nn_r;

		Xc[1] = r;  Xc[2] = r;
		alpha_beg = 0;   alpha_end = M_PI;
		circular_segment(Nodes, Nn_max + Nn_, 2 * Nn_r0, Xc, r0, alpha_beg, alpha_end);
		Nn_ += 2 * Nn_r0;

		Xbeg[1] = r - r0;  Xbeg[2] = r;
		Xend[1] = r - r0;  Xend[2] = r0;
		line_segment(Nodes, Nn_max + Nn_, Nn_r, Xbeg, Xend);
		Nn_ += Nn_r;

		Xbeg[1] = r - r0;  Xbeg[2] = r0;
		Xend[1] =     r0;  Xend[2] = r0;
		line_segment(Nodes, Nn_max + Nn_, Nn_r, Xbeg, Xend);
		Nn_ += Nn_r;

		// U
		copy_solid_mesh(Nodes, Nn_max + 0 * Nn_, Nn_max + 1 * Nn_, Nn_);

		// L
		copy_solid_mesh(Nodes, Nn_max + 1 * Nn_, Nn_max + 2 * Nn_, Nn_);

		// D
		copy_solid_mesh(Nodes, Nn_max + 2 * Nn_, Nn_max + 3 * Nn_, Nn_);

	}
	else if (shape == 16) {  // propeller
		int Nn_ = 0;
		GeomVec Xbeg = ZeroVec();
		GeomVec Xend = ZeroVec();
		GeomVec Xc   = ZeroVec();
		double alpha_beg, alpha_end;

		// R
		Xbeg[1] = r + r0;  Xbeg[2] = 0;
		Xend[1] = r + r0;  Xend[2] = r;
		line_segment(Nodes, Nn_max + Nn_, Nn_r, Xbeg, Xend);
		Nn_ += Nn_r;

		Xc[1] = r;  Xc[2] = r;
		alpha_beg = 0;   alpha_end = M_PI;
		circular_segment(Nodes, Nn_max + Nn_, 2*Nn_r0, Xc, r0, alpha_beg, alpha_end);
		Nn_ += 2 * Nn_r0;

		Xbeg[1] = r - r0;  Xbeg[2] = r;
		Xend[1] = r - r0;  Xend[2] = 0;
		line_segment(Nodes, Nn_max + Nn_, Nn_r, Xbeg, Xend);
		Nn_ += Nn_r;

		Xc[1] = r;  Xc[2] = 0;
		alpha_beg = M_PI;   alpha_end = 2*M_PI;
		circular_segment(Nodes, Nn_max + Nn_, 2*Nn_r0, Xc, r0, alpha_beg, alpha_end);
		Nn_ += 2 * Nn_r0;

		// U
		copy_solid_mesh(Nodes, Nn_max + 0 * Nn_, Nn_max + 1 * Nn_, Nn_);

		// L
		copy_solid_mesh(Nodes, Nn_max + 1 * Nn_, Nn_max + 2 * Nn_, Nn_);

		// D
		copy_solid_mesh(Nodes, Nn_max + 2 * Nn_, Nn_max + 3 * Nn_, Nn_);

	}
	else {
		std::cout << "fill_solid_coordinates: unknown solid shape" << std::endl;
	}

	GeomVec o = ZeroVec();
	o[3] = alpha;
	for (size_t i = 0; i < Nn; ++i) {
		int Ind = Nn_max + i;
		Nodes[Ind].x_n = rotate_Vector_around_vector(Nodes[Ind].x_n, o);
		Nodes[Ind].x = Nodes[Ind].x_n;
	}

	fill_solid_ds(Nodes, Nn_max, Nn, shape, dxy);
}

void line_segment(std::vector<Node> &Nodes, const int N_start, const int Nn, GeomVec Xbeg, GeomVec Xend) {
	for (size_t i = 0; i < Nn; ++i) {
		Nodes[N_start+i].x_n = Xbeg + (Xend - Xbeg) * i / Nn;
	}
}

void circular_segment(std::vector<Node> &Nodes, const int N_start, const int Nn, GeomVec Xc, double e, double alpha_beg, double alpha_end) {
	double alpha;
	for (size_t i = 0; i < Nn; ++i) {
		alpha = alpha_beg + (alpha_end - alpha_beg) * i / Nn;
		Nodes[N_start + i].x_n[1] = Xc[1] + e*cos(alpha);
		Nodes[N_start + i].x_n[2] = Xc[2] + e*sin(alpha);
	}
}

void copy_solid_mesh(std::vector<Node> &Nodes, const int N_beg_from, const int N_beg_to, const int Nn) {
	GeomVec ooo = ZeroVec();
	ooo[3] = M_PI/2;
	for (size_t i = 0; i < Nn; ++i) {
		Nodes[N_beg_to + i] = Nodes[N_beg_from + i];
		Nodes[N_beg_to + i].x_n = rotate_Vector_around_vector(Nodes[N_beg_to + i].x_n, ooo);
	}
}


void fill_solid_ds(std::vector<Node> &Nodes, const int Nn_max, const int Nn, const int shape, const double dxy) {

	for (size_t i = 1; i < Nn - 1; ++i) {
		int Ind = Nn_max + i;
		Nodes[Ind].ds = 0.5 * length(Nodes[Ind + 1].x_n - Nodes[Ind - 1].x_n) * dxy;
	}
	if (shape == 0 || shape == 2 || shape == 4 || shape == 8) {    // ellipse or cross
		Nodes[Nn_max + 0     ].ds = 0.5 * length(Nodes[Nn_max + 1].x_n - Nodes[Nn_max + Nn - 1].x_n) * dxy;
		Nodes[Nn_max + Nn - 1].ds = 0.5 * length(Nodes[Nn_max + 0].x_n - Nodes[Nn_max + Nn - 2].x_n) * dxy;
	}
	else if (shape == 16) {  // propeller
		int N04 = Nn_max;
		int N14 = Nn_max + 1 * Nn / 4;
		int N24 = Nn_max + 2 * Nn / 4;
		int N34 = Nn_max + 3 * Nn / 4;
		int N44 = Nn_max + 4 * Nn / 4;
		Nodes[N04    ].ds = 0.5 * length(Nodes[N04 + 1].x_n - Nodes[N14 - 1].x_n) * dxy;
		Nodes[N14 - 1].ds = 0.5 * length(Nodes[N04 + 0].x_n - Nodes[N14 - 2].x_n) * dxy;
		Nodes[N14    ].ds = 0.5 * length(Nodes[N14 + 1].x_n - Nodes[N24 - 1].x_n) * dxy;
		Nodes[N24 - 1].ds = 0.5 * length(Nodes[N14 + 0].x_n - Nodes[N24 - 2].x_n) * dxy;
		Nodes[N24    ].ds = 0.5 * length(Nodes[N24 + 1].x_n - Nodes[N34 - 1].x_n) * dxy;
		Nodes[N34 - 1].ds = 0.5 * length(Nodes[N24 + 0].x_n - Nodes[N34 - 2].x_n) * dxy;
		Nodes[N34    ].ds = 0.5 * length(Nodes[N34 + 1].x_n - Nodes[N44 - 1].x_n) * dxy;
		Nodes[N44 - 1].ds = 0.5 * length(Nodes[N34 + 0].x_n - Nodes[N44 - 2].x_n) * dxy;
	}

}

void velocities(std::vector<Solid>::iterator &Solid, std::vector<Node> &Nodes) {
	for (size_t k = 0; k < Solid->Nn; ++k) {
		int Ind = Solid->IndNodes[k];
		Nodes[Ind].us = 0.5 * (Solid->u + Solid->u + x_product(Solid->omega+ Solid->omega, Nodes[Ind].x));
	}
}

void coordinates(std::vector<Solid>::iterator &Solid, std::vector<Node> &Nodes) {
	for (size_t k = 0; k < Solid->Nn; ++k) {
		int Ind = Solid->IndNodes[k];
		Nodes[Ind].x_s = Solid->x + Nodes[Ind].x;
	}
}

void Solid::log_init(std::string WorkDir) {
	std::ofstream output;
	std::string filename = WorkDir + "Solids/" + std::to_string(name) + ".plt";
	output.open(filename);

	output << "title = " << '"' << "Solid" << name << '"' << std::endl;
	output << "Variables = n x y u v fx fy omega tau alpha" << std::endl;
}

void Solid::log(std::string WorkDir, int n) {
	std::ofstream output;
	std::string filename = WorkDir + "Solids/" + std::to_string(name) + ".plt";
	output.open(filename, std::ios::app);

	output << std::setprecision(8);
	output << n << "   " 
		   << x_n_plt[1] << "   "
	       << x_n_plt[2] << "   "
	       << u_n[1] << "   "
	       << u_n[2] << "   "
	       << f_new[1] << "   "
	       << f_new[2] << "   "
	       << omega_n[3] << "   "
	       << tau_new[3]   << "   "
		   << alpha[3] << "   "
	       << std::endl;
}

void Read_Solids(std::string filename, std::vector<Solid>& Solids, std::vector<Node> &Nodes, Param &par) {
	std::ifstream input;
	std::string line;

	input.open(filename.c_str());
	if (input.is_open()) {
		while (getline(input, line)) { // read line from file to string $line$
			if (line == "circle{" || line == "line{" || line == "cross{") {
				double x = par.L*0.1;
				double y = par.H*0.5;
				double ux = 0;
				double uy = 0;
				double omega = 0;
				double rho = par.rho;
				int Nn_ = par.Nn_;
				int moving = 1;
				int shape = 0;
				double r0 = par.r0;
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
						else if (PAR == "Nn_")        Nn_          = stoi(VALUE);
						else if (PAR == "moving")     moving       = stoi(VALUE);
						else if (PAR == "Poiseuille") Poiseuille   = bool(stoi(VALUE));
						else if (PAR == "shape")      shape        = stoi(VALUE);
						else if (PAR == "r0")         r0           = stod(VALUE);
						else if (PAR == "r")          r            = stod(VALUE);
						else if (PAR == "e")          e            = stod(VALUE);
						else if (PAR == "alpha")      alpha        = stod(VALUE) * M_PI / 180;
						else    std::cout << "Read_Solids: unknown parameter " << PAR << std::endl;
					}
					else {
						std::cout << "Read_Solids: no value inputed" << std::endl;
					}
				}
				if (Poiseuille) {
					ux = par.u_in * ux_Poiseuille(y, par.H);
					uy = 0;
					omega = - dux_dy_Poiseuille(y, par.H);
				}

				Solid c(x, y, ux, uy, alpha, omega, rho, Nn_, moving, par.SolidName_max+1, shape, r0, r, e);
				//std::cout << c.I << std::endl;
				c.add_Nodes(Nodes, par.Nn_max);
				fill_solid_coordinates(Nodes, par.Nn_max, c.Nn, c.Nn_r0, c.Nn_r, c.shape, c.r0, c.r, c.e, alpha, 0.5*(par.d_x+par.d_y));
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

void Add_Solids(std::vector<Solid>& Solids, std::vector<Node>& Nodes, Param &par) {
	if (   (par.N_step - par.AddSolids_start) % par.AddSolids_interval == 0   &&   (par.N_step >= par.AddSolids_start)) { //create new solids starting from $AddSolids_start$ iteration with interval of $AddSolids_interval$ iterations
		for (int i = 0; i < par.AddSolids_N; i++) { // add $AddSolids_N$ solids
			GeomVec x = ZeroVec();
			x[1] = (par.L)  * (0.5 + 1*(1 - 3 * par.r / par.L) * (double(rand()) - RAND_MAX / 2) / RAND_MAX);
			x[2] = (par.H)  * (0.5 + 1*(1 - 3 * par.r / par.H) * (double(rand()) - RAND_MAX / 2) / RAND_MAX);

			// check if new Solid does not cross other Solids
			bool add = true;
			for (auto solid = Solids.begin(); solid != Solids.end(); solid++) {
				if (length(x - solid->x) < 1.5 * (par.r + solid->r)) {
				//if (length(x - solid->x) < 1.25 * (par.r + solid->r) && solid->moving == 1) { // flow inside Solid
					add = false;
					i--;
					break;
				}
			}
			if (add) {
				Solid c(x[1], x[2], par);
				c.add_Nodes(Nodes, par.Nn_max);
				fill_solid_coordinates(Nodes, par.Nn_max, c.Nn, c.Nn_r0, c.Nn_r, c.shape, par.r0, par.r, par.e, 0., 0.5*(par.d_x + par.d_y));
				par.Nn_max += c.Nn;
				if (c.name > par.SolidName_max) par.SolidName_max = c.name;
				Solids.push_back(c);
				c.log_init(par.WorkDir);
			}
		}
	}

}

double Distance_2Solids(Solid& s1, Solid& s2, Param& par, GeomVec& r) {
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
	dist = rrr - (s1.r + s2.r);
	r = r / rrr;
	return dist;
}

double Distance_2Solids_(Solid& s1, Solid& s2, std::vector<Node> Nodes, Param& par, GeomVec& r_out, GeomVec& x1, GeomVec& x2, GeomVec& u1, GeomVec& u2) {
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
				if (length(r_plus) < rrr ) { r = r_plus ; rrr = length(r_plus) ; };
				if (length(r_minus) < rrr) { r = r_minus; rrr = length(r_minus); };
			}
			if (rrr < dist) {
				dist = rrr;
				x1 = Nodes[Ind1].x_s;
				x2 = Nodes[Ind2].x_s;
				u1 = s1.u_n + x_product(s1.omega_n, x1);
				u2 = s2.u_n + x_product(s2.omega_n, x2);
				r_out = r / rrr;
			}
		}
	}
	//std::cout << x1 << std::endl;
	//std::cout << x2 << std::endl;

	return dist;
}

void Distance_2Walls(Solid& s1, std::vector<Node> Nodes, Param& par,
	double& d_up, double& d_down, double& d_left, double& d_right,
	GeomVec x_up, GeomVec x_down, GeomVec x_left, GeomVec x_right ) {
	d_up = +1.e99;
	d_down = +1.e99;
	d_left = +1.e99;
	d_right = +1.e99;

	for (size_t k1 = 0; k1 < s1.Nn; ++k1) {
		int Ind1 = s1.IndNodes[k1];
			double d_up_    = 2 * abs(Nodes[Ind1].x_s[2] - par.H);
			double d_down_  = 2 * abs(Nodes[Ind1].x_s[2]);
			double d_left_  = 2 * abs(Nodes[Ind1].x_s[1]);
			double d_right_ = 2 * abs(Nodes[Ind1].x_s[1] - par.L);

			if (d_up_    < d_up   ) { d_up    = d_up_   ; x_up    = Nodes[Ind1].x;}
			if (d_down_  < d_down ) { d_down  = d_down_ ; x_down  = Nodes[Ind1].x;}
			if (d_left_  < d_left ) { d_left  = d_left_ ; x_left  = Nodes[Ind1].x;}
			if (d_right_ < d_right) { d_right = d_right_; x_right = Nodes[Ind1].x;}
	}
}

GeomVec F_collide(GeomVec norm, double dist, GeomVec u1, GeomVec u2, double dist_u, double dist_r, double alpha, double beta, double friction) {
	GeomVec F_collide = ZeroVec();
	if (dist <= dist_u) {
		GeomVec d_u12 = u1 - u2;
		GeomVec d_u12_norm = dot_product(d_u12, norm)*norm;
		double d_u12_proj = dot_product(d_u12, norm);
		GeomVec d_u12_tau = d_u12 - d_u12_norm;
		if (d_u12_proj < 0.0) { // if Solids move to each other
			F_collide -=            alpha * 2 * (dist_u - dist) * (dist_u - dist) / dist_u / dist_u * d_u12_norm;
			F_collide -= friction * alpha * 2 * (dist_u - dist) * (dist_u - dist) / dist_u / dist_u * d_u12_tau;
			std::cout << "Collision u:  F_collide = " << F_collide[1] << "   " << F_collide[2] << std::endl;
		}
		if (dist <= dist_r) {
			F_collide += beta * (dist_r - dist)*(dist_r - dist) / dist_r / dist_r * norm;
			std::cout << "Collision r:  F_collide = " << F_collide[1] << "   " << F_collide[2] << std::endl;
			std::getchar();
		}
	}
	return F_collide;
}

void Collide(Solid& s1, Solid& s2, std::vector<Node> &Nodes, Param par, double dist_u, double dist_r, double alpha, double beta, double friction) {

	GeomVec r = s1.x - s2.x;
	GeomVec x1 = s1.x;
	GeomVec x2 = s2.x;
	GeomVec u1 = s1.u_n;
	GeomVec u2 = s2.u_n;
	double dist = Distance_2Solids(s1, s2, par, r);
	if (dist < s1.r + s2.r) dist = Distance_2Solids_(s1, s2, Nodes, par, r, x1, x2, u1, u2);

	double m1 = s1.rho * s1.V;
	double m2 = s2.rho * s2.V;
	double M = m1 * m2 / (m1 + m2);
	double I1 = s1.rho * s1.I;
	double I2 = s2.rho * s2.I;
	double I = I1 * I2 / (I1 + I2);
	GeomVec F = F_collide(r, dist, u1, u2, dist_u, dist_r, alpha, beta, friction);

	if (s1.moving == 1) s1.a_collide += F * M / m1;
	if (s2.moving == 1) s2.a_collide -= F * M / m2;
	if (s1.moving == 1) s1.d_omega_collide += x_product(x1, F) * I / I1;
	if (s2.moving == 1) s2.d_omega_collide -= x_product(x2, F) * I / I2;

}

void Solids_collide(std::vector<Solid> &solidList, std::vector<Node> &Nodes, Param par) {

	double kr = 0.25;    // fraction of distance where distance force switches on
	double dist_u = par.k_dist*par.d_x;
	double dist_r = par.k_dist*par.d_x * kr;

	double alpha = 20. * sqrt(par.Gravity_module / dist_u); // coefficient for the collision force based on velocity value
	double beta  = 100. * std::fmax(par.Gravity_module, 1000 * std::fmax(abs(par.u_up - par.u_down), abs(par.u_in)) / par.H / par.H / par.Re); // coefficient for the collision force based on distance between particles value
	double friction = 0.04; // coefficient for the friction force based on velocity value

	double  d_up, d_down, d_left, d_right;
	GeomVec x_up = ZeroVec(), x_down = ZeroVec(), x_left = ZeroVec(), x_right = ZeroVec();
	GeomVec u_up = ZeroVec(), u_down = ZeroVec(), u_left = ZeroVec(), u_right = ZeroVec();
	GeomVec n_up = ZeroVec(), n_down = ZeroVec(), n_left = ZeroVec(), n_right = ZeroVec();
	
	u_up[1]   = par.u_up;
	u_down[1] = par.u_down;

	n_left[1]  =  1.;
	n_right[1] = -1.;
	n_up[2]    = -1.;
	n_down[2]   = 1.;

	for (auto one = solidList.begin(); one != solidList.end(); one++) {
		if (one->moving == 1){
			for (auto two = next(one); two != solidList.end(); two++) {
				Collide(*one, *two, Nodes, par, dist_u, dist_r, alpha, beta, friction);
			}

			Distance_2Walls(*one, Nodes, par, d_up, d_down, d_left, d_right, x_up, x_down, x_left, x_right);

			if (one->moving == 1) {
				GeomVec F_up    = F_collide(n_up   , d_up   , one->u_n, u_up   , dist_u, dist_r, alpha, beta, friction);
				GeomVec F_down  = F_collide(n_down , d_down , one->u_n, u_down , dist_u, dist_r, alpha, beta, friction);
				one->a_collide += F_up + F_down;
				one->d_omega_collide += x_product(x_up, F_up) + x_product(x_down, F_down);
			if (par.BC == box) {
				GeomVec F_left  = F_collide(n_left , d_left , one->u_n, u_left , dist_u, dist_r, alpha, beta, friction);
				GeomVec F_right = F_collide(n_right, d_right, one->u_n, u_right, dist_u, dist_r, alpha, beta, friction);
				one->a_collide += F_left + F_right;
				one->d_omega_collide += x_product(x_left, F_left) + x_product(x_right, F_right);
			}
			}
		}
	}

}



void Solids_move(std::vector<Solid> &solidList, std::vector<Node> &Nodes, Param par) {
	for (auto it = solidList.begin(); it != solidList.end();) {
			if (it->moving > 0) {
				it->u_n    = it->u;
				it->omega_n = it->omega;
				it->x_n    = it->x;
				it->x_n_plt = it->x_plt;
				for (size_t k = 0; k < it->Nn; ++k) {
					int Ind = it->IndNodes[k];
					Nodes[Ind].x_n = Nodes[Ind].x;
				}
			}
			it->log(par.WorkDir, par.N_step);

			//Right boundary conditions for Solids
			if (it->x_n[1] < 0) {
				if (par.BC == periodical) {
					it->x_n[1] += par.L;
					it++;
				}
				else
					it = solidList.erase(it);
			}
			else if(it->x_n[1] < par.L) {
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


void h_average_of_Solids_Layer(std::vector<Solid> &solidList, Param par, double& h_average) {
	//solidList.sort(std::greater<Solid>());
	sort(solidList.begin(), solidList.end(), std::greater<Solid>());
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


void Solids_zero_force(std::vector<Solid>& Solids, std::vector<Node>& Nodes, int Nn_max) {
	for (auto& it : Solids) {
		it.f = ZeroVec();
		it.tau = ZeroVec();
		it.a_collide = ZeroVec();
		it.d_omega_collide = ZeroVec();
		
	}
	for (size_t k = 0; k < Nn_max; ++k) {
		Nodes[k].f = ZeroVec();
	}
}

void Solids_velocity_new(std::vector<Solid>& Solids, Param par) {

	for (auto& it : Solids) {
		if (it.moving == 0) {
			it.u = ZeroVec();
			it.omega[3] = 0.;
		}
		else if (it.moving == 1) {
			it.u_s     = it.u;
			it.omega_s = it.omega;

			GeomVec a   = (it.integralV_du_dt  - it.f_new  ) / it.V / it.rho * 1 + par.Gravity * (it.rho - 1.) / it.rho + it.a_collide;  // fluid density equals 1
			GeomVec tau = (it.integralV_dur_dt - it.tau_new) / it.I / it.rho * 1 + it.d_omega_collide;  // angular moment I is normalized with density

			it.u     = it.u_n     + a   * par.d_t;
			it.omega = it.omega_n + tau * par.d_t;
		}
		else if (it.moving == 2) {
			it.u = ZeroVec();
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

void Solids_position_new(std::vector<Solid>& Solids, std::vector<Node>& Nodes, Param par) {
	for (auto& it : Solids) {
		if (it.moving>0) {
			it.x     = it.x_n     + 0.5 * (it.u_n + it.u) * par.d_t;
			it.x_plt = it.x_n_plt + 0.5 * (it.u_n + it.u) * par.d_t;
			it.alpha = it.alpha + 0.5 * (it.omega_n + it.omega) * par.d_t;
			for (size_t k = 0; k < it.Nn; ++k) {
				int Ind = it.IndNodes[k];
				Nodes[Ind].x = rotate_Vector_around_vector(Nodes[Ind].x_n, 0.5 * (it.omega_n + it.omega) * par.d_t); //rotate
			}
		}
	}
}

void Solid::integrals(Matrix U_n, Matrix V_n, Matrix U_new, Matrix V_new, Param par) {
	GeomVec integralV_un = ZeroVec();
	GeomVec integralV_unew = ZeroVec();
	GeomVec integralV_un_r = ZeroVec();
	GeomVec integralV_unew_r = ZeroVec();

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
			GeomVec un = ZeroVec();
			un[1] = U_n[i_real][j];
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
			GeomVec unew = ZeroVec();
			unew[1] = U_new[i_real][j];
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
			GeomVec vn = ZeroVec();
			vn[2] = V_n[i_real][j];
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
			GeomVec vnew = ZeroVec();
			vnew[2] = V_new[i_real][j];
			integralV_unew_r += Frac * x_product(xv - x, vnew);
		}
	}

	//std::cout << integralV_unew_r[3] << "   " << integralV_un_r[3] << std::endl;

	integralV_du_dt  = (integralV_unew   - integralV_un  ) / par.d_t;
	integralV_dur_dt = (integralV_unew_r - integralV_un_r) / par.d_t;

}


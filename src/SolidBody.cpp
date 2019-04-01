#include "stdafx.h"
#include "SolidBody.h"


SolidBody::SolidBody(double x, double y, double ux, double uy, double omega, double rho, int Nn, bool moving, int name)
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
	Nodes.resize(Nn);
	this->moving = moving;
	this->name = name;
	this->x     = this->x_n;
	this->u     = this->u_n;
	this->omega = this->omega_n;
	this->Fr = 0.;
	this->n_moving = 0;
}


SolidBody::~SolidBody()
{

}

GeomVec Circle_Equation(GeomVec xc, double r, double theta) {
	GeomVec xy;
	xy[1] = xc[1] + cos(theta) * r;
	xy[2] = xc[2] + sin(theta) * r;
	xy[3] = 0;
	return xy;
}

Circle::Circle(double x, double y, double ux, double uy, double omega, double rho, int Nn, bool moving, int name, double r) :
     SolidBody(       x,        y,        ux,        uy,        omega,        rho,     Nn,      moving,     name) {
	this->r = r;
	int N_Solid_Layers = 1;
	int Nn2 = Nn / N_Solid_Layers;
	for (int i = 0; i < Nn2; ++i){
		Nodes[i].x_n[1] = cos(i * 2.0 * M_PI / Nn2) * r;
		Nodes[i].x_n[2] = sin(i * 2.0 * M_PI / Nn2) * r;
		Nodes[i].n[1]  = cos(i * 2.0 * M_PI / Nn2);
		Nodes[i].n[2]  = sin(i * 2.0 * M_PI / Nn2);
		Nodes[i].x = Nodes[i].x_n;
	}
	/*for (int i = Nn2; i < 2*Nn2; ++i) {
		Nodes[i].x_n[1] = cos(i * 2.0 * M_PI / Nn2) * r * 0.95;
		Nodes[i].x_n[2] = sin(i * 2.0 * M_PI / Nn2) * r * 0.95;
		Nodes[i].n[1] = cos(i * 2.0 * M_PI / Nn2);
		Nodes[i].n[2] = sin(i * 2.0 * M_PI / Nn2);
		Nodes[i].x = Nodes[i].x_n;
	}*/
	/*for (int i = 2*Nn2; i < 3*Nn2; ++i) {
		Nodes[i].x_n[1] = cos(i * 2.0 * M_PI / Nn2) * r * 0.5;
		Nodes[i].x_n[2] = sin(i * 2.0 * M_PI / Nn2) * r * 0.5;
		Nodes[i].n[1] = cos(i * 2.0 * M_PI / Nn2);
		Nodes[i].n[2] = sin(i * 2.0 * M_PI / Nn2);
		Nodes[i].x = Nodes[i].x_n;
	}

	for (int i = 3*Nn2; i < Nn; ++i) {
		Nodes[i].x_n[1] = cos(i * 2.0 * M_PI / Nn2) * r * 0.25;
		Nodes[i].x_n[2] = sin(i * 2.0 * M_PI / Nn2) * r * 0.25;
		Nodes[i].n[1]  = cos(i * 2.0 * M_PI / Nn2);
		Nodes[i].n[2]  = sin(i * 2.0 * M_PI / Nn2);
		Nodes[i].x = Nodes[i].x_n;
	}*/
	V = M_PI * r * r;
	I =  V * r * r / 2.0; // angular momentum for unit density
}

Circle::Circle(double x, double y, Param &par):
	    Circle(       x,        y, 0.0, 0.0, 0.0, par.rho, par.Nn, true, par.SolidName_max+1, par.r) {
	this->u_n[1]    =  ux_Poiseuille(y, par.H);
	this->omega_n[3] = -dux_dy_Poiseuille(y, par.H);
	this->omega     = this->omega_n;
	this->u    = this->u_n;
}
 
Circle::~Circle()
{
}

void SolidBody::velocities() {
	for (int i = 0; i < Nn; ++i) {
		Nodes[i].us = u_n + x_product(omega_n, Nodes[i].x_n);
	}
}


double SolidBody::ds(size_t i) {
	GeomVec xL, xR;
	if (i > 0)    xL = Nodes[i-1].x;
	else          xL = Nodes[Nn-1].x;
	if (i < Nn-1) xR = Nodes[i+1].x;
	else          xR = Nodes[0].x;
	double result = length(xL - xR) / 2;
	return result;
}

void SolidBody::log_init(std::string WorkDir) {
	std::ofstream output;
	std::string filename = WorkDir + "Solids/" + std::to_string(name) + ".plt";
	output.open(filename);

	output << "title = " << '"' << "Solid" << name << '"' << std::endl;
	output << "Variables = x y u v fx fy omega tau n Fx_hd Fy_hd Fx_L Fy_L IntTau tau_hd" << std::endl;
}

void SolidBody::log(std::string WorkDir, int n) {
	std::ofstream output;
	std::string filename = WorkDir + "Solids/" + std::to_string(name) + ".plt";
	output.open(filename, std::ios::app);

	output << std::setprecision(8);
	output << x_n[1] << "   "
	       << x_n[2] << "   "
	       << u_n[1] << "   "
	       << u_n[2] << "   "
	       << f_new[1] << "   "
	       << f_new[2] << "   "
	       << omega_n[3] << "   "
	       << tau[3]   << "   "
	       <<  n        << "   "
	       << -F_hd[1] << "   "
	       << -F_hd[2] << "   "
	       << f_L[1] << "   "
	       << f_L[2] << "   "
	       << integralV_dur_dt[3] << "   "
	       << tau_hd[3] << "   "
	       << std::endl;
}

void Read_Solids(std::string filename, std::list<Circle>& Solids, Param &par) {
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
				bool moving = true;
				double r = par.r;
				int n_moving = 0;
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
						else if (PAR == "moving")     moving       = bool(stoi(VALUE));
						else if (PAR == "Poiseuille") Poiseuille   = bool(stoi(VALUE));
						else if (PAR == "r")          r            = stod(VALUE);
						else if (PAR == "n_moving")	  n_moving	   = stoi(VALUE);
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

void Add_Solids(std::list<Circle>& Solids, Param &par) {
	if (   (par.N_step - par.AddSolids_start) % par.AddSolids_interval == 0   &&   (par.N_step >= par.AddSolids_start)) { //create new solids starting from $AddSolids_start$ iteration with interval of $AddSolids_interval$ iterations
		for (int i = 0; i < par.AddSolids_N; i++) { // add $AddSolids_N$ solids
			GeomVec x;
			x[0] = 0;
			x[1] = (par.L)  * (0.5 + (1 - 4 * par.r / par.L) * (double(rand()) - RAND_MAX / 2) / RAND_MAX);
			x[2] = (par.H)  * (0.5 + (1 - 3 * par.r / par.H) * (double(rand()) - RAND_MAX / 2) / RAND_MAX);
			x[3] = 0;
			Circle c(x[1], x[2], par);
			if (c.name > par.SolidName_max) par.SolidName_max = c.name;
			// check if new Solid does not cross other Solids
			bool add = true;
			for (auto solid = Solids.begin(); solid != Solids.end(); solid++) {
				if (length(x - solid->x) < par.k_dist * (par.r + solid->r)) {
					add = false;
					break;
				}
			}
			if (add) {
				Solids.push_back(c);
				c.log_init(par.WorkDir);
			}
		}
	}

}

bool Collide(Circle& s1, Circle& s2, Param par) {
	bool result = false;
	GeomVec r = s1.x - s2.x;
	double distance = length(r); //<----distance between two particles
	if (par.BC == periodical) {
		GeomVec r_plus  = r;
		GeomVec r_minus = r;
		r_plus[1]  += par.L;
		r_minus[1] -= par.L;
		if (length(r_plus ) < distance) { r = r_plus ; distance = length(r_plus ); };
		if (length(r_minus) < distance) { r = r_minus; distance = length(r_minus); };
	}
	if (distance <= par.k_dist*(s1.r + s2.r)) {
		result = true;
		if (par.InelasticCollision) { //Perfectly inelastic collision
			s1.u = (s1.u + s2.u) / 2;
			s2.u = s1.u;
		}
		else { //Perfectly elastic collision
			r = r / distance;
			double u1_before = dot_product(s1.u, r);
			double u2_before = dot_product(s2.u, r);
			if (u1_before - u2_before < 0.0) { // if Solids move to each other
				double m1 = s1.rho * s1.V;
				double m2 = s2.rho * s2.V;
				double u1_after = (u1_before * (m1 - m2) + 2 * m2 * u2_before) / (m1 + m2);
				double u2_after = (u2_before * (m2 - m1) + 2 * m1 * u1_before) / (m1 + m2);
				s1.u += (u1_after - u1_before) * r;
				s2.u += (u2_after - u2_before) * r;
			}
		}
	}
	return result;
}

void Solids_move(std::list<Circle> &solidList, Param par) {

	///--------------collisions between particles-----------------
	for (auto one = solidList.begin(); one != solidList.end(); one++) {
		for (auto two = next(one); two != solidList.end(); two++) {
			if (Collide(*one, *two, par)) if (Debug) std::cout << "Collision detected" << std::endl;
		}
	}
	///-------------collision with walls--------------------------
	for (auto one = solidList.begin(); one != solidList.end(); one++) {
		double DistUpper = par.H - one->x[2];//<----distance to upper wall
		double DistLower = one->x[2];//<-------distance to lower wall
		if ((DistUpper < par.k_dist * one->r && one->u[2] > 0) ||
		    (DistLower < par.k_dist * one->r && one->u[2] < 0)) {
			if (Debug) std::cout << "Collision with wall detected" << std::endl;
			one->u[2] = -one->u[2];
		}
	}

	///
	for (auto it = solidList.begin(); it != solidList.end();) {

		if (it->n_moving <= par.N_step) {
			if (it->moving) {
				it->u_n    = it->u;
				it->omega_n = it->omega;
				it->x_n    = it->x;
				for (size_t k = 0; k < it->Nn; ++k) {
					it->Nodes[k].x_n = it->Nodes[k].x;
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
}

void Solids_zero_force(std::list<Circle>& Solids) {
	for (auto& it : Solids) {
		std::fill(it.f.begin(), it.f.end(), 0.0);
		std::fill(it.tau.begin(), it.tau.end(), 0.0);
		it.Fr_all = 0.;
		std::fill(it.f_L.begin(), it.f_L.end(), 0.0);
		for (size_t k = 0; k < it.Nn; ++k) {
			std::fill(it.Nodes[k].f.begin(), it.Nodes[k].f.end(), 0.0);
		}
	}
}

void Solids_velocity_new(std::list<Circle>& Solids, Param par) {
	for (auto& it : Solids) {
		if (it.moving) {
			it.u_s     = it.u;
			it.omega_s = it.omega;

			//GeomVec d_uc    = par.d_t * (it.integralV_du_dt  - it.f  ) / it.V / it.rho * 1;  // fluid density equals 1
			//GeomVec d_omega = par.d_t * (it.integralV_dur_dt - it.tau) / it.I / it.rho * 1;  // angular moment I is normalized with density

			//double alpha = 0.1;
			//it.u     = (1 - alpha) * it.u_s     + alpha * (it.u_n    + d_uc);
			//it.omega = (1 - alpha) * it.omega_s + alpha * (it.omega_n + d_omega);

			it.u     = it.u_n     - it.f_new   * par.d_t / (it.rho - 1) / it.V;  // fluid density equals 1
			it.omega = it.omega_n - it.tau     * par.d_t / (it.rho - 1) / it.I;  // angular moment I is normalized with density

			//it.u     = it.u_n     + it.F_hd   / it.rho / it.V * par.d_t;  // fluid density equals 1
			//it.omega = it.omega_n + it.tau_hd / it.rho / it.I * par.d_t;  // angular moment I is normalized with density

			it.x = it.x_n + 0.5 * (it.u_n + it.u_n) * par.d_t;
			for (size_t k = 0; k < it.Nn; ++k) {
				it.Nodes[k].x = rotate_Vector_around_vector(it.Nodes[k].x_n, 0.5 * (it.omega_n + it.omega_n) * par.d_t); //rotate
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


#include "stdafx.h"
#include "Output.h"

#pragma warning(disable : 4996)//for using <chrono>

void OutputPressure(Matrix p, int n, std::list<Circle> iList, Param par) {

	std::ofstream output;
	std::string filename = par.WorkDir + "solution_pressure" + std::to_string(n) + ".plt";

	output.open(filename);

	output << "title = " << '"' << "sample mesh" << '"' << std::endl;
	output << "Variables = x y p" << std::endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << par.N1 + 1 << ", j=" << par.N2 + 1 << ", f=point" << std::endl;
	output << "SolutionTime = " << n << std::endl;

	for (int j = 0; j <= par.N2; ++j) {
		for (int i = 0; i <= par.N1; ++i) {
			GeomVec xp = x_p(i, j, par);
			output << xp[1] << ' ' << xp[2] << ' ' << p[i][j] << std::endl;
		}
	}

	for (auto& solid : iList) {
		output << "zone T = circle" << ",  i=" << solid.Nn << ", f=point" << std::endl;
		output << "SolutionTime = " << n << std::endl;
		for (int i = 0; i < solid.Nn; ++i) {
			output << solid.Nodes[i].x[1] << " "
				<< solid.Nodes[i].x[2] << " "
				<< 0 << std::endl;
		}
	}

	output.close();
}


void OutputVelocity_U(Matrix u, int n, std::list<Circle> iList, Param par) {

	std::ofstream output;
	std::string filename = par.WorkDir + "solution_velocity_u" + std::to_string(n) + ".plt";

	output.open(filename);

	output << "title = " << '"' << "sample mesh" << '"' << std::endl;
	output << "Variables = x y u v" << std::endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << par.N1 << ", j=" << par.N2 + 1 << ", f=point" << std::endl;
	output << "SolutionTime = " << n << std::endl;

	for (int j = 0; j <= par.N2; ++j) {
		for (int i = 0; i < par.N1; ++i) {
			GeomVec xu = x_u(i, j, par);
			output << xu[1] << ' ' << xu[2] << ' ' << u[i][j] << ' ' << 0 << std::endl;
		}
	}

	for (auto& solid : iList) {
		output << "zone T = circle" << ",  i=" << solid.Nn << ", f=point" << std::endl;
		output << "SolutionTime = " << n << std::endl;
		for (int i = 0; i < solid.Nn; ++i) {
			output << solid.Nodes[i].x[1] << " "
				<< solid.Nodes[i].x[2] << " "
				<< solid.Nodes[i].f[1] << " "
				<< solid.Nodes[i].f[2] << " " << std::endl;
		}
	}

	output.close();
}


void OutputVelocity_V(Matrix v, int n, std::list<Circle> iList, Param par) {

	std::ofstream output;

	// Solution v

	std::string filename = par.WorkDir + "solution_velocity_v" + std::to_string(n) + ".plt";

	output.open(filename);

	output << "title = " << '"' << "sample mesh" << '"' << std::endl;
	output << "Variables = x y u v" << std::endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << par.N1 + 1 << ", j=" << par.N2 << ", f=point" << std::endl;
	output << "SolutionTime = " << n << std::endl;

	for (int j = 0; j < par.N2; ++j) {
		for (int i = 0; i <= par.N1; ++i) {
			GeomVec xv = x_v(i, j, par);
			output << xv[1] << ' ' << xv[2] << ' ' << 0 << ' ' << v[i][j] << std::endl;
		}
	}

	for (auto& solid : iList) {
		output << "zone T = circle" << ",  i=" << solid.Nn << ", f=point" << std::endl;
		output << "SolutionTime = " << n << std::endl;
		for (int i = 0; i < solid.Nn; ++i) {
			output << solid.Nodes[i].x[1] << " "
				<< solid.Nodes[i].x[2] << " "
				<< solid.Nodes[i].f[1] << " "
				<< solid.Nodes[i].f[2] << " " << std::endl;
		}
	}

	output.close();
}

void Output(Matrix p, Matrix u, Matrix v, Matrix Fx, Matrix Fy, int n, std::list<Circle> iList, Param par) {

	std::ofstream output;
	std::string filename = par.WorkDir + "step" + std::to_string(n) + ".plt";

	output.open(filename);

	output << "title = " << '"' << "sample mesh" << '"' << std::endl;
	output << "Variables = x y p u v fx fy tx ty" << std::endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << par.N1 + 1 << ", j=" << par.N2 + 1 << ", f=point" << std::endl;
	output << "SolutionTime = " << n << std::endl;

	for (int j = 0; j <= par.N2; ++j) {
		for (int i = 0; i <= par.N1; ++i) {
			GeomVec xp = x_p(i, j, par);

			double u_, v_, Fx_, Fy_;
			if (i == 0) { u_ = u[0][j];	Fx_ = Fx[0][j]; }
			else if (i == par.N1) { u_ = u[par.N1 - 1][j];	Fx_ = Fx[par.N1 - 1][j]; }
			else { u_ = (u[i - 1][j] + u[i][j]) * 0.5;	Fx_ = (Fx[i - 1][j] + Fx[i][j]) * 0.5; }
			if (j == 0) { v_ = v[i][0];	Fy_ = Fy[i][0]; }
			else if (j == par.N2) { v_ = v[i][par.N2 - 1]; ;	Fy_ = Fy[i][par.N2 - 1]; }
			else { v_ = (v[i][j - 1] + v[i][j]) * 0.5;	Fy_ = (Fy[i][j - 1] + Fy[i][j]) * 0.5; }
			output << xp[1] << " " << xp[2] << " " << p[i][j] << " " << u_ << " " << v_ << " " << Fx_ << " " << Fy_ << " " << 0 << " " << 0 << std::endl;
		}
	}

	for (auto& solid : iList) {
		output << "zone T = circle" << ",  i=" << solid.Nn + 2 << ", f=point" << std::endl;
		output << "SolutionTime = " << n << std::endl;

		output << solid.xc[1] << " "
		       << solid.xc[2] << " "
		       << 0 << " "
		       << solid.uc[1] << " "
		       << solid.uc[2] << " "
		       << solid.f[1] << " "
		       << solid.f[2] << " "
		       << solid.F_hd[1] << " "
		       << solid.F_hd[2] << " "
		       << std::endl;
		for (int i = 0; i < solid.Nn; ++i) {

			output << solid.Nodes[i].x[1] << " "
			       << solid.Nodes[i].x[2] << " "
			       << solid.Nodes[i].p << " "
			       << solid.Nodes[i].uf[1] << " "
			       << solid.Nodes[i].uf[2] << " "
			       << solid.Nodes[i].f[1] << " "
			       << solid.Nodes[i].f[2] << " "
			       << solid.Nodes[i].t[1] << " "
			       << solid.Nodes[i].t[2] << " "
			       << std::endl;

		}
		output << solid.Nodes[0].x[1] << " "
		       << solid.Nodes[0].x[2] << " "
		       << solid.Nodes[0].p << " "
		       << solid.Nodes[0].uf[1] << " "
		       << solid.Nodes[0].uf[2] << " "
		       << solid.Nodes[0].f[1] << " "
		       << solid.Nodes[0].f[2] << " "
		       << solid.Nodes[0].t[1] << " "
		       << solid.Nodes[0].t[2] << " "
		       << std::endl;
	}

	output.close();
}


void CreateDirectory(fs::path directory) {


	try
	{
		if (exists(directory)) {
			//std::cout << "This folder is already exist! \nRemove all files inside? (Y/N)" << std::endl;
			//char answer = std::cin.get();
			//if ( (answer == 'Y') ||  (answer == 'y') ) {
			//	fs::remove_all(directory);
			//	fs::create_directory(directory);
			//}
			//else {
			//	std::cout << "Start program? (Y/N)" << std::endl;
			//	std::cin >> answer;
			//	if ((answer != 'Y') && (answer != 'y'))exit(0);
			//}
		}
		else {
			try {
				fs::create_directory(directory);
			}
			catch (const fs::filesystem_error& ex) {
				CreateDirectory(directory.parent_path());
				fs::create_directory(directory);
			}
		}
	}
	catch (const fs::filesystem_error& ex)
	{
		std::cout << ex.what() << std::endl << "Enter any key";
		std::cin.get();
	}
}
void SetLog(std::ostream& log, Param par) {

	log << "The IBM program starts.		";
	time_t t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());   // get time now
	log << ctime(&t) << std::endl;
	log << "The parameters are the following:" << std::endl;
	log << "Reynolds number               : Re  = " << par.Re << std::endl;
	log << "Channel length                : L   = " << par.L << std::endl;
	log << "Channel width                 : W   = " << par.H << std::endl;
	log << "Number of nodes on            : N1  = " << par.N1 << std::endl;
	log << "Number of nodes on            : N2  = " << par.N2 << std::endl;
	log << "Number of nodes for a particle: Nn  = " << par.Nn << std::endl;
	log << "Time step                     : tau = " << par.d_t << std::endl;
	/*log << "Force parameter alpha         : alpha = " << par.alpha_f << std::endl;
	log << "Force parameter beta          : beta  = " << par.beta_f << std::endl;*/
	log << "Tolerance for Zeidel method   : tol = " << par.Zeidel_eps << std::endl;
	log << "Step number for start of Solids adding  : AddSolids_start    = " << par.AddSolids_start << std::endl;
	log << "Step interval for Solids adding         : AddSolids_interval = " << par.AddSolids_interval << std::endl;
	log << "Number of added Solids                  : AddSolids_N        = " << par.AddSolids_N << std::endl;
	log << std::endl;

}

void PushLog(std::ostream& log, int n, double eps_u, double eps_v) {
	time_t t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()); // get time now
	std::string s_time = ctime(&t);
	s_time.erase(7, 1);
	s_time.erase(0, 4);
	s_time.erase(s_time.size() - 6, 5);
	log << "n = " << std::setw(6) << n << "\t eps_u = " << std::fixed << eps_u << "\t eps_v = " << std::fixed << eps_v << "\t" << s_time;
	std::cout << "n = " << std::setw(6) << n << "\t eps_u = " << std::fixed << eps_u << "\t eps_v = " << std::fixed << eps_v << "\t" << s_time;
}


void Output_dp(Matrix dp, int n, Param par) {


	std::ofstream output;
	std::string filename = par.WorkDir + "dp" + std::to_string(n) + ".plt";

	output.open(filename);

	output << "title = " << '"' << "sample mesh" << '"' << std::endl;
	output << "Variables = x y dp" << std::endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << dp.size() << ", j=" << dp[0].size() << ", f=point" << std::endl;
	output << "SolutionTime = " << n << std::endl;

	for (int j = 0; j <= par.N2; ++j) {
		for (int i = 0; i <= par.N1; ++i) {
			GeomVec xp = x_p(i, j, par);
			output << xp[1] << " " << xp[2] << " " << dp[i][j] << std::endl;
		}
	}

	output.close();
}

void Output_c(Matrix c, int n, Param par) {

	std::ofstream output;
	std::string filename = par.WorkDir + "c" + std::to_string(n) + ".plt";

	output.open(filename);

	output << "title = " << '"' << "sample mesh" << '"' << std::endl;
	output << "Variables = x y c" << std::endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << c.size() << ", j=" << c[0].size() << ", f=point" << std::endl;
	output << "SolutionTime = " << n << std::endl;

	for (int j = 0; j < par.N2; ++j) {
		for (int i = 0; i < par.N1; ++i) {
			GeomVec xc = x_c(i, j, par);
			output << xc[1] << " " << xc[2] << " " << c[i][j] << std::endl;
		}
	}

	output.close();
}

void MakeHibernationFile(int n, Param& par, std::list<Circle>& solidList, Matrix& U_n, Matrix& V_n, Matrix& P) {
	std::ofstream output;
	std::string filename = par.WorkDir + "hibernation.txt";
	output.open(filename);
	output << "n = " << n << std::endl;
	output << "par{" << std::endl;
	output << "Re = " << par.Re << std::endl;
	output << "L = " << par.L << std::endl;
	output << "H = " << par.H << std::endl;
	output << "N1 = " << par.N1 << std::endl;
	output << "N2 = " << par.N2 << std::endl;
	output << "d_t = " << par.d_t << std::endl;
	output << "Nn = " << par.Nn << std::endl;
	output << "rho = " << par.rho << std::endl;
	output << "r = " << par.r << std::endl;
	output << "output_step = " << par.output_step << std::endl;
	output << "N_max = " << par.N_max << std::endl;
	output << "N_Zeidel = " << par.N_Zeidel << std::endl;
	output << "Zeidel_eps = " << par.Zeidel_eps << std::endl;
	output << "InelasticCollision = " << par.InelasticCollision << std::endl;
	output << "k_dist = " << par.k_dist << std::endl;
	output << "AddSolids_N = " << par.AddSolids_N << std::endl;
	output << "AddSolids_start = " << par.AddSolids_start << std::endl;
	output << "AddSolids_interval = " << par.AddSolids_interval << std::endl;
	output << "BC = " << par.BC << std::endl;
	output << "output_step = " << par.output_step << std::endl;
	//output << "WorkDir = " << par.WorkDir << std::endl;
	output << "SolidName_max = " << par.SolidName_max << std::endl;
	output << "d_x = " << par.d_x << std::endl;
	output << "d_y = " << par.d_y << std::endl;
	output << "}" << std::endl;

	output << "U_n{" << std::endl;
	for (int i = 0; i < par.N1; i++) {
		for (int j = 0; j < par.N2 + 1; j++) output << U_n[i][j] << ' ';
		output << std::endl;
	}
	output << "}" << std::endl;

	output << "V_n{" << std::endl;
	for (int i = 0; i < par.N1 + 1; i++) {
		for (int j = 0; j < par.N2; j++) output << V_n[i][j] << ' ';
		output << std::endl;
	}
	output << "}" << std::endl;

	output << "P{" << std::endl;
	for (int i = 0; i < par.N1 + 1; i++) {
		for (int j = 0; j < par.N2 + 1; j++) output << P[i][j] << ' ';
		output << std::endl;
	}
	output << "}" << std::endl;

	output << "<Solidlist>" << std::endl;
	for (auto one = solidList.begin(); one != solidList.end(); one++) {
		output << "<Solid>" << std::endl;
		output << "moving = " << one->moving << std::endl;
		output << "xc = " << std::endl << one->xc << std::endl;
		output << "uc = " << std::endl << one->uc << std::endl;
		output << "uc_n = " << std::endl << one->uc_n << std::endl;
		output << "omega = " << std::endl << one->omega << std::endl;
		output << "omega_n = " << std::endl << one->omega_n << std::endl;
		output << "f = " << std::endl << one->f << std::endl;
		output << "Fr = " << one->Fr << std::endl;
		output << "Fr_all = " << one->Fr_all << std::endl;
		output << "F_hd = " << std::endl << one->F_hd << std::endl;
		output << "tau_hd = " << std::endl << one->tau_hd << std::endl;
		output << "S = " << one->S << std::endl;
		output << "tau = " << std::endl << one->tau << std::endl;
		output << "I = " << one->I << std::endl;
		output << "rho = " << one->rho << std::endl;
		output << "V = " << one->V << std::endl;
		output << "Nn = " << one->Nn << std::endl;
		output << "name = " << one->name << std::endl;
		output << "r = " << one->r << std::endl;
		output << "<Nodes>" << std::endl;
		for (int j = 0; j < par.Nn; j++) {
			output << "Node{" << std::endl;
			output << "x = " << std::endl << one->Nodes[j].x << std::endl;
			output << "uf = " << std::endl << one->Nodes[j].uf << std::endl;
			output << "f = " << std::endl << one->Nodes[j].f << std::endl;
			output << "f_tmp = " << std::endl << one->Nodes[j].f_tmp << std::endl;
			output << "n = " << std::endl << one->Nodes[j].n << std::endl;
			output << "Eps = " << std::endl << one->Nodes[j].Eps << std::endl;
			output << "p = " << one->Nodes[j].p << std::endl;
			output << "}" << std::endl;
		}
		output << "<\\Nodes>" << std::endl;
		output << "<\\Solid>" << std::endl;
	}
	output << "<\\Solidlist>" << std::endl;
	output.close();
}

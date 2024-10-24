#include "stdafx.h"
#include "Output.h"

#pragma warning(disable : 4996)//for using <chrono>

void Output_U(Matrix u, std::string filename, int n, Param par) {

	std::ofstream output;
	filename = par.WorkDir + "/" + filename + std::to_string(n) + ".plt";

	output.open(filename);

	output << "title = " << '"' << filename << '"' << std::endl;
	output << "Variables = i j x y u" << std::endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << u.size() << ", j=" << u[0].size() << ", f=point" << std::endl;
	output << "SolutionTime = " << n << std::endl;

	for (int j = 0; j < u[0].size(); ++j) {
	//for (int j = u[0].size()-1 ; j >=0; --j) {
		for (int i = 0; i < u.size(); ++i) {
			GeomVec xu = x_u(i, j, par);
			double U = u[i][j];
			if (par.BC == Taylor_Green) {
				if (j == 0) {
					xu = (x_u(i, 0, par) + x_u(i, 1, par))*0.5;
					U = (u[i][0] + u[i][1])*0.5;
				}
				if (j == par.N2_u - 1) {
					xu = (x_u(i, par.N2_u - 1, par) + x_u(i, par.N2_u - 2, par))*0.5;
					U = (u[i][par.N2_u - 1] + u[i][par.N2_u - 2])*0.5;
				}
			}
			output << i << ' ' << j << ' ' << xu[1] << ' ' << xu[2] << ' ' << U << std::endl;
		}
	}

	output.close();
}


void Output_V(Matrix v, std::string filename, int n, Param par) {

	std::ofstream output;

	// Solution v

	filename = par.WorkDir + "/" + filename + std::to_string(n) + ".plt";

	output.open(filename);

	output << "title = " << '"' << filename << '"' << std::endl;
	output << "Variables = i j x y v" << std::endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << v.size() << ", j=" << v[0].size() << ", f=point" << std::endl;
	output << "SolutionTime = " << n << std::endl;

	for (int j = 0; j < v[0].size(); ++j) {
	//for (int j = v[0].size() - 1; j >= 0; --j) {
		for (int i = 0; i < v.size(); ++i) {
			GeomVec xv = x_v(i, j, par);
			double V = v[i][j];
			if (par.BC == Taylor_Green) {
				if (i == 0) {
					xv = (x_v(0, j, par) + x_v(1, j, par))*0.5;
					V = (v[0][j] + v[1][j])*0.5;
				}
				if (i == par.N1_v - 1) {
					xv = (x_v(par.N1_v - 1, j, par) + x_v(par.N1_v - 2, j, par))*0.5;
					V = (v[par.N1_v - 1][j] + v[par.N1_v - 2][j])*0.5;
				}
			}
			output << i << ' ' << j << ' ' << xv[1] << ' ' << xv[2] << ' '  << V << std::endl;
		}
	}

	output.close();
}

void Output(Matrix p, Matrix u, Matrix v, Matrix Fx, Matrix Fy, int n, std::vector<Solid> iList, std::vector<Node> Nodes, Param par) {

	std::ofstream output;
	std::string filename = par.WorkDir + "/" + "step" + std::to_string(n) + ".plt";

	output.open(filename);
	output << std::setprecision(15);

	output << "title = " << '"' << filename << '"' << std::endl;
	output << "Variables = x y p u v fx fy" << std::endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << par.N1_p << ", j=" << par.N2_p << ", f=point" << std::endl;
	output << "SolutionTime = " << par.time << std::endl;

	for (int j = 0; j < par.N2_p; ++j) {
	//for (int j = par.N2_p - 1 ; j >= 0; --j) {
		for (int i = 0; i < par.N1_p; ++i) {
			GeomVec xp = x_p(i, j, par);

			double u_, v_, Fx_, Fy_;
			u_ = (u[i][j] + u[i + 1][j]) * 0.5;		Fx_ = (Fx[i][j] + Fx[i + 1][j]) * 0.5;
			v_ = (v[i][j] + v[i][j + 1]) * 0.5;		Fy_ = (Fy[i][j] + Fy[i][j + 1]) * 0.5;
			output << xp[1] << " " << xp[2] << " " << p[i][j] + par.grad_p_x * (par.L - xp[1]) << " " << u_ << " " << v_ << " " << Fx_ << " " << Fy_  << std::endl;
		}
	}

	for (auto& solid : iList) {
		output << "zone T = circle" << ",  i=" << solid.Nn + 2 << ", f=point" << std::endl;
		output << "SolutionTime = " << par.time << std::endl;

		output << solid.x_n[1] << " "
		       << solid.x_n[2] << " "
		       << 0 << " "
		       << solid.u[1] << " "
		       << solid.u[2] << " "
		       //<< solid.a_collide[1] << " "
		       //<< solid.a_collide[2] << " "
			   << 0 << " "
			   << 0 << " "
		       << std::endl;
		for (int k = 0; k < solid.Nn; ++k) {
			int Ind = solid.IndNodes[k];

			output << Nodes[Ind].x_n[1] + solid.x_n[1] << " "
			       << Nodes[Ind].x_n[2] + solid.x_n[2] << " "
			       << Nodes[Ind].p << " "
			       << Nodes[Ind].us[1] << " "
			       << Nodes[Ind].us[2] << " "
			       << Nodes[Ind].f_r_collide[1] << " "
			       << Nodes[Ind].f_r_collide[2] << " "
			       << std::endl;

		}
		int Ind = solid.IndNodes[0];
		output << Nodes[Ind].x_n[1] + solid.x_n[1] << " "
		       << Nodes[Ind].x_n[2] + solid.x_n[2] << " "
		       << Nodes[Ind].p << " "
		       << Nodes[Ind].us[1] << " "
		       << Nodes[Ind].us[2] << " "
		       << Nodes[Ind].f_r_collide[1] << " "
		       << Nodes[Ind].f_r_collide[2] << " "
		       << std::endl;
	}

	output.close();
}

bool Read_plt(std::string filename, Param &par, std::vector<Solid>& solidList) {

	std::ifstream input;
	std::string line;
	
	input.open(filename);

	if (input.is_open()) {
		while (getline(input, line)) { // read line from file to string $line$
			if (line.substr(0, 15) == "zone T = circle") {
				std::vector<Node> Nodes;
				Solid c(par); //Create circle
				
				//Get the number of nodes in particle
				//auto strings = split_string(line, ",");  // split line into parts by delimeter ","
				//int i = 1;
				//for (auto itr = strings.begin(); itr != strings.end(); itr++) {
				//	std::string PAR, VALUE;
				//	GetParValue(*itr, PAR, VALUE);
				//	if (boost::trim_copy(PAR) == "i") {
				//		c.Nn = stoi(VALUE);
				//		//std::cout << "Nn = " << c.Nn << std::endl;
				//	}
				//}

				getline(input, line); // SolutionTime =
				std::string PAR, VALUE;
				GetParValue(line, PAR, VALUE);
				par.N_step = stoi(VALUE);

				getline(input, line); // Center of the particle
				auto strings = split_string(line, " "); // split line into parts by delimeter " "
				int i = 1;
				for (auto itr = strings.begin(); itr != strings.end(); itr++) {
					//std::cout << *itr << std::endl;
					if (i == 1) c.x_n[1]  = stod(*itr);
					if (i == 2) c.x_n[2]  = stod(*itr);
					if (i == 4) c.u_n[1]  = stod(*itr);
					if (i == 5) c.u_n[2]  = stod(*itr);
					if (i == 6) c.f[1]    = stod(*itr);
					if (i == 7) c.f[2]    = stod(*itr);
					i++;
				}

				solidList.push_back(c);

				//std::cin.ignore();
			}
		}
		input.close();
		return true;
	}
	else {
		return false;
	}

}

void CreateDir(fs::path directory) {


	try
	{
		if (exists(directory)) {

		}
		else {
			try {
				fs::create_directory(directory);
			}
			catch (const fs::filesystem_error& ex) {
				CreateDir(directory.parent_path());
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


void history_init(std::string WorkDir, std::string file, boundary_conditions BC) {
	std::ofstream output;
	std::string filename = WorkDir + "/" + file + ".plt";

	output << "title = " << '"' << "history" << '"' << std::endl;
	if (BC == box) {
		output.open(filename);
		output << "Variables = t h_average zero zero" << std::endl;
	}
	else if (BC == Taylor_Green || BC == Lamb_Oseen || BC == Line_Vortex) {
		output.open(filename);
		output << "Variables = t error_P error_U error_V" << std::endl;
	}
}

void history_log(std::string WorkDir, std::string file, double t, double var1, double var2, double var3) {
	std::ofstream output;
	std::string filename = WorkDir + "/" + file + ".plt";
	output.open(filename, std::ios::app);

	output << std::setprecision(8);
	output << t << "   "
	       << var1 << "   "
	       << var2 << "   "
	       << var3 << "   "
	       << std::endl;
}


void Output_P(Matrix P, std::string filename, int n, Param par) {


	std::ofstream output;
	filename = par.WorkDir + "/" + filename + std::to_string(n) + ".plt";

	output.open(filename);

	output << std::setprecision(15);
	output << "title = " << '"' << filename << '"' << std::endl;
	output << "Variables = i j x y p" << std::endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << P.size() << ", j=" << P[0].size() << ", f=point" << std::endl;
	output << "SolutionTime = " << n << std::endl;

	for (int j = 0; j < par.N2_p; ++j) {
	//for (int j = par.N2_p-1; j >=0; --j) {
		for (int i = 0; i < par.N1_p; ++i) {
			GeomVec xp = x_p(i, j, par);
			output << i << " " << j << " " << xp[1] << " " << xp[2] << " " << P[i][j] << std::endl;
		}
	}

	output.close();
}

void Output_c(Matrix c, std::string filename, int n, Param par) {

	std::ofstream output;
	filename = par.WorkDir + "/" + filename + std::to_string(n) + ".plt";

	output.open(filename);

	output << "title = " << '"' << filename << '"' << std::endl;
	output << "Variables = x y c" << std::endl;
	output << "zone T=" << '"' << n << '"' << ",  i=" << c.size() << ", j=" << c[0].size() << ", f=point" << std::endl;
	output << "SolutionTime = " << n << std::endl;

	for (int j = 0; j < par.N2 + 1; ++j) {
		for (int i = 0; i < par.N1 + 1; ++i) {
			GeomVec xc = x_c(i, j, par);
			output << xc[1] << " " << xc[2] << " " << c[i][j] << std::endl;
		}
	}

	output.close();
}

void Output_Matrix(Matrix A, std::string WorkDir, std::string Variable, int n) {
	std::ofstream output;
	std::string filename = WorkDir + "/" + Variable + std::to_string(n) + ".plt";

	output.open(filename);

	output << std::setprecision(15);
	output << "Variables = i j " << Variable << std::endl;
	output << "Zone T=" << '"' << "Flow" << '"' << ",  I =  " << A[0].size() << ", J =  " << A.size() << ", Datapacking = Point" << std::endl;

	for (int i = 0; i < A.size(); ++i) {
		for (int j = 0; j < A[0].size(); ++j) {
		//for (int j = A[0].size()-1; j >= 0; --j) {
			output << i << ' ' << j << ' ' << A[i][j] << std::endl;
		}
	}

	output.close();
}

void Output_Matrix_mid(Matrix A, std::string WorkDir, std::string Variable, int n) {
	std::ofstream output;
	std::string filename = WorkDir + "/" + Variable + std::to_string(n) + ".plt";

	output.open(filename);

	output << std::setprecision(15);
	output << "Variables = i j " << Variable << std::endl;
	output << "Zone T=" << '"' << "Flow" << '"' << ",  I =  " << A[0].size() - 2 << ", J =  " << A.size() - 2 << ", Datapacking = Point" << std::endl;

	for (int i = 1; i < A.size() - 1; ++i) {
		for (int j = 1; j < A[0].size() - 1; ++j) {
			output << i << ' ' << j << ' ' << A[i][j] << std::endl;
		}
	}

	output.close();
}

void Output_2DArray(double* A, int Nx, int Ny, std::string WorkDir, std::string Variable, int n) {
	std::ofstream output;
	std::string filename = WorkDir + "/" + Variable + std::to_string(n) + ".plt";

	output.open(filename);

	output << std::setprecision(15);
	output << "Variables = i j " << Variable << std::endl;
	output << "Zone T=" << '"' << "Array" << '"' << ",  I =  " << Ny << ", J =  " << Nx << ", Datapacking = Point" << std::endl;

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
		//for (int j = Ny - 1; j >= 0; --j) {
			//for (int j = Ny-1; j >= 0; --j) {
			output << i << ' ' << j << ' ' << A[i + j*Nx] << std::endl;
		}
	}

	output.close();
}


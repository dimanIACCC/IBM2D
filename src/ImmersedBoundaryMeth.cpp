#include "stdafx.h"
#include "CalculateForce.h"
#include "Calculate_u_p.h"
#include "Output.h"
#include "PredictVel.h"
#include "Testing.h"
#include "BodyOfProgram.h"





#pragma warning(disable : 4996)//for using <chrono>
#pragma warning(disable : 4244)//for GetInfluenceArea

void Awake(int& n0, Param& par, std::list<Circle>& solidList, Matrix& U_n, Matrix& V_n, Matrix& P);

int main(int argc, char *argv[]) {

	fs::path WorkDir = L"\Result\\";
	std::list<Circle> solidList; // list of immersed solids

	for (int i = 1; i < argc; i++) {
		std::string line, PAR, VALUE;
		line = (std::string)argv[i];
		GetParValue(line, PAR, VALUE);
		if (PAR == "-d") {
			DoTesting();
			std::cout << "Start main program? (Y/N)" << std::endl;
			char ch = std::cin.get();
			if ((ch != 'Y') || (ch != 'y')) return 0;

		}
		else if (PAR == "-hibernation"|| PAR == "-h") {
			if (VALUE.size() > 0) WorkDir = VALUE + '\\';
			Matrix U_n, V_n, P;
			Param par;
			int n;
			par.WorkDir = WorkDir.string();
			Awake(n, par, solidList, U_n, V_n, P);
			BodyOfProgram(par, solidList, U_n, V_n, P, n+1,false,false);
			return 0;
		}
		else if (PAR == "-dir") if (VALUE.size() > 0) WorkDir = VALUE + '/';

	}

	CreateDirectory(WorkDir);
	CreateDirectory(WorkDir.string() + "/Solids");


	Param par(WorkDir.string(), "input.txt");					// create the variable which contains parameters according to input data

	Read_Solids(par.WorkDir + "Solids.txt", solidList, par);	// read Solids from file and write them into list of solids

	CreateMatrix(U_n, par.N1, par.N2 + 1);						// creation matrices for velocity
	CreateMatrix(V_n, par.N1 + 1, par.N2);						// and pressure
	CreateMatrix(P, par.N1 + 1, par.N2 + 1);					//
	ApplyInitialData(U_n, P, par,solidList);					// Applying initial conditions for velocity and pressure according to 
																// boundary conditions


	BodyOfProgram(par, solidList, U_n, V_n, P);					// start solver


	std::cout << "The End" << std::endl;
	getchar();

	return 0;
}

void Awake(int& n0, Param& par, std::list<Circle>& solidList, Matrix& U_n, Matrix& V_n, Matrix& P) {
	std::ifstream hibernation_source;
	std::string line;
	std::string PAR, VALUE;
	char ch;
	
	hibernation_source.open(par.WorkDir+"hibernation.txt");
	if (hibernation_source.is_open()) {
		while (getline(hibernation_source, line)) { // read line from file to string $line$
			GetParValue(line, PAR, VALUE);
			if (PAR == "n")				n0 = stoi(VALUE);
			else if (line == "par{") {
				while (line != "}") {
					getline(hibernation_source, line);
					if (line == "}") break;
					std::string PAR, VALUE;
					GetParValue(line, PAR, VALUE);
					if (VALUE.size() > 0) {
						if (PAR == "Re")				par.Re = stod(VALUE);
						else if (PAR == "L")            par.L = stod(VALUE);
						else if (PAR == "H")            par.H = stod(VALUE);
						else if (PAR == "N1")           par.N1 = stoi(VALUE);
						else if (PAR == "N2")           par.N2 = stoi(VALUE);
						else if (PAR == "d_t")          par.d_t = stod(VALUE);
						else if (PAR == "Nn")           par.Nn = stoi(VALUE);
						else if (PAR == "rho")          par.rho = stod(VALUE);
						else if (PAR == "r")            par.r = stod(VALUE);
						else if (PAR == "output_step")  par.output_step = stoi(VALUE);
						else if (PAR == "N_max")        par.N_max = stoi(VALUE);
						else if (PAR == "N_Zeidel")     par.N_Zeidel = stoi(VALUE);
						else if (PAR == "Zeidel_eps")   par.Zeidel_eps = stod(VALUE);
						else if (PAR == "InelasticCollision")   par.InelasticCollision = bool(stoi(VALUE));
						else if (PAR == "k_dist")				par.k_dist = (stod(VALUE));
						else if (PAR == "AddSolids_N")          par.AddSolids_N = stoi(VALUE);
						else if (PAR == "AddSolids_start")      par.AddSolids_start = stoi(VALUE);
						else if (PAR == "AddSolids_interval")   par.AddSolids_interval = stoi(VALUE);
						else if (PAR == "BC")                   par.BC = string_to_BC(VALUE);
						else if (PAR == "output_step")          par.output_step = stod(VALUE);
						else if (PAR == "SolidName_max")        par.SolidName_max = stoi(VALUE);
						else if (PAR == "d_x")                  par.d_x = stod(VALUE);
						else if (PAR == "d_y")                  par.d_y = stod(VALUE);
						else    std::cout << "unknown parameter into par." << PAR << std::endl;
					}
					else {
						std::cout << "par: no value inputed" << std::endl;
					}
				}
			}
		    else if (line == "U_n{") {
				U_n.resize(par.N1);
				for (int i = 0; i < par.N1; i++)U_n[i].resize(par.N2 + 1);
				for (int i = 0; i < par.N1; i++)
					for (int j = 0; j < par.N2 + 1; j++) hibernation_source >> U_n[i][j];


				if (hibernation_source >> ch && ch != '}') std::cout << "Some troubles in U_n" << std::endl;
			}
			else if (line == "V_n{") {
				V_n.resize(par.N1 + 1);
				for (int i = 0; i < par.N1 + 1; i++)V_n[i].resize(par.N2);
				for (int i = 0; i < par.N1 + 1; i++)
					for (int j = 0; j < par.N2; j++) hibernation_source >> V_n[i][j];

				if (hibernation_source >> ch && ch != '}') std::cout << "Some troubles in V_n" << std::endl;
			}
			else if (line == "P{") {
				P.resize(par.N1 + 1);
				for (int i = 0; i < par.N1 + 1; i++)P[i].resize(par.N2 + 1);
				for (int i = 0; i < par.N1 + 1; i++)
					for (int j = 0; j < par.N2 + 1; j++) hibernation_source >> P[i][j];

				if (hibernation_source >> ch && ch != '}') std::cout << "Some troubles in P" << std::endl;
			}
			else if (line == "<Solidlist>")
				while (line != "<\\Solidlist>") {
					getline(hibernation_source, line);
					Circle c(0, 0, par);
					if (line == "<\\Solidlist>") break;
					if (line == "<Solid>") {
						while (line != "<\\Solid>") {
							getline(hibernation_source, line);
							if (line == "<\\Solid>") {
								solidList.push_back(c);
								break;
							}
							double x = par.L*0.1;
							double y = par.H*0.5;
							double ux = 0;
							double uy = 0;
							double omega = 0;
							double rho = par.rho;
							int Nn = par.Nn;
							bool moving = true;
							double r = par.r;

							GetParValue(line, PAR, VALUE);
							
								if (PAR == "moving")			 c.moving = bool(stoi(VALUE));
								else if (PAR == "xc")			 hibernation_source >> c.xc;//
								else if (PAR == "uc")			 hibernation_source >> c.uc;//
								else if (PAR == "uc_n")			 hibernation_source >> c.uc_n;//
								else if (PAR == "omega")		 hibernation_source >> c.omega;//
								else if (PAR == "omega_n")       hibernation_source >> c.omega_n;//
								else if (PAR == "f")			 hibernation_source >> c.f;//
								else if (PAR == "Fr")			 c.Fr = stod(VALUE);
								else if (PAR == "Fr_all")        c.Fr_all = stod(VALUE);
								else if (PAR == "F_hd")          hibernation_source >> c.F_hd;//
								else if (PAR == "tau_hd")        hibernation_source >> c.tau_hd;//
								else if (PAR == "S")			 c.S = stod(VALUE);
								else if (PAR == "tau")			 hibernation_source >> c.tau;//
								else if (PAR == "I")			 c.I = stod(VALUE);
								else if (PAR == "rho")			 c.rho = stod(VALUE);
								else if (PAR == "V")			 c.V = stod(VALUE);
								else if (PAR == "Nn")			 c.Nn = stoi(VALUE);
								else if (PAR == "name")			 c.name = stoi(VALUE);
								else if (PAR == "r")			 c.r = stod(VALUE);
								else if (PAR == "n_moving")			 c.n_moving = stoi(VALUE);
								else if (PAR == "<Nodes>") {
									while (line != "<\\Nodes>") {
										getline(hibernation_source, line);
										if (line == "<\\Nodes>") break;
										for (int j = 0; j < par.Nn; j++) {
											if (line == "Node{") {
												while (line != "}") {
													getline(hibernation_source, line);
													if (line == "}") { 
														getline(hibernation_source, line);
														break; 
													}
													GetParValue(line, PAR, VALUE);
													if (PAR == "x")					hibernation_source >> c.Nodes[j].x;
													else if (PAR == "uf")			hibernation_source >> c.Nodes[j].uf;
													else if (PAR == "f")			hibernation_source >> c.Nodes[j].f;
													else if (PAR == "f_tmp")		hibernation_source >> c.Nodes[j].f_tmp;
													else if (PAR == "n")			hibernation_source >> c.Nodes[j].n;
													else if (PAR == "Eps")			hibernation_source >> c.Nodes[j].Eps;
													else if (PAR == "p")			c.Nodes[j].p = stod(VALUE);
												}
											}
										}
									}
								}
						}
					}
				}
			}
		}
	else std::cout << "Hibernation file is not found. Please check the path\n The correct command is like -h=ResultFolder" << std::endl;
	std::cout << "End of awaking" << std::endl;

	
}
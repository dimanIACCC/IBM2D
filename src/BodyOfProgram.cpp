#include "stdafx.h"
#include "BodyOfProgram.h"

void BodyOfProgram(std::string WorkDir, int Re, bool TEST) {
	Param par;
	if (!TEST) {
		Param par(WorkDir + "input.txt"); // Construct Parameters using file input.txt
	}
	else
	{
		par.BC = u_infinity;
		par.Re = Re;
		par.alpha_f = -8e-5;
		par.beta_f = -2e3;
		if (Re > 43) par.N_max = 200e3;
	}

#pragma region SetMatrices
	CreateMatrix(U_n, par.N1, par.N2 + 1);
	CreateMatrix(U_new, par.N1, par.N2 + 1);
	CreateMatrix(Fx, par.N1, par.N2 + 1);
	CreateMatrix(V_n, par.N1 + 1, par.N2);
	CreateMatrix(V_new, par.N1 + 1, par.N2);
	CreateMatrix(B_v, par.N1 + 1, par.N2);
	CreateMatrix(Fy, par.N1 + 1, par.N2);
	CreateMatrix(P, par.N1 + 1, par.N2 + 1);
	ublas::matrix<Template> A_u(par.N1, par.N2 + 1);
	ublas::matrix<Template> A_v(par.N1 + 1, par.N2);
	Calculate_A(A_u, par, par.Re, Du);
	Calculate_A(A_v, par, par.Re, Dv);
#pragma endregion SetMatrices

	std::ofstream log;
	log.open(WorkDir + "log.txt", std::ios::out);
	SetLog(log, par);

	std::ofstream force;
	force.open(WorkDir + "force.plt", std::ios::out);
	force << "Variables = n, Fx, Fy" << std::endl;

	ApplyInitialData(U_n, P, par); // Applying initial data to velocity

	std::list<Circle> solidList; // list of immersed solids
	if (!TEST) {
		Read_Solids(WorkDir + "Solids.txt", solidList, par); // read Solids from file
	}
	else
	{
		Circle c(10, 3.5, 0, 0, 0, par.rho, par.Nn, false, 1);
		solidList.push_back(c);
	}
	Output(P, U_n, V_n, -1, solidList, par, WorkDir);

	for (int n = 0; n <= par.N_max; ++n) {
		if(!TEST) Add_Solids(solidList, n, par);

		CalculateForce(Fx, Fy, solidList, U_n, V_n, par);
		force << n << " " << Summ(Fx) << " " << Summ(Fy) << std::endl;
		Calculate_u_p(U_n, V_n, U_new, V_new, P, Fx, Fy, A_u, A_v, solidList, par, WorkDir);

		double eps_u = diff(U_n, U_new);
		double eps_v = diff(V_n, V_new);

		U_n = U_new;
		V_n = V_new;

		Solids_move(solidList, par);

		PushLog(log, n, eps_u, eps_v);
		log.flush();

		if (n % par.output_step == 0) {

			Output(P, U_new, V_new, n, solidList, par, WorkDir);
		}


		const double epsilon = 1e-7;
		if (eps_u < epsilon && eps_v < epsilon) {

			Output(P, U_new, V_new, n, solidList, par, WorkDir);

			if (TEST&& Re < 43) {
				 
					//the line is drawn with two points (Cd1,Nx1) & (Cd2,Nx2)
					double Cd1 = -3.15387;
					double Cd2 = -2.13; // Russel and Wang, 2003 

					int Nx1 = 101;
					int Nx2 = 301;
					if (par.N1 / (double)Nx1 > 1.5) {
						if (abs((solidList.front().f[1] - Cd1) / (Cd2 - Cd1) - (par.N1 - Nx1) / (Nx2 - Nx1)) < 0.5) {
							std::cout << "OK!" << std::endl;
							log << "OK!" << std::endl;
						}
						else {
							std::cout << "Not OK" << std::endl;
							log << "Not OK" << std::endl;
						}
						
					}
					else {
						double Cd_expected = Cd1 * Nx1 / (double)par.N1;
						if (abs(Cd_expected - solidList.front().f[1]) < 0.8) {
							std::cout << "OK!" << std::endl;
							log << "OK!" << std::endl;
						}
						else {
							std::cout << "Not OK" << std::endl;
							log << "Not OK" << std::endl;
						}
					}
				}
				else {
					std::cout << "That`s strange" << "Re =" << par.Re << std::endl;
					log << "That`s strange" << "Re =" << std::endl;
				} 

				break;
			}
			
		}

	
	log.close();
}
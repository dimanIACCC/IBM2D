#include "stdafx.h"
#include "BodyOfProgram.h"

void BodyOfProgram(Param par, std::list<Circle> solidList, bool TEST) {

#pragma region SetMatrices
	CreateMatrix(U_n, par.N1, par.N2 + 1);
	CreateMatrix(U_new, par.N1, par.N2 + 1);
	CreateMatrix(V_n, par.N1 + 1, par.N2);
	CreateMatrix(V_new, par.N1 + 1, par.N2);
	CreateMatrix(B_v, par.N1 + 1, par.N2);
	CreateMatrix(P, par.N1 + 1, par.N2 + 1);
	CreateMatrix(Fx, par.N1, par.N2 + 1);
	CreateMatrix(Fy, par.N1 + 1, par.N2);
	ublas::matrix<Template> A_u(par.N1, par.N2 + 1);
	ublas::matrix<Template> A_v(par.N1 + 1, par.N2);
	Calculate_A(A_u, par, par.Re, Du);
	Calculate_A(A_v, par, par.Re, Dv);
#pragma endregion SetMatrices

	std::ofstream log;
	log.open(par.WorkDir + "log.txt", std::ios::out);
	SetLog(log, par);

	std::ofstream force;
	force.open(par.WorkDir + "force.plt", std::ios::out);
	force << "Variables = n, Fx, Fy" << std::endl;

	ApplyInitialData(U_n, P, par); // Applying initial data to velocity

	Output(P, U_n, V_n, Fx, Fy, -1, solidList, par);

	for (int n = 0; n <= par.N_max; ++n) {
		Calculate_u_p(U_n, V_n, U_new, V_new, P, Fx, Fy, A_u, A_v, solidList, par);

		double eps_u = diff(U_n, U_new);
		double eps_v = diff(V_n, V_new);

		U_n = U_new;
		V_n = V_new;

		Solids_move(solidList, par,n);

		PushLog(log, n, eps_u, eps_v);
		log.flush();

		if (n % par.output_step == 0) {

			Output(P, U_new, V_new, Fx, Fy, n, solidList, par);
		}


		const double epsilon = 1e-7;
		if (eps_u < epsilon && eps_v < epsilon) {

			Output(P, U_new, V_new, Fx, Fy, n, solidList, par);

			if (TEST&& par.Re < 43) {
				 
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
				else if(TEST) {
					std::cout << "That`s strange" << "Re =" << par.Re << std::endl;
					log << "That`s strange" << "Re =" << std::endl;
				} 

				break;
			}
			
		}

	
	log.close();
}
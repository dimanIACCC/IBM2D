#include "stdafx.h"
#include "BodyOfProgram.h"

void BodyOfProgram(Param par, std::list<Circle> solidList, Matrix U_n, Matrix V_n, Matrix P, int n0, bool TEST, bool NeedNewLog) {
																		 
#pragma region SetMatrices 
	CreateMatrix(U_new, par.N1, par.N2 + 1);							// service region responsible for 
	CreateMatrix(V_new, par.N1 + 1, par.N2);							// creation matrices and one`s initial 
	CreateMatrix(B_v, par.N1 + 1, par.N2);								// filling
	CreateMatrix(Fx, par.N1, par.N2 + 1);								//	
	CreateMatrix(Fy, par.N1 + 1, par.N2);
	ublas::matrix<Template> A_u(par.N1, par.N2 + 1);
	ublas::matrix<Template> A_v(par.N1 + 1, par.N2);
	Calculate_A(A_u, par, par.Re, Du);
	Calculate_A(A_v, par, par.Re, Dv);
#pragma endregion SetMatrices


	std::ofstream log;														// creation output stream for writing log.
	if (NeedNewLog) {														// condition for creation new log file or 
		log.open(par.WorkDir + "log.txt", std::ios::out);					// continue writing if it`s after hibernation 
		SetLog(log, par);													// 
	}else
		log.open(par.WorkDir + "log.txt", std::ios::app);
	


	Output(P, U_n, V_n, Fx, Fy, -1, solidList, par);

	for (int n = n0; n <= par.N_max; ++n) {												// main cycle of time iterations 
		MakeHibernationFile(n-1, par, solidList, U_n, V_n, P);							// writting hibernation file for prior time step
		Calculate_u_p(U_n, V_n, U_new, V_new, P, Fx, Fy, A_u, A_v, solidList, par);		// calculation all hydrodynamics: velocities, pressure and forces
			

		double eps_u = diff(U_n, U_new);												// calculate maximum difference between current and prior velocity fields
		double eps_v = diff(V_n, V_new);												//

		U_n = U_new;
		V_n = V_new;
		
		Solids_move(solidList, par,n);												// moving solids if it necessary (checking it up inside) 
																					// and detection of collisions 

		PushLog(log, n, eps_u, eps_v);												// writting log
		log.flush();

		if (n % par.output_step == 0|| n<20) {

			Output(P, U_new, V_new, Fx, Fy, n, solidList, par);
		}


		const double epsilon = 1e-7;
		if (eps_u < epsilon && eps_v < epsilon && n>20) {
			
			Output(P, U_new, V_new, Fx, Fy, n, solidList, par);

				break;
			}
			
		}

	
	log.close();
}
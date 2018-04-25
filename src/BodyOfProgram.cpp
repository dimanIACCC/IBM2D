#include "stdafx.h"
#include "BodyOfProgram.h"

void BodyOfProgram(Param par, std::list<Circle> solidList, Matrix U_n, Matrix V_n, Matrix P, int n0, bool TEST, bool NeedNewLog) {
																		 
#pragma region SetMatrices 
	CreateMatrix(U_new, par.N1, par.N2 + 1);
	CreateMatrix(V_new, par.N1 + 1, par.N2);
	CreateMatrix(P_new, par.N1 + 1, par.N2 + 1);
	CreateMatrix(Fx_n  , par.N1, par.N2 + 1);
	CreateMatrix(Fy_n  , par.N1 + 1, par.N2);
	CreateMatrix(Fx_new, par.N1, par.N2 + 1);
	CreateMatrix(Fy_new, par.N1 + 1, par.N2);
	Template A_u;
	Template A_v;
	Calculate_A(A_u, par, par.Re, Du);
	Calculate_A(A_v, par, par.Re, Dv);
#pragma endregion SetMatrices


	std::ofstream log;														// creation output stream for writing log.
	if (NeedNewLog) {														// condition for creation new log file or 
		log.open(par.WorkDir + "log.txt", std::ios::out);					// continue writing if it`s after hibernation 
		SetLog(log, par);													// 
	}else
		log.open(par.WorkDir + "log.txt", std::ios::app);
	

	Multidirect_Forcing_Method(Fx_n, Fy_n, U_n, V_n, solidList, par);

	Output(P, U_n, V_n, Fx_n, Fy_n, -1, solidList, par);

	for (int n = n0; n <= par.N_max; ++n) {                                 // main cycle of time iterations
		MakeHibernationFile(n-1, par, solidList, U_n, V_n, P);              // writting hibernation file for prior time step
		Add_Solids(solidList, n, par);										// add solids if the conditions are fulfilled
		Calculate_u_p(U_n, V_n, U_new, V_new, P, P_new, Fx_n, Fy_n, A_u, A_v, solidList, par, n);  // calculate velocity and pressure at the new time step

		double eps_u = diff(U_n, U_new);												// calculate maximum difference between current and prior velocity fields
		double eps_v = diff(V_n, V_new);

		U_n = U_new;
		V_n = V_new;
		P   = P_new;
		
		//workaround for moving walls
		//for (auto& solid : solidList) {
		//	for (size_t i = 0; i < U_n.size(); i++)
		//		for (size_t j = 0; j < U_n[0].size(); j++)
		//			U_n[i][j] = U_n[i][j] - solid.uc[1];
		//	par.u_wall -= solid.uc[1];
		//	solid.uc[1] = 0;
		//}

		Solids_move(solidList, par, n);												// moving solids if it necessary (checking it up inside)
																					// and detection of collisions

		PushLog(log, n, eps_u, eps_v);                                              // writting log into log-file
		log.flush();

		if (n % par.output_step == 0|| n<1000) {

			Output(P, U_new, V_new, Fx_n, Fy_n, n, solidList, par);
		}


		const double epsilon = 1e-7;
		if (eps_u < epsilon && eps_v < epsilon && n>20) {
			
			Output(P, U_new, V_new, Fx_n, Fy_n, n, solidList, par);

				break;
			}
			
		}

	
	log.close();
}
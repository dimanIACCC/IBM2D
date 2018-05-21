#include "stdafx.h"
#include "BodyOfProgram.h"

void BodyOfProgram(Param par, std::list<Circle> solidList, Matrix U_n, Matrix V_n, Matrix P_n, int n0, bool NeedNewLog) {
																		 
#pragma region SetMatrices 
	CreateMatrix(U_new, par.N1_u, par.N2_u);
	CreateMatrix(V_new, par.N1_v, par.N2_v);
	CreateMatrix(P_new, par.N1 + 1, par.N2 + 1);
	CreateMatrix(Fx_n  , par.N1_u, par.N2_u);
	CreateMatrix(Fy_n  , par.N1_v, par.N2_v);
	CreateMatrix(Fx_new, par.N1_u, par.N2_u);
	CreateMatrix(Fy_new, par.N1_v, par.N2_v);
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

	Output(P_n, U_n, V_n, Fx_n, Fy_n, -1, solidList, par);

	for (par.N_step = n0; par.N_step <= par.N_max; ++par.N_step) {                                 // main cycle of time iterations
		MakeHibernationFile(par.N_step-1, par, solidList, U_n, V_n, P_n);              // writting hibernation file for prior time step
		Add_Solids(solidList, par);                                                     // add solids if the conditions are fulfilled

		Calculate_u_p(U_n, U_new, V_n, V_new, P_n, P_new, Fx_n, Fx_new, Fy_n, Fy_new, A_u, A_v, solidList, par);  // calculate velocity and pressure at the new time step

		double eps_u = diff(U_n, U_new);												// calculate maximum difference between current and prior velocity fields
		double eps_v = diff(V_n, V_new);

		U_n = U_new;
		V_n = V_new;
		P_n = P_new;
		
		//workaround for moving walls
		//for (auto& solid : solidList) {
		//	for (size_t i = 0; i < U_n.size(); i++)
		//		for (size_t j = 0; j < U_n[0].size(); j++)
		//			U_n[i][j] = U_n[i][j] - solid.uc_new[1];
		//	par.u_wall -= solid.uc_new[1];
		//	solid.uc_new[1] = 0;
		//}

		Solids_move(solidList, par);												// moving solids if it is necessary (checking it up inside)
																					// and detection of collisions

		PushLog(log, par.N_step, eps_u, eps_v);                                              // writting log into log-file
		log.flush();

		if (par.N_step % par.output_step == 0|| par.N_step < 1000) {
			Output(P_n, U_new, V_new, Fx_new, Fy_new, par.N_step, solidList, par);
			//OutputVelocity_U(U_new, par.N_step, par);
			//OutputVelocity_V(V_new, par.N_step, par);
		}

		const double epsilon = 1e-7;
		if (eps_u < epsilon && eps_v < epsilon && par.N_step > 1000) {
			break;
		}
	}

	
	log.close();
}
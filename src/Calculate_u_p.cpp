#include "Calculate_u_p.h"

#pragma warning(disable : 4996)//for using <chrono>
#pragma warning(disable : 4244)//for GetInfluenceArea

// calculate velocity $U_new, V_new$ and pressure $P$ at the new time step
// Solve the Navier-Stokes equation
// (U_new-U_n) / d_t + 0.5 * ( U_n \nabla U_n + U_s \nabla U_s ) = - nabla P + 1/Re * 0.5 * (\Delta U_n + \Delta U_new)
// Move U_new to left side and the rest terms to right side
// U_new / d_t - 1/Re * 0.5 \Delta U_new  = U_n / d_t - 0.5 * ( U_n \nabla U_n + U_s \nabla U_s ) + 1/Re * 0.5 \Delta U_n - nabla P
// Left-hand side operator is calculated in Template A_u and A_v (subroutines Calculate_A(...) and Operator_Ax(...))
// Right-hand side term is calculated in B_u and B_v (subroutine CalculateB() )
void Calculate_u_p(Matrix &U_n  , Matrix &V_n,
                   Matrix &U_new, Matrix &V_new,
                   Matrix &P    ,
                   Matrix &Fx   , Matrix &Fy,
                   Template A_u, Template A_v,
                   std::list<Circle> &solidList, Param par, int N_step) {

	CreateMatrix(U_s, par.N1, par.N2 + 1);
	CreateMatrix(V_s, par.N1 + 1, par.N2);
	CreateMatrix(P_Right, par.N1 + 1, par.N2 + 1);
	CreateMatrix(Delta_P, par.N1 + 1, par.N2 + 1);

	CreateMatrix(Exx, par.N1 + 1, par.N2 + 1);
	CreateMatrix(Eyy, par.N1 + 1, par.N2 + 1);
	CreateMatrix(Exy, par.N1, par.N2);

	int N_BiCGStab_u, N_BiCGStab_v, N_DeltaP;

	U_s = U_n;
	V_s = V_n;

	std::clock_t begin, end;

	// output of iterations information
	/*std::ofstream output;
	std::string filename = par.WorkDir + "/iterations" + std::to_string(N_step) + ".plt";
	output.open(filename);

	output << "title = iterations_step" << N_step << std::endl;
	output << "Variables = s ux uy omega f IntU tau IntUr" << std::endl;
	output << "Variables = s  N_BiCGStab_u  N_BiCGStab_v  N_DeltaP" << std::endl;*/

	int s_max = 10000;
	// start iterations for pressure and velocity to fulfill the conituity equation
	for (int s = 0; s <= s_max; ++s) {

		begin = std::clock();
		#pragma region Velocity
			CreateMatrix(B_u, par.N1, par.N2 + 1);
			CreateMatrix(B_v, par.N1 + 1, par.N2);

			B_u = CalculateB(U_n, V_n, U_s, V_s, P, par, Du);
			B_v = CalculateB(V_n, U_n, V_s, U_s, P, par, Dv);

			U_new = U_s;
			V_new = V_s;
			#pragma omp parallel sections num_threads(2)
			{
				#pragma omp section
				{
					BiCGStab(U_new, A_u, B_u, par, Du, N_BiCGStab_u);
				}
				#pragma omp section
				{
					BiCGStab(V_new, A_v, B_v, par, Dv, N_BiCGStab_v);
				}
			}

			end = std::clock();
			//std::cout << "time of Velocity     : " << end - begin << " " << std::endl;
		#pragma endregion Velocity

			//Output(P, U_new, V_new, Fx, Fy, s, solidList, par);

		#pragma region Force
			begin = std::clock();

			// apply force from immersed particles for several times to fulfill no-slip BC
			Multidirect_Forcing_Method(Fx, Fy, U_new, V_new, solidList, par);

			end = std::clock();
			//std::cout << "time of Force        : " << end - begin << " " << std::endl;
		#pragma endregion Force

			for (auto& it : solidList) {
				it.integrals(U_n, V_n, U_new, V_new, par);
			}


		#pragma region Pressure
			begin = std::clock();

			P_Right = Calculate_Press_Right(U_new, V_new, par);
			double Delta_P_max = Calculate_Press_correction(Delta_P, P_Right, par, N_DeltaP);
			double P_max = std::max(max(P), 1.e-4);
			double relax = 0.02 * std::max(pow(P_max / Delta_P_max, 0.5), 1.);

			std::cout << "s = " << s << ", delta_P / P = " << Delta_P_max / P_max << std::endl;

			end = std::clock();
			//std::cout << "time of Pressure     : " << end - begin << " " << std::endl;
		#pragma endregion Pressure

		#pragma region New P and U
			P += Delta_P * relax;

			for (size_t i = 0; i < U_new.size(); ++i) {
				for (size_t j = 0; j < U_new[0].size(); ++j) {
					U_new[i][j] -= relax * par.d_t * (Delta_P[i + 1][j] - Delta_P[i][j]) / par.d_x;
				}
			}

			for (size_t i = 0; i < V_new.size(); ++i) {
				for (size_t j = 0; j < V_new[0].size(); ++j) {
					V_new[i][j] -= relax * par.d_t * (Delta_P[i][j + 1] - Delta_P[i][j]) / par.d_y;
				}
			}

		#pragma endregion New P and U

		//Output_dp(Delta_P, s, par);
		//Output(P, U_new, V_new, Fx, Fy, s, solidList, par);
		


		// code for iterations solid u and omega iterations
		/*bool key_solid = false;
		for (auto& it : solidList) {

			double eps_uc = length(it.uc - it.uc_s) / (length(it.uc) + 1e-4);
			double eps_omega = length(it.omega - it.omega_s) / std::max( length(it.omega), length(it.uc)/it.r );

			double eps_max = 5e-6;
			if (eps_uc < eps_max && eps_omega < eps_max) key_solid = true;

			//output << s << "   " << it.uc[1] << "   " << it.uc[2] << "   " << it.omega[3] << "   " << it.f[1] << "   " << it.integralV_du_dt[1] << "   " << it.tau[3] << "   " << it.integralV_dur_dt[3] << "   " << std::endl;
		}*/


		// output of the iterations number
		// output << s << "   " << N_BiCGStab_u  << "   " << N_BiCGStab_v << "   " << N_DeltaP << std::endl;


		if (Delta_P_max / P_max < par.eps_P) {
			Solids_velocity_new(solidList, par);
			std::cout << "s iterations: " << s << std::endl;
			// if (key_solid == true)
			break;
		}

		U_s = U_new;
		V_s = V_new;
	}

	// Calculation of the strain rate, pressure and HydroDynamic (HD) force in Lagrange mesh
	deformation_velocity(U_new, V_new, Exx, Eyy, Exy, par);
	//Output_c(Exy, -1, par);
	Solids_deformation_velocity_pressure(solidList, Exx, Eyy, Exy, P, par);
	Solids_Force(solidList, par.Re);

}


// Apply initial data for velocity
void ApplyInitialData(Matrix &u, Matrix &p, Param par) {

	for (size_t i = 0; i < u.size(); ++i) {
		for (size_t j = 0; j < u[0].size(); ++j) {
			GeomVec xu = x_u(i, j, par);
			switch (par.BC) {
				case u_infinity: u[i][j] = 1.0; break;
				case u_inflow:   u[i][j] = 1.0 * ux_Poiseuille(xu[2], par.H); break;
				case periodical: u[i][j] = 1.0 * ux_Poiseuille(xu[2], par.H); break;
				default: std::cout << "ApplyInitialData: unknown BC" << std::endl;
			}
		}
	}

	for (size_t i = 0; i < p.size(); ++i) {
		for (size_t j = 0; j < p[0].size(); ++j) {
			GeomVec xp = x_p(i, j, par);
			p[i][j] = (par.L - xp[1]) * dpdx_Poiseuille(par.H, par.Re);
		}
	}

	size_t Nx = p.size();
	if (par.BC == periodical) {
		for (size_t j = 0; j < p[0].size(); ++j) {
			//p[0][j]      = par.L * dpdx_Poiseuille(par.H, par.Re);
			//p[1][j]      = par.L * dpdx_Poiseuille(par.H, par.Re);
			//p[Nx-2][j]   = 0;
			//p[Nx-1][j]   = 0;
		}
	}

}

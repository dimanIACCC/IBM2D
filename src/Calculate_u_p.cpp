#include "Calculate_u_p.h"

#pragma warning(disable : 4996)//for using <chrono>
#pragma warning(disable : 4244)//for GetInfluenceArea

// calculate velocity $U_new, V_new$ and pressure $P_n$ at the new time step
// Solve the Navier-Stokes equation
// (U_new-U_n) / d_t + 0.5 * ( U_n \nabla U_n + U_s \nabla U_s ) = - nabla P_n + 1/Re * 0.5 * (\Delta U_n + \Delta U_new)
// Move U_new to left side and the rest terms to right side
// U_new / d_t - 1/Re * 0.5 \Delta U_new  = U_n / d_t - 0.5 * ( U_n \nabla U_n + U_s \nabla U_s ) + 1/Re * 0.5 \Delta U_n - nabla P_n
// Left-hand side operator is calculated in Template A_u and A_v (subroutines Calculate_A(...) and Operator_Ax(...))
// Right-hand side term is calculated in B_u and B_v (subroutine CalculateB() )
void Calculate_u_p(Matrix &U_n   , Matrix &U_new,
                   Matrix &V_n   , Matrix &V_new,
                   Matrix &P_n   , Matrix &P_new,
                   Matrix &Fx_n  , Matrix &Fx_new,
                   Matrix &Fy_n  , Matrix &Fy_new,
                   Template A_u  , Template A_v,
                   std::list<Circle> &solidList, Param par) {

	CreateMatrix(U_s, par.N1_u, par.N2_u);
	CreateMatrix(V_s, par.N1_v, par.N2_v);
	CreateMatrix(P_Right, par.N1 + 1, par.N2 + 1);
	CreateMatrix(Delta_P, par.N1 + 1, par.N2 + 1);

	CreateMatrix(Exx, par.N1 + 1, par.N2 + 1);
	CreateMatrix(Eyy, par.N1 + 1, par.N2 + 1);
	CreateMatrix(Exy, par.N1, par.N2);

	int N_BiCGStab_u, N_BiCGStab_v, N_DeltaP;

	// if (N_step == 0)
	for (auto& solid : solidList) {
		solid.xc = solid.xc_n;
		solid.uc = solid.uc_n;
		solid.omega = solid.omega_n;
	}
	Multidirect_Forcing_Method(Fx_n, Fy_n, U_n, V_n, solidList, par);
	for (auto& solid : solidList) {
		solid.f_n = solid.f;
		solid.tau_n = solid.tau;
	}

	U_s = U_n;
	V_s = V_n;
	P_new = P_n;

	std::clock_t begin, end;

	// output of iterations information
	/*std::ofstream output;
	std::string filename = par.WorkDir + "/iterations" + std::to_string(N_step) + ".plt";
	output.open(filename);

	output << "title = iterations_step" << N_step << std::endl;
	output << "Variables = s ux uy omega f IntU tau IntUr" << std::endl;
	output << "Variables = s  N_BiCGStab_u  N_BiCGStab_v  N_DeltaP" << std::endl;*/

	int s_max = 1000;
	// start iterations for pressure and velocity to fulfill the conituity equation
	for (int s = 0; s <= s_max; ++s) {													// cycle while (delta_P / P > eps_P)
		begin = std::clock();
		#pragma region Velocity
			CreateMatrix(B_u, par.N1_u, par.N2_u);										// create matrix filled by 0
			CreateMatrix(B_v, par.N1_v, par.N2_v);										//

			B_u = CalculateB(U_n, V_n, U_s, V_s, P_n, P_new, par, Du);                           // RHS for Navier-Stokes non-linear equation
			B_v = CalculateB(V_n, U_n, V_s, U_s, P_n, P_new, par, Dv);

			U_new = U_s;
			V_new = V_s;
			#pragma omp parallel sections num_threads(2)
			{
				#pragma omp section
				{
					BiCGStab(U_new, A_u, B_u, par, Du, N_BiCGStab_u);                   // solving A_u * U_new = B_u
				}
				#pragma omp section
				{
					BiCGStab(V_new, A_v, B_v, par, Dv, N_BiCGStab_v);                   // solving A_v * V_new = B_v
				}
			}

			end = std::clock();
			//std::cout << "time of Velocity     : " << end - begin << " " << std::endl;
		#pragma endregion Velocity

			//Output(P_new, U_new, V_new, Fx, Fy, s, solidList, par);

		#pragma region Force
			begin = std::clock();

			// apply force from immersed particles for several times to fulfill no-slip BC
			for (auto& solid : solidList) {
				solid.xc = solid.xc_new;
				solid.uc = solid.uc_new;
				solid.omega = solid.omega_new;
			}
			Multidirect_Forcing_Method(Fx_new, Fy_new, U_new, V_new, solidList, par);
			for (auto& solid : solidList) {
				solid.f_new = solid.f;
				solid.tau_new = solid.tau;
			}

			end = std::clock();
			//std::cout << "time of Force        : " << end - begin << " " << std::endl;
		#pragma endregion Force

			for (auto& it : solidList) {
				it.integrals(U_n, V_n, U_new, V_new, par);
			}


		#pragma region Pressure
			begin = std::clock();

			P_Right = Calculate_Press_Right(U_new, V_new, par);                                     // calculating P_Right = 1/dt ( {U_i,j - U_i-1,j}/{h_x} + {V_i,j - V_i,j-1}/{h_x} )

			double Delta_P_max = Calculate_Press_correction(Delta_P, P_Right, par, N_DeltaP);       // Zeidel method for solving Poisson equation

			double P_max = std::max(max(P_new), 1.e-4);
			double relax = 0.05 * std::max(pow(P_max / Delta_P_max, 0.5), 1.);						// coefficient of relaxation


			std::cout << "s = " << s << ", delta_P / P = " << Delta_P_max / P_max << std::endl;

			end = std::clock();
			//std::cout << "time of Pressure     : " << end - begin << " " << std::endl;
		#pragma endregion Pressure

		#pragma region New P and U																
																										// area of correction pressure and velocity fields
			P_new += Delta_P * relax;                                                                   // correction of pressure

			if (par.BC == Taylor_Green) {
				//for (size_t i = 0; i < par.N1 + 1; ++i) {
				//	GeomVec x_D = x_p(i, 0     , par);
				//	GeomVec x_U = x_p(i, par.N2, par);
				//	P_new[i][0]      = 0.5 * (pow(cos(M_PI* x_D[2]),2) - pow(sin(M_PI * x_D[1]),2)) * exp(-4 * (M_PI*M_PI) / par.Re * par.d_t * par.N_step);      // D
				//	P_new[i][par.N2] = 0.5 * (pow(cos(M_PI* x_U[2]),2) - pow(sin(M_PI * x_U[1]),2)) * exp(-4 * (M_PI*M_PI) / par.Re * par.d_t * par.N_step);      // U
				//}
				//for (size_t j = 0; j < par.N2 + 1; ++j) {
				//	GeomVec x_L = x_p(0     , j, par);
				//	GeomVec x_R = x_p(par.N1, j, par);
				//	P_new[0][j]      = 0.5 * (pow(cos(M_PI* x_L[2]), 2) - pow(sin(M_PI * x_L[1]), 2)) * exp(-4 * (M_PI*M_PI) / par.Re * par.d_t * par.N_step);    // L
				//	P_new[par.N1][j] = 0.5 * (pow(cos(M_PI* x_R[2]), 2) - pow(sin(M_PI * x_R[1]), 2)) * exp(-4 * (M_PI*M_PI) / par.Re * par.d_t * par.N_step);    // R
				//}
				//P_new[par.N1 / 2][par.N2 / 2] = 0;
			}

			if (par.N_step < 2)
			P_n = P_new;

			for (size_t i = 1; i < U_new.size() - 1; ++i) {
				for (size_t j = 1; j < U_new[0].size() - 1; ++j) {
					U_new[i][j] -= relax * par.d_t * (Delta_P[i][j] - Delta_P[i - 1][j]) / par.d_x;		// correction of predicted velocity U_new
				}
			}

			for (size_t i = 1; i < V_new.size() - 1; ++i) {
				for (size_t j = 1; j < V_new[0].size() - 1; ++j) {
					V_new[i][j] -= relax * par.d_t * (Delta_P[i][j] - Delta_P[i][j - 1]) / par.d_y;		// correction of predicted velocity V_new
				}
			}

			//OutputVelocity_U(U_new, s, par);
			//OutputVelocity_V(V_new, s, par);

		#pragma endregion New P and U

		//Output_dp(Delta_P, s, par);
		//Output(P_new, U_new, V_new, Fx_n, Fy_n, s, solidList, par);
		


		// code for solid u and omega_new iterations
		/*bool key_solid = false;
		for (auto& it : solidList) {

			double eps_uc = length(it.uc_new - it.uc_s) / (length(it.uc_new) + 1e-4);
			double eps_omega = length(it.omega_new - it.omega_s) / std::max( length(it.omega_new), length(it.uc_new)/it.r );

			double eps_max = 5e-6;
			if (eps_uc < eps_max && eps_omega < eps_max) key_solid = true;

			//output << s << "   " << it.uc_new[1] << "   " << it.uc_new[2] << "   " << it.omega_new[3] << "   " << it.f[1] << "   " << it.integralV_du_dt[1] << "   " << it.tau[3] << "   " << it.integralV_dur_dt[3] << "   " << std::endl;
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
	Solids_deformation_velocity_pressure(solidList, Exx, Eyy, Exy, P_new, par);
	Solids_Force(solidList, par.Re);

}


// Apply initial data for velocity
void ApplyInitialData(Matrix &u, Matrix &v, Matrix &p, Param par) {

	for (size_t i = 0; i < u.size(); ++i) {
		for (size_t j = 0; j < u[0].size(); ++j) {
			GeomVec xu = x_u(i, j, par);
			switch (par.BC) {
				case u_infinity:   u[i][j] = 1.0; break;
				case u_inflow:     u[i][j] = 1.0 * ux_Poiseuille(xu[2], par.H); break;
				case periodical:   u[i][j] = 1.0 * ux_Poiseuille(xu[2], par.H); break;
				case Taylor_Green: u[i][j] = sin(M_PI * xu[1] / par.L) * cos(M_PI * xu[2] / par.H);  break;
				default: std::cout << "ApplyInitialData: unknown BC" << std::endl;
			}
		}
	}

	for (size_t i = 0; i < v.size(); ++i) {
		for (size_t j = 0; j < v[0].size(); ++j) {
			GeomVec xv = x_v(i, j, par);
			if (par.BC == Taylor_Green)
				v[i][j] = -par.H / par.L *sin(M_PI * xv[2] / par.H) * cos(M_PI * xv[1] / par.L);
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

		if (par.BC == Taylor_Green) {
			for (size_t i = 0; i < par.N1 + 1; ++i) {
				for (size_t j = 0; j < par.N2 + 1; ++j) {
					GeomVec x = x_p(i, j, par);
					double k1 = M_PI / par.L;
					double k2 = M_PI / par.H;
					p[i][j] = 0.5 * (pow(cos(k2 * x[2]) * k1 / k2, 2) - pow(sin(k1 * x[1]), 2)) * exp(-2 * (k1*k1 + k2*k2) / par.Re * par.d_t * par.N_step);
				}
			}
		}

	}

void Zero_velocity_in_Solids(Matrix &u, Param par, std::list<Circle> iList) {

	for (auto solid : iList) {
		for (int i = ((solid.xc_new[1] - solid.r) / par.d_x + 1); i < (solid.xc_new[1] + solid.r) / par.d_x + 1; i++)
			for (int j = (solid.xc_new[2] - solid.r) / par.d_y + 1; j < (solid.xc_new[2] + solid.r) / par.d_y + 1; j++) {
				double distance = sqrt(pow(par.d_x*i - solid.xc_new[1], 2) + pow(par.d_y*j - solid.xc_new[2], 2));
				if (distance <= solid.r) u[i][j] = 0;
				//p[i][j] /= solid.GetSweetness(i, j, par);
			}
	}

}

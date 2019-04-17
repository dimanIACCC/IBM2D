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
	CreateMatrix(P_RHS  , par.N1_p, par.N2_p);
	CreateMatrix(Delta_P, par.N1_p, par.N2_p);
	CreateMatrix(P      , par.N1_p, par.N2_p);

	CreateMatrix(Exx, par.N1_p, par.N2_p);
	CreateMatrix(Eyy, par.N1_p, par.N2_p);
	CreateMatrix(Exy, par.N1+1, par.N2+1);

	int N_BiCGStab_u, N_BiCGStab_v, N_DeltaP;

	for (auto& solid : solidList) {
		solid.x = solid.x_n;
		solid.u = solid.u_n;
		solid.omega = solid.omega_n;
	}

	U_s = U_n;
	V_s = V_n;
	P_new = P_n;

	std::clock_t begin, end, time_velocity, time_pressure, time_force;

	// output of iterations information
	std::ofstream output;
	std::string filename = par.WorkDir + "/iterations" + std::to_string(par.N_step) + ".plt";
	//output.open(filename);

	//output << "title = iterations_step" << par.N_step << std::endl;
	//output << "Variables = s ux uy omega f IntU tau IntUr" << std::endl;
	//output << "Variables = s  N_BiCGStab_u  N_BiCGStab_v  N_DeltaP time_velocity time_pressure time_force" << std::endl;

	/*if (par.BC == Taylor_Green) {
		int i = par.N1 / 2;
		int j = par.N2 / 2;
		double p_fix = 2 * exact_p(x_p(i, j, par), par, par.d_t * (par.N_step + 0.5)) - 2 * P_n[i][j];

		for (i = 0; i < P_new.size(); ++i)
		for (j = 0; j < P_new[0].size(); ++j)
			P_new[i][j] = P_n[i][j] + p_fix;
	}*/

	int s_max = 1000;
	// start iterations for pressure and velocity to fulfill the conituity equation
	for (int s = 0; s <= s_max; ++s) {													// cycle while (delta_P / P > eps_P)

		#pragma region Velocity
		begin = std::clock();
			CreateMatrix(B_u, par.N1_u, par.N2_u);										// create matrix filled by 0
			CreateMatrix(B_v, par.N1_v, par.N2_v);										//

			B_u = CalculateB(U_n, V_n, U_s, V_s, P_n, P_new, Fx_new, par, Du);                           // RHS for Navier-Stokes non-linear equation
			B_v = CalculateB(V_n, U_n, V_s, U_s, P_n, P_new, Fy_new, par, Dv);

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
		time_velocity = end - begin;
		#pragma endregion Velocity

		//Output(P_new, U_new, V_new, Fx_new, Fy_new, s, solidList, par);

		#pragma region Force
		begin = std::clock();

		// apply force from immersed particles for several times to fulfill no-slip BC
			Multidirect_Forcing_Method(Fx_new, Fy_new, U_new, V_new, solidList, par);
			for (auto& solid : solidList) {
				solid.f_new   = solid.f;
			}

		end = std::clock();
		time_force = end - begin;

		//Output_U(Fx_new, "U", par.N_step, par);
		//Output_V(Fy_new, "V", par.N_step, par);
		//std::cin.get();

		#pragma endregion Force

			for (auto& it : solidList) {
				it.integrals(U_n, V_n, U_new, V_new, par);
			}


		#pragma region Pressure
		begin = std::clock();

			P_RHS = Pressure_RHS(U_new, V_new, par);                                     // calculating P_RHS = 1/dt ( {U_i,j - U_i-1,j}/{h_x} + {V_i,j - V_i,j-1}/{h_x} )

			double Delta_P_max = Pressure_correction_solve(Delta_P, P_RHS, par, N_DeltaP);     // solve pressure correction equation  -Laplace Delta_P = P_RHS

			//Output_P(Delta_P, "helmholtz", s, par);
			//std::cin.get();

			double P_max = std::max(max(P_new), 1.e-14);
			double relax = std::min(0.05 * std::max(pow(P_max / Delta_P_max, 1), 1.), 1.0);						// coefficient of relaxation
			relax = 1.; // workaround

			std::cout  << "s = "
				<< std::setw(3) << s << ", delta_P / P = "
				<< std::setprecision(6) << std::scientific << Delta_P_max / P_max << ", relax = "
				<< std::fixed << relax << std::endl;

		end = std::clock();
		time_pressure = end - begin;
		#pragma endregion Pressure

		/*if (s < 3) {
			Output_Matrix(U_new  , par.WorkDir, "u_predict", s);
			Output_Matrix(Delta_P, par.WorkDir, "Delta_P"  , s);
		}*/

		//if (s == 0)
		//	Output(Delta_P, U_new, V_new, Fx_new, Fy_new, par.N_step, solidList, par);

		#pragma region New P and U																

			if ((par.BC == periodical || par.BC == u_inflow || par.BC == u_infinity && par.N_step < 10) || s > 50) {
				P_new += Delta_P * (relax / 2.);   // 1st order approximation for pressure
				P_n = P_new;                       //
			}
			else {
				P_new += Delta_P * relax;          // 2nd order approximation for pressure
			}

			for (size_t i = 1; i < U_new.size() - 1; ++i) {
				for (size_t j = 1; j < U_new[0].size() - 1; ++j) {
					U_new[i][j] -= relax * par.d_t * (Delta_P[i][j] - Delta_P[i - 1][j]) / par.d_x / 2.;		// correction of predicted velocity U_new
				}
			}

			for (size_t i = 1; i < V_new.size() - 1; ++i) {
				for (size_t j = 1; j < V_new[0].size() - 1; ++j) {
					V_new[i][j] -= relax * par.d_t * (Delta_P[i][j] - Delta_P[i][j - 1]) / par.d_y / 2.;		// correction of predicted velocity V_new
				}
			}

			/*if (s < 3) {
				Output_Matrix(U_new, par.WorkDir, "u_correct", s);
				Output_Matrix(P_new, par.WorkDir, "p_correct", s);
			}*/

			if (par.BC == Lamb_Oseen || par.BC == Line_Vortex || par.BC == Taylor_Green) {
				BC_exact_p(P_new, par, par.d_t * (par.N_step + 1));
			}

		#pragma endregion New P and U

		// code for solid u and omega_new iterations
		/*bool key_solid = false;
		for (auto& it : solidList) {

			double eps_uc = length(it.u - it.u_s) / (length(it.u) + 1e-4);
			double eps_omega = length(it.omega - it.omega_s) / std::max( length(it.omega), length(it.u)/it.r );

			double eps_max = 5e-6;
			if (eps_uc < eps_max && eps_omega < eps_max) key_solid = true;

			//output << s << "   " << it.u[1] << "   " << it.u[2] << "   " << it.omega[3] << "   " << it.f[1] << "   " << it.integralV_du_dt[1] << "   " << it.tau[3] << "   " << it.integralV_dur_dt[3] << "   " << std::endl;
		}*/


		//output of the iterations number
		output << s << "   " << N_BiCGStab_u  << "   " << N_BiCGStab_v << "   " << N_DeltaP << "  " << time_velocity << "  " << time_pressure << "  " << time_force << std::endl;

		if (Delta_P_max / P_max < par.eps_P) {
			Solids_velocity_new(solidList, par);
			if (par.BC == Lamb_Oseen) {
				for (auto& it : solidList) {
					it.omega[3] = Lamb_Oseen_omega(it.r, par.Re, par.d_t*(par.N_step + 1), par.Lamb_Oseen_r0);
				}
			}
			// std::cout << "s iterations: " << s << std::endl;
			//if (key_solid == true)
			break;

			//Output_Matrix(U_new, par.WorkDir, "u_finish", s);
			//Output_Matrix(V_new, par.WorkDir, "v_finish", s);
			//Output_Matrix(P_new, par.WorkDir, "p_finish", s);
		}

		U_s = U_new;
		V_s = V_new;
	}

	//std::cin.get();

	Fx_n = Fx_new;
	Fy_n = Fy_new;

	// Calculation of the strain rate, pressure and HydroDynamic (HD) force in Lagrange mesh
	deformation_velocity(U_new, V_new, Exx, Eyy, Exy, par);
	//Output_P(Exx, "Exx", -1, par);
	//Output_P(Eyy, "Eyy", -1, par);
	//Output_c(Exy, "Exy", -1, par);
	P = (P_n + P_new) * 0.5;
	Solids_deformation_velocity_pressure(solidList, Exx, Eyy, Exy, P, par);
	Solids_Force(solidList, par.Re);

}

void Zero_velocity_in_Solids(Matrix &u, Param par, std::list<Circle> iList) {

	for (auto solid : iList) {
		for (int i = ((solid.x[1] - solid.r) / par.d_x + 1); i < (solid.x[1] + solid.r) / par.d_x + 1; i++)
			for (int j = (solid.x[2] - solid.r) / par.d_y + 1; j < (solid.x[2] + solid.r) / par.d_y + 1; j++) {
				double distance = sqrt(pow(par.d_x*i - solid.x[1], 2) + pow(par.d_y*j - solid.x[2], 2));
				if (distance <= solid.r) u[i][j] = 0;
			}
	}

}

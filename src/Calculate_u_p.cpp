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
                   Matrix &P,
                   Matrix &Fx,
                   Matrix &Fy,
                   std::list<Circle> &solidList, Param par) {

	CreateMatrix(U_s, par.N1_u, par.N2_u);
	CreateMatrix(V_s, par.N1_v, par.N2_v);
	CreateMatrix(dFx, par.N1_u, par.N2_u);
	CreateMatrix(dFy, par.N1_v, par.N2_v);
	CreateMatrix(P_RHS  , par.N1_p, par.N2_p);
	CreateMatrix(Delta_P, par.N1_p, par.N2_p);

	CreateMatrix(Exx, par.N1_p, par.N2_p);
	CreateMatrix(Eyy, par.N1_p, par.N2_p);
	CreateMatrix(Exy, par.N1+1, par.N2+1);

	int N_BiCGStab_u=0, N_BiCGStab_v=0, N_DeltaP=0;

	for (auto& solid : solidList) {
		solid.x = solid.x_n;
		solid.u = solid.u_n;
		solid.omega = solid.omega_n;
	}

	U_s = U_n;
	V_s = V_n;

	std::clock_t begin, end, time_velocity, time_pressure, time_force;

	// output of iterations information
	std::ofstream output;
	std::string filename = par.WorkDir + "/iterations" + std::to_string(par.N_step) + ".plt";
	//output.open(filename);

	output << "title = iterations_step" << par.N_step << std::endl;
	//output << "Variables = s ux uy omega f IntU tau IntUr" << std::endl;
	output << "Variables = s  N_BiCGStab_u  N_BiCGStab_v  N_DeltaP time_velocity time_pressure time_force" << std::endl;


	CreateMatrix(RHS_u, par.N1_u, par.N2_u);
	CreateMatrix(RHS_v, par.N1_v, par.N2_v);

	Fx = Fx * 0.;
	Fy = Fy * 0.;

	for (auto& it : solidList) {
		std::fill(it.f_new.begin(), it.f_new.end(), 0.0);
		std::fill(it.tau_new.begin(), it.tau_new.end(), 0.0);
	}

	double Delta_P_max = 1.;
	double P_max = 1.;

	// start iterations for pressure and velocity to fulfill the continuity equation
	for (int s = 1; s <= par.s_max; ++s) {													// cycle while (delta_P / P > eps_P)

		#pragma region Velocity
		begin = std::clock();

			U_new = U_s;
			V_new = V_s;

			make_uv_RHS(RHS_u, RHS_v, U_n, V_n, U_s, V_s, P, Fx, Fy, par);
			predict_uv(U_new, V_new, RHS_u, RHS_v, par);


		end = std::clock();
		time_velocity = end - begin;
		#pragma endregion Velocity

		#pragma region Force
		begin = std::clock();

			CreateMatrix(U_f, par.N1_u, par.N2_u);
			CreateMatrix(V_f, par.N1_v, par.N2_v);

			U_f = U_new;
			V_f = V_new;

			// apply force from immersed particles for several times to fulfill no-slip BC
			Multidirect_Forcing_Method(dFx, dFy, U_f, V_f, solidList, par);

			if (par.IBM == 0) {
				U_new += dFx * (par.d_t);
				V_new += dFy * (par.d_t);
			}
			else if (par.IBM == 1) {
				double coef = std::max(1., 0.01*par.d_t * (par.ldxdx + par.ldydy) / par.Re);
				Fx += dFx * coef;
				Fy += dFy * coef;
			}

			for (auto& solid : solidList) {
				if (par.IBM == 0) {
					solid.f_new = solid.f;
					solid.tau_new = solid.tau;
				}
				else if (par.IBM == 1) {
					solid.f_new += solid.f;
					solid.tau_new += solid.tau;
				}
			}

		end = std::clock();
		time_force = end - begin;

		#pragma endregion Force

			for (auto& it : solidList) {
				it.integrals(U_n, V_n, U_new, V_new, par);
			}


		#pragma region Pressure
		begin = std::clock();

			P_RHS = Pressure_RHS(U_new, V_new, par);                                     // calculating P_RHS = 1/dt ( {U_i,j - U_i-1,j}/{h_x} + {V_i,j - V_i,j-1}/{h_x} )

			Delta_P_max = Pressure_correction_solve(Delta_P, P_RHS, par, N_DeltaP);     // solve pressure correction equation  -Laplace Delta_P = P_RHS

			P_max = std::max(max(P), 1.e-14);
			
			std::cout  << "s = "
				<< std::setw(3) << s << ", delta_P / P = "
				<< std::setprecision(6) << std::scientific << Delta_P_max / P_max << std::endl;

		end = std::clock();
		time_pressure = end - begin;
		#pragma endregion Pressure

		#pragma region New P and U																

			P += Delta_P;

			if (par.BC == Lamb_Oseen || par.BC == Line_Vortex) {
				BC_exact_p(P, par, par.d_t * (par.N_step + 0.5));
			}

			if (par.BC == Taylor_Green) {

				//Dirichlet BC for pressure
				//BC_exact_p(P, par, par.d_t * (par.N_step + 0.5));

				// Correct pressure for Neumann BC in Taylor-Green problem
				int i = par.N1 / 2;
				int j = par.N2 / 2;
				double p_fix = exact_p(x_p(i, j, par), par, par.d_t * (par.N_step + 0.5)) - P[i][j];

				for (i = 0; i < P.size(); ++i)
					for (j = 0; j < P[0].size(); ++j)
						P[i][j] = P[i][j] + p_fix;
			}

			if (par.BC == u_infinity || par.BC == u_inflow || par.BC == periodical || par.BC == box || par.BC == Taylor_Green) {

				// Up-Down BC
				for (size_t i = 0; i <= par.N1 + 1; ++i) {
					P[i][0         ] = P[i][1     ];       // D
					P[i][par.N2 + 1] = P[i][par.N2];       // U
				}

				// Left-Right BC
				for (size_t j = 0; j <= par.N2 + 1; ++j) {
					P[0         ][j] = P[1     ][j];       // L
					P[par.N1 + 1][j] = P[par.N1][j];       // R
					if (par.BC == u_inflow || par.BC == u_infinity) {
						P[par.N1 + 1][j] = -P[par.N1][j];
					}
				}

			}

			// Periodical Left-Right BC
			if (par.BC == periodical) {
				for (size_t j = 0; j <= par.N2 + 1; ++j) {
					P[0         ][j] = P[par.N1][j];        // L
					P[par.N1 + 1][j] = P[1     ][j];        // R
				}
			}


			for (size_t i = 1; i < P.size() - 1; ++i) {
				for (size_t j = 1; j < P[0].size() - 1; ++j) {
					P[i][j] -= 0.5 * par.d_t / par.Re * (par.ldxdx * (Delta_P[i + 1][j] - 2.*Delta_P[i][j] + Delta_P[i - 1][j])
						                                +par.ldydy * (Delta_P[i][j + 1] - 2.*Delta_P[i][j] + Delta_P[i][j - 1]));
				}
			}

			for (size_t i = 1; i < U_new.size() - 1; ++i) {
				for (size_t j = 1; j < U_new[0].size() - 1; ++j) {
					U_new[i][j] -= par.d_t * (Delta_P[i][j] - Delta_P[i - 1][j]) / par.d_x;		// correction of predicted velocity U_new
				}
			}

			for (size_t i = 1; i < V_new.size() - 1; ++i) {
				for (size_t j = 1; j < V_new[0].size() - 1; ++j) {
					V_new[i][j] -= par.d_t * (Delta_P[i][j] - Delta_P[i][j - 1]) / par.d_y;		// correction of predicted velocity V_new
				}
			}

			U_s = U_new;
			V_s = V_new;


		#pragma endregion New P and U

		//output of the iterations number
		output << s << "   " << N_BiCGStab_u  << "   " << N_BiCGStab_v << "   " << N_DeltaP << "  " << time_velocity << "  " << time_pressure << "  " << time_force << std::endl;
		//std::cin.get();

		Solids_collide(solidList, par);
		Solids_velocity_new(solidList, par);

		if (par.BC == Lamb_Oseen) {
			for (auto& it : solidList) {
				it.omega[3] = Lamb_Oseen_omega(it.r, par.Re, par.d_t*(par.N_step + 1), par.Lamb_Oseen_r0);
			}
		}

		if (Delta_P_max / P_max < par.eps_P && s>2) break;

	}

	// Calculation of the strain rate, pressure and HydroDynamic (HD) force in Lagrange mesh
	deformation_velocity(U_new, V_new, Exx, Eyy, Exy, par);
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

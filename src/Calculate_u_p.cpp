#include "Calculate_u_p.h"

#pragma warning(disable : 4996)//for using <chrono>
#pragma warning(disable : 4244)//for GetInfluenceArea


void Calculate_u_p(Matrix &U_n  , Matrix &V_n,
                   Matrix &U_new, Matrix &V_new,
                   Matrix &P    ,
                   Matrix &Fx   , Matrix &Fy,
                   ublas::matrix<Template> A_u,
                   ublas::matrix<Template> A_v, std::list<Circle> &solidList, Param par) {

	CreateMatrix(U_s, par.N1, par.N2 + 1);
	CreateMatrix(V_s, par.N1 + 1, par.N2);
	CreateMatrix(Fx_tmp, par.N1, par.N2 + 1);
	CreateMatrix(Fy_tmp, par.N1 + 1, par.N2);
	CreateMatrix(P_Right, par.N1 + 1, par.N2 + 1);
	CreateMatrix(Delta_P, par.N1 + 1, par.N2 + 1);

	CreateMatrix(Exx, par.N1 + 1, par.N2 + 1);
	CreateMatrix(Eyy, par.N1 + 1, par.N2 + 1);
	CreateMatrix(Exy, par.N1, par.N2);

	U_s = U_n;
	V_s = V_n;

	std::clock_t begin, end;

	int f_max = 20;
	int s_max = 1000;
	for (int s = 0; s <= s_max; ++s) {

		begin = std::clock();
		#pragma region Velocity
			CreateMatrix(B_u, par.N1, par.N2 + 1);
			CreateMatrix(B_v, par.N1 + 1, par.N2);

			B_u = CalculateB(U_n, V_n, U_s, V_s, P, Fx, par, Du);
			B_v = CalculateB(V_n, U_n, V_s, U_s, P, Fy, par, Dv);

			U_new = U_s;
			V_new = V_s;
			#pragma omp parallel sections num_threads(2)
			{
				#pragma omp section
				{
					BiCGStab(U_new, A_u, B_u, par, Du);
				}
				#pragma omp section
				{
					BiCGStab(V_new, A_v, B_v, par, Dv);
				}
			}

			end = std::clock();
			//std::cout << "time of Velocity     : " << end - begin << " " << std::endl;
		#pragma endregion Velocity

		#pragma region Force
			begin = std::clock();

			Solids_zero_force(solidList);
			Fx = Fx * 0.;
			Fy = Fy * 0.;
			for (int f = 0; f <= f_max; ++f) {

				for (auto& it : solidList) {
					it.Fr = 0.;
					it.S = 0.;
				}

				//deformation_velocity(U_new, V_new, Exx, Eyy, Exy, par);
				//Output_c(Exy, s, par);

				CalculateForce(Fx_tmp, Fy_tmp, solidList, U_new, V_new, par);
				U_new += Fx_tmp * (par.d_t);
				V_new += Fy_tmp * (par.d_t);
				Fx += Fx_tmp;
				Fy += Fy_tmp;

				//Output(P, U_new, V_new, Fx, Fy, f, solidList, par);

			}

			Solids_velocity_new(solidList, par);

			Solids_Force(solidList, par.Re);

			end = std::clock();
			//std::cout << "time of Force        : " << end - begin << " " << std::endl;
		#pragma endregion Force

		#pragma region Pressure
			begin = std::clock();

			P_Right = Calculate_Press_Right(U_new, V_new, par);
			double Delta_P_max = Calculate_Press_correction(Delta_P, P_Right, par);
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

		if (s % par.output_step == 0) {
			//Output_dp(Delta_P, s, par);
			//Output(P, U_new, V_new, Fx, Fy, s, solidList, par);
			//OutputVelocity_V(V_new, s, solidList, par);
		}

		double eps_P = 1.e-1;
		if (Delta_P_max / P_max < eps_P) {

			std::cout << "s iterations: " << s << std::endl;
			break;
		}

		U_s = U_new;
		V_s = V_new;
	}

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

void ApplyInitialData(Matrix &u, Matrix &p, Param par, std::list<Circle> iList) {

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
	for (auto solid : iList) {
		for (int i = ((solid.xc[1] - solid.r) / par.d_x + 1); i < (solid.xc[1] + solid.r) / par.d_x + 1; i++)
			for (int j = (solid.xc[2] - solid.r) / par.d_y + 1; j < (solid.xc[2] + solid.r) / par.d_y + 1; j++) {
				double distance = sqrt(pow(par.d_x*i - solid.xc[1], 2) + pow(par.d_y*j - solid.xc[2], 2));
				if (distance <= solid.r) u[i][j] = 0;
				//p[i][j] /= solid.GetSweetness(i, j, par);
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

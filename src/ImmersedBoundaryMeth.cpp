#include "CalculateForce.h"
#include "Output.h"
#include "BiCGStab.h"
#include "Calculate_press.h"
#include "PredictVel.h"
#include "Testing.h"



#pragma warning(disable : 4996)//for using <chrono>
#pragma warning(disable : 4244)//for GetInfluenceArea

// functions

void SetLog(std::ostream &log, Param par);
void PushLog(std::ostream &log, int n, double eps_u, double eps_v);
void ApplyInitialData(Matrix &u, Matrix &p, Param par);

int sgn(double x);



int main(int argc, char *argv[]) {

	std::string WorkDir = "";
	for (int i = 1; i < argc; i++) {
		std::string line, PAR, VALUE;
		line = (std::string)argv[i];
		GetParValue(line, PAR, VALUE);
		if      (PAR == "-d") DoSomeTest();
		else if (PAR == "-dir") if (VALUE.size() > 0) WorkDir = VALUE + '/';
	}

	Param par(WorkDir + "input.txt"); // Construct Parameters using file input.txt

	#pragma region SetMatrices
	CreateMatrix(U_n, par.N1, par.N2 + 1);
	CreateMatrix(U_new, par.N1, par.N2 + 1);
	CreateMatrix(U_prev, par.N1, par.N2 + 1);
	CreateMatrix(B_u, par.N1, par.N2 + 1);
	CreateMatrix(Force_x, par.N1, par.N2 + 1);
	CreateMatrix(V_n, par.N1 + 1, par.N2);
	CreateMatrix(V_new, par.N1 + 1, par.N2);
	CreateMatrix(V_prev, par.N1 + 1, par.N2);
	CreateMatrix(B_v, par.N1 + 1, par.N2);
	CreateMatrix(Force_y, par.N1 + 1, par.N2);
	CreateMatrix(P, par.N1 + 1, par.N2 + 1);
	CreateMatrix(Delta_P, par.N1 + 1, par.N2 + 1);
	CreateMatrix(P_Right, par.N1 + 1, par.N2 + 1);
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

	ApplyInitialData(U_new, P, par); // Applying initial data to velocity
	U_n = U_new;
	U_prev = U_new;

	std::list<Circle> solidList; // list of immersed solids
	Read_Solids(WorkDir + "Solids.txt", solidList, par); // read Solids from file

	Output(P, U_new, V_new, -1, solidList, par, WorkDir);

	int n = 0; // iteration counter
	while (n <= par.N_max) {

		Add_Solids(solidList, n, par);

		CalculateForce(Force_x, Force_y, solidList, U_new, V_new, par);
		force << n << " " << Summ(Force_x) << " " << Summ(Force_y) << std::endl;

		#pragma region Prediction of velocity
		B_u = CalculateB(U_n, V_n, U_prev, V_prev, P, Force_x, par, Du);
		B_v = CalculateB(V_n, U_n, V_prev, U_prev, P, Force_y, par, Dv);

		#pragma omp parallel sections num_threads(2)
		{

			#pragma omp section
			{
				BiCGStab(U_new, par.N1, par.N2 + 1, A_u, B_u, par, Du);
			}
			#pragma omp section
			{
				BiCGStab(V_new, par.N1 + 1, par.N2, A_v, B_v, par, Dv);
			}

		}

		#pragma endregion Prediction of velocity

		P_Right = Calculate_Press_Right(U_n, V_n, par);

		double eps_p = Calculate_Press_correction(Delta_P, P_Right, par);
		// if (n % par.output_step == 0) Output_dp(Delta_P, n, par, WorkDir);


		double relax = 0.1;
		for (int i = 0; i < par.N1 + 1; ++i) {
			for (int j = 0; j < par.N2 + 1; ++j) {
				P[i][j] = P[i][j] + relax * Delta_P[i][j];
			}
		}
		if (par.BC == periodical) {
			size_t n1 = P.size();
			size_t n2 = P[0].size();
			for (size_t j = 0; j < n2; ++j) {
				P[0][j] = par.L * dpdx_Poiseuille(par.H, par.Re);
				P[1][j] = par.L * dpdx_Poiseuille(par.H, par.Re);
			}
		}

		for (int i = 1; i < par.N1 - 1; ++i) {
			for (int j = 1; j < par.N2; ++j) {
				U_new[i][j] -= par.d_t * (Delta_P[i + 1][j] - Delta_P[i][j]) / par.d_x;
			}
		}

		for (int j = 1; j < par.N2; ++j) {
			int i = par.N1 - 1;
			U_new[i][j] = U_new[i - 1][j];
		}

		for (int i = 1; i < par.N1 + 1; ++i) {
			for (int j = 1; j < par.N2 - 1; ++j) {
				V_new[i][j] -= par.d_t * (Delta_P[i][j + 1] - Delta_P[i][j]) / par.d_y;
			}
		}

		#pragma region calculating eps_u and eps_v
		double eps_u = 0.0;
		for (int i = 0; i < par.N1; ++i) {
			for (int j = 0; j < par.N2 + 1; ++j) {
				if (fabs(U_n[i][j] - U_new[i][j]) > eps_u) {
					eps_u = fabs(U_n[i][j] - U_new[i][j]);
				}

				U_prev[i][j] = U_n[i][j];
				U_n[i][j] = U_new[i][j];
			}
		}
		double eps_v = 0.0;
		for (int i = 0; i < par.N1 + 1; ++i) {
			for (int j = 0; j < par.N2; ++j) {
				if (fabs(V_n[i][j] - V_new[i][j]) > eps_v) {
					eps_v = fabs(V_n[i][j] - V_new[i][j]);
				}
				V_prev[i][j] = V_n[i][j];
				V_n[i][j] = V_new[i][j];
			}
		}
		#pragma endregion calculating eps_u and eps_v


		#pragma region collisions check
		///--------------collisions between particles-----------------
		for (auto one = solidList.begin(); one != solidList.end(); one++) {
			for (auto two = next(one); two != solidList.end(); two++) {
				if (Collide(*one, *two, par)) if (Debug) std::cout << "Collision detected" << std::endl;
			}
		}
		///-------------collision with walls--------------------------
		for (auto one = solidList.begin(); one != solidList.end(); one++) {
			double DistUpper = par.H - one->xc[2];//<----distance to upper wall
			double DistLower = one->xc[2];//<-------distance to lower wall
			if (DistUpper < par.k_dist * one->r || DistLower < par.k_dist * one->r) {
				if (Debug) std::cout << "Collision with wall detected" << std::endl;
				one->uc[2] = -one->uc[2];
			}
		}
		#pragma endregion collisions check


		for (auto it = solidList.begin(); it != solidList.end();) {

			it->move(par.d_t);

			//delete bodies which move 95% of length
			if (it->xc[1] > par.L*0.95) {
				solidList.erase(it++);
			}
			else {
				++it;
			}
		}

		PushLog(log, n, eps_u, eps_v);
		log.flush();

		if (n % par.output_step == 0) {
			Output(P, U_new, V_new, n, solidList, par, WorkDir);
		}


		const double epsilon = 1e-5;
		if (eps_u < epsilon && eps_v < epsilon) {
			Output(P, U_new, V_new, n, solidList, par, WorkDir);
			break;
		}


		++n;
	}
	log.close();
	std::cout << "Over" << std::endl;
	getchar();

	return 0;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SetLog(std::ostream& log, Param par) {

	log << "The IBM program starts.		";
	time_t t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());   // get time now
	log << ctime(&t) << std::endl;
	log << "The parameters are the following:" << std::endl;
	log << "Reynolds number               : Re  = " << par.Re << std::endl;
	log << "Channel length                : L   = " << par.L << std::endl;
	log << "Channel width                 : W   = " << par.H << std::endl;
	log << "Number of nodes on            : N1  = " << par.N1 << std::endl;
	log << "Number of nodes on            : N2  = " << par.N2 << std::endl;
	log << "Number of nodes for a particle: Nn  = " << par.Nn << std::endl;
	log << "Time step                     : tau = " << par.d_t << std::endl;
	log << "Force parameter alpha         : alpha = " << par.alpha_f << std::endl;
	log << "Force parameter beta          : beta  = " << par.beta_f << std::endl;
	log << "Tolerance for Zeidel method   : tol = " << par.Zeidel_eps << std::endl;
	log << std::endl;

}

void PushLog(std::ostream& log, int n, double eps_u, double eps_v) {
	time_t t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()); // get time now
	std::string s_time = ctime(&t);
	s_time.erase(7, 1);
	s_time.erase(0, 4);
	s_time.erase(s_time.size() - 6, 5);
	log       << "n = " << std::setw(6) << n << "\t eps_u = " << std::fixed << eps_u << "\t eps_v = " << std::fixed << eps_v << "\t" << s_time;
	std::cout << "n = " << std::setw(6) << n << "\t eps_u = " << std::fixed << eps_u << "\t eps_v = " << std::fixed << eps_v << "\t" << s_time;
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
		for (size_t j = 0; j < p[1].size(); ++j) {
			GeomVec xp = x_p(i, j, par);
			p[i][j] = 1. * (par.L - xp[1]) * dpdx_Poiseuille(par.H, par.Re);
		}
	}
}



int sgn(double x)
{
	(x >= 0) ? x = 1 : x = -1;
	return x;
}
